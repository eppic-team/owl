package proteinstructure;

import gnu.getopt.Getopt;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;

import tools.MySQLConnection;

public class AtomTypeScorer extends TypeScorer {

	private static final File DEFAULT_LIST_FILE = new File("/project/StruPPi/jose/emp_potential/cullpdb_pc20_res1.6_R0.25_d090728_chains1627.list");
	private static final double DEFAULT_CUTOFF = 4.0;
	private static final int DEFAULT_MIN_SEQ_SEP = 3;
	
	private static final int NUM_ATOM_TYPES = 167;

	/**
	 * Constructs an AtomTypeScorer by taking a list of structure ids (pdbCodes+pdbChainCodes),
	 * and parameters distance cutoff and minimum sequence separation used for calculating 
	 * the scoring matrix. 
	 * @param listFile the file with the list of PDB ids (pdbCode+pdbChainCode)
	 * @param cutoff the distance cutoff to be used as definition of contacts
	 * @param minSeqSep the minimum sequence separation to be used when counting type pairs
	 * @throws SQLException if can't establish connection to db server
	 */
	public AtomTypeScorer(File listFile, double cutoff, int minSeqSep) throws SQLException {
		
		this.scoringMethod = ScoringMethod.ATOMTYPE;
		
		this.structureIds = new ArrayList<String>();
		this.listFile = listFile;
		this.ct = null;
		this.cutoff = cutoff;
		this.minSeqSep = minSeqSep;
		
		this.numEntities = NUM_ATOM_TYPES;
		entityCounts = new int[numEntities];
		pairCounts = new int[numEntities][numEntities];
		
		this.conn = new MySQLConnection();
		
	}

	/**
	 * Constructs an AtomTypeScorer given a file with a scoring matrix. All values of the matrix 
	 * and some other necessary fields will be read from the file.
	 * @param scMatFile
	 * @throws IOException
	 * @throws FileFormatError
	 */
	public AtomTypeScorer(File scMatFile) throws IOException, FileFormatError  {
		this.scoringMethod = ScoringMethod.ATOMTYPE;
		readScMatFromFile(scMatFile);
	}
	
	private void initAtomMap() {
		types2indices = new HashMap<String, Integer>();
		int i=0;
		for (String resType:AAinfo.getAAs()) {
			for (String atom:AAinfo.getAtomsForCTAndRes("ALL", resType)) {
				types2indices.put(resType+atom,i);
				i++;				
			}
		}
		indices2types = new HashMap<Integer, String>();
		for (String resType:types2indices.keySet()) {
			indices2types.put(types2indices.get(resType), resType);
		}
	}
	
	@Override
	public void countPairs() throws SQLException, IOException {
		this.initAtomMap();

		for (String id:TemplateList.readIdsListFile(listFile)) {
			String pdbCode = id.substring(0,4);
			String pdbChainCode = id.substring(4,5);
			Pdb pdb = null;
			try {
				pdb = new PdbasePdb(pdbCode,DB,conn);
				pdb.load(pdbChainCode);
				if (!isValidPdb(pdb)) {
					System.err.println(id+" didn't pass the quality checks to be included in training set");
					continue;
				}
				System.out.println(id);
				this.structureIds.add(id);
				
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Couldn't find pdb "+pdbCode);
				continue;
			} catch (PdbLoadError e) {
				System.err.println("Couldn't load pdb "+pdbCode);
				continue;
			}
			AIGraph graph = pdb.getAllAtomGraph(cutoff);
			graph.restrictContactsToMinRange(minSeqSep);
			for (AIGNode node:graph.getVertices()) {
				// we have to avoid OXT atoms, they are not in the map
				if (node.getAtomName().equals("OXT")) continue;
				countEntity(types2indices.get(node.getParent().getResidueType()+node.getAtomName()));
			}
			for (AIGEdge edge:graph.getEdges()) {
				AIGNode inode = graph.getEndpoints(edge).getFirst();
				AIGNode jnode = graph.getEndpoints(edge).getSecond();
				// we have to avoid OXT atoms, they are not in the map
				if (inode.getAtomName().equals("OXT") || jnode.getAtomName().equals("OXT")) continue;
				countPair(types2indices.get(inode.getParent().getResidueType()+inode.getAtomName()),
						types2indices.get(jnode.getParent().getResidueType()+jnode.getAtomName()));
			}
		}
		this.totalStructures = structureIds.size();
	}

	@Override
	public double scoreIt(Pdb pdb) {
		AIGraph graph = pdb.getAllAtomGraph(this.cutoff);
		graph.restrictContactsToMinRange(minSeqSep);
		
		double totalScore = 0;
		for (AIGEdge edge:graph.getEdges()) {
			AIGNode inode = graph.getEndpoints(edge).getFirst();
			AIGNode jnode = graph.getEndpoints(edge).getSecond();
			if (inode.getAtomName().equals("OXT") || jnode.getAtomName().equals("OXT")) continue;
			int i = types2indices.get(inode.getParent().getResidueType()+inode.getAtomName());
			int j = types2indices.get(jnode.getParent().getResidueType()+jnode.getAtomName());
			if (j>=i) {
				totalScore+=scoringMat[i][j];
			} else {
				totalScore+=scoringMat[j][i];
			}
		}

		return (totalScore/graph.getEdgeCount());
	}
	
	public static void main(String[] args) throws Exception {
		String help = 
		"\nCompiles a scoring matrix based on atom types from a given file with a list of \n" +
		"pdb structures.\n" +
		"Usage:\n" +
		"AtomTypeScorer -o <output_matrix_file> [-l <list_file> -c <cutoff> -m <out>]\n"+
		"  -o <file>     : file to write the scoring matrix to\n" +
		"  -l <file>     : file with list of pdbCodes+pdbChainCodes to use as training set \n" +
		"                  the scoring matrix. Default is "+DEFAULT_LIST_FILE+"\n" +
		"  -c <float>    : distance cutoff for the atom contacts. Default: "+DEFAULT_CUTOFF+"\n" +
		"  -m <int>      : minimum sequence separation to consider a contact. Default: "+DEFAULT_MIN_SEQ_SEP+"\n";

		File listFile = DEFAULT_LIST_FILE;
		File scMatFile = null;
		double cutoff = DEFAULT_CUTOFF;
		int minSeqSep = DEFAULT_MIN_SEQ_SEP;

		Getopt g = new Getopt("AtomTypeScorer", args, "l:c:o:m:h?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'l':
				listFile = new File(g.getOptarg());
				break;
			case 'c':
				cutoff = Double.parseDouble(g.getOptarg());
				break;				
			case 'o':
				scMatFile = new File(g.getOptarg());
				break;
			case 'm':
				minSeqSep = Integer.parseInt(g.getOptarg());
				break;				
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		if (scMatFile==null) {
			System.err.println("Missing argument output matrix file (-o)");
			System.exit(1);
		}

		AtomTypeScorer sc = new AtomTypeScorer(listFile,cutoff,minSeqSep);
		sc.countPairs();
		sc.calcScoringMat();
		sc.writeScMatToFile(scMatFile,false);
		
		
	}


}

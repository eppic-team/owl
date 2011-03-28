package owl.decoyScoring;

import gnu.getopt.Getopt;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;

import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbasePdb;
import owl.core.structure.TemplateList;
import owl.core.structure.graphs.AIGNode;
import owl.core.structure.graphs.AIGraph;
import owl.core.util.FileFormatException;
import owl.core.util.MySQLConnection;



public class AtomCountScorer extends CountScorer {

	
	private static final File DEFAULT_LIST_FILE = new File("/project/StruPPi/jose/emp_potential/cullpdb_pc20_res1.6_R0.25_d090728_chains1627.list");
	private static final double DEFAULT_CUTOFF = 4.0;
	private static final int DEFAULT_MIN_SEQ_SEP = 0;
	
	private static final int NUM_ATOM_TYPES = 167;
	
	/**
	 * Constructs an AtomCountScorer by taking a list of structure ids (pdbCodes+pdbChainCodes),
	 * and parameters distance cutoff and minimum sequence separation used for 
	 * calculating the scoring matrix. 
	 * @param listFile the file with the list of PDB ids (pdbCode+pdbChainCode)
	 * @param cutoff the distance cutoff to be used as definition of contacts
	 * @param minSeqSep the minimum sequence separation to be used when counting type pairs
	 * @throws SQLException if can't establish connection to db server
	 */
	public AtomCountScorer(File listFile, double cutoff, int minSeqSep) throws SQLException {
		
		this.scoringMethod = ScoringMethod.ATOMCOUNT;

		this.types2indices = new HashMap<String, Integer>();
		this.indices2types = new HashMap<Integer, String>();

		this.listFile = listFile;
		this.structureIds = new ArrayList<String>();
		this.cutoff = cutoff;
		this.minSeqSep = minSeqSep;
		
		this.numCountBins = NUM_COUNT_BINS;
		this.numTypes = NUM_ATOM_TYPES;
		binCountsPerType = new int[numCountBins][numTypes];
		
		this.conn = new MySQLConnection();

	}

	/**
	 * Constructs an AtomCountScorer given a file with a scoring matrix. All values of the matrix 
	 * and some other necessary fields will be read from the file.
	 * @param scMatFile
	 * @throws IOException
	 * @throws FileFormatException
	 */
	public AtomCountScorer(File scMatFile) throws IOException, FileFormatException  {
		this.scoringMethod = ScoringMethod.ATOMCOUNT;
		readScMatFromFile(scMatFile);
	}
	
	@Override
	public void countNodes() throws SQLException, IOException {
		Scorer.initAtomMap(types2indices, indices2types);

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
				structureIds.add(id);
				
			} catch (PdbCodeNotFoundException e) {
				System.err.println("Couldn't find pdb "+pdbCode);
				continue;
			} catch (PdbLoadException e) {
				System.err.println("Couldn't load pdb "+pdbCode);
				continue;
			}
			AIGraph graph = pdb.getAllAtomGraph(cutoff);
			graph.restrictContactsToMinRange(minSeqSep);
			for (AIGNode node:graph.getVertices()) {
				// we have to avoid OXT atoms, they are not in the map
				if (node.getAtomName().equals("OXT")) continue;
				int nbrCount = graph.getNeighborCount(node);
				count(nbrCount,types2indices.get(node.getParent().getResidueType()+node.getAtomName()));
			}
		}
		this.totalStructures = structureIds.size();

	}

	@Override
	public double scoreIt(Pdb pdb) {
		AIGraph graph = pdb.getAllAtomGraph(this.cutoff);
		graph.restrictContactsToMinRange(minSeqSep);
		return scoreIt(graph);	
	}
	
	protected double scoreIt(AIGraph graph) {
		double totalScore = 0;
		for (AIGNode node:graph.getVertices()) {
			if (node.getAtomName().equals("OXT")) continue;
			int nbrCount = graph.getNeighborCount(node);
			int typeIdx = types2indices.get(node.getParent().getResidueType()+node.getAtomName());
			if (nbrCount<scoringMat.length) {
				totalScore+=scoringMat[nbrCount][typeIdx];
			} else {
				// it can occur that a certain bin is not in the matrix, 
				// that basically means that in the background data there was not
				// a single instance of that neighbor-count bin. Thus we consider is
				// as a TOO_FEW_COUNTS_SCORE
				totalScore+=TOO_FEW_COUNTS_SCORE;
			}
		}

		return (totalScore/graph.getVertexCount());		
	}
	
	public static void main(String[] args) throws Exception {
		String help = 
			"\nCompiles a scoring matrix based on residue counts from a given file with a list of \n" +
			"pdb structures.\n" +
			"Usage:\n" +
			"AtomCountScorer -o <output_matrix_file> [-l <list_file> -c <cutoff> -m <min_seq_sep>]\n"+
			"  -o <file>     : file to write the scoring matrix to\n" +
			"  -l <file>     : file with list of pdbCodes+pdbChainCodes to use as training set \n" +
			"                  the scoring matrix. Default is "+DEFAULT_LIST_FILE+"\n" +
			"  -c <float>    : distance cutoff for the contacts. Default: "+DEFAULT_CUTOFF+"\n" +
			"  -m <int>      : minimum sequence separation to consider a contact. Default: "+DEFAULT_MIN_SEQ_SEP+"\n";

			File listFile = DEFAULT_LIST_FILE;
			File scMatFile = null;
			double cutoff = DEFAULT_CUTOFF;
			int minSeqSep = DEFAULT_MIN_SEQ_SEP;

			Getopt g = new Getopt("AtomCountScorer", args, "l:c:o:m:h?");
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

		AtomCountScorer sc = new AtomCountScorer(listFile,cutoff,minSeqSep);
		sc.countNodes();
		sc.calcScoringMat();
		sc.writeScMatToFile(scMatFile,false);

		// testing reading of score matrix file
		//AtomCountScorer newSc = new AtomCountScorer(scMatFile);
		//newSc.writeScMatToFile(new File(scMatFile+".tmp"), false);
		
	}



}

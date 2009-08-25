package proteinstructure;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;

import tools.MySQLConnection;

public class AtomTypeScorer extends TypeScorer {

	private static final int NUM_ATOM_TYPES = 167;

	/**
	 * Constructs an AtomTypeScorer by taking a list of structure ids (pdbCodes+pdbChainCodes),
	 * and parameters distance cutoff and minimum sequence separation used for calculating 
	 * the scoring matrix. 
	 * @param structureIds the list of PDB ids (pdbCode+pdbChainCode)
	 * @param cutoff the distance cutoff to be used as definition of contacts
	 * @param minSeqSep the minimum sequence separation to be used when counting type pairs
	 * @throws SQLException if can't establish connection to db server
	 */
	public AtomTypeScorer(String[] structureIds, double cutoff, int minSeqSep) throws SQLException {
		
		this.structureIds = structureIds;
		this.ct = null;
		this.cutoff = cutoff;
		this.minSeqSep = minSeqSep;
		
		this.numEntities = NUM_ATOM_TYPES;
		entityCounts = new int[numEntities];
		pairCounts = new int[numEntities][numEntities];
		totalStructures = 0;
		
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
		readFromFile(scMatFile);
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
	
	/**
	 * Performs the counts of the atom type pairs and stores them in the internal arrays.
	 * Use subsequently {@link #calcScoringMat()} to compute the scoring matrix from counts arrays.
	 * @throws SQLException
	 */
	public void countAtomTypes() throws SQLException {
		this.initAtomMap();

		for (String id:structureIds) {
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
				totalStructures++;
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
	}

	@Override
	public double scoreIt(Pdb pdb, int minSeqSep) {
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
		double CUTOFF_ALL = 4;
		File listFile = new File("/project/StruPPi/jose/emp_potential/cullpdb_pc20_res1.6_R0.25_d090728_chains1627.list");
		File scMatFile = new File("/project/StruPPi/jose/emp_potential/scoremat.atom.cullpdb20");
		//File outScMatFile = new File("/project/StruPPi/jose/emp_potential/scoremat.atom.cullpdb20.tmp");
		String[] ids = TemplateList.readIdsListFile(listFile);
		AtomTypeScorer sc = new AtomTypeScorer(ids,CUTOFF_ALL,3);
		sc.countAtomTypes();
		sc.calcScoringMat();
		sc.writeScMatToFile(scMatFile,false);
		
//		AtomScorer scFromFile = new AtomScorer(scMatFile);


		
	}


}

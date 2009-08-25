package proteinstructure;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;

import tools.MySQLConnection;

public class ResTypeScorer extends TypeScorer {
	
	private static final int NUM_RES_TYPES = 20;
	

	/**
	 * Constructs a ResTypeScorer by taking a list of structure ids (pdbCodes+pdbChainCodes),
	 * and parameters contact type, distance cutoff and minimum sequence separation used for 
	 * calculating the scoring matrix. 
	 * @param structureIds the list of PDB ids (pdbCode+pdbChainCode)
	 * @param ct the contact type to be used as definition of contacts
	 * @param cutoff the distance cutoff to be used as definition of contacts
	 * @param minSeqSep the minimum sequence separation to be used when counting type pairs
	 * @throws SQLException if can't establish connection to db server
	 */
	public ResTypeScorer(String[] structureIds, String ct, double cutoff, int minSeqSep) throws SQLException {
		
		this.structureIds = structureIds;
		this.ct = ct;
		this.cutoff = cutoff;
		this.minSeqSep = minSeqSep;
		
		this.numEntities = NUM_RES_TYPES;
		entityCounts = new int[numEntities];
		pairCounts = new int[numEntities][numEntities];
		totalStructures = 0;
		
		this.conn = new MySQLConnection();
		
	}

	/**
	 * Constructs a ResTypeScorer given a file with a scoring matrix. All values of the matrix 
	 * and some other necessary fields will be read from the file.
	 * @param scMatFile
	 * @throws IOException
	 * @throws FileFormatError
	 */
	public ResTypeScorer(File scMatFile) throws IOException, FileFormatError  {
		readFromFile(scMatFile);
	}

	
	private void initResMap() {
		types2indices = new HashMap<String, Integer>();
		int i=0;
		for (String resType:AAinfo.getAAs()) {
			types2indices.put(resType,i);
			i++;
		}
		indices2types = new HashMap<Integer, String>();
		for (String resType:types2indices.keySet()) {
			indices2types.put(types2indices.get(resType), resType);
		}
	}

	/**
	 * Performs the counts of the residue type pairs and stores them in the internal arrays.
	 * Use subsequently {@link #calcScoringMat()} to compute the scoring matrix from counts arrays.
	 * @throws SQLException
	 */
	public void countResTypes() throws SQLException {
		this.initResMap();

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

			RIGraph graph = pdb.get_graph(ct, cutoff);
			graph.restrictContactsToMinRange(minSeqSep);
			for (RIGNode node:graph.getVertices()) {
				countEntity(types2indices.get(node.getResidueType()));
			}
			for (RIGEdge edge:graph.getEdges()) {
				String iRes = graph.getEndpoints(edge).getFirst().getResidueType();
				String jRes = graph.getEndpoints(edge).getSecond().getResidueType();
				countPair(types2indices.get(iRes),types2indices.get(jRes));
			}
		}
	}
	
	@Override
	public double scoreIt(Pdb pdb, int minSeqSep) {
		RIGraph graph = pdb.get_graph(this.ct, this.cutoff);
		graph.restrictContactsToMinRange(minSeqSep);
		
		double totalScore = 0;
		for (RIGEdge edge:graph.getEdges()) {
			String iRes = graph.getEndpoints(edge).getFirst().getResidueType();
			String jRes = graph.getEndpoints(edge).getSecond().getResidueType();
			int i = types2indices.get(iRes);
			int j = types2indices.get(jRes);
			if (j>=i) {
				totalScore+= scoringMat[i][j];
			} else {
				totalScore+= scoringMat[j][i];
			}
		}
		
		return (totalScore/graph.getEdgeCount());
	}

	public static void main(String[] args) throws Exception {
		String CT = "Cb";
		double CUTOFF = 8;
		File listFile = new File("/project/StruPPi/jose/emp_potential/cullpdb_pc20_res1.6_R0.25_d090728_chains1627.list");
		File scMatFile = new File("/project/StruPPi/jose/emp_potential/scoremat.res.cullpdb20");
		//File outScMatFile = new File("/project/StruPPi/jose/emp_potential/scoremat.res.cullpdb20.tmp");
		String[] ids = TemplateList.readIdsListFile(listFile);
		ResTypeScorer sc = new ResTypeScorer(ids,CT,CUTOFF,3);
		sc.countResTypes();
		sc.calcScoringMat();
		sc.writeScMatToFile(scMatFile,false);
		
//		ResScorer scFromFile = new ResScorer(scMatFile);
		
	}


}

package owl.decoyScoring;

import gnu.getopt.Getopt;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;

import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbasePdb;
import owl.core.structure.TemplateList;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatError;
import owl.core.util.MySQLConnection;



public class ResCountScorer extends CountScorer {

	private static final File DEFAULT_LIST_FILE = new File("/project/StruPPi/jose/emp_potential/cullpdb_pc20_res1.6_R0.25_d090728_chains1627.list");
	private static final double DEFAULT_CUTOFF = 8.0;
	private static final int DEFAULT_MIN_SEQ_SEP = 0;
	private static final String DEFAULT_CT = "Cb";

	private static final int NUM_RES_TYPES = 20;

	/**
	 * Constructs a ResCountScorer by taking a list of structure ids (pdbCodes+pdbChainCodes),
	 * and parameters contact type, distance cutoff and minimum sequence separation used for 
	 * calculating the scoring matrix. 
	 * @param listFile the file with the list of PDB ids (pdbCode+pdbChainCode)
	 * @param ct the contact type to be used as definition of contacts
	 * @param cutoff the distance cutoff to be used as definition of contacts
	 * @param minSeqSep the minimum sequence separation to be used when counting type pairs
	 * @throws SQLException if can't establish connection to db server
	 */
	public ResCountScorer(File listFile, String ct, double cutoff, int minSeqSep) throws SQLException {
		
		this.scoringMethod = ScoringMethod.RESCOUNT;
		
		this.types2indices = new HashMap<String, Integer>();
		this.indices2types = new HashMap<Integer, String>();

		this.listFile = listFile;
		this.structureIds = new ArrayList<String>();
		this.ct = ct;
		this.cutoff = cutoff;
		this.minSeqSep = minSeqSep;
		
		this.numCountBins = NUM_COUNT_BINS;
		this.numTypes = NUM_RES_TYPES;
		binCountsPerType = new int[numCountBins][numTypes];
		
		this.conn = new MySQLConnection();
		
	}

	/**
	 * Constructs a ResCountScorer given a file with a scoring matrix. All values of the matrix 
	 * and some other necessary fields will be read from the file.
	 * @param scMatFile
	 * @throws IOException
	 * @throws FileFormatError
	 */
	public ResCountScorer(File scMatFile) throws IOException, FileFormatError  {
		this.scoringMethod = ScoringMethod.RESCOUNT;
		readScMatFromFile(scMatFile);
	}

	@Override
	public void countNodes() throws SQLException, IOException {
		Scorer.initResMap(types2indices, indices2types);

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
				
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Couldn't find pdb "+pdbCode);
				continue;
			} catch (PdbLoadError e) {
				System.err.println("Couldn't load pdb "+pdbCode);
				continue;
			}

			RIGraph graph = pdb.getRIGraph(ct, cutoff);
			graph.restrictContactsToMinRange(minSeqSep);
			for (RIGNode node:graph.getVertices()) {
				int nbrCount = graph.getNeighborCount(node);
				count(nbrCount,types2indices.get(node.getResidueType()));
			}
		}
		this.totalStructures = structureIds.size();
	}

	@Override
	public double scoreIt(Pdb pdb) {
		RIGraph graph = pdb.getRIGraph(this.ct, this.cutoff);
		graph.restrictContactsToMinRange(minSeqSep);
		return scoreIt(graph);
	}
	
	protected double scoreIt(RIGraph graph) {
		double totalScore = 0;
		for (RIGNode node:graph.getVertices()) {
			String res = node.getResidueType();
			int nbrCount = graph.getNeighborCount(node);
			int typeIdx = types2indices.get(res);
			if (nbrCount<scoringMat.length) {
				totalScore+= scoringMat[nbrCount][typeIdx];
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
			"ResCountScorer -o <output_matrix_file> [-l <list_file> -c <cutoff> -m <min_seq_sep>]\n"+
			"  -o <file>     : file to write the scoring matrix to\n" +
			"  -l <file>     : file with list of pdbCodes+pdbChainCodes to use as training set \n" +
			"                  the scoring matrix. Default is "+DEFAULT_LIST_FILE+"\n" +
			"  -t <string>   : contact type. Default: "+DEFAULT_CT+"\n"+
			"  -c <float>    : distance cutoff for the contacts. Default: "+DEFAULT_CUTOFF+"\n" +
			"  -m <int>      : minimum sequence separation to consider a contact. Default: "+DEFAULT_MIN_SEQ_SEP+"\n";

			File listFile = DEFAULT_LIST_FILE;
			File scMatFile = null;
			double cutoff = DEFAULT_CUTOFF;
			int minSeqSep = DEFAULT_MIN_SEQ_SEP;
			String ct = DEFAULT_CT;

			Getopt g = new Getopt("ResCountScorer", args, "l:t:c:o:m:h?");
			int c;
			while ((c = g.getopt()) != -1) {
				switch(c){
				case 'l':
					listFile = new File(g.getOptarg());
					break;
				case 't':
					ct = g.getOptarg();
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

		ResCountScorer sc = new ResCountScorer(listFile,ct,cutoff,minSeqSep);
		sc.countNodes();
		sc.calcScoringMat();
		sc.writeScMatToFile(scMatFile,false);

		// testing reading of score matrix file
		//ResCountScorer newSc = new ResCountScorer(scMatFile);
		//newSc.writeScMatToFile(new File(scMatFile+".tmp"), false);
		
	}

}

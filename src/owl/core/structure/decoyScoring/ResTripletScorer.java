package owl.core.structure.decoyScoring;

import gnu.getopt.Getopt;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;

import owl.core.structure.FileFormatError;
import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbasePdb;
import owl.core.structure.RIGNode;
import owl.core.structure.RIGraph;
import owl.core.structure.TemplateList;
import owl.core.util.MySQLConnection;
import owl.core.util.Triplet;



public class ResTripletScorer extends TripletScorer {
	
	private static final File DEFAULT_LIST_FILE = new File("/project/StruPPi/jose/emp_potential/cullpdb_pc20_res1.6_R0.25_d090728_chains1627.list");
	private static final double DEFAULT_CUTOFF = 8.0;
	private static final int DEFAULT_MIN_SEQ_SEP = 3;
	private static final String DEFAULT_CT = "Cb";

	private static final int NUM_RES_TYPES = 20;
	

	/**
	 * Constructs a ResTripletScorer by taking a list of structure ids (pdbCodes+pdbChainCodes),
	 * and parameters contact type, distance cutoff and minimum sequence separation used for 
	 * calculating the scoring matrix. 
	 * @param listFile the file with a list of PDB ids (pdbCode+pdbChainCode)
	 * @param ct the contact type to be used as definition of contacts
	 * @param cutoff the distance cutoff to be used as definition of contacts
	 * @param minSeqSep the minimum sequence separation for a contact to be considered
	 * @throws SQLException if can't establish connection to db server
	 */
	public ResTripletScorer(File listFile, String ct, double cutoff, int minSeqSep) throws SQLException {
		
		this.scoringMethod = ScoringMethod.RESTRIPLET;
		
		this.types2indices = new HashMap<String, Integer>();
		this.indices2types = new HashMap<Integer, String>();

		this.structureIds = new ArrayList<String>();
		this.listFile = listFile;
		this.ct = ct;
		this.cutoff = cutoff;
		this.minSeqSep = minSeqSep;
		
		this.numEntities = NUM_RES_TYPES;
		entityCounts = new int[numEntities];
		tripletCounts = new int[numEntities][numEntities][numEntities];
		
		this.conn = new MySQLConnection();
		
	}

	/**
	 * Constructs a ResTripletScorer given a file with a scoring matrix. All values of the matrix 
	 * and some other necessary fields will be read from the file.
	 * @param scMatFile
	 * @throws IOException
	 * @throws FileFormatError
	 */
	public ResTripletScorer(File scMatFile) throws IOException, FileFormatError  {
		this.scoringMethod = ScoringMethod.RESTRIPLET;
		readScMatFromFile(scMatFile);
	}

	
	@Override
	public void countTriplets() throws SQLException, IOException {
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
				this.structureIds.add(id);
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
				countEntity(types2indices.get(node.getResidueType()));
			}
			//for (Triplet<RIGNode> trip:graph.getTriangles()) {
			for (Triplet<RIGNode> trip:graph.getTriplets()) {
				String iRes = trip.getFirst().getResidueType();
				String jRes = trip.getSecond().getResidueType();
				String kRes = trip.getThird().getResidueType();
				countTriplet(types2indices.get(iRes),types2indices.get(jRes),types2indices.get(kRes));
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
		//Collection<Triplet<RIGNode>> triplets = graph.getTriangles();
		Collection<Triplet<RIGNode>> triplets = graph.getTriplets();
 		for (Triplet<RIGNode> trip:triplets) {
			String iRes = trip.getFirst().getResidueType();
			String jRes = trip.getSecond().getResidueType();
			String kRes = trip.getSecond().getResidueType();
			int[] ind = {types2indices.get(iRes),types2indices.get(jRes),types2indices.get(kRes)};
			Arrays.sort(ind);
			totalScore+= scoringMat[ind[0]][ind[1]][ind[2]];
		}
		
		return (totalScore/triplets.size());
	}
	
	
	public static void main(String[] args) throws Exception {
		String help = 
			"\nCompiles a scoring matrix based on residue types from a given file with a list of \n" +
			"pdb structures.\n" +
			"Usage:\n" +
			"ResTripletScorer -o <output_matrix_file> [-l <list_file> -c <cutoff> -m <min_seq_sep>]\n"+
			"  -o <file>     : file to write the scoring matrix to\n" +
			"  -l <file>     : file with list of pdbCodes+pdbChainCodes to use as training set \n" +
			"                  the scoring matrix. Default is "+DEFAULT_LIST_FILE+"\n" +
			"  -t <string>   : contact type. Default: "+DEFAULT_CT+"\n"+
			"  -c <float>    : distance cutoff for the contacts. Default: "+DEFAULT_CUTOFF+"\n" +
			"  -m <int>      : minimum sequence separation to consider a contact. Default: "+DEFAULT_MIN_SEQ_SEP+"\n" +
			"  -w            : write counts matrix as well as scoring matrix \n";

			File listFile = DEFAULT_LIST_FILE;
			File scMatFile = null;
			double cutoff = DEFAULT_CUTOFF;
			int minSeqSep = DEFAULT_MIN_SEQ_SEP;
			String ct = DEFAULT_CT;
			boolean writeCounts = false;

			Getopt g = new Getopt("ResTripletScorer", args, "l:t:c:o:m:wh?");
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
				case 'w':
					writeCounts = true;
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


		ResTripletScorer sc = new ResTripletScorer(listFile,ct,cutoff,minSeqSep);
		sc.countTriplets();
		sc.calcScoringMat();
		sc.writeScMatToFile(scMatFile,writeCounts);

	}

}

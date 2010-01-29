package proteinstructure.decoyScoring;

import gnu.getopt.Getopt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.AIGEdge;
import proteinstructure.AIGNode;
import proteinstructure.AIGraph;
import proteinstructure.FileFormatError;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbasePdb;
import proteinstructure.TemplateList;
import tools.MySQLConnection;

public class DistanceScorer extends Scorer {
	
	private static final File DEFAULT_LIST_FILE = new File("/project/StruPPi/jose/emp_potential/cullpdb_pc20_res1.6_R0.25_d090728_chains1627.list");
	private static final int DEFAULT_MIN_SEQ_SEP = 3;


	private static final int TOO_FEW_COUNTS_THRESHOLD = 20;
	private static final double TOO_FEW_COUNTS_SCORE = 0.0;
	
	private static final int NUM_ATOM_TYPES = 167;

	// our convention is for the index k to be for distances between DISTANCE_BINS[k] and DISTANCE_BINS[k-1]. k=0 will have all distances below DISTANCE_BINS[0]
	private static final double[] DISTANCE_BINS = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0}; 
	
	private static final String ATOM_DISTANCE_COUNTS_TABLE = "atom_distance_counts";
	
	
	private double[] distanceBins;
	
	private int[][][] pairCounts;		// the counts of the pairs, indices: 1st i type, 2nd j type, 3rd distance bin (for types see types2indices and indices2types maps)
	private double[][][] scoringMat;	// the scoring matrix
	
	private int[][] sumOverDistanceBins;// counts of pairs of type a,b (irrespective of distance bin)
	private int[] sumOverPairTypes;
	private long totalPairsCount;
		
	private HashMap<String,Integer> types2indices;	// map of types to indices of the above arrays
	private HashMap<Integer,String> indices2types;	// map of indices of the above arrays to types

	
	public DistanceScorer(File listFile, int minSeqSep) throws SQLException {
		this.scoringMethod = ScoringMethod.ATOMDISTANCEDEP;
		
		this.types2indices = new HashMap<String, Integer>();
		this.indices2types = new HashMap<Integer, String>();
		
		this.distanceBins = DISTANCE_BINS;
		this.structureIds = new ArrayList<String>();
		this.listFile = listFile;
		this.ct = null;
		this.cutoff = distanceBins[distanceBins.length-1];
		this.minSeqSep = minSeqSep;
		
		this.pairCounts = new int[NUM_ATOM_TYPES][NUM_ATOM_TYPES][distanceBins.length];
		
		this.conn = new MySQLConnection();
	}
	
	public DistanceScorer(File countsFile) throws IOException, FileFormatError {
		this.scoringMethod = ScoringMethod.ATOMDISTANCEDEP;

		this.types2indices = new HashMap<String, Integer>();
		this.indices2types = new HashMap<Integer, String>();
		
		readCountsFromFile(countsFile);
		calcScoringMat();
	}
	
	/**
	 * Performs the counts of the type pairs and stores them in the internal arrays.
	 * Use subsequently {@link #calcScoringMat()} to compute the scoring matrix from counts arrays.
	 * @throws SQLException if database server can't be accessed to get PDB data
	 * @throws IOException if list file can't be read  
	 */
	public void countPairs() throws SQLException, IOException {
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
			//for (AIGNode node:graph.getVertices()) {
			//	// we have to avoid OXT atoms, they are not in the map
			//	if (node.getAtomName().equals("OXT")) continue;
			//	countType(types2indices.get(node.getParent().getResidueType()+node.getAtomName()));
			//}
			for (AIGEdge edge:graph.getEdges()) {
				AIGNode inode = graph.getEndpoints(edge).getFirst();
				AIGNode jnode = graph.getEndpoints(edge).getSecond();
				int distBin = getDistanceBin(edge.getDistance());
				// we have to avoid OXT atoms, they are not in the map
				if (inode.getAtomName().equals("OXT") || jnode.getAtomName().equals("OXT")) continue;
				countPair(types2indices.get(inode.getParent().getResidueType()+inode.getAtomName()),
						types2indices.get(jnode.getParent().getResidueType()+jnode.getAtomName()),
						distBin);
			}
		}
		this.totalStructures = structureIds.size();

	}
	
	
	protected void countPair(int i, int j, int k) {
		if (j>i){
			pairCounts[i][j][k]++;
		} else {
			pairCounts[j][i][k]++;
		}
	}
	
	/**
	 * Computes the scoring matrix from the counts arrays. The score is log2(p_obs/p_exp)
	 * Whenever the pair count of a certain pair is below the threshold ({@value #TOO_FEW_COUNTS_THRESHOLD}
	 * the score assigned for that pair is {@value #TOO_FEW_COUNTS_SCORE}
	 */
	public void calcScoringMat() {
		countTotals();
		this.scoringMat = new double[NUM_ATOM_TYPES][NUM_ATOM_TYPES][distanceBins.length];
		double logof2 = Math.log(2);
		for(int i=0;i<scoringMat.length;i++) {
			for(int j=i;j<scoringMat[i].length;j++) {
				for (int k=0;k<distanceBins.length;k++) {
					if (pairCounts[i][j][k]<TOO_FEW_COUNTS_THRESHOLD) { // When counts are too small, we can't get significant statistics for them
						scoringMat[i][j][k] = TOO_FEW_COUNTS_SCORE;     // We use a score of 0 in this case, i.e. no information
					} else {
						scoringMat[i][j][k] = 
							Math.log(
									((double)pairCounts[i][j][k]/(double)sumOverDistanceBins[i][j])/
									((double)sumOverPairTypes[k]/(double)totalPairsCount)
							)/logof2; // we explicitely cast every int/long to double to avoid overflows, java doesn't check for them!
					}
				}
			}
		}
	}	
	
	public String getTypeFromIndex(int i) {
		return indices2types.get(i);
	}
	
	public int getIndexFromType(String type) {
		return types2indices.get(type);
	}
	
	private int getDistanceBin(double distance) {
		for (int k=0;k<distanceBins.length;k++) {
			int bin = k;
			if (distance<=distanceBins[k]) {
				return bin;
			}
		}
		return -1; //shouldn't even happen: the cutoff for the AIGraph limits this
	}
	
	private void countTotals() {
		sumOverDistanceBins = new int[NUM_ATOM_TYPES][NUM_ATOM_TYPES];
		sumOverPairTypes = new int[distanceBins.length];
		totalPairsCount=0;
		
		for (int i=0;i<NUM_ATOM_TYPES;i++) {
			for (int j=i;j<NUM_ATOM_TYPES;j++) {
				for (int k=0;k<distanceBins.length;k++) {
					totalPairsCount+=pairCounts[i][j][k];
					sumOverDistanceBins[i][j]+=pairCounts[i][j][k];
					sumOverPairTypes[k]+=pairCounts[i][j][k];
				}
			}
		}
	}
	
	@Override
	public double scoreIt(Pdb pdb) {
		AIGraph graph = pdb.getAllAtomGraph(this.cutoff);
		graph.restrictContactsToMinRange(minSeqSep);
		return scoreIt(graph);
	}

	protected double scoreIt(AIGraph graph) {
		double totalScore = 0;
		for (AIGEdge edge:graph.getEdges()) {
			AIGNode inode = graph.getEndpoints(edge).getFirst();
			AIGNode jnode = graph.getEndpoints(edge).getSecond();
			if (inode.getAtomName().equals("OXT") || jnode.getAtomName().equals("OXT")) continue;
			int i = types2indices.get(inode.getParent().getResidueType()+inode.getAtomName());
			int j = types2indices.get(jnode.getParent().getResidueType()+jnode.getAtomName());
			int k = getDistanceBin(edge.getDistance());
			if (j>=i) {
				totalScore+=scoringMat[i][j][k];
			} else {
				totalScore+=scoringMat[j][i][k];
			}
		}

		return (totalScore/graph.getEdgeCount());
	}

	public void writeCountsToFile (File file) throws FileNotFoundException {
		PrintWriter pw = new PrintWriter(file);

		pw.println("# SCORE METHOD: "+getScoringMethod().getDescription());
		pw.println("# contact type: "+ct);
		pw.println("# cutoff: "+cutoff);
		pw.println("# min sequence separation: "+minSeqSep);
		pw.println("# structures: "+sizeOfTrainingSet());
		pw.println("# list: "+getListFile().toString());
		pw.println("# pairs: "+totalPairsCount);
		pw.print("# distance bins: ");
		for (int k=0;k<distanceBins.length;k++) {
			pw.printf("%3.1f ",distanceBins[k]);
		}
		pw.println();

		for (int i=0;i<NUM_ATOM_TYPES;i++) {
			for (int j=i;j<NUM_ATOM_TYPES;j++) {
				for (int k=0;k<distanceBins.length;k++) {
					pw.println(indices2types.get(i)+"\t"+indices2types.get(j)+"\t"+k+"\t"+pairCounts[i][j][k]);
				}
			}
		}
		pw.close();
	}
	
	public void readCountsFromFile(File file) throws IOException, FileFormatError {
		Scorer.initAtomMap(types2indices, indices2types);
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		int lineCount = 0;
		String line;
		while ((line=br.readLine())!=null) {
			lineCount++;
			Pattern p = Pattern.compile("^# contact type: (.*)$");
			Matcher m = p.matcher(line);
			if (m.matches()) {
				this.ct = m.group(1);
			}
			p = Pattern.compile("^# cutoff: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.cutoff = Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^# min sequence separation: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.minSeqSep = Integer.parseInt(m.group(1));
			}			
			p = Pattern.compile("^# structures: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.totalStructures = Integer.parseInt(m.group(1));
			}
			p = Pattern.compile("^# list: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.listFile = new File(m.group(1));
			}			
			p = Pattern.compile("^# pairs: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.totalPairsCount = Long.parseLong(m.group(1));
			}
			p = Pattern.compile("^# distance bins: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				String[] tokens = m.group(1).split("\\s");
				distanceBins = new double[tokens.length];
				for (int k=0;k<distanceBins.length;k++) {
					distanceBins[k] = Double.parseDouble(tokens[k]);
				}
				this.pairCounts = new int[NUM_ATOM_TYPES][NUM_ATOM_TYPES][distanceBins.length];
			}
			if (!line.startsWith("#")) {
				String[] tokens = line.split("\\s");
				if (tokens.length!=4) throw new FileFormatError("Counts file "+file+" has incorrect number of columns in line "+lineCount);
				int i = types2indices.get(tokens[0]);
				int j = types2indices.get(tokens[1]);
				int k = Integer.parseInt(tokens[2]);
				int count = Integer.parseInt(tokens[3]);
				pairCounts[i][j][k] = count;
			}
		}
		
		if (cutoff!=distanceBins[distanceBins.length-1]) {
			throw new FileFormatError("Cutoff line doesn't coincide with last value in distance bins line in counts file "+file);
		}
		br.close();
	}
	
	public void writeCountsToDb(String database) throws SQLException {
		//TODO this is a (working) stub implementation, should write metadata with parameters like distanceBins and so on to separate table
		Statement stmt = conn.createStatement();
		for (int i=0;i<NUM_ATOM_TYPES;i++) {
			for (int j=i;j<NUM_ATOM_TYPES;j++) {
				for (int k=0;k<distanceBins.length;k++) {

					String sql = "INSERT INTO "+database+"."+ATOM_DISTANCE_COUNTS_TABLE+" (i_res, j_res, dist_bin, count) " +
							" VALUES ('"+indices2types.get(i)+"', '"+indices2types.get(j)+"', "+k+", "+pairCounts[i][j][k]+")";
					stmt.executeUpdate(sql);
				}
			}
		}
		stmt.close();
	}
	
	public void readCountsFromDb(String database) throws SQLException {
		//TODO this is a (working) stub implementation, missing many things: distanceBins and other parameters must be read from db (a metadata table)
		Scorer.initAtomMap(types2indices, indices2types);
		pairCounts = new int[NUM_ATOM_TYPES][NUM_ATOM_TYPES][distanceBins.length];
		Statement stmt = conn.createStatement();
		String sql = "SELECT i_res, j_res, dist_bin, count FROM "+database+"."+ATOM_DISTANCE_COUNTS_TABLE;
		ResultSet rsst = stmt.executeQuery(sql);
		while (rsst.next()) {
			int i = types2indices.get(rsst.getString(1));
			int j = types2indices.get(rsst.getString(2));
			int k = rsst.getInt(3);
			int count = rsst.getInt(4);
			if (i>j) {
				System.err.println("Warning: indices for atom types in wrong order in database!");
			}
			pairCounts[i][j][k]=count;
		}
		rsst.close();
		stmt.close();
		this.listFile = new File("unknown"); 
	}
	
	public static void main(String[] args) throws Exception {
		String help = 
			"\nCompiles a scoring matrix for a distance based atom types potential from a given file \n" +
			"with a list of pdb structures.\n" +
			"Usage:\n" +
			"DistanceScorer -o <output_counts_file> [-l <list_file> -m <out>]\n"+
			"  -o <file>     : file to write the counts matrix to\n" +
			"  -l <file>     : file with list of pdbCodes+pdbChainCodes to use as training set \n" +
			"                  the scoring matrix. Default is "+DEFAULT_LIST_FILE+"\n" +
			"  -m <int>      : minimum sequence separation to consider a contact. Default: "+DEFAULT_MIN_SEQ_SEP+"\n";

			File listFile = DEFAULT_LIST_FILE;
			File outFile = null;
			int minSeqSep = DEFAULT_MIN_SEQ_SEP;

			Getopt g = new Getopt("AtomTypeScorer", args, "l:c:o:m:h?");
			int c;
			while ((c = g.getopt()) != -1) {
				switch(c){
				case 'l':
					listFile = new File(g.getOptarg());
					break;
				case 'o':
					outFile = new File(g.getOptarg());
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

			if (outFile==null) {
				System.err.println("Missing argument output file (-o)");
				System.exit(1);
			}

			DistanceScorer sc = new DistanceScorer(listFile,minSeqSep);
			sc.countPairs();
			//sc.calcScoringMat();
			sc.writeCountsToFile(outFile);
			

	}
}

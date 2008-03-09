package proteinstructure;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.uci.ics.jung.graph.util.Pair;

/**
 * This class encapsulates calling the Maxcluster command line application and parsing its output.
 * @author stehr
 *
 */
public class MaxClusterRunner {

	private static final String TMP_MATRIX_FILE = "matrix.tmp."+System.currentTimeMillis();
	
	public enum ScoreType {RMSD,GDT};
	
	public class MaxClusterRow {
		
		private String fileName;
		private int index;
		private double score;
		private int rank;
		
		public MaxClusterRow(String fileName, int index, double score) {
			this.fileName = fileName;
			this.index = index;
			this.score = score;
			this.rank = -1;			// i.e. not ranked yet
		}

		public String getFileName() {
			return fileName;
		}

		public int getRank() {
			return rank;
		}

		public double getScore() {
			return score;
		}
		
		public int getIndex() {
			return index;
		}
		
		public void setRank(int rank) {
			this.rank = rank;
		}
		
		public String toString() {
			return String.format("%3d %s %3d %6.3f", index, fileName, rank, score);
		}
		
	}
	
	/*--------------------------- member variables --------------------------*/
	
	String maxClusterExecutable;
	
	/*----------------------------- constructors ----------------------------*/
	
	/** Create a new MaxClusterRunner object */
	public MaxClusterRunner(String maxClusterExecutable) throws FileNotFoundException {
		this.maxClusterExecutable = maxClusterExecutable;
		File file = new File(maxClusterExecutable);
		if(!file.canRead()) {
			throw new FileNotFoundException("Could not read file " + maxClusterExecutable);
		}
	}
	
	/**
	 * Compares two pdb files and returns the GDT-TS/RMSD score.
	 * @param prediction
	 * @param experiment
	 * @param scoreType
	 * @return the gdt/rmsd score or -1 if something went wrong
	 */
	public double calculatePairwiseScore(String prediction, String experiment, ScoreType scoreType) throws IOException {
		double score = -1;
		String scoreTypeStr = "";
		if (scoreType==ScoreType.GDT) scoreTypeStr = "gdt";
		if (scoreType==ScoreType.RMSD) scoreTypeStr = "rmsd";
		String cmdLine = String.format("%s %s %s -%s", maxClusterExecutable, prediction, experiment, scoreTypeStr);
		Process maxClusterProcess = Runtime.getRuntime().exec(cmdLine);
		BufferedReader maxClusterOutput = new BufferedReader(new InputStreamReader(maxClusterProcess.getInputStream()));
		String line;
		Pattern p = Pattern.compile("^(GDT|RMSD)=\\s*(\\d+\\.\\d+).*");
		while((line = maxClusterOutput.readLine()) != null) {
			Matcher m = p.matcher(line);
			if(m.matches()) {
				score = Double.parseDouble(m.group(2));
			}
		}
		return score;
	}

	
	/**
	 * Compares all files from a predictionList with experiment file returning an ArrayList of scores
	 * @param predictionList
	 * @param experiment
	 * @param scoreType 
	 * @return list of the rankings
	 */
	public ArrayList<MaxClusterRow> calculateRanking (String predictionList, String experiment, ScoreType scoreType) throws IOException {
		String scoreTypeStr = "";
		if (scoreType==ScoreType.GDT) scoreTypeStr = "gdt";
		if (scoreType==ScoreType.RMSD) scoreTypeStr = "rmsd";
		String cmdLine = String.format("%s -l %s -e %s -%s -nosort", maxClusterExecutable, predictionList, experiment, scoreTypeStr);
		Process maxClusterProcess = Runtime.getRuntime().exec(cmdLine);
		BufferedReader maxClusterOutput = new BufferedReader(new InputStreamReader(maxClusterProcess.getInputStream()));
		return readMaxClusterRanking(maxClusterOutput, scoreType);
	}

	/**
	 * Performs all against all comparison of files from a predictionList returning a matrix of pairwise scores
	 * @param predictionList
	 * @param scoreType 
	 * @return list of the rankings or null if something goes wrong
	 */
	public HashMap<Pair<Integer>, Double> calculateMatrix (String predictionList, ScoreType scoreType) throws IOException {
		String scoreTypeStr = "";
		File outFile = new File(System.getProperty("java.io.tmpdir"),TMP_MATRIX_FILE);
		outFile.deleteOnExit();
		if (scoreType==ScoreType.GDT) scoreTypeStr = "gdt";
		if (scoreType==ScoreType.RMSD) scoreTypeStr = "rmsd";
		String cmdLine = String.format("%s -l %s -C 0 -R %s -%s", maxClusterExecutable, predictionList, outFile, scoreTypeStr);
		Process maxClusterProcess = Runtime.getRuntime().exec(cmdLine);
		try {
			maxClusterProcess.waitFor();
		} catch (InterruptedException e) {
			return null;
		}
		return readMaxclusterMatrix(outFile.getAbsolutePath());
	}
	
	/**
	 * Reads a maxCluster ranking BufferedReader and returns a list of MaxClusterRows
	 * @param in
	 * @return list of the rankings or null if something goes wrong
	 */
	public ArrayList<MaxClusterRow> readMaxClusterRanking(BufferedReader in, ScoreType scoreType) throws IOException{
		ArrayList<MaxClusterRow> table = new ArrayList<MaxClusterRow>();

		String line;
		int lineNum = 0;
		while((line = in.readLine()) != null) {
			Pattern p = Pattern.compile("INFO  : +(\\d+). (.+) vs. (.+)   (GDT|RMSD)= *(\\d+\\.\\d+)");
			Matcher m = p.matcher(line);
			if(m.find()) {
				lineNum++;
				//String rank = m.group(1);
				//String targetName = m.group(2);
				String fileName = m.group(3);
				//String scoreType = m.group(4);
				String score = m.group(5);
				MaxClusterRow row = new MaxClusterRow(fileName,lineNum,Double.parseDouble(score));
				table.add(row);
			}
		}
		// calculate ranks
		
		// order by score
		if(scoreType == ScoreType.GDT) {
		Collections.sort(table, new Comparator<MaxClusterRow>() {
			public int compare(MaxClusterRow arg0, MaxClusterRow arg1) {
				return Double.compare(arg1.getScore(), arg0.getScore());
			}		
		});
		} else
			if(scoreType == ScoreType.RMSD){
			Collections.sort(table, new Comparator<MaxClusterRow>() {
				public int compare(MaxClusterRow arg0, MaxClusterRow arg1) {
					return Double.compare(arg0.getScore(), arg1.getScore());
				}		
			});			
		} else return null;	// unknown score type
		// add ranks
		int c = 1;
		for(MaxClusterRow r:table) {
			r.setRank(c);
			c++;
		}
		// order by index
		Collections.sort(table, new Comparator<MaxClusterRow>() {
			public int compare(MaxClusterRow arg0, MaxClusterRow arg1) {
				return new Integer(arg0.getIndex()).compareTo(arg1.getIndex());
			}		
		});		
		
		return table;
	}
	
	/**
	 * Read a distance matrix output file of the Maxcluster program
	 * @param fileName
	 */
	public HashMap<Pair<Integer>,Double> readMaxclusterMatrix(String fileName) throws IOException{
		HashMap<Pair<Integer>,Double> matrix = new HashMap<Pair<Integer>, Double>();
		BufferedReader in = new BufferedReader(new FileReader(fileName));
		String line;
		while((line = in.readLine()) != null) {
			if(line.startsWith("DIST")) {
				String[] t = line.split("\\s+");
				int i = Integer.parseInt(t[2]);
				int j = Integer.parseInt(t[3]);
				double dist = Double.parseDouble(t[4]);
				matrix.put(new Pair<Integer>(i,j),dist);
			}
		}
		return matrix;
	}
	
	/* for testing */
	public static void main(String[] args) throws IOException {
		String nativePdbFileName = "/project/StruPPi/projects/tinker_runs/run75x20_Cb9_weighted/out/3ezm_A_Cb_9_r1.native.pdb";
		String maxClusterExecutable = "/project/StruPPi/bin/maxcluster";
		int numModels = 20;
		String modelDir = "/project/StruPPi/projects/tinker_runs/run75x20_Cb9_weighted/out";
		String baseName = "3ezm_A_Cb_9_r1";
		String listFile = "/project/StruPPi/projects/tinker_runs/test.list";
		ScoreType scoreType = ScoreType.GDT;
		double[] gdtScores = new double[numModels];
		
		// TODO: check whether maxClusterExecutable and nativePdbFileName exist
		
		MaxClusterRunner maxCluster = new MaxClusterRunner(maxClusterExecutable);
		
		// check pairwise
		System.out.println("Pairwise:");
		for(int i=1; i <= numModels; i++) {
			String ext = String.format(".%03d",i); // 001, 002, 003, ...
			File resultPdbFile = new File(modelDir, baseName+ext+".pdb");
			// TODO: check whether resultPdbFile exists
			String modelFileName = resultPdbFile.getAbsolutePath();
			gdtScores[i-1] = maxCluster.calculatePairwiseScore(modelFileName, nativePdbFileName,scoreType);
		}		
		// print results
		for(int i=1; i <= numModels; i++) {
			System.out.println(i + "\t" + gdtScores[i-1]);
		}
		System.out.println();
		
		// check ranking
		System.out.println("Ranking:");
		ArrayList<MaxClusterRow> table = maxCluster.calculateRanking(listFile, nativePdbFileName, scoreType);
		for (MaxClusterRow row:table){
			System.out.println(row);
		}
		System.out.println();
		
		// check matrix
		System.out.println("Matrix:");
		HashMap<Pair<Integer>,Double> matrix = maxCluster.calculateMatrix(listFile, scoreType);
		for (int i=1;i<=numModels;i++){
			for (int j=i+1;j<=numModels;j++) {
				System.out.printf("%5.2f ",matrix.get(new Pair<Integer>(i,j)));
			}
			System.out.println();
		}
		System.out.println();
	}

}

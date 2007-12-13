package proteinstructure;

import java.io.*;
import java.util.ArrayList;
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
		
		String fileName;
		int rank;
		double score;
		
		public MaxClusterRow(String fileName, int rank, double score) {
			this.fileName = fileName;
			this.rank = rank;
			this.score = score;
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
		while((line = maxClusterOutput.readLine()) != null) {
			if(line.startsWith(scoreTypeStr.toUpperCase()+"=")) {
				String gdtStr = line.substring(5);
				score = Double.parseDouble(gdtStr);
				break;
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
		String cmdLine = String.format("%s -l %s -e %s -%s", maxClusterExecutable, predictionList, experiment, scoreTypeStr);
		Process maxClusterProcess = Runtime.getRuntime().exec(cmdLine);
		BufferedReader maxClusterOutput = new BufferedReader(new InputStreamReader(maxClusterProcess.getInputStream()));
		return readMaxClusterRanking(maxClusterOutput);
	}

	/**
	 * Compares all files from a predictionList with experiment file returning an ArrayList of scores
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
		return readFromMaxclusterMatrix(outFile.getAbsolutePath());
	}
	
	/**
	 * Reads a maxCluster ranking BufferedReader and returns a list of MaxClusterRows
	 * @param in
	 * @return list of the rankings
	 */
	public ArrayList<MaxClusterRow> readMaxClusterRanking(BufferedReader in) throws IOException{
		ArrayList<MaxClusterRow> table = new ArrayList<MaxClusterRow>();

		String line;
		while((line = in.readLine()) != null) {
			Pattern p = Pattern.compile("INFO  : +(\\d+). (.+) vs. (.+)   (GDT|RMSD)= *(\\d+\\.\\d+)");
			Matcher m = p.matcher(line);
			if(m.find()) {
				String rank = m.group(1);
				//String targetName = m.group(2);
				String fileName = m.group(3);
				//String scoreType = m.group(4);
				String score = m.group(5);
				MaxClusterRow row = new MaxClusterRow(fileName,Integer.parseInt(rank),Double.parseDouble(score));
				table.add(row);

			}
		}
		return table;
	}
	
	/**
	 * Read a distance matrix output file of the Maxcluster program
	 * @param fileName
	 */
	public HashMap<Pair<Integer>,Double> readFromMaxclusterMatrix(String fileName) throws IOException{
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
		double[] gdtScores = new double[numModels];
		
		// TODO: check whether maxClusterExecutable and nativePdbFileName exist
		
		MaxClusterRunner maxCluster = new MaxClusterRunner(maxClusterExecutable);
		
		for(int i=1; i <= numModels; i++) {
			String ext = String.format(".%03d",i); // 001, 002, 003, ...
			File resultPdbFile = new File(modelDir, baseName+ext+".pdb");
			// TODO: check whether resultPdbFile exists
			String modelFileName = resultPdbFile.getAbsolutePath();
			gdtScores[i-1] = maxCluster.calculatePairwiseScore(modelFileName, nativePdbFileName,ScoreType.GDT);
		}
		
		// print results
		for(int i=1; i <= numModels; i++) {
			System.out.println(i + "\t" + gdtScores[i-1]);
		}
		
		ArrayList<MaxClusterRow> table = maxCluster.calculateRanking(listFile, nativePdbFileName, ScoreType.GDT);
		for (MaxClusterRow row:table){
			System.out.println(row.fileName+" "+row.rank+" "+row.score);
		}
		
		HashMap<Pair<Integer>,Double> matrix = maxCluster.calculateMatrix(listFile, ScoreType.GDT);
		for (int i=1;i<=numModels;i++){
			for (int j=i+1;j<=numModels;j++) {
				System.out.printf("%5.2f ",matrix.get(new Pair<Integer>(i,j)));
			}
			System.out.println();
		}
	}

}

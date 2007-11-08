package proteinstructure;

import java.io.*;

/**
 * This class encapsulates calling the Maxcluster command line application and parsing its output.
 * @author stehr
 *
 */
public class MaxClusterRunner {

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
	 * Compares two pdb files and returns the GDT-TS score.
	 * @return the gdt score or -1 if something went wrong
	 */
	public double calculateGdt(String prediction, String experiment) throws IOException {
		double gdtScore = -1;
		String cmdLine = String.format("%s %s %s -gdt", maxClusterExecutable, prediction, experiment);
		Process maxClusterProcess = Runtime.getRuntime().exec(cmdLine);
		BufferedReader maxClusterOutput = new BufferedReader(new InputStreamReader(maxClusterProcess.getInputStream()));
		String line;
		while((line = maxClusterOutput.readLine()) != null) {
			if(line.startsWith("GDT=")) {
				String gdtStr = line.substring(5);
				gdtScore = Double.parseDouble(gdtStr);
				break;
			}
		}
		return gdtScore;
	}
	
	/* for testing */
	public static void main(String[] args) throws IOException {
		String nativePdbFileName = "/project/StruPPi/projects/tinker_runs/run75x20_Cb9_weighted/out/3ezm_A_Cb_9_r1.native.pdb";
		String maxClusterExecutable = "/project/StruPPi/bin/maxcluster";
		int numModels = 20;
		String modelDir = "/project/StruPPi/projects/tinker_runs/run75x20_Cb9_weighted/out";
		String baseName = "3ezm_A_Cb_9_r1";
		double[] gdtScores = new double[numModels];
		
		// TODO: check whether maxClusterExecutable and nativePdbFileName exist
		
		MaxClusterRunner maxCluster = new MaxClusterRunner(maxClusterExecutable);
		
		for(int i=1; i <= numModels; i++) {
			String ext = String.format(".%03d",i); // 001, 002, 003, ...
			File resultPdbFile = new File(modelDir, baseName+ext+".pdb");
			// TODO: check whether resultPdbFile exists
			String modelFileName = resultPdbFile.getAbsolutePath();
			gdtScores[i-1] = maxCluster.calculateGdt(modelFileName, nativePdbFileName);
		}
		
		// print results
		for(int i=1; i <= numModels; i++) {
			System.out.println(i + "\t" + gdtScores[i-1]);
		}
		
	}

}

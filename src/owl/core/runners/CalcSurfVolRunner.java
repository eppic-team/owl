package owl.core.runners;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;

import owl.core.structure.PdbChain;

/**
 * Class to run the calc-volume and calc-surface programs:
 * NR Voss, MB Gerstein. Calculation of standard atomic volumes for RNA and comparison
 * with proteins: RNA is packed more tightly. J Mol. Biol. v346(2): 2005, pp 477-492.
 * doi:10.1016/j.jmb.2004.11.072 
 * 
 * Software downloadable from http://geometry.molmovdb.org/
 *
 */
public class CalcSurfVolRunner {

	
	private static final String TMP_DIR = System.getProperty("java.io.tmpdir");
	
	
	/** 
	 * Runs an external calc-surface executable and returns the total ASA 
	 * @param calcExecutable
	 * @param calcParameters
	 * @return
	 */
	public static double calcSurface(PdbChain pdb, String calcExecutable, String calcParameters) throws IOException {
		String line;
		double surface = 0;
		
		File test = new File(calcExecutable);
		if(!test.canRead()) throw new IOException("calc-surface Executable is not readable");
		Process myCalc = Runtime.getRuntime().exec(calcExecutable + " -i - " + calcParameters);
		PrintStream calcInput = new PrintStream(myCalc.getOutputStream());
		BufferedReader calcOutput = new BufferedReader(new InputStreamReader(myCalc.getInputStream()));
		BufferedReader calcError = new BufferedReader(new InputStreamReader(myCalc.getErrorStream()));
		pdb.writeAtomLines(calcInput);	// pipe atom lines to calc-surface
		calcInput.close();
		while((line = calcOutput.readLine()) != null) {
			surface += Double.valueOf(line.substring(66,73).trim());
		}
		calcOutput.close();
		calcError.close();
		
		return surface;
	}

	/** 
	 * Runs an external calc-volume executable and returns the volume 
	 * @param calcExecutable
	 * @param calcParameters
	 * @return
	 */
	public static double calcVolume(PdbChain pdb, String calcExecutable, String calcParameters) throws IOException {
		String line;
		double vol = 0;
		File pdbFile = new File(TMP_DIR,pdb.getPdbCode()+pdb.getChainCode()+"_"+System.currentTimeMillis()+".pdb");
		pdbFile.deleteOnExit();
		String pdbFileName = pdbFile.getAbsolutePath();
		pdb.writeToPDBFile(pdbFileName);
		
		File test = new File(calcExecutable);
		if(!test.canRead()) throw new IOException("calc-volume Executable is not readable");
		Process myCalc = Runtime.getRuntime().exec(calcExecutable + " -i "+pdbFileName+" "+ calcParameters);
		BufferedReader calcOutput = new BufferedReader(new InputStreamReader(myCalc.getInputStream()));
		BufferedReader calcError = new BufferedReader(new InputStreamReader(myCalc.getErrorStream()));
		while((line = calcOutput.readLine()) != null) {
			if (line.substring(80,82).trim().equals("0")) {
				vol += Double.valueOf(line.substring(66,79).trim());
			}
		}
		calcOutput.close();
		calcError.close();
		
		return vol;
	}
	

}

package owl.core.runners;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import owl.core.runners.blast.BlastError;
import owl.core.runners.blast.BlastRunner;

/**
 * A psipred runner, adapted from shell script provided by David Jones.
 * See also {@link SecondaryStructure.readPsiPredHorizFile} and {@link JPredConnection}.
 *
 */
public class PsipredRunner {
	
	private static final String PASS1_PROG = "psipred";
	private static final String PASS2_PROG = "psipass2";

	private String pass1Prog;
	private String pass2Prog;
	
	private String weights1Dat;
	private String weights2Dat;
	private String weights3Dat;
	private String weights4Dat;
	private String weightsP2Dat;
	
	public PsipredRunner(String psipredHomeDir) {
		String binDir = new File(psipredHomeDir,"bin").getAbsolutePath();
		this.pass1Prog = new File(binDir, PASS1_PROG).getAbsolutePath();
		this.pass2Prog = new File(binDir, PASS2_PROG).getAbsolutePath();
		String dataDir = new File(psipredHomeDir,"data").getAbsolutePath();
		this.weights1Dat = new File(dataDir, "weights.dat").getAbsolutePath();
		this.weights2Dat = new File(dataDir, "weights.dat2").getAbsolutePath();
		this.weights3Dat = new File(dataDir, "weights.dat3").getAbsolutePath();
		this.weights4Dat = new File(dataDir, "weights.dat4").getAbsolutePath();
		this.weightsP2Dat = new File(dataDir, "weights_p2.dat").getAbsolutePath();
		
	}
	
	private void runPass1(File inMtxFile, File outSsFile) throws IOException, PsipredError {
		String cmdLine = pass1Prog + " "+ inMtxFile.getAbsolutePath() +" "+weights1Dat+" "+weights2Dat+" "+weights3Dat+" "+weights4Dat;
		Process pass1Proc = Runtime.getRuntime().exec(cmdLine);
		BufferedReader pass1Output = new BufferedReader(new InputStreamReader(pass1Proc.getInputStream()));
		PrintWriter out = new PrintWriter(new FileWriter(outSsFile));
		String line;
		while ((line=pass1Output.readLine())!=null) {
			out.println(line);
		}
		out.close();
		
		try {
			int exitValue = pass1Proc.waitFor();
			if (exitValue>0) {
				throw new PsipredError(PASS1_PROG + " exited with error value " + exitValue);
			}
		} catch (InterruptedException e) {
			System.err.println("Unexpected error while running psipred: "+e.getMessage());
		}
	}
	
	private void runPass2(File inSsFile, File outSs2File, File outHorizFile) throws IOException, PsipredError {
		String cmdLine = pass2Prog + " " + weightsP2Dat + " 1 1.0 1.0 "+ outSs2File.getAbsolutePath()+" "+ inSsFile.getAbsolutePath();
		Process pass2Proc = Runtime.getRuntime().exec(cmdLine);
		BufferedReader pass2Output = new BufferedReader(new InputStreamReader(pass2Proc.getInputStream()));
		PrintWriter out = new PrintWriter(new FileWriter(outHorizFile));
		String line;
		while ((line=pass2Output.readLine())!=null) {
			out.println(line);
		}
		out.close();
		try {
			int exitValue = pass2Proc.waitFor();
			if (exitValue>0) {
				throw new PsipredError(PASS2_PROG + " exited with error value " + exitValue);
			}
		} catch (InterruptedException e) {
			System.err.println("Unexpected error while running psipred: "+e.getMessage());
		}		
	}
	
	private void copyFile(File inFile, File outFile) throws FileNotFoundException, IOException {
		FileInputStream in = new FileInputStream(inFile);
		FileOutputStream out = new FileOutputStream(outFile);
		int b;
		while ((b=in.read())!=-1) {
			out.write(b);
		}
		in.close();
		out.close();
	}
	
	/**
	 * Gets a file basename, stripped of its path and its extension if there is one
	 * @param file
	 * @return
	 */
	private String getFileBasename(File file) {
		String name = file.getName();
		if (name.contains(".")) {
			name = name.substring(0,name.lastIndexOf("."));
		}
		return name;
	}
	
	private void writeStringToFile(File file, String string) throws FileNotFoundException, IOException{
		PrintWriter pw = new PrintWriter(new FileWriter(file));
		pw.println(string);
		pw.close();
	}
	
	/**
	 * Runs a secondary structure prediction with psipred using a psi-blast chk file as input profile
	 * @param inSeqFile
	 * @param outSs2File
	 * @param outHorizFile
	 * @param blastChkFile
	 * @param blastBinDir
	 */
	public void run(File inSeqFile, File outSs2File, File outHorizFile, File blastChkFile, String blastBinDir) throws IOException, PsipredError {
		// copy seq file and chk profile file to tmp dir and get basename of chk file 
		String tmpDir = System.getProperty("java.io.tmpdir");
		File tmpChkFile = new File(tmpDir,blastChkFile.getName());
		tmpChkFile.deleteOnExit();
		File tmpSeqFile = new File(tmpDir,inSeqFile.getName());
		tmpSeqFile.deleteOnExit();
		copyFile(blastChkFile,tmpChkFile);
		copyFile(inSeqFile, tmpSeqFile);
		String basename = getFileBasename(tmpChkFile);
		// create the .cn, .sn file required by makemat, containing the names of the seq file and profile chk file
		File pnFile = new File(tmpDir,basename+".pn");
		pnFile.deleteOnExit();
		writeStringToFile(pnFile, tmpChkFile.getAbsolutePath());
		File snFile = new File(tmpDir,basename+".sn");
		snFile.deleteOnExit();
		writeStringToFile(snFile, tmpSeqFile.getAbsolutePath());
		
		// running makemat (NOTE: we won't blast with this BlastRunner, we do just makemat, thus we don't need a blastDbDir)
		BlastRunner blastRunner = new BlastRunner(blastBinDir, null);
		try {
			blastRunner.runMakemat(tmpDir, basename);
		} catch (BlastError e) {
			throw new PsipredError("Makemat step of psipred failed to run, error: "+e.getMessage());
		}
		// running psipred pass 1 and pass 2 with output of makemat (.mtx file)
		File inMtxFile = new File(tmpDir,basename+".mtx");
		inMtxFile.deleteOnExit();
		File outSsFile = new File(tmpDir,basename+".ss");
		outSsFile.deleteOnExit();
		runPass1(inMtxFile, outSsFile);
		runPass2(outSsFile, outSs2File, outHorizFile);
		
		// cleaning up other temp files generated by makemat/psipred
		new File(tmpDir, basename+".aux").deleteOnExit();
		new File(tmpDir, basename+".mn").deleteOnExit();
	}
	
	
}

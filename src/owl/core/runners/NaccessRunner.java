package owl.core.runners;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;

import owl.core.structure.Pdb;
import owl.core.structure.Residue;



/**
 * A class to run naccess and parse its output
 * 
 * @author duarte
 *
 */
public class NaccessRunner {

	
	private static final boolean DEBUG = false; // set to true to keep output files and write out naccess command line
	
	private File naccessExecutable;
	private String naccessParameters;
	
	/**
	 * Constructs a NaccessRunner by passing the executable location and parameters
	 * @param naccessExecutable
	 * @param naccessParameters
	 * @throws IOException if naccessExecutable not readable
	 */
	public NaccessRunner(File naccessExecutable, String naccessParameters) throws IOException {
		this.naccessExecutable = naccessExecutable;
		this.naccessParameters = naccessParameters;
		if(!naccessExecutable.canRead()) 
			throw new IOException("Naccess executable "+naccessExecutable+" is not readable");
	}
	
	/**
	 * Runs naccess updating the given pdb object with rsa and sc rsa values.
	 * All temporary files are removed unless the {@link #DEBUG} flag is set.
	 * @throws IOException if I/O problems running naccess or if it finishes with an 
	 * error exit status
	 */
	public void runNaccess(Pdb pdb) throws IOException {

		File tempDir = new File(System.getProperty("java.io.tmpdir"));
		String prefix = ((pdb.getPdbCode()+pdb.getChainCode()).length() < 3)?"1xxx":(pdb.getPdbCode()+pdb.getChainCode());
		File pdbFile = File.createTempFile(prefix, ".pdb",tempDir);
		String baseName = pdbFile.getName().substring(0, pdbFile.getName().lastIndexOf('.'));

		pdb.writeToPDBFile(pdbFile.getAbsolutePath());
		String line;
		int errorLineCount = 0;

		String cmd = naccessExecutable + " " + pdbFile.getName() + " " + naccessParameters;
		if (DEBUG) System.out.println(cmd);
		// we have to run with working dir the tempDir, naccess always writes output to current dir
		Process myNaccess = Runtime.getRuntime().exec(cmd,null,tempDir);
		BufferedReader naccessOutput = new BufferedReader(new InputStreamReader(myNaccess.getInputStream()));
		BufferedReader naccessError = new BufferedReader(new InputStreamReader(myNaccess.getErrorStream()));
		while((line = naccessOutput.readLine()) != null) {
		}
		while((line = naccessError.readLine()) != null) {
			errorLineCount++;
		}
		naccessOutput.close();
		naccessError.close();
		int exitVal = 1;
		try {
			exitVal = myNaccess.waitFor();
		} catch (InterruptedException e) {
			System.err.println("Unexpected error while waiting for naccess to finish");
		}
		if ((exitVal == 1) || (errorLineCount > 0)) {
			throw new IOException("Naccess error: wrong arguments or pdb file format!");
		}

		File rsa = new File(tempDir,baseName+".rsa");
		if (rsa.exists()) {
			BufferedReader rsaInput = new BufferedReader(new FileReader(rsa));
			while ((line = rsaInput.readLine()) != null) {
				if (line.startsWith("RES")) {
					int resser = Integer.valueOf(line.substring(9,13).trim());
					double allrsa = Double.valueOf(line.substring(22,28).trim());
					double scrsa = Double.valueOf(line.substring(35,41).trim());
					if (pdb.containsResidue(resser)) {
						Residue residue = pdb.getResidue(resser);
						residue.setRsa(allrsa);
						residue.setScRsa(scrsa);
					}
				}
			}
			rsaInput.close();
		} else {
			throw new IOException("Naccess output file "+rsa+" wasn't found");
		}
		if (!DEBUG) {
			String[] filesToDelete = {".pdb", ".rsa", ".asa", ".log" };
			for (int i=0; i < filesToDelete.length; i++) {
				File fileToDelete = new File(tempDir,baseName+filesToDelete[i]);
				if (fileToDelete.exists()) {
					fileToDelete.deleteOnExit();
				}
			}
		}
		pdb.setHasASA(true);
	}

}

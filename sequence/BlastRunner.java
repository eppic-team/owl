package sequence;

import java.io.File;
import java.io.IOException;

/**
 * A class to run blast programs
 *
 */
public class BlastRunner {

	// constants
	private static final String BLASTALL_PROG = "blastall";
	private static final String BLASTPGP_PROG = "blastpgp";
	
	private static final String BLASTP_PROGRAM_NAME = "blastp";
	
	// members
	private String blastDbDir;
	private String blastallProg;
	private String blastpgpProg;
	
	/**
	 * Constructs a BlastRunner object given a blast bin directory and a blast 
	 * db directory
	 * @param blastBinDir
	 * @param blastDbDir
	 */
	public BlastRunner(String blastBinDir, String blastDbDir) {
		this.blastDbDir = blastDbDir;
		this.blastallProg = new File(blastBinDir,BLASTALL_PROG).getAbsolutePath();
		this.blastpgpProg = new File(blastBinDir,BLASTPGP_PROG).getAbsolutePath();
	}
	
	private String getCommonOptionsStr(File queryFile, String db, File outFile, int outputType) { 
		String dbFullPath = new File(blastDbDir,db).getAbsolutePath();
		String options = " -m "+outputType+
						 " -i "+queryFile.getAbsolutePath()+
						 " -d "+dbFullPath+
						 " -o "+outFile.getAbsolutePath() + " ";
		return options;
	}
	
	private void checkIO(File queryFile, String db) throws IOException{
		if (!queryFile.canRead()) throw new IOException("Can't read query file "+queryFile);
		String dbFullPath = new File(blastDbDir,db).getAbsolutePath();
		if (! (new File(dbFullPath+".pal").canRead() || new File(dbFullPath+".psq").canRead()))
			throw new IOException("Can't read the database files for database in "+dbFullPath);
	}
	
	/**
	 * Runs psi-blast for given query file against given db
	 * @param queryFile
	 * @param db the identifier of the blast database to run against
	 * @param outFile the blast output file
	 * @param maxIter
	 * @param outProfileFile the output profile file for the -C option, null if 
	 * not required. Must have .chk extension
	 * @param inProfileFile the input profile file for the -R option, null if a 
	 * psi-blast from scratch wanted. Must have .chk extension
	 * @param outputType 0 classic, 7 XML, 8 tabular, 9 tabular with comments 
	 * @throws IOException
	 * @throws BlastError if exit status of the program is not 0
	 */
	public void runPsiBlast(File queryFile, String db, File outFile, int maxIter, File outProfileFile, File inProfileFile, int outputType) 
	throws IOException, BlastError {
		
		checkIO(queryFile, db);
		
		String outProfileOpt = "";
		String inProfileOpt = "";
		if (outProfileFile!=null) outProfileOpt = " -C "+outProfileFile.getAbsolutePath();
		if (inProfileFile!=null) inProfileOpt = " -R "+inProfileFile.getAbsolutePath();
		String cmdLine = blastpgpProg + getCommonOptionsStr(queryFile, db, outFile, outputType) +
						" -j " + maxIter+
						outProfileOpt + inProfileOpt;
		Process blastpgpProc = Runtime.getRuntime().exec(cmdLine);

		try {
			int exitValue = blastpgpProc.waitFor();
			if (exitValue>0) {
				throw new BlastError(BLASTALL_PROG + " exited with error value " + exitValue);
			}
		} catch (InterruptedException e) {
			System.err.println("Unexpected error while running blast: "+e.getMessage());
		}
	}
	
	/**
	 * Runs the blastall program given by prog, for given query file against database db
	 * @param queryFile
	 * @param db the identifier of the blast database to run against
	 * @param outFile the blast output file
	 * @param prog the blast program: one of blastp, blastn, blastx, tblastn, tblastx 
	 * @param outputType 0 classic, 7 XML, 8 tabular, 9 tabular with comments
	 * @throws IOException
	 * @throws BlastError if exit statis of the program is not 0
	 */
	public void runBlast(File queryFile, String db, File outFile, String prog, int outputType) throws IOException, BlastError {
		
		checkIO(queryFile, db);
		
		String cmdLine = blastallProg + " -p " + prog + getCommonOptionsStr(queryFile, db, outFile, outputType);
		Process blastallProc = Runtime.getRuntime().exec(cmdLine);
		
		try {
			int exitValue = blastallProc.waitFor();
			if (exitValue>0) {
				throw new BlastError(BLASTALL_PROG + " exited with error value " + exitValue);
			}
		} catch (InterruptedException e) {
			System.err.println("Unexpected error while running blast: "+e.getMessage());
		}
	}
	
	/**
	 * Runs protein blast against given db for given input query file 
	 * @param queryFile
	 * @param db
	 * @param outFile
	 * @param outputType 0 classic, 7 XML, 8 tabular, 9 tabular with comments 
	 * @throws IOException
	 * @throws BlastError
	 */
	public void runBlastp(File queryFile, String db, File outFile, int outputType) throws IOException, BlastError {
		runBlast(queryFile, db, outFile, BLASTP_PROGRAM_NAME, outputType);
	}
	
}

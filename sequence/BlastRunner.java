package sequence;

import java.io.File;
import java.io.IOException;

import proteinstructure.TemplateList;
import tools.MySQLConnection;

/**
 * A class to run blast programs
 *
 */
public class BlastRunner {

	// constants
	private static final String BLASTALL_PROG = "blastall";
	private static final String BLASTPGP_PROG = "blastpgp";
	private static final String MAKEMAT_PROG  = "makemat_withcd"; // customised makemat script that cds to working dir before running makemat
	
	private static final String BLASTP_PROGRAM_NAME = "blastp";
	
	// members
	private String blastDbDir;
	private String blastallProg;
	private String blastpgpProg;
	private String makematProg;
	
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
		this.makematProg = new File(blastBinDir,MAKEMAT_PROG).getAbsolutePath();
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
	
	/**
	 * 
	 * @param workDir
	 * @param basename
	 */
	public void runMakemat(String workDir, String basename) throws IOException, BlastError {
		String cmdLine = makematProg + " "+workDir+" -P " + basename;
		Process makematProc = Runtime.getRuntime().exec(cmdLine);
		
		try {
			int exitValue = makematProc.waitFor();
			if (exitValue>0) {
				throw new BlastError(MAKEMAT_PROG + " exited with error value " + exitValue);
			}
		} catch (InterruptedException e) {
			System.err.println("Unexpected error while running blast: "+e.getMessage());
		}
	}
	
	/**
	 * Test blast runner
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		String BLASTBIN_DIR = "/project/StruPPi/bin";
		String BLASTDB_DIR = "/project/StruPPi/CASP8/blast_dbs";
		int OUTPUT_TYPE = 9;
		double eValCutoff = 1e-30;
		File queryFile = new File("/project/StruPPi/CASP8/example_files/2uczA_ubc7.fasta");
		
		// read sequence
		
		Sequence query = new Sequence();
		try {
			query.readFromFastaFile(queryFile);
		} catch(IOException e) {
			System.out.println("Error reading Fasta file: " + e.getMessage());
		}
		int queryLength = query.getLength();

		// run blast
		
		BlastRunner br = new BlastRunner(BLASTBIN_DIR, BLASTDB_DIR);
		File outFile = new File("out.blast");
		File outFile2 = new File("out2.blast");
		File outFile3 = new File("out3.blast");
		String nrdb = "nr";
		String pdbdb = "fasta_file_from_pdbase.fix.reps.fa";
		int maxIter = 2;
		File outProfileFile = new File("out.chk");
		
		// BlastP
		System.out.println("Running blast against PDB...");
		br.runBlastp(queryFile, pdbdb, outFile, OUTPUT_TYPE);
		BlastTabularParser blastParser1 = new BlastTabularParser(outFile);
		BlastHitList hits1 = blastParser1.getHits();
		hits1.setQueryLength(queryLength);
		System.out.println("Number of hits: "+hits1.size());
		if(hits1.size() > 0) System.out.println("Best E-value: "+hits1.getBestHit().getEValue());
		hits1.printSome(5);
		
		// PsiBlast
		System.out.println("Create psi-blast profile with " + maxIter + " iterations against " + nrdb + "...");
		br.runPsiBlast(queryFile, nrdb, outFile2, maxIter, outProfileFile, null, OUTPUT_TYPE);
		System.out.println("Running psi-blast against PDB...");
		br.runPsiBlast(queryFile, pdbdb, outFile3, 1, null, outProfileFile, OUTPUT_TYPE);
		
		BlastTabularParser blastParser2 = new BlastTabularParser(outFile);
		BlastHitList hits2 = blastParser2.getHits();
		hits2.setQueryLength(queryLength);
		System.out.println("Number of hits: "+hits2.size());
		if(hits2.size() > 0) {
			System.out.println("Best E-value: "+hits2.getBestHit().getEValue());
			System.out.println("Filtering with cutoff " + eValCutoff);
			hits2.applyCutoff(eValCutoff);
			System.out.println("Number of hits: "+hits2.size());
			hits2.printSome(5);
			System.out.println("Generating cluster graph...");
			String[] ids = hits2.getTemplateIds();
			TemplateList templates = new TemplateList(ids);
			BlastUtils.writeClusterGraph(templates, new MySQLConnection("white","duarte", "nieve"), "pdbase",new File("out.gdl"), new File("out.matrix"), 6.0);
		}
	}
}

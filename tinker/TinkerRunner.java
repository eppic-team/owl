package tinker;

import edu.uci.ics.jung.graph.util.Pair;
import graphAveraging.ConsensusSquare;
import graphAveraging.GraphAverager;
import graphAveraging.GraphAveragerError;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.channels.FileChannel;
import java.util.Formatter;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Point3d;

import org.ggf.drmaa.DrmaaException;
import org.ggf.drmaa.JobInfo;
import org.ggf.drmaa.JobTemplate;
import org.ggf.drmaa.Session;
import org.ggf.drmaa.SessionFactory;

import proteinstructure.AAinfo;
import proteinstructure.IntPairSet;
import proteinstructure.Pdb;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbfilePdb;
import proteinstructure.MaxClusterRunner;
import proteinstructure.RIGEdge;
import proteinstructure.RIGEnsemble;
import proteinstructure.RIGNode;
import proteinstructure.RIGraph;

public class TinkerRunner {
	
	/*------------------------------ constants ------------------------------*/

	private static final String PROTEIN_PROG = "protein";
	private static final String DISTGEOM_PROG = "distgeom";
	private static final String PDBXYZ_PROG = "pdbxyz";
	private static final String XYZPDB_PROG = "xyzpdb";
	private static final String MINIMIZE_PROG = "minimize";
	private static final String ANALYZE_PROG = "analyze";
	
	private static final String CYCLISE_PROTEIN_STR = "N";
	private static final String DGEOM_DEFAULT_PARAMS = "Y N Y Y N N ";
	private static final String REFINE_VIA_ANNEALING = "A";
	private static final String REFINE_VIA_MINIMIZATION = "M";
	
	public static final double DEFAULT_FORCECONSTANT_DISTANCE = 100.0;
	public static final double DEFAULT_FORCECONSTANT_TORSION = 1.0;
	
	private static final String TINKER_ERROR_STR = " TINKER is Unable to Continue";
	private static final String CHECKXYZ_WARNING = " CHKXYZ";
	
	private static final PRMInfo.PRMType DEFAULT_FF_FILE_TYPE = PRMInfo.PRMType.amber;
	
	public static final String DEFAULT_RECONSTR_CHAIN_CODE = Pdb.NULL_CHAIN_CODE;
	
	private static final String ANALYZE_ENERGY_MODE = "E"; // energy mode of analyze program, other valid modes are: A, L, D, M, P (see tinker's docs)
	
	private static final int MAXSEED = 2000000000; // max seed that tinker supports
	
	private static final boolean DEBUG = false;              // if true temp files for parallel distgeom run are kept
	private static final long PARALLEL_JOBS_TIMEOUT = 36000; // (36000s = 10h) timeout for SGE jobs to finish (in seconds)
	private static final String SGE_JOBS_PREFIX = "RC_";     // prefix for SGE jobs
	private static final String SGE_QUEUE = "-q all.q";      // queue were SGE jobs will run
	private static final int MAX_RETRIES_FIND_OUTPUT = 10;   // max number of retries for checking output files of a parallel distgeom run
	private static final long RETRY_TIME_FIND_OUTPUT = 2000; // time between retries for checking output files of a parallel disgeom run
	
	/*--------------------------- member variables --------------------------*/
	// input parameters
	private String tinkerBinDir;
	private String forceFieldFileName;
	private double forceConstant;
	
	private String dgeomParams;
	
	// derived parameters
	private String proteinProg;
	private String distgeomProg;
	private String pdbxyzProg;
	private String xyzpdbProg;
	private String minimizeProg;
	private String analyzeProg;
	
	// variables for storing distgeom output data
	private double[] errorFunctionVal; 
	private int[] numUpperBoundViol;
	private int[] numLowerBoundViol;
	private double[] maxUpperBoundViol;
	private double[] maxLowerBoundViol;
	private double[] rmsBoundViol;
	private int[] numUpperViol;
	private int[] numLowerViol;
	private double[] maxUpperViol;
	private double[] maxLowerViol;
	private double[] rmsRestViol;

	// information about last reconstruct run
	String lastOutputDir;       // output directory of last reconstruction run
	String lastBaseName;		// basename of last reconstruction run
	int lastNumberOfModels;		// number of models in last reconstruction run
	
	// for parallel runs
	Session session;			// the DRMAA session to communicate with the sge system
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Constructs a TinkerRunner object by passing initial parameters
	 * @param tinkerBinDir The directory where the tinker executables are
	 * @param forceFieldFileName The force field file
	 * @throws FileNotFoundException
	 */
	public TinkerRunner(String tinkerBinDir, String forceFieldFileName) throws FileNotFoundException {
		this.tinkerBinDir = tinkerBinDir;
		File tinkerbindir = new File(tinkerBinDir);
		if (!(tinkerbindir.isDirectory() && tinkerbindir.canRead())) throw new FileNotFoundException("Can't read tinker bin directory "+tinkerBinDir);
		this.forceFieldFileName = forceFieldFileName;
		this.proteinProg = new File(this.tinkerBinDir,PROTEIN_PROG).getAbsolutePath();
		this.distgeomProg = new File(this.tinkerBinDir,DISTGEOM_PROG).getAbsolutePath();
		this.pdbxyzProg = new File(this.tinkerBinDir,PDBXYZ_PROG).getAbsolutePath();
		this.xyzpdbProg = new File(this.tinkerBinDir,XYZPDB_PROG).getAbsolutePath();
		this.minimizeProg = new File(this.tinkerBinDir,MINIMIZE_PROG).getAbsolutePath();
		this.analyzeProg = new File(this.tinkerBinDir,ANALYZE_PROG).getAbsolutePath();
		this.dgeomParams = DGEOM_DEFAULT_PARAMS+REFINE_VIA_ANNEALING; // default: Annealing refinement 
		
		this.forceConstant = -1;
		this.lastOutputDir = null;
		this.lastBaseName = null;
		this.lastNumberOfModels = 0;
	}
		
	/*---------------------------- private methods --------------------------*/
	
	/**
	 * To get the expected File that a tinker program will output given an input 
	 * file and an extension for the output files
	 * The directory where the input file is will be scanned to see if it contains 
	 * files of the form basename.ext, basename.ext_2, basename.ext_3 etc.
	 * @param file
	 * @param ext
	 * @return
	 */
	private File getTinkerOutputFileName(File file, String ext){
		String basename = getBasename(file);
		String dirname = file.getParent();
		
		String tinkerOutFileName = basename + "." + ext;
		
		if (new File(dirname,tinkerOutFileName).exists()) {
			int i = 2;
			tinkerOutFileName = basename + "." + ext + "_" + i;
			while (new File(dirname,tinkerOutFileName).exists()) {
				i++;
				tinkerOutFileName = basename + "." + ext + "_" + i;
			}
		}
		return new File(dirname,tinkerOutFileName);
	}
	
	/**
	 * Runs tinker's protein program to generate an elongated protein structure 
	 * given a sequence
	 * @param sequence
	 * @param outPath The directory where output files will be written
	 * @param outBasename The base name for the output files
	 * @param log A PrintWriter for logging output
	 * @throws IOException
	 * @throws TinkerError
	 */
	private void runProtein(String sequence, String outPath, String outBasename, PrintWriter log) throws IOException, TinkerError {
		boolean tinkerError = false; // to store the exit state of the tinker program
		
		if (!new File(outPath).exists()) {
			throw new FileNotFoundException("Specified directory "+outPath+" does not exist");
		}
		File tinkerxyzout = getTinkerOutputFileName(new File(outPath,outBasename+".xyz"),"xyz");
		File tinkerintout = getTinkerOutputFileName(new File(outPath,outBasename+".int"),"int");
		tinkerintout.deleteOnExit();
		File tinkerseqout = getTinkerOutputFileName(new File(outPath,outBasename+".seq"),"seq");
		
		// running protein program in outPath dir (so that output files are written to outPath)
		Process protProc = Runtime.getRuntime().exec(proteinProg, null, new File(outPath));
		// piping input
		PrintWriter protInput = new PrintWriter(protProc.getOutputStream());
		protInput.println(outBasename);
		protInput.println("Unfolded chain created by tinker's protein program");
		protInput.println(forceFieldFileName);
		for (int i=0;i<sequence.length();i++) {
			// we've got to use 3 letter code for CYS, otherwise tinker takes the default CYX (which means cystein with disulfide bridge)
			if (sequence.charAt(i)=='C') {
				protInput.println("CYS");
			} else {
				protInput.println(sequence.charAt(i));
			}
		}
		protInput.println();
		protInput.println(CYCLISE_PROTEIN_STR);
		protInput.close();
		
		// logging output
		BufferedReader protOutput = new BufferedReader(new InputStreamReader(protProc.getInputStream()));
		String line;
		while((line = protOutput.readLine()) != null) {
			log.println(line);
			if (line.startsWith(TINKER_ERROR_STR)) {
				tinkerError = true;
			}
		}
		
		tinkerxyzout.renameTo(new File(outPath,outBasename+".xyz"));
		tinkerseqout.renameTo(new File(outPath,outBasename+".seq"));
		
		if (tinkerError) {
			log.flush();
			throw new TinkerError("Tinker error, revise log file. ");
		}
		
		log.flush();
	}

	/**
	 * Runs tinker's distgeom program capturing output with restrain violation statistics into member variable arrays
	 * that can be retrieved using the getters: getMaxLowerBoundViol, getMaxUpperBoundViol, getMaxLowerViol etc... 
	 * Two files are needed as input for distgeom: an xyz file and a key file, the latter is not passed but instead implicitely 
	 * defined by xyzFile: must be in same directory and must have same basename with extension .key 
	 * @param xyzFile
	 * @param outPath directory where output files will be written
	 * @param outBasename base name of the output files
	 * @param n number of models that we want distgeom to produce
	 * @param log a PrintWriter for logging output
	 * @throws TinkerError if an error seen in tinker's output
	 * @throws IOException
	 */
	private void runDistgeom(File xyzFile, String outPath, String outBasename, int n, PrintWriter log) throws TinkerError, IOException {

		checkDistgeomInput(xyzFile, outPath);
		// running distgeom program
		boolean tinkerError = false; // to store the exit state of the tinker program
		String cmdLine = distgeomProg+" "+xyzFile.getAbsolutePath()+" "+n+" "+dgeomParams;
		Process dgeomProc = Runtime.getRuntime().exec(cmdLine);
		// logging and capturing output
		BufferedReader dgeomOutput = new BufferedReader(new InputStreamReader(dgeomProc.getInputStream()));
		log.println("#cmd: "+cmdLine);
		tinkerError = parseDistgeomOutput(dgeomOutput, n, log);

		// getting names of tinker output files
		File[] tinkerout = new File[n+1];
		for (int i=1;i<=n;i++) {
			tinkerout[i] = getTinkerOutputFileName(xyzFile, String.format("%03d", i));
		}
		//renaming files to our chosen outBasename+ext
		for (int i=1;i<=n;i++) {
			tinkerout[i].renameTo(new File(outPath,outBasename+"."+String.format("%03d", i)));
		}
		// throwing exception if error string was caught in output
		if (tinkerError) {
			log.flush();
			throw new TinkerError("Tinker error, revise log file. ");
		}
		int exitValue = 1;
		try {
			exitValue = dgeomProc.waitFor();
		} catch (InterruptedException e) {
			throw new TinkerError("Unexpected error when waiting for distgeom to finish. Error: "+e.getMessage());
		}
		if (exitValue==137 || exitValue==139) {
			log.flush();
			throw new TinkerError("Distgeom was killed with exit code "+exitValue+". Not enough memory.");

		}
		// this is to catch all other possible errors not caught already by the parse of the error string in output
		else if (exitValue!=0) { 
			log.flush();
			throw new TinkerError("Distgeom exited with a non 0 exit code: "+exitValue+". Unknown error.");
		}
		
		log.flush();
	}
	
	/**
	 * Creates a shell script with a distgeom command to submit to the SGE system through qsub
	 * Requires java6 to be able to set the execute permission
	 * @param submitScript the script file
	 * @param xyzFile distgeom input xyz file
	 * @param n number of models for distgeom
	 * @throws FileNotFoundException
	 */
	private void createDistgeomSubmitScript (File submitScript, File xyzFile, int n) throws FileNotFoundException {
		PrintStream out = new PrintStream(submitScript);
		out.println("#!/bin/sh");
		out.println(distgeomProg+" "+xyzFile.getAbsolutePath()+" "+n+" "+dgeomParams);
		out.close();
		submitScript.setExecutable(true, true); // if script not executable sge can't submit it 
	}
	
	/**
	 * Runs distgeom by submitting it to a DRMAA system. Only SGE (qsub) supported at the moment
	 * @param xyzFile
	 * @param outPath
	 * @param outBasename
	 * @param n
	 * @return
	 * @throws TinkerError
	 * @throws IOException
	 * @throws DrmaaException
	 */
	private String runDistgeomDRMAA(File xyzFile, String outPath, String outBasename, int n) throws TinkerError, IOException, DrmaaException {
		checkDistgeomInput(xyzFile, outPath);		
		// running distgeom program
		File submitScript = new File(outPath,outBasename+".sh");
		if (!DEBUG) submitScript.deleteOnExit();
		createDistgeomSubmitScript(submitScript, xyzFile, n);
		
		JobTemplate jt = this.session.createJobTemplate();
		jt.setRemoteCommand(submitScript.getAbsolutePath());
		jt.setJobName(SGE_JOBS_PREFIX+outBasename);
		jt.setOutputPath(":"+outPath); // for some reason out/error paths require a ":" to start with
		jt.setErrorPath(":"+outPath); 
		jt.setNativeSpecification(SGE_QUEUE);
		
		String jobId = this.session.runJob(jt);
		
		this.session.deleteJobTemplate(jt);

		return jobId;
	}
	
	/**
	 * Checks that distgeom input paths and files are there
	 * @param xyzFile
	 * @param outPath
	 * @throws FileNotFoundException
	 */
	private void checkDistgeomInput(File xyzFile, String outPath) throws FileNotFoundException {
		if (!new File(outPath).exists()) {
			throw new FileNotFoundException("Specified directory "+outPath+" does not exist");
		}
		if (!xyzFile.exists()){
			throw new FileNotFoundException("Specified xyz file "+xyzFile.getAbsolutePath()+" does not exist");
		}
		String basename = getBasename(xyzFile);
		File keyFile = new File(xyzFile.getParent(),basename+".key");
		if (! keyFile.exists()) {
			throw new FileNotFoundException("Key file "+keyFile.getAbsolutePath()+" not present in input directory "+xyzFile.getParent());
		}
	}
	
	/**
	 * Initialises the drmaa session to communicate with the SGE queuing system
	 * @throws DrmaaException
	 */
	private void initDRMAASession() throws DrmaaException {
		SessionFactory factory = SessionFactory.getFactory();
		session = factory.getSession();

		session.init("");
	}
	
	/**
	 * Finalises the drmaa session
	 */
	private void finaliseDRMAASession()  {
		try {
			session.exit();
		} catch (DrmaaException e) {
			System.err.println("Couldn't finalise SGE session, some cleanup may have not been done. Error "+e.getMessage());
		}
	}
	
	/**
	 * Parses distgeom output of restrain violation statistics storing it into member variable arrays
	 * that can be retrieved using the getters: getMaxLowerBoundViol, getMaxUpperBoundViol, getMaxLowerViol etc...   
	 * @param dgeomOutput
	 * @param n
	 * @param log if not null the dgeomOutput is also written to log
	 * @return
	 * @throws IOException
	 */
	private boolean parseDistgeomOutput (BufferedReader dgeomOutput, int n, PrintWriter log) throws IOException {
		// initialising arrays were we store captured output data
		errorFunctionVal = new double[n+1];
		numUpperBoundViol = new int[n+1];
		numLowerBoundViol = new int[n+1];
		maxUpperBoundViol = new double[n+1];
		maxLowerBoundViol = new double[n+1];
		rmsBoundViol = new double[n+1];
		numUpperViol = new int[n+1];
		numLowerViol = new int[n+1];
		maxUpperViol = new double[n+1];
		maxLowerViol = new double[n+1];
		rmsRestViol = new double[n+1];

		boolean tinkerError = false;
		String line;
		int i=1;
		while((line = dgeomOutput.readLine()) != null) {
			if (log!=null) log.println(line);
			if (line.startsWith(TINKER_ERROR_STR)) {
				tinkerError = true;
			}
			Pattern p = Pattern.compile("^ Final Error Function Value :\\s+(\\d+\\.\\d+)");
			Matcher m = p.matcher(line);
			if (m.find()) {
				errorFunctionVal[i]=Double.parseDouble(m.group(1));
			}						
			p = Pattern.compile("^ Num Upper Bound Violations :\\s+(\\d+)");
			m = p.matcher(line);
			if (m.find()) {
				numUpperBoundViol[i]=Integer.parseInt(m.group(1));
			}
			p = Pattern.compile("^ Num Lower Bound Violations :\\s+(\\d+)");
			m = p.matcher(line);
			if (m.find()) {
				numLowerBoundViol[i]=Integer.parseInt(m.group(1));
			}
			p = Pattern.compile("^ Max Upper Bound Violation :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				maxUpperBoundViol[i]=Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^ Max Lower Bound Violation :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				maxLowerBoundViol[i]=Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^ RMS Deviation from Bounds :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				rmsBoundViol[i]=Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^ Num Upper Restraint Violations :\\s+(\\d+)");
			m = p.matcher(line);
			if (m.find()) {
				numUpperViol[i]=Integer.parseInt(m.group(1));
			}
			p = Pattern.compile("^ Num Lower Restraint Violations :\\s+(\\d+)");
			m = p.matcher(line);
			if (m.find()) {
				numLowerViol[i]=Integer.parseInt(m.group(1));
			}
			p = Pattern.compile("^ Max Upper Restraint Violation :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				maxUpperViol[i]=Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^ Max Lower Restraint Violation :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				maxLowerViol[i]=Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^ RMS Restraint Dist Violation :\\s+(\\d+\\.\\d\\d\\d\\d)");
			m = p.matcher(line);
			if (m.find()) {
				rmsRestViol[i]=Double.parseDouble(m.group(1));
				//System.out.println("Done model "+i+". Violations: "+numUpperViol[i]+" upper, "+numLowerViol[i]+" lower");
				i++;
				
			}
		}
		return tinkerError;
	}
	
	/**
	 * Copies inFile to outFile
	 * @param inFile
	 * @param outFile
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
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
	 * Copies inKeyFile to outKeyFile adding a RANDOMSEED tinker directive 
	 * @param inKeyFile
	 * @param outKeyFile
	 * @param seed the random seed to be written, must be below {@value #MAXSEED}
	 * @throws IOException
	 */
	private void addSeedToKeyFile(File inKeyFile, File outKeyFile, int seed) throws IOException {
		FileInputStream in = new FileInputStream(inKeyFile);
		PrintStream out = new PrintStream(outKeyFile);
		out.println("RANDOMSEED "+seed);
		int b;
		while ((b=in.read())!=-1) {
			out.write(b);
		}
		in.close();
		out.close();		
	}
	
	/**
	 * Kills all jobs submitted in this drmaa session.
	 * If some job doesn't exist or if it can't be killed for any other reason, 
	 * then it prints an error message.
	 */
	private void killJobs() {		
		try {
			this.session.control(Session.JOB_IDS_SESSION_ALL, Session.TERMINATE);
		} catch (DrmaaException e) {
			System.err.println("Couldn't kill job while cleaning up. Error: "+e.getMessage());
		}
	}
	
	/**
	 * Runs distgeom in parallel using one instance of distgeom per model.
	 * Each of the distgeom instances are run through a DRMAA system (only SGE supported at the moment)
	 * @param xyzFile
	 * @param outPath
	 * @param outBasename
	 * @param n
	 * @param log
	 * @throws TinkerError
	 * @throws IOException
	 */
	private void runParallelDistgeom(File xyzFile, String outPath, String outBasename, int n, PrintWriter log) throws TinkerError, IOException {
		Random ran = new Random();
		File[] nXyzOutFiles = new File[n+1]; 
		File[] nOLogFiles = new File[n+1];
		File[] nELogFiles = new File[n+1];
		try {
			initDRMAASession();
		} catch (DrmaaException e) {
			throw new TinkerError("Couldn't contact the SunGridEngine system. Error "+e.getMessage());
		}
		
		String[] jobIds = new String[n+1]; 
		for (int i=1;i<=n;i++) {
			String nBaseName = getBasename(xyzFile)+"_"+i;
			File keyFile = new File(xyzFile.getParent(),getBasename(xyzFile)+".key"); // input key file matching input xyz file
			File nXyzFile = new File(xyzFile.getParent(),nBaseName+".xyz");
			if (!DEBUG) nXyzFile.deleteOnExit();
			File nKeyFile = new File(xyzFile.getParent(),nBaseName+".key");
			if (!DEBUG) nKeyFile.deleteOnExit();
			copyFile(xyzFile,nXyzFile);
			addSeedToKeyFile(keyFile, nKeyFile, ran.nextInt(MAXSEED));
			
			// we predict what's going to be the xyz output file name (tinker will add a _2, _3, ... if file already exists with .001)
			nXyzOutFiles[i] = getTinkerOutputFileName(nXyzFile, "001");
			
			// submit jobs
			try {
				jobIds[i] = runDistgeomDRMAA(nXyzFile, outPath, nBaseName, 1);
			} catch (DrmaaException e) {
				// we are going for the moment for a conservative error handling. We don't allow any submission to fail
				// if only one of them fails we abort by throwing exception
				// we have to clean up the ones that did already start,
				killJobs();
				throw new TinkerError("SGE job "+nBaseName+" failed to run. "+e.getMessage());
			}
			
			// we now have to define the log file names from the name that sge gives the out/err files
			// we use SGE_JOBS_PREFIX+nBaseName because is the value to which we set the name of the job with jt.setName() in runDistgeomDRMAA
			nOLogFiles[i] = new File(outPath,SGE_JOBS_PREFIX+nBaseName+".o"+jobIds[i]);
			if (!DEBUG) nOLogFiles[i].deleteOnExit();
			nELogFiles[i] = new File(outPath,SGE_JOBS_PREFIX+nBaseName+".e"+jobIds[i]);
			if (!DEBUG) nELogFiles[i].deleteOnExit();
						
		}		
		
		// 1 first check jobs are finished and successful
		try {
			// synchronize (waits for all jobs) seems to be broken, it crashes the JVM! don't use it!
			//List<String> jobIdsList = Arrays.asList(jobIds);
			//this.session.synchronize(jobIdsList, PARALLEL_JOBS_TIMEOUT, true);
			//int status = this.session.getJobProgramStatus(jobIds[i]);
			
			// we wait for one job at a time, so one not finishing will block all the others
			// TODO implement some better approach, ideas: 
			// a) variable timeout: estimate it from the size of the protein or the runtime of the fastest job + 20% margin
			// b) multithreaded waiting and then allow for some failure rate (like 10%)
			for (int i=1;i<=n;i++) {
				JobInfo info = this.session.wait(jobIds[i],PARALLEL_JOBS_TIMEOUT);
				int exitStatus = 999;
				if (info.hasExited()) {
					exitStatus = info.getExitStatus();
				} else {
					killJobs(); // we cleanup
					throw new TinkerError("Job "+jobIds[i]+" exited abnormally or didn't exit yet. Couldn't produce output file "+nXyzOutFiles[i]);
				}
				if (exitStatus!=0) {
					killJobs();
					throw new TinkerError("Job "+jobIds[i]+" finished with exit status "+exitStatus+". Couldn't produce output file "+nXyzOutFiles[i]);
				} 
			}
		} catch (DrmaaException e) {
			throw new TinkerError(e);
		}
		
		// 2 if jobs successful check output files are there
		// we've got to retry a few times with waiting of 2s in between
		// to leave a bit of time between jobs finish and check of output, 
		// file system seems to have problems to detect the files inmediately after jobs finished
		boolean allFilesFound = false;
		int retries = 0;
		while (!allFilesFound && retries<MAX_RETRIES_FIND_OUTPUT) {
			retries++;
			try {
				Thread.sleep(RETRY_TIME_FIND_OUTPUT);
			} catch (InterruptedException e) {
				System.err.println("Unexpected error: couldn't sleep to wait for output files");
				System.exit(1);
			}
			for (int i=1;i<=n;i++) {
				if (!nXyzOutFiles[i].exists()) {
					//throw new TinkerError("All jobs finished but output file "+nXyzOutFiles[i]+" not present");
					break;
				}
				allFilesFound = true;
			}
		}
		if (!allFilesFound) 
			throw new TinkerError("All jobs finished but some output files were not found");
		
		// now we need to convert the split output into a normal (non-parallel output)
		double[] errorFunctionVal = new double[n+1];
		int[] numUpperBoundViol = new int[n+1];
		int[] numLowerBoundViol = new int[n+1];
		double[] maxUpperBoundViol = new double[n+1];
		double[] maxLowerBoundViol = new double[n+1];
		double[] rmsBoundViol = new double[n+1];
		int[] numUpperViol = new int[n+1];
		int[] numLowerViol = new int[n+1];
		double[] maxUpperViol = new double[n+1];
		double[] maxLowerViol = new double[n+1];
		double[] rmsRestViol = new double[n+1];
		
		for (int i=1;i<=n;i++) {
			// 1 rename all bn_i.001 files to bn.iii
			nXyzOutFiles[i].renameTo(new File(outPath,outBasename+String.format(".%03d",i)));
			// 2 parse log to get distgeom statistics and write individual logs to main log file
			BufferedReader l = new BufferedReader(new FileReader(nOLogFiles[i]));
			parseDistgeomOutput(l, 1, log);
			// 3 get the individual violation statistics and put them into a new array 
			errorFunctionVal[i] = this.getErrorFunctionVal()[1];
			numUpperBoundViol[i] = this.getNumUpperBoundViol()[1];
			numLowerBoundViol[i] = this.getNumLowerBoundViol()[1];
			maxUpperBoundViol[i] = this.getMaxUpperBoundViol()[1];
			maxLowerBoundViol[i] = this.getMaxLowerBoundViol()[1];
			rmsBoundViol[i] = this.getRmsBoundViol()[1];
			numUpperViol[i] = this.getNumUpperViol()[1];
			numLowerViol[i] = this.getNumLowerViol()[1];
			maxUpperViol[i] = this.getMaxUpperViol()[1];
			maxLowerViol[i] = this.getMaxLowerViol()[1];
			rmsRestViol[i] = this.getRmsRestViol()[1];
		}
		this.errorFunctionVal = errorFunctionVal;
		this.numUpperBoundViol = numUpperBoundViol;
		this.numLowerBoundViol = numLowerBoundViol;
		this.maxUpperBoundViol = maxUpperBoundViol;
		this.maxLowerBoundViol = maxLowerBoundViol;
		this.rmsBoundViol = rmsBoundViol;
		this.numUpperViol = numUpperViol;
		this.numLowerViol = numLowerViol;
		this.maxUpperViol = maxUpperViol;
		this.maxLowerViol = maxLowerViol;
		this.rmsRestViol = rmsRestViol;

		finaliseDRMAASession();
	}
	
	/**
	 * Runs tinker's xyzpdb program to convert a given xyzFile (needing also a seqFile) to a pdbFile
	 * @param xyzFile
	 * @param seqFile
	 * @param pdbFile
	 * @param log A PrintWriter for logging output
	 * @throws IOException 
	 * @throws TinkerError If an error seen in tinker's output
	 */
	private void runXyzpdb(File xyzFile, File seqFile, File pdbFile, PrintWriter log) throws IOException, TinkerError {
		boolean tinkerError = false; // to store the exit state of the tinker program
		if (!xyzFile.exists()){
			throw new FileNotFoundException("Specified xyz file "+xyzFile.getAbsolutePath()+" does not exist");
		}
		if (!seqFile.exists()){
			throw new FileNotFoundException("Specified seq file "+seqFile.getAbsolutePath()+" does not exist");
		}
		
		String basename = getBasename(xyzFile);
		File tmpSeqFile = new File(seqFile.getParent(),basename+".seq");
				
		// if seqFile doesn't follow the naming convention (basename of xyzFile+seq extension) that tinker expects, we copy it to tmpSeqFile (which has right name) 
		if (!tmpSeqFile.equals(seqFile)) {
			FileChannel srcChannel = new FileInputStream(seqFile).getChannel();
			FileChannel dstChannel = new FileOutputStream(tmpSeqFile).getChannel();
			dstChannel.transferFrom(srcChannel, 0, srcChannel.size());
			srcChannel.close();
			dstChannel.close();
			// if we copied then that means tmpSeqFile is different from seqFile and thus we want to delete the tmp file on exit
			tmpSeqFile.deleteOnExit(); 
		}
        
		File tinkerpdbout = getTinkerOutputFileName(xyzFile, "pdb");
		
		// running tinker's xyzpdb
		// beware: it takes as a silent input the seq file seqFile (or tmpSeqFile if the input seqFile didn't have the right name)
		String cmdLine = xyzpdbProg+" "+xyzFile.getAbsolutePath()+" "+forceFieldFileName;
		Process xyzpdbProc = Runtime.getRuntime().exec(cmdLine);

		// logging output
		BufferedReader xyzpdbOutput = new BufferedReader(new InputStreamReader(xyzpdbProc.getInputStream()));
		String line;
		log.println("#cmd: "+cmdLine);
		while((line = xyzpdbOutput.readLine()) != null) {
			log.println(line);
			if (line.startsWith(TINKER_ERROR_STR)) {
				tinkerError = true;
			}
		}
		
		tinkerpdbout.renameTo(pdbFile);
		
		if (tinkerError) {
			log.flush();
			throw new TinkerError("Tinker error while running xyzpdb, revise log file.");
		}
		
		log.flush();
	}
	
	/**
	 * Runs tinker's pdbxyz program to convert a pdbFile to a xyzFile
	 * @param pdbFile
	 * @param xyzFile
	 * @param log A PrintWriter for logging output	
	 * @throws IOException
	 * @throws TinkerError If an error seen in tinker's output
	 */
	private void runPdbxyz(File pdbFile, File xyzFile, PrintWriter log) throws IOException, TinkerError{
		boolean tinkerError = false; // to store the exit state of the tinker program
		if (!pdbFile.exists()){
			throw new FileNotFoundException("Specified pdb file "+pdbFile.getAbsolutePath()+" does not exist");
		}
		File tinkerxyzout = getTinkerOutputFileName(pdbFile, "xyz");
		File tinkerseqout = getTinkerOutputFileName(pdbFile, "seq");
		// running tinker's pdbxyz
		String cmdLine = pdbxyzProg+" "+pdbFile.getAbsolutePath()+" "+forceFieldFileName;
		Process pdbxyzProc = Runtime.getRuntime().exec(cmdLine);

		// logging output
		BufferedReader pdbxyzOutput = new BufferedReader(new InputStreamReader(pdbxyzProc.getInputStream()));
		String line;
		log.println("#cmd: "+cmdLine);
		while((line = pdbxyzOutput.readLine()) != null) {
			log.println(line);
			if (line.startsWith(TINKER_ERROR_STR)) {
				tinkerError = true;
			}
			// TODO: we want here to catch cases where there is missing data in input and tinker prompts for something else (e.g. if PDB file has alt locs, tinker asks the alt loc code we want to choose)
			// PROBLEM: (that's why following code doesn't work) line hasn't been written fully so readLine() didn't read it yet, how to get around it? reading character by character?
			//if (line.startsWith(" Enter")) { // i.e. tinker program is waiting for a prompt
			//	log.flush();
			//	throw new TinkerError("Some missing parameter to pdbxyz. Tinker is waiting for data: \""+line+"\"");
			//}
		}
		tinkerxyzout.renameTo(xyzFile);
		// for the seq file we rename to a file with seq extension with the same path and basename as xyzFile
		tinkerseqout.renameTo(new File(xyzFile.getAbsolutePath().substring(0, xyzFile.getAbsolutePath().lastIndexOf("."))+".seq"));
		if (tinkerError) {
			// we flush so that we are sure the log file gets written, we don't want to close it though, because the calling program might continue using the PrintWriter
			log.flush(); 
			throw new TinkerError("Tinker error while running pdbxyz, revise log file.");
		}
		
		log.flush();
	}
	
	/**
	 * Gets the basename of given file: no path, chopped off bit after last dot.
	 * @param file
	 * @return
	 */
	private String getBasename(File file) {
		String baseName = file.getName();
		// for files with no extension we don't have to do more, if file has a extension we take everything up to the dot
		if (baseName.contains(".")) { 
			baseName = baseName.substring(0, baseName.lastIndexOf("."));
		}
		return baseName;
	}
	
	/**
	 * Run tinker's minimize program for the given xyz file. The final minimized structure is written back to the same file.
	 * @param xyzFile
	 * @param rmsGradient  target rms gradient value for which minimization will terminate (in kcal/mole/A)
	 * @param log
	 * @return  the final energy of the minimized structure
	 * @throws IOException  if temporary/result files can't be accessed
	 * @throws TinkerError  if tinker finishes in error status or prints an error in the output
	 */
	private double runMinimize(File xyzFile, double rmsGradient, PrintWriter log) throws IOException, TinkerError {
		double energy = 0.0;
		boolean tinkerError = false;
		// minimize and get final energy
		File xyzOutFile = getTinkerOutputFileName(xyzFile, "xyz");
		String cmdLine = minimizeProg + " " + xyzFile.getAbsolutePath() + " " + forceFieldFileName + " " + rmsGradient;
		Process minimizeProc = Runtime.getRuntime().exec(cmdLine);
		// logging and capturing output
		BufferedReader minimizeOutput = new BufferedReader(new InputStreamReader(minimizeProc.getInputStream()));
		String line;
		log.println("#cmd: "+cmdLine);
		while((line = minimizeOutput.readLine()) != null) {
			log.println(line);
			if (line.startsWith(TINKER_ERROR_STR)) {
				tinkerError = true;
			}
			if (line.startsWith(CHECKXYZ_WARNING)) {
				log.flush();
				throw new TinkerError("minimize gave a warning about the input xyz file. Revise log file");
			}
			Pattern p = Pattern.compile("^ Final Function Value :\\s+(-?\\d+\\.\\d+)");
			Matcher m = p.matcher(line);
			if (m.find()) {
				energy = Double.parseDouble(m.group(1));
			}			
		}
		if (tinkerError) {
			log.flush();
			throw new TinkerError("Tinker error while running minimize. Revise log file. ");
		}
		try {
			int exitValue = minimizeProc.waitFor(); 
			if (exitValue!=0) {
				log.flush();
				throw new TinkerError("Tinker exited with "+exitValue+" status. Some error occurred while running minimize. Revise log file.");
			}
		} catch (InterruptedException e) {
			throw new TinkerError(e);
		}
		
		// rename the output file 
		xyzOutFile.renameTo(xyzFile);
		
		log.flush();
		
		return energy;
	}
	
	/**
	 * Run tinker's analyze program to get the potential energy of the given xyz file.
	 * @param xyzFile
	 * @param mode
	 * @param log
	 * @return
	 * @throws IOException
	 * @throws TinkerError
	 */
	private double runAnalyze(File xyzFile, PrintWriter log) throws IOException, TinkerError {
		double energy = 0.0;
		boolean tinkerError = false;
		// get energy by running analyze

		String cmdLine = analyzeProg + " " + xyzFile.getAbsolutePath() + " " + forceFieldFileName + " " + ANALYZE_ENERGY_MODE;
		Process analyzeProc = Runtime.getRuntime().exec(cmdLine);
		// logging and capturing output
		BufferedReader analyzeOutput = new BufferedReader(new InputStreamReader(analyzeProc.getInputStream()));
		String line;
		log.println("#cmd: "+cmdLine);
		while((line = analyzeOutput.readLine()) != null) {
			log.println(line);
			if (line.startsWith(TINKER_ERROR_STR)) {
				tinkerError = true;
			}
			Pattern p = Pattern.compile("^ Total Potential Energy :\\s+(-?\\d+\\.\\d+)");
			Matcher m = p.matcher(line);
			if (m.find()) {
				energy = Double.parseDouble(m.group(1));
			}			
		}
		if (tinkerError) {
			log.flush();
			throw new TinkerError("Tinker error while running analyze. Revise log file. ");
		}
		try {
			int exitValue = analyzeProc.waitFor(); 
			if (exitValue!=0) {
				log.flush();
				throw new TinkerError("Tinker exited with "+exitValue+" status. Some error occurred while running analyze. Revise log file.");
			}
		} catch (InterruptedException e) {
			throw new TinkerError(e);
		}
		
		log.flush();
		
		return energy;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/** 
	 * Returns the index of the structure with the least bound violations.
	 * Note: The structure can then be retrieved using getStructure(i);
	 * @returns
	 */ 
	public int pickByLeastBoundViols() {
		int minIdx = 1;
		int minViols = numLowerBoundViol[1] + numUpperBoundViol[1];
		for (int i = 2; i <= lastNumberOfModels; i++) {
			int score = numLowerBoundViol[i] + numUpperBoundViol[i];
			if(score < minViols) {
				minIdx = i;
				minViols = score;
			}
		}
		return minIdx;
	}
	
	// reconstruction
	
	/** 
	 * Reconstructs the given graph(s) and returns a putative good model as a Pdb object
	 * using the default forceConstant.
	 * The model is picked based on the number of bound violations reported by Tinker.
	 * See reconstruct() for more details.
	 * @param sequence sequence of the structure to be generated
	 * @param graphs array of graph objects containing constraints
	 * @param phiPsiConsensus Map containing consensus phi/psi angle for each position of the sequence, 
	 * if null no phi/psi restraints will be applied
	 * @param forceTransOmega to force omega angles into the trans conformation (180) 
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param parallel whether to run the parallel version or not (needs a sun grid engine cluster to work) 
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 * @returns A pdb object containg the generated structure
	 */
	public Pdb reconstruct(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, boolean forceTransOmega, int numberOfModels, boolean parallel) 
	throws TinkerError, IOException {
		
		return reconstruct(sequence, graphs, phiPsiConsensus, forceTransOmega, numberOfModels, 
				DEFAULT_FORCECONSTANT_DISTANCE, DEFAULT_FORCECONSTANT_TORSION, parallel);
	}	
	
	/** 
	 * Reconstructs the given graph(s) and returns a putative good model as a Pdb object
	 * using the default forceConstant.
	 * The model is picked based on the number of bound violations reported by Tinker.
	 * See reconstruct() for more details.
	 * @param sequence sequence of the structure to be generated
	 * @param graphs array of graph objects containing constraints
	 * @param phiPsiConsensus Map containing consensus phi/psi angle for each position of the sequence, 
	 * if null no phi/psi restraints will be applied
	 * @param forceTransOmega to force omega angles into the trans conformation (180) 
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param parallel whether to run the parallel version or not (needs a sun grid engine cluster to work) 
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 * @returns A pdb object containg the generated structure
	 */
	public Pdb reconstructFast(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, boolean forceTransOmega, int numberOfModels, boolean parallel) 
	throws TinkerError, IOException {
		
		return reconstructFast(sequence, graphs, phiPsiConsensus, forceTransOmega,
				numberOfModels, DEFAULT_FORCECONSTANT_DISTANCE, DEFAULT_FORCECONSTANT_TORSION, parallel);
	}

	/** 
	 * Reconstructs the given graph(s) and returns a putative good model as a Pdb object.
	 * The model is picked based on the number of bound violations reported by Tinker.
	 * See reconstruct() for more details.
	 * @param sequence sequence of the structure to be generated
	 * @param graphs array of graph objects containing constraints
	 * @param phiPsiConsensus Map containing consensus phi/psi angle for each position of the sequence, 
	 * if null no phi/psi restraints will be applied
	 * @param forceTransOmega to force omega angles into the trans conformation (180) 
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param forceConstantDist the force constant to be used for all our given distance restraints
	 * @param forceConstantTorsion the force constant to be used for all our given torsion restraints
	 * @param parallel whether to run the parallel version or not (needs a sun grid engine cluster to work) 
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 * @returns a Pdb object containg the generated structure
	 */
	public Pdb reconstruct(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, boolean forceTransOmega,
			int numberOfModels, double forceConstantDist, double forceConstantTorsion, boolean parallel) 
	throws TinkerError, IOException {
		
		Pdb resultPdb = null;
		
		String outputDir = System.getProperty("java.io.tmpdir");
		String baseName = Long.toString(System.currentTimeMillis());	// some hopefully unique basename
		boolean cleanUp = true;						
		
		reconstruct(sequence, graphs, phiPsiConsensus, forceTransOmega, numberOfModels, 
				forceConstantDist, forceConstantTorsion, true, outputDir, baseName, cleanUp, parallel);
		int pickedIdx = pickByLeastBoundViols();
		resultPdb = getStructure(pickedIdx);
		
		return resultPdb;
	}
	
	/** 
	 * Reconstructs the given graph(s) and returns a putative good model as a Pdb object.
	 * The model is picked based on the number of bound violations reported by Tinker.
	 * See reconstruct() for more details.
	 * @param sequence sequence of the structure to be generated
	 * @param graphs array of graph objects containing constraints
	 * @param phiPsiConsensus Map containing consensus phi/psi angle for each position of the sequence, 
	 * if null no phi/psi restraints will be applied
	 * @param forceTransOmega to force omega angles into the trans conformation (180) 
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param forceConstantDist the force constant to be used for all our given distance restraints
	 * @param forceConstantTorsion the force constant to be used for all our given torsion restraints
	 * @param parallel whether to run the parallel version or not (needs a sun grid engine cluster to work) 
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 * @returns a Pdb object containg the generated structure
	 */
	public Pdb reconstructFast(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, boolean forceTransOmega, 
			int numberOfModels, double forceConstantDist, double forceConstantTorsion, boolean parallel) 
	throws TinkerError, IOException {
		
		Pdb resultPdb = null;
		
		String outputDir = System.getProperty("java.io.tmpdir");
		String baseName = Long.toString(System.currentTimeMillis());	// some hopefully unique basename
		boolean cleanUp = true;						
		
		reconstruct(sequence, graphs, phiPsiConsensus, forceTransOmega, numberOfModels, forceConstantDist, forceConstantTorsion, false, outputDir, baseName, cleanUp, parallel);
		int pickedIdx = pickByLeastBoundViols();
		resultPdb = getStructure(pickedIdx);
		
		return resultPdb;
	}

	/**
	 * Reconstruct the given graph(s). Edge types and distance cutoffs are taken from the graph object(s).
	 * Simulated annealing is used for refinement
	 * If multiple graphs are specified, tinker simply takes the union of constraints for the reconstruction. 
	 * Output files starting with baseName will be written to the given outputDir.
	 * If cleanUp is true, all output files are marked to be deleted on shutdown of the virtual machine.
	 * The results of the last call to a reconstruct... method can be retrieved by the
	 * getMax..., getNum..., getRms... and getErrorFunctionVal() methods.  
	 * @param sequence sequence of the structure to be generated
	 * @param graphs array of graph objects containing constraints
	 * @param phiPsiConsensus Map containing consensus phi/psi angle for each position of the sequence, 
	 * if null no phi/psi restraints will be applied
	 * @param forceTransOmega to force omega angles into the trans conformation (180) 
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param forceConstantDist the force constant to be used for all our given distance restraints
	 * @param forceConsantTorsion the force constant to be used for all our given torsion restraints
	 * @param outputDir the directory where the temporary and result files will be written to
	 * @param baseName the basename of the temporary and result files
	 * @param cleanUp whether to mark all created files to be deleted on shutdown
	 * @param parallel whether to run the parallel version or not (needs a sun grid engine cluster to work) 
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 */
	public void reconstruct(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, boolean forceTransOmega, 
			int numberOfModels, double forceConstantDist, double forceConsantTorsion, 
			String outputDir, String baseName, boolean cleanUp, boolean parallel) 
	throws TinkerError, IOException {
		
		reconstruct(sequence, graphs, phiPsiConsensus, forceTransOmega, numberOfModels, forceConstantDist, forceConsantTorsion, true, outputDir, baseName, cleanUp, parallel);
	}
	
	/**
	 * Reconstruct the given graph(s). Edge types and distance cutoffs are taken from the graph object(s).
	 * Minimization is used for refinement (a lot faster than simulated annealing)
	 * If multiple graphs are specified, tinker simply takes the union of constraints for the reconstruction. 
	 * Output files starting with baseName will be written to the given outputDir.
	 * If cleanUp is true, all output files are marked to be deleted on shutdown of the virtual machine.
	 * The results of the last call to a reconstruct... method can be retrieved by the
	 * getMax..., getNum..., getRms... and getErrorFunctionVal() methods.  
	 * @param sequence sequence of the structure to be generated
	 * @param graphs array of graph objects containing constraints
	 * @param phiPsiConsensus Map containing consensus phi/psi angle for each position of the sequence, 
	 * if null no phi/psi restraints will be applied
	 * @param forceTransOmega to force omega angles into the trans conformation (180) 
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param forceConstantDist the force constant to be used for all our given distance restraints
	 * @param forceConstantTorsion the force constant to be used for all our given torsion restraints
	 * @param outputDir the directory where the temporary and result files will be written to
	 * @param baseName the basename of the temporary and result files
	 * @param cleanUp whether to mark all created files to be deleted on shutdown
	 * @param parallel whether to run the parallel version or not (needs a sun grid engine cluster to work) 
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 */	
	public void reconstructFast(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, boolean forceTransOmega,  
			int numberOfModels, double forceConstantDist, double forceConstantTorsion, 
			String outputDir, String baseName, boolean cleanUp, boolean parallel) 
	throws TinkerError, IOException {
		
		reconstruct(sequence, graphs, phiPsiConsensus, forceTransOmega, 
				numberOfModels, 
				forceConstantDist, forceConstantTorsion, 
				false, outputDir, baseName, cleanUp, parallel);
	}
	
	/** 
	 * Reconstruct the given graph(s). Edge types and distance cutoffs are taken from the graph object(s).
	 * If multiple graphs are specified, tinker simply takes the union of constraints for the reconstruction. 
	 * Output files starting with baseName will be written to the given outputDir.
	 * If cleanUp is true, all output files are marked to be deleted on shutdown of the virtual machine.
	 * The results of the last call to a reconstruct... method can be retrieved by the
	 * getMax..., getNum..., getRms... and getErrorFunctionVal() methods. 
	 * @param sequence sequence of the structure to be generated
	 * @param graphs array of graph objects containing constraints
	 * @param phiPsiConsensus Map containing consensus phi/psi angle for each position of the sequence, 
	 * if null no phi/psi restraints will be applied
	 * @param forceTransOmega to force omega angles into the trans conformation (180)
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param forceConstantDist the force constant to be used for all our given distance restraints
	 * @param forceConstantTorsion the force constant to be used for all our given torsion restraints
	 * @param annealing if true simulated annealing is used for refinement (slow), if false minimization is used for refinement (fast) 
	 * @param outputDir the directory where the temporary and result files will be written to
	 * @param baseName the basename of the temporary and result files
	 * @param cleanUp whether to mark all created files to be deleted on shutdown
	 * @param parallel whether to run the parallel version or not (needs a sun grid engine cluster to work)
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 */
	private void reconstruct(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, boolean forceTransOmega,
			int numberOfModels, double forceConstantDist, double forceConstantTorsion, 
			boolean annealing, String outputDir, String baseName, boolean cleanUp, 
			boolean parallel) 
	throws TinkerError, IOException {
		
		if (!annealing) dgeomParams = DGEOM_DEFAULT_PARAMS+REFINE_VIA_MINIMIZATION;
		this.forceConstant = forceConstantDist; 
		
		// defining files
		File prmFile = new File(this.forceFieldFileName);					// don't delete this one
		File xyzFile = new File(outputDir,baseName+".xyz");
		File seqFile = new File(outputDir,baseName+".seq");
		File pdbFile = new File(outputDir,baseName+".pdb");
		File keyFile = new File(outputDir,baseName+".key");
		//File reportFile = new File(outputDir,baseName+".report");
		File logFile = new File(outputDir,baseName+".tinker.log");		
		PrintWriter log = new PrintWriter(new FileOutputStream(logFile));
		
		
		// make sure that files are deleted on exit
		if(cleanUp) {
			//logFile.deleteOnExit();
			xyzFile.deleteOnExit();
			seqFile.deleteOnExit();
			pdbFile.deleteOnExit();
			keyFile.deleteOnExit();
			//reportFile.deleteOnExit();
		}
		
		// 1. run tinker's protein program (create unfolded protein chain)	
		runProtein(sequence, outputDir, baseName, log);

		// 1a. convert xyz file to pdb to be able to map atom serials after
		runXyzpdb(xyzFile, seqFile, pdbFile, log);

		// 2. creating constraints into key file
		ConstraintsMaker cm = null;
		try {
			cm = new ConstraintsMaker(pdbFile,xyzFile,prmFile,DEFAULT_FF_FILE_TYPE,keyFile);
		} catch(PdbLoadError e) {
			throw new TinkerError(e);
		}
		for(RIGraph graph:graphs) {
			cm.createDistanceConstraints(graph,forceConstantDist);
		}
		
		if (phiPsiConsensus!=null) {
			cm.createPhiPsiConstraints(phiPsiConsensus, forceConstantTorsion);
		}
		
		if (forceTransOmega) {
			cm.createOmegaConstraints(forceConstantTorsion);
		}
		
		cm.closeKeyFile();

		// 3. run tinker's distgeom
		if (parallel) {
			runParallelDistgeom(xyzFile, outputDir, baseName, numberOfModels, log);
		} else {
			runDistgeom(xyzFile, outputDir, baseName, numberOfModels, log);
		}

		// 4. converting xyz output files to pdb files and calculating rmsds
		for (int i = 1; i<=numberOfModels; i++) {
			String ext = new Formatter().format(".%03d",i).toString(); // 001, 002, 003, ...
			File outputXyzFile = new File(outputDir, baseName+ext);
			File outputPdbFile = new File(outputDir, baseName+ext+".pdb");
			if(cleanUp) {
				outputXyzFile.deleteOnExit();
				outputPdbFile.deleteOnExit();
			}
			runXyzpdb(outputXyzFile, seqFile, outputPdbFile, log);
		}					
		log.close();
		
		this.lastOutputDir = outputDir;
		this.lastBaseName = baseName;
		this.lastNumberOfModels = numberOfModels;
	}
	
	/**
	 * Reconstructs the given contact map with the given number of models returning
	 * an average graph containing a union of edges of all reconstructed models. The edges are weighted 
	 * by the fraction of occurrence in the graphs derived from the reconstructed models.
	 * The idea is that this method geometrizes a given graph which is not necessarily embeddable. From the
	 * ensemble output one would see which edges were physically possible compare to the input contact map.
	 * We are not sure if this works it's just an idea, but if it does would be great!   
	 * @param graph
	 * @param numberOfModels
	 * @param fastReconstruct
	 * @param outputDir the output directory, if left null then a temp directory is used and all output removed
	 * @param baseName the base name of the output files, if left null then a random name is used and all output removed
	 * @param parallel whether to run the parallel version or not (needs a sun grid engine cluster to work) 
	 * @return
	 * @throws IOException 
	 * @throws TinkerError 
	 */
	public RIGraph geometrizeContactMap(RIGraph graph, int numberOfModels, boolean fastReconstruct, String outputDir, String baseName, boolean parallel) 
	throws TinkerError, IOException {
		boolean cleanUp = false;
		if (outputDir == null || baseName == null) {
			cleanUp = true;
		}
		RIGraph[] graphs = {graph};
		if (fastReconstruct) {
			reconstructFast(graph.getSequence(), graphs, null, false, 
					numberOfModels, DEFAULT_FORCECONSTANT_DISTANCE, DEFAULT_FORCECONSTANT_TORSION, 
					outputDir, baseName, cleanUp, parallel);
		} else {
			reconstruct(graph.getSequence(), graphs, null, false, 
					numberOfModels, DEFAULT_FORCECONSTANT_DISTANCE, DEFAULT_FORCECONSTANT_TORSION,
					outputDir, baseName, cleanUp, parallel);
		}
		
		RIGEnsemble ensemble = new RIGEnsemble(graph.getContactType(), graph.getCutoff());
		for (File file:getLastModelFiles()) {
			PdbfilePdb pdb = new PdbfilePdb(file.getAbsolutePath());
			try {
				pdb.load(Pdb.NULL_CHAIN_CODE);
			} catch (PdbLoadError e) {
				throw new TinkerError(e);
			} 
			ensemble.addRIG(pdb.get_graph(graph.getContactType(), graph.getCutoff()));

		}
		
		GraphAverager ga = null;
		try {
			ga = new GraphAverager(ensemble);
		} catch (GraphAveragerError e) {
			throw new TinkerError(e);
		}
		return ga.getAverageGraph();

	}
	
	// retrieving results
	
	/**
	 * Returns the number of models generated in the last reconstruction run.
	 * @return the number of models or 0 if no run was completed yet
	 */
	public int getLastNumberOfModels() {
		return lastNumberOfModels;
	}

	/**
	 * Returns the directory where temporary and output files of the last reconstruction run were written to.
	 * @return the output directory or null if no run was completed yet
	 */
	public String getLastOutputDir() {
		return lastOutputDir;
	}
	
	/**
	 * Returns the basename pf the temporary and output files of the last reconstruction run.
	 * @return the base name or null if no run was completed yet
	 */
	public String getLastBaseName() {
		return lastBaseName;
	}
	
	/**
	 * Returns the force constant used in the last reconstruction run.
	 * @return the force constant or -1 if no run was completed yet
	 */
	public double getLastForceConstant() {
		return forceConstant;
	}	
	
	/**
	 * Returns one of the structures generated in the last reconstruction run as a pdb object.
	 * Structures are numbered from 1 to getLastNumberOfModels().
	 * @param i the number of the model to be returned.
	 * @return the generated structure as a pdb object
	 */
	public Pdb getStructure(int i) throws TinkerError {
		Pdb resultPdb = null;
		
		String ext = String.format(".%03d",i); // 001, 002, 003, ...
		File resultPdbFile = new File(lastOutputDir, lastBaseName+ext+".pdb");
		try {
			resultPdb = new PdbfilePdb(resultPdbFile.getAbsolutePath());
			resultPdb.load(DEFAULT_RECONSTR_CHAIN_CODE); 
		} catch(PdbLoadError e) {
			throw new TinkerError(e);
		}
		
		return resultPdb;
	}
	
	public double[] getErrorFunctionVal() {
		return errorFunctionVal;
	}

	public double[] getMaxLowerBoundViol() {
		return maxLowerBoundViol;
	}

	public double[] getMaxLowerViol() {
		return maxLowerViol;
	}

	public double[] getMaxUpperBoundViol() {
		return maxUpperBoundViol;
	}

	public double[] getMaxUpperViol() {
		return maxUpperViol;
	}

	public int[] getNumLowerBoundViol() {
		return numLowerBoundViol;
	}

	public int[] getNumLowerViol() {
		return numLowerViol;
	}

	public int[] getNumUpperBoundViol() {
		return numUpperBoundViol;
	}

	public int[] getNumUpperViol() {
		return numUpperViol;
	}

	public double[] getRmsBoundViol() {
		return rmsBoundViol;
	}

	public double[] getRmsRestViol() {
		return rmsRestViol;
	}
	
	/**
	 * Returns a vector of model file objects for the last reconstruction run.
	 * Note: To get a particular structures as a Pdb object, use getStructure();
	 * @return a vector of file objects
	 */
	public Vector<File> getLastModelFiles() {
		Vector<File> ret = new Vector<File>(lastNumberOfModels);
		for(int i=1; i <= lastNumberOfModels; i++) {
			String ext = String.format(".%03d",i); // 001, 002, 003, ...
			File resultPdbFile = new File(lastOutputDir, lastBaseName + ext + ".pdb");
			ret.add(resultPdbFile);
		}
		return ret;
	}
	
	/**
	 * Calculates a vector of GDT_TS scores for the models of the last run versus the given structure using maxcluster.
	 * @param nativePdbFileName the name of the pdb file with the native structure
	 * @param maxClusterExecutable the path to the maxcluster executable
	 * @return an array of gdt scores
	 * @throws IOException if anything goes wrong
	 */
	public double[] getGdtsToNative(String nativePdbFileName, String maxClusterExecutable) throws IOException {
		double[] gdtScores = new double[lastNumberOfModels];
		
		File file = new File(nativePdbFileName);
		if(!file.canRead()) {
			throw new FileNotFoundException("Could not read from file " + nativePdbFileName);
		}
		
		MaxClusterRunner maxCluster = new MaxClusterRunner(maxClusterExecutable);
		Vector<File> modelFiles = getLastModelFiles();
		
		for(int i=1; i <= lastNumberOfModels; i++) {
			String modelFileName = modelFiles.get(i-1).getAbsolutePath(); 
			gdtScores[i-1] = maxCluster.calculatePairwiseScore(modelFileName, nativePdbFileName, MaxClusterRunner.ScoreType.GDT);
		}
		return gdtScores;
	}
	
	/**
	 * Given a pdb file computes its energy for the currently set forceField
	 * Energy minimization of the structure will be performed previously to 
	 * calculating the energy.
	 * The minimized pdb file is written to same directory as input pdb file with 
	 * same base name and extension ".min.pdb"
	 * Temporary files .xyz and .seq will be removed on exit.
	 * @param pdbFile
	 * @param log
	 * @return
	 * @throws TinkerError
	 * @throws IOException 
	 */
	public double minimize(File pdbFile, double rmsGradient, PrintWriter log) throws TinkerError, IOException {

		boolean cleanUp = true; // for debugging
		
		String baseName = getBasename(pdbFile);
		String outputDir = pdbFile.getParent();
		File xyzFile = new File(outputDir,baseName+".xyz");
		File seqFile = new File(outputDir,baseName+".seq"); // this file is created implicitely by pdbxyz
		if (cleanUp) {
			xyzFile.deleteOnExit();
			seqFile.deleteOnExit();
		}
		// convert pdb file to xyz file
		runPdbxyz(pdbFile, xyzFile, log);

		double energy = runMinimize(xyzFile, rmsGradient, log);
		
		// convert final structure to pdb file with extension "min.pdb"
		runXyzpdb(xyzFile, seqFile, new File(outputDir,baseName+".min.pdb"), log);
		
		return energy;
	}
	
	/**
	 * Given a pdb file computes its energy for the currently set forceField
	 * Temporary files .xyz and .seq will be removed on exit.
	 * @param pdbFile
	 * @param log
	 * @return
	 * @throws TinkerError
	 * @throws IOException
	 */
	public double computeEnergy(File pdbFile, PrintWriter log) throws TinkerError, IOException {
		
		boolean cleanUp = true; // for debugging
		
		String baseName = getBasename(pdbFile);
		String outputDir = pdbFile.getParent();
		File xyzFile = new File(outputDir,baseName+".xyz");
		File seqFile = new File(outputDir,baseName+".seq"); // this file is created implicitely by pdbxyz
		if (cleanUp) {
			xyzFile.deleteOnExit();
			seqFile.deleteOnExit();
		}
		// convert pdb file to xyz file
		runPdbxyz(pdbFile, xyzFile, log);

		return runAnalyze(xyzFile, log);
	}	
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Returns a list of violated edges given an input graph to a reconstruction 
	 * and the output reconstructed structure.
	 * @return
	 */
	public static IntPairSet getViolatedEdges(RIGraph graph, Pdb reconstructedPdb) {
		double cutoff = graph.getCutoff();
		String ct = graph.getContactType();
		
		if (!graph.getSequence().equals(reconstructedPdb.getSequence())) {
			throw new IllegalArgumentException("Given graph and reconstructedPdb don't have identical sequence.");
		}
		if (!AAinfo.isValidSingleAtomContactType(ct)) {
			throw new IllegalArgumentException("Given contact type "+ct+" is invalid for getting violations. Only single atom contact types are supported.");
		}
		IntPairSet edgeSet = new IntPairSet();

		for (RIGEdge edge:graph.getEdges()) {
			Pair<RIGNode> nodePair = graph.getEndpoints(edge);
			int i = nodePair.getFirst().getResidueSerial();
			int j = nodePair.getSecond().getResidueSerial();
			Set<String> iatoms = AAinfo.getAtomsForCTAndRes(ct, nodePair.getFirst().getResidueType());
			Set<String> jatoms = AAinfo.getAtomsForCTAndRes(ct, nodePair.getSecond().getResidueType());
			// iatoms and jatoms have only 1 member (we are forcing ct to be single atom) so this is not an iteration
			for (String iatom:iatoms) {  
				for (String jatom:jatoms) {
					Point3d iCoord = reconstructedPdb.getAtomCoord(i, iatom);
					Point3d jCoord = reconstructedPdb.getAtomCoord(j, jatom);
					if (iCoord.distance(jCoord)>cutoff) {
						edgeSet.add(new Pair<Integer>(i,j));
					}
				}
			}
		}
		return edgeSet;
	}
}



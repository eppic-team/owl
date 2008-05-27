package tinker;

import graphAveraging.ConsensusSquare;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.channels.FileChannel;
import java.util.Formatter;
import java.util.TreeMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import proteinstructure.Pdb;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbfilePdb;
import proteinstructure.MaxClusterRunner;
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
	
	public static final double DEFAULT_FORCECONSTANT_DISTANCE_CONST = 10.0;
	public static final double DEFAULT_FORCECONSTANT_TORSION_CONST = 1.0;
	
	private static final String TINKER_ERROR_STR = " TINKER is Unable to Continue";
	private static final String CHECKXYZ_WARNING = " CHKXYZ";
	
	private static final PRMInfo.PRMType DEFAULT_FF_FILE_TYPE = PRMInfo.PRMType.amber;
	
	public static final String DEFAULT_RECONSTR_CHAIN_CODE = Pdb.NULL_CHAIN_CODE;
	
	private static final String ANALYZE_ENERGY_MODE = "E"; // energy mode of analyze program, other valid modes are: A, L, D, M, P (see tinker's docs)
	
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
	 * @param outPath Directory where output files will be written
	 * @param outBasename Base name of the output files
	 * @param n Number of models that we want distgeom to produce
	 * @param log A PrintWriter for logging output
	 * @throws TinkerError If an error seen in tinker's output
	 * @throws IOException
	 */
	private void runDistgeom(File xyzFile, String outPath, String outBasename, int n, PrintWriter log) throws TinkerError, IOException, InterruptedException {
		boolean tinkerError = false; // to store the exit state of the tinker program
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
		// getting names of tinker output files
		File[] tinkerout = new File[n+1];
		for (int i=1;i<=n;i++) {
			String ext = new Formatter().format("%03d", i).toString();
			tinkerout[i] = getTinkerOutputFileName(xyzFile, ext);
		}
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
		// running distgeom program
		String cmdLine = distgeomProg+" "+xyzFile.getAbsolutePath()+" "+n+" "+dgeomParams;
		Process dgeomProc = Runtime.getRuntime().exec(cmdLine);
		// logging and capturing output
		BufferedReader dgeomOutput = new BufferedReader(new InputStreamReader(dgeomProc.getInputStream()));
		String line;
		int i=1;
		log.println("#cmd: "+cmdLine);
		while((line = dgeomOutput.readLine()) != null) {
			log.println(line);
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
		//renaming files to our chosen outBasename+ext
		for (i=1;i<=n;i++) {
			String ext = new Formatter().format("%03d", i).toString();
			tinkerout[i].renameTo(new File(outPath,outBasename+"."+ext));
		}
		// throwing exception if error string was caught in output
		if (tinkerError) {
			log.flush();
			throw new TinkerError("Tinker error, revise log file. ");
		}
		int exitValue = dgeomProc.waitFor();
		// throwing exception if exit state is 137: happens in Linux when another instance of distgeom is running in same machine, the OS kills it with exit state 137 
		if (exitValue==137) {
			log.flush();
			throw new TinkerError("Distgeom was killed by OS. There may be another instance of distgeom running in this computer" +
					" or Tinker could not allocate enough memory.");
		}
		else if (exitValue==139) {
			log.flush();
			throw new TinkerError("Distgeom was killed with exit code 139. Not enough memory.");

		}
		// this is to catch all other possible errors not caught already by the parse of the error string in output
		else if (exitValue!=0) { 
			log.flush();
			throw new TinkerError("Distgeom exited with a non 0 exit code: "+exitValue+". Unknown error.");
		}
		
		log.flush();
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
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 * @returns A pdb object containg the generated structure
	 */
	public Pdb reconstruct(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, int numberOfModels) 
	throws TinkerError, IOException {
		
		return reconstruct(sequence, graphs, phiPsiConsensus, numberOfModels, 
				DEFAULT_FORCECONSTANT_DISTANCE_CONST, DEFAULT_FORCECONSTANT_TORSION_CONST);
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
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 * @returns A pdb object containg the generated structure
	 */
	public Pdb reconstructFast(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, int numberOfModels) 
	throws TinkerError, IOException {
		
		return reconstructFast(sequence, graphs, phiPsiConsensus, 
				numberOfModels, DEFAULT_FORCECONSTANT_DISTANCE_CONST, DEFAULT_FORCECONSTANT_TORSION_CONST);
	}

	/** 
	 * Reconstructs the given graph(s) and returns a putative good model as a Pdb object.
	 * The model is picked based on the number of bound violations reported by Tinker.
	 * See reconstruct() for more details.
	 * @param sequence sequence of the structure to be generated
	 * @param graphs array of graph objects containing constraints
	 * @param phiPsiConsensus Map containing consensus phi/psi angle for each position of the sequence, 
	 * if null no phi/psi restraints will be applied
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param forceConstantDist the force constant to be used for all our given distance restraints
	 * @param forceConstantTorsion the force constant to be used for all our given torsion restraints
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 * @returns a Pdb object containg the generated structure
	 */
	public Pdb reconstruct(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus,
			int numberOfModels, double forceConstantDist, double forceConstantTorsion) 
	throws TinkerError, IOException {
		
		Pdb resultPdb = null;
		
		String outputDir = System.getProperty("java.io.tmpdir");
		String baseName = Long.toString(System.currentTimeMillis());	// some hopefully unique basename
		boolean cleanUp = true;						
		
		reconstruct(sequence, graphs, phiPsiConsensus, numberOfModels, 
				forceConstantDist, forceConstantTorsion, true, outputDir, baseName, cleanUp);
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
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param forceConstantDist the force constant to be used for all our given distance restraints
	 * @param forceConstantTorsion the force constant to be used for all our given torsion restraints
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 * @returns a Pdb object containg the generated structure
	 */
	public Pdb reconstructFast(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, 
			int numberOfModels, double forceConstantDist, double forceConstantTorsion) 
	throws TinkerError, IOException {
		
		Pdb resultPdb = null;
		
		String outputDir = System.getProperty("java.io.tmpdir");
		String baseName = Long.toString(System.currentTimeMillis());	// some hopefully unique basename
		boolean cleanUp = true;						
		
		reconstruct(sequence, graphs, phiPsiConsensus, numberOfModels, forceConstantDist, forceConstantTorsion, false, outputDir, baseName, cleanUp);
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
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param forceConstantDist the force constant to be used for all our given distance restraints
	 * @param forceConsantTorsion the force constant to be used for all our given torsion restraints
	 * @param outputDir the directory where the temporary and result files will be written to
	 * @param baseName the basename of the temporary and result files
	 * @param cleanUp whether to mark all created files to be deleted on shutdown
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 */
	public void reconstruct(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, 
			int numberOfModels, double forceConstantDist, double forceConsantTorsion, 
			String outputDir, String baseName, boolean cleanUp) 
	throws TinkerError, IOException {
		
		reconstruct(sequence, graphs, phiPsiConsensus, numberOfModels, forceConstantDist, forceConsantTorsion, true, outputDir, baseName, cleanUp);
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
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param forceConstantDist the force constant to be used for all our given distance restraints
	 * @param forceConstantTorsion the force constant to be used for all our given torsion restraints
	 * @param outputDir the directory where the temporary and result files will be written to
	 * @param baseName the basename of the temporary and result files
	 * @param cleanUp whether to mark all created files to be deleted on shutdown
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 */	
	public void reconstructFast(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, 
			int numberOfModels, double forceConstantDist, double forceConstantTorsion, 
			String outputDir, String baseName, boolean cleanUp) 
	throws TinkerError, IOException {
		
		reconstruct(sequence, graphs, phiPsiConsensus, 
				numberOfModels, 
				forceConstantDist, forceConstantTorsion, 
				false, outputDir, baseName, cleanUp);
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
	 * @param numberOfModels number of reconstructions to be done by tinker
	 * @param forceConstantDist the force constant to be used for all our given distance restraints
	 * @param forceConstantTorsion the force constant to be used for all our given torsion restraints
	 * @param annealing if true simulated annealing is used for refinement (slow), if false minimization is used for refinement (fast) 
	 * @param outputDir the directory where the temporary and result files will be written to
	 * @param baseName the basename of the temporary and result files
	 * @param cleanUp whether to mark all created files to be deleted on shutdown
	 * @throws TinkerError  if reconstruction fails because of problems with Tinker
	 * @throws IOException  if some temporary or result file could not be accessed
	 */
	private void reconstruct(String sequence, RIGraph[] graphs, TreeMap<Integer, ConsensusSquare> phiPsiConsensus, 
			int numberOfModels, double forceConstantDist, double forceConstantTorsion, 
			boolean annealing, String outputDir, String baseName, boolean cleanUp) 
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
		
		cm.closeKeyFile();

		// 3. run tinker's distgeom
		try {
			runDistgeom(xyzFile, outputDir, baseName, numberOfModels, log);
		} catch(InterruptedException e) {
			throw new TinkerError(e);
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
	
}

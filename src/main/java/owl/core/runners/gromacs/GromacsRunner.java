package owl.core.runners.gromacs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GromacsRunner {

	private static final String PDB2GMX_PROG = "pdb2gmx";
	private static final String EDITCONF_PROG = "editconf";
	private static final String GENBOX_PROG = "genbox";
	private static final String GENION_PROG = "genion";
	private static final String GROMPP_PROG = "grompp";
	private static final String MDRUN_PROG = "mdrun";
	private static final String MDRUNMPI_PROG = "mdrun_mpi";
	private static final String MAKENDX_PROG = "make_ndx";
	
	private static final String MPIRUN_PROG = "mpirun";
	
	private static final String DEFAULT_SOLVENT_MODEL = "spc216.gro";
	private static final String DEFAULT_POS_ION = "NA+";
	private static final String DEFAULT_NEG_ION = "CL-";
	private static final String DEFAULT_BOXTYPE = "dodecahedron";
	private static final double DEFAULT_BOXSIDES = 0.5;
	
	private static final String SOLVENT_GROUP = "SOL";
	
	private static final String GROMACS_ERROR_STR = "Fatal error";
	
	private static final String TMP_DIR = System.getProperty("java.io.tmpdir");
	private static final long currentTimeMS = System.currentTimeMillis(); // for a unique name of temp files
	
	private static final boolean debug = false;
	
	// force fields
	public static final String GROMOS96_43a1 = "G43a1"; //GROMOS96 43a1 (option 0)
	
	// members
	private File gmxBinDir;
	private File pdb2gmxProg;
	private File editconfProg;
	private File genboxProg;
	private File genionProg;
	private File gromppProg;
	private File mdrunProg;
	private File mdrunmpiProg;
	private File make_ndxProg;
	
	private int numProc; // number of processes (for parallel gromacs using MPI)
	
	private File logFile;
	private PrintWriter gmxLog;
	
	/**
	 * Constructs a new GromacsRunner with 1 as number of processes, i.e. will not use parallel gromacs
	 * @param gmxBinDir
	 * @param logFile
	 * @throws FileNotFoundException if given logFile does not denote a writable file path
	 */
	public GromacsRunner(File gmxBinDir, File logFile) throws FileNotFoundException {
		new GromacsRunner(gmxBinDir, logFile, 1);
	}
	
	/**
	 * Constructs a new GromacsRunner
	 * @param gmxBinDir
	 * @param logFile
	 * @param numProc
	 * @throws FileNotFoundException if given logFile does not denote a writable file path
	 */
	public GromacsRunner(File gmxBinDir, File logFile, int numProc) throws FileNotFoundException {
		this.gmxBinDir = gmxBinDir;
		this.numProc = numProc;
		
		this.pdb2gmxProg = new File(this.gmxBinDir, PDB2GMX_PROG);
		this.editconfProg = new File(this.gmxBinDir, EDITCONF_PROG);
		this.genboxProg = new File(this.gmxBinDir, GENBOX_PROG);
		this.genionProg = new File(this.gmxBinDir, GENION_PROG);
		this.gromppProg = new File(this.gmxBinDir, GROMPP_PROG);
		this.mdrunProg = new File(this.gmxBinDir, MDRUN_PROG);
		this.mdrunmpiProg = new File(this.gmxBinDir, MDRUNMPI_PROG);
		this.make_ndxProg = new File(this.gmxBinDir, MAKENDX_PROG);
		this.logFile = logFile;
		this.gmxLog = new PrintWriter(this.logFile);
	}
	
	/**
	 * Runs pdb2gmx to convert a pdb file to a gro file and parses 
	 * the charge value from the output.
	 * @param inPdb
	 * @param outGro
	 * @param outTop
	 * @param forceField
	 * @return the total charge that gromacs report for the given pdb file
	 * @throws GromacsError
	 */
	private int runPdb2gmx(File inPdb, File outGro, File outTop, File outPosreItp, String forceField) throws GromacsError {
		int charge = 0;
		String cmdLine = pdb2gmxProg+" -f "+inPdb+" -o "+outGro+" -p "+outTop+ " -i "+outPosreItp+" -ff "+forceField;
		
		try {
			Process proc = Runtime.getRuntime().exec(cmdLine);
			// logging and capturing output
			BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));
			BufferedReader stderr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

			boolean gromacsError = false;
			String line;

			gmxLog.println("#cmd: "+cmdLine);
			gmxLog.println("#################");
			gmxLog.println("# "+PDB2GMX_PROG+" stdout ");
			gmxLog.println("#################");
			while((line = stdout.readLine()) != null) {
				gmxLog.println(line);
				if (line.startsWith(GROMACS_ERROR_STR)) {
					gromacsError = true;
				}
			}
			gmxLog.println("#################");
			gmxLog.println("# "+PDB2GMX_PROG+" stderr ");
			gmxLog.println("#################");
			while((line = stderr.readLine()) != null) {
				gmxLog.println(line);
				if (line.startsWith(GROMACS_ERROR_STR)) {
					gromacsError = true;
				}				
				Pattern p = Pattern.compile("^Total charge (-?\\d+)\\.\\d+ e");
				Matcher m = p.matcher(line);
				if (m.find()) {
					charge = Integer.parseInt(m.group(1));
				}						
			}


			// throwing exception if error string was caught in output
			if (gromacsError) {
				gmxLog.flush();
				throw new GromacsError("Error caught in output of "+PDB2GMX_PROG+". Revise log file "+logFile);
			}
			try {
				int exitValue = proc.waitFor();
				// throwing exception if exit state is not 0 
				if (exitValue!=0) {
					gmxLog.flush();
					throw new GromacsError(PDB2GMX_PROG + " exited with value "+exitValue+". Revise log file "+logFile);
				}
			} catch (InterruptedException e) {
				System.err.println("Unexpected error while waiting for "+PDB2GMX_PROG+" to exit. Error: "+e.getMessage());
				System.exit(1);
			}

		} catch (IOException e) {
			throw new GromacsError("IO error while trying to run "+PDB2GMX_PROG+": "+e.getMessage());
		}
		
		return charge;
	}
	
	private void createBox(File inGro, File outGro, String boxType, double boxSide) throws GromacsError {
		String cmdLine = editconfProg+" -f "+inGro+" -o "+outGro+" -c -d "+boxSide+" -bt "+boxType;
		runProg(cmdLine, EDITCONF_PROG, null);
	}
	
	/**
	 * Runs genbox to add solvent molecules to the given gro file
	 * @param inGro
	 * @param outGro
	 * @param ioTop both the input and output topology file
	 * @throws GromacsError
	 */
	private void runGenbox(File inGro, File outGro, File ioTop) throws GromacsError {
		//-cp is for solute, -cs for solvent
		String cmdLine = genboxProg+" -cp "+inGro+" -cs "+DEFAULT_SOLVENT_MODEL+" -o "+outGro+" -p "+ioTop;
		runProg(cmdLine, GENBOX_PROG, null);
	}
	
	/**
	 * Runs genion to add the specified number of negative (CL-) and positive (NA+) ions 
	 * writing output to outGro
	 * @param inTpr
	 * @param outGro
	 * @param ioTop both the input and output topology file
	 * @param log
	 * @param negCharge
	 * @param posCharge
	 * @throws GromacsError
	 */
	private void runGenion(File inTpr, File outGro, File ioTop, File log, int negCharge, int posCharge) throws GromacsError {
		if (negCharge == 0 || posCharge==0) {
			throw new GromacsError("Invalid input for runGenion, one of the given charges was 0: positive "+posCharge+", negative "+negCharge);
		}
		String cmdLine = genionProg + " -s "+inTpr+" -o "+outGro+" -p "+ioTop+" -g "+log+
		" -nname "+DEFAULT_NEG_ION+" -nn "+negCharge+" -pname "+DEFAULT_POS_ION+" -np "+posCharge;
		runProg(cmdLine, GENION_PROG, SOLVENT_GROUP);
	}
	
	private void runGrompp(File inMdp, File inGro, File inTop, File outTpr, int numProc) throws GromacsError {
		// by default grompp writes a mdout.mdp which is just a confirmation of what it read from inMdp, pretty useless except for debugging
		File mdoutMdp = new File(TMP_DIR,"mdout_"+currentTimeMS+".mdp");  
		
		String cmdLine = gromppProg+" -f "+inMdp+" -c "+inGro+" -p "+inTop+" -o "+outTpr + " -po " + mdoutMdp + " -np "+numProc;
		runProg(cmdLine, GROMPP_PROG, null);
		
		// deleteOnExit is not enough, because if runnining several mdruns within one java program then mdoutMdp won't be deleted from one run to the next
		if (!debug) mdoutMdp.delete();
	}
	
	private void runMdrun(File inTpr, File outTrr, File outGro, File outEdr, File outXtc, File log) throws GromacsError {
		String mdRunCmd = mdrunProg.toString();
		File mdRunProg = mdrunProg;
		if (numProc>1) {
			mdRunProg = mdrunmpiProg;
			mdRunCmd = MPIRUN_PROG + " -np " + numProc + " " + mdRunProg;
		}
		String cmdLine = mdRunCmd+" -s "+inTpr+" -o "+outTrr+" -c "+outGro+" -e "+outEdr+" -x "+outXtc+" -g "+log + " -np "+numProc;
		runProg(cmdLine, mdRunProg.getName(), null);
	}
	
	private void runProg(String cmdLine, String progName, String stringToPipe) throws GromacsError{
		
		try {
			Process proc = Runtime.getRuntime().exec(cmdLine);
			if (stringToPipe!=null) {
				// piping input
				PrintWriter stdin = new PrintWriter(proc.getOutputStream());
				stdin.println(stringToPipe);
				stdin.close();
			}
			
			// logging and capturing stdout/stderr
			BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));
			BufferedReader stderr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

			boolean gromacsError = false;
			String line;

			gmxLog.println("#cmd: "+cmdLine);
			gmxLog.println("#################");
			gmxLog.println("# "+progName+" stdout ");
			gmxLog.println("#################");
			while((line = stdout.readLine()) != null) {
				gmxLog.println(line);
				if (line.startsWith(GROMACS_ERROR_STR)) {
					gromacsError = true;
				}
			}
			gmxLog.println("#################");
			gmxLog.println("# "+progName+" stderr ");
			gmxLog.println("#################");
			while((line = stderr.readLine()) != null) {
				gmxLog.println(line);
				if (line.startsWith(GROMACS_ERROR_STR)) {
					gromacsError = true;
				}				
			}


			// throwing exception if error string was caught in output
			if (gromacsError) {
				gmxLog.flush();
				throw new GromacsError("Error caught in output of "+progName+". Revise log file "+logFile);
			}
			try {
				int exitValue = proc.waitFor();
				// throwing exception if exit state is not 0 
				if (exitValue!=0) {
					gmxLog.flush();
					throw new GromacsError(progName + " exited with value "+exitValue+". Revise log file "+logFile);
				}
			} catch (InterruptedException e) {
				System.err.println("Unexpected error while waiting for "+progName+" to exit. Error: "+e.getMessage());
				System.exit(1);
			}

		} catch (IOException e) {
			throw new GromacsError("IO error while trying to run "+progName+": "+e.getMessage());
		}

	}
	
	/*--------------------------------- public methods ---------------------------------------------*/
	
	public void closeLog() {
		this.gmxLog.close();
	} 

	/**
	 * Converts a gro file into a pdb file
	 * @param inGro
	 * @param outPdb 
	 * @param whichGroup either "System": everything in gro file will be in pdb, 
	 * "Protein": only protein atoms go in pdb file,
	 * "Protein-H": only non-Hydrogen protein atoms go in pdb file
	 */
	public void convertGro2Pdb(File inGro, File outPdb, String whichGroup) throws GromacsError {
		File indexFile = new File(TMP_DIR,"index_"+currentTimeMS+".ndx");
		if (!debug) indexFile.deleteOnExit();
		String cmdLine = make_ndxProg+" -f "+inGro+" -o "+indexFile;
		runProg(cmdLine, MAKENDX_PROG, "q");
		cmdLine = editconfProg+" -f "+inGro+" -o "+outPdb+" -n "+indexFile;
		runProg(cmdLine, EDITCONF_PROG, whichGroup);
	}
	
	/**
	 * Converts a pdb file into a gro file
	 * @param inPdb
	 * @param outGro
	 * @param outTop
	 * @param forceField
	 * @return the total charge that gromacs report for the given pdb file
	 * @throws GromacsError
	 */
	public int convertPdb2Gro(File inPdb, File outGro, File outTop, File outPosreItp, String forceField) throws GromacsError {
		return runPdb2gmx(inPdb, outGro, outTop, outPosreItp, forceField);
	}
	
	/**
	 * Adds solvent (water) box 
	 * @param inGro
	 * @param outGro
	 * @param ioTop both the input and output topology file
	 * @throws GromacsError
	 */
	public void addSolvent(File inGro, File outGro, File ioTop) throws GromacsError {
		File boxGro = new File(TMP_DIR,"box_"+currentTimeMS+".gro");
		if (!debug) boxGro.deleteOnExit();
		createBox(inGro, boxGro, DEFAULT_BOXTYPE, DEFAULT_BOXSIDES);
		runGenbox(boxGro, outGro, ioTop);
	}
	
	/**
	 * Adds CL- and NA+ ions to compensate the given charge. 
	 * It will always add both type of ions resulting in a total charge that compensates the given one,
	 * e.g. for charge = +2: 3 CL- and 1 NA+ are added.
	 * @param inGro
	 * @param ioTop both the input and output topology file
	 * @param outGro
	 * @param mdp
	 * @param log
	 * @param charge
	 * @throws GromacsError
	 */
	public void addIons(File inGro, File ioTop, File outGro, File mdp, File log, int charge) throws GromacsError{
		File ionTpr = new File(TMP_DIR,"ions_"+currentTimeMS+".tpr");
		if (!debug) ionTpr.deleteOnExit();
		runGrompp(mdp, inGro, ioTop, ionTpr, 1);
		int negCharge = 1;
		int posCharge = 1;
		if (charge>0) {
			negCharge = charge + 1;
			posCharge = 1;
		} else if (charge<0) {
			negCharge = 1;
			posCharge = Math.abs(charge)+1;
		}
		System.out.println("Adding "+negCharge+" CL- and "+posCharge+" NA+");
		runGenion(ionTpr, outGro, ioTop, log, negCharge, posCharge);
	}

	/**
	 * Runs a molecular dynamics simulationg first using grompp (preprocessor) and then mdrun.
	 * The parameters for the simulation are specified through the inMdp file.  
	 * @param inGro
	 * @param inTop
	 * @param inMdp
	 * @param tpr
	 * @param outGro
	 * @param outTrr
	 * @param outEdr
	 * @param outXtc
	 * @param log
	 * @throws GromacsError
	 */
	public void doSimulation(File inGro, File inTop, File inMdp, File tpr, File outGro, File outTrr, File outEdr, File outXtc, File log) 
	throws GromacsError {
		runGrompp(inMdp, inGro, inTop, tpr, this.numProc);
		runMdrun(tpr, outTrr, outGro, outEdr, outXtc, log);		
	}	
	
	/*--------------------------------- static methods ---------------------------------------------*/
	
	/**
	 * Returns "i386" (for our 32 bit machines) or "amd64" (for our 64 bit machines) 
	 * depending on architecture of machine where the code is executed
	 * @return
	 */
	public static String getArchitecture() {
		return System.getProperty("os.arch");
	}
}

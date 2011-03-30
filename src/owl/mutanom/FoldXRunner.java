package owl.mutanom;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.net.UnknownHostException;
import java.util.Collection;
import java.util.LinkedList;


import owl.core.structure.AminoAcid;
import owl.core.structure.PdbChain;
import owl.mutanom.core.Gene;
import owl.mutanom.core.Mutation;
import owl.mutanom.core.Substructure;
import owl.mutanom.core.TargetList;


/**
 * Interface class to run the FoldX algorithms for the calculation of energy changes upon mutation.
 * Provides methods to run FoldX and to retrieve the latest results which are stored in this object.
 * See: Francois Stricher, Tom Lenaerts, Joost Schymkowitz, Frederic Rousseau and Luis Serrano (2008). FoldX 3.0
 * @author stehr
 */
public class FoldXRunner {

	/*------------------------------ constants ------------------------------*/
	public static final String FOLDX_EXE = "FoldX.linux64";
	public static final String ROTABASE = "rotabase.txt";
	public static final String INDIVIDUAL_LIST_FILE = "individual_list.txt";
	
	// for running a local job (files on nfs)
	public static final String FOLDX_PATH = "/project/StruPPi/Software/FoldX/";
	
	// for running a job on the cluster (these need to be copied using cp2cluster.sh)
	public static final String FOLDX_PATH_CL = "/scratch/local/foldx/copy2nodes/";	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Creates a new FoldXRunner
	 */
	public FoldXRunner() {
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Calculates the energy change of the given mutation on the given structure.
	 */
	public double getEnergyChangeUponMutation(PdbChain pdb, Mutation m) {
		double result = 0.0;
		return result;
		
		/*		
  		<TITLE>FOLDX_runscript;
		<JOBSTART>#;
		<PDBS>#;
		<BATCH>list.txt;
		<COMMANDS>FOLDX_commandfile;
		<BuildModel>#,individual_list.txt;
		<END>#;
		<OPTIONS>FOLDX_optionfile;
		<Temperature>298;
		<R>#;
		<pH>7;
		<IonStrength>0.050;
		<water>-CRYSTAL;
		<metal>-CRYSTAL;
		<VdWDesign>2;
		<OutPDB>true;
		<pdb_hydrogens>false;
		<complex_with_DNA> true;
		<END>#;
		<JOBEND>#;
		<ENDFILE>#;
		*/
	}
	
	/**
	 * Returns the free energy of the given structure as caculated by FoldX. Note that
	 * absolute energies should not be trusted. FoldX only claims to give useful results
	 * for energy changes.
	 * @param pdb
	 * @return
	 */
	public double getEnergy(PdbChain pdb) {
		return 0;
		
		/*
		<TITLE>FOLDX_runscript;
		<JOBSTART>#;
		<PDBS>#;
		<BATCH>list.txt;
		<COMMANDS>FOLDX_commandfile;
		<Stability>Stability.txt;
		<END>#;
		<OPTIONS>FOLDX_optionfile;
		<Temperature>298;
		<R>#;
		<pH>7;
		<IonStrength>0.050;
		<water>-CRYSTAL;
		<metal>-CRYSTAL;
		<VdWDesign>2;
		<OutPDB>false;
		<pdb_hydrogens>false;
		<END>#;
		<JOBEND>#;
		<ENDFILE>#;
		*/
	}

	/**
	 * Runs the repair function of FoldX and returns the resulting structure
	 * @param pdb
	 * @return
	 */
	public PdbChain repairStructure(PdbChain pdb) {
		return null;
		
		/*
		<TITLE>FOLDX_runscript;
		<JOBSTART>#;
		<PDBS>#;
		<BATCH>list.txt;
		<COMMANDS>FOLDX_commandfile;
		<RepairPDB>#;
		<END>#;
		<OPTIONS>FOLDX_optionfile;
		<Temperature>298;
		<R>#;
		<pH>7;
		<IonStrength>0.050;
		<water>-CRYSTAL;
		<metal>-CRYSTAL;
		<VdWDesign>2;
		<OutPDB>true;
		<pdb_hydrogens>false;
		<END>#;
		<JOBEND>#;
		<ENDFILE>#; 
		*/
	}
	
	/**
	 * Returns a list of the interface residues in the given complex. This is actually a side-product
	 * of calculating the interaction energy of the complex. For the moment, we ignore this energy.
	 * @return
	 */
	public Collection<Integer> getInterfaceResidues(String pdbCode, String chainCode) {
		LinkedList<Integer> interfaceResidues = null;
		return interfaceResidues;
		
		/*
		<TITLE>FOLDX_runscript;
		<JOBSTART>#;
		<PDBS>#;
		<BATCH>list.txt;
		<COMMANDS>FOLDX_commandfile;
		<AnalyseComplex>#,A;
		<END>#;
		<OPTIONS>FOLDX_optionfile;
		<Temperature>298;
		<R>#;
		<pH>7;
		<IonStrength>0.050;
		<water>-CRYSTAL;
		<metal>-CRYSTAL;
		<VdWDesign>2;
		<OutPDB>false;
		<pdb_hydrogens>false;
		<END>#;
		<JOBEND>#;
		<ENDFILE>#;		 
		*/
	}
	
	/**
	 * Calculates the mutant structure based on the given pdb and set of mutations.
	 * To get useful results, repairStructure() should be run on the input previously.
	 * @param pdb
	 * @return
	 */
	public PdbChain getMutatedStructure(PdbChain pdb, Collection<Mutation> mutations) {
		return null;
	}
	
	/**
	 * Writes a shell script which can be called to run a comprehensive mutational
	 * screening (all positions and all possible target amino acids) on the cluster.
	 * The results will be written to the defined outDir.
	 * Substructures and Pdbs from the target genes need to be loaded.
	 * The input structures need to be called RepairPDB_XXXX.pdb where XXXX is the PDB code.
	 * Requires a script runit.sh in the output dir which looks like this:
	 * #!/bin/sh
	 * CLASSPATH=...
	 * java $@
	 * @throws FileNotFoundException 
	 */
	public static void writeMutScanMasterShellScript(File scriptFile, TargetList targets, File inDir, File outDir, boolean useSingleChainPdbFiles) throws FileNotFoundException {
		
		// directory with pre-calculated pdb files which have been "repaired" using FoldX
		File repairedPdbFilesDir = inDir;
		
		// if true, single chain files from Pdbase are used, otherwise downloaded pdb files with all chains
		//boolean useSingleChainPdbFiles = false;	
		int numJobs = 0;
		PrintWriter out = new PrintWriter(scriptFile);
		out.println("#!/bin/sh");
		out.println("# This script was automatically created by mutanom.FoldXRunner.writeMutScanMasterShellScript()");
		out.println("# Run to submit the specified jobs to the cluster. Output goes to " + outDir);		
		for(Gene g:targets.getTargets()) {
			if(!g.areSubstructuresLoaded()) {
				System.err.println("Skipping gene " + g.getGeneName() + ": substructures not loaded.");
			}
			for(Substructure ss:g.getSubstructures()) {
				if(!ss.isPdbLoaded()) {
					System.err.println("Skipping " + g.getGeneName() + " " + ss.getRange() + " " + ss.getPdbCode()+ss.getChainCode() + ": Structure not loaded.");
				} else {
					for(int p: ss.getPdb().getAllSortedResSerials()) {	// only observed residues
						// assuming cif residue numbers are fine
						String pdbCode = ss.getPdbCode();
						String chain = ss.getChainCode();
						File pdbFile = null;
						if(useSingleChainPdbFiles) {
							pdbFile = new File(repairedPdbFilesDir, "RepairPDB_" + pdbCode + chain + ".pdb");
						} else {
							pdbFile = new File(repairedPdbFilesDir, "RepairPDB_" + pdbCode + ".pdb");
						}
						if(!pdbFile.canRead()) {
							System.err.println("Error, repaired pdb file " + pdbFile + " not found.");
							continue;
						}
						String pdbFileName = pdbFile.getAbsolutePath();
						String aaBefore = Character.toString(ss.getPdb().getResidue(p).getAaType().getOneLetterCode());
						String position = "?";
						if(useSingleChainPdbFiles) {
							position = Integer.toString(p);
						} else {
							position = ss.getPdb().getResidue(p).getPdbSerial();	// cif -> pdb serial
						}
						String outDirPath = outDir.getAbsolutePath();
						String jobName = getJobName(pdbCode, chain, position);
						out.printf("qsub -V -q all.q -e %s -o %s -N %s runit.sh mutanom.FoldXRunner %s %s %s %s %s %s\n", 
								outDirPath, outDirPath, jobName, pdbFileName, position, pdbCode, chain, aaBefore, outDirPath);
						numJobs++;
					}					
				}
			}
		}
		out.println("echo \""+ numJobs + " jobs submitted on `date`\"");		
		out.close();
	}
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Creates a unique temporary directory for putting FoldX output files or null if something went wrong.
	 */
	private static File createTempDir() {
		File tempDir;
		try {
			tempDir = File.createTempFile("FoldX", "");
			tempDir.delete();
			if(!tempDir.mkdir()) return null;
			tempDir.deleteOnExit();
		} catch (IOException e) {
			return null;
		}
		return tempDir;
	}
	
	private static void writeOptionFile(File outFile, boolean outPdb) throws FileNotFoundException {
		PrintWriter out = new PrintWriter(outFile);
		out.println("<TITLE>FOLDX_optionfile;");
		out.println("<Temperature>298;");
		out.println("<R>#;");
		out.println("<pH>7;");
		out.println("<IonStrength>0.050;");
		out.println("<water>-CRYSTAL;");
		out.println("<metal>-CRYSTAL;");
		out.println("<VdWDesign>2;");
		out.printf("<OutPDB>%b;\n", outPdb);
		out.println("<pdb_hydrogens>false;");
		out.close();
	}
	
//	private static void writeBuildModelCommandFile(File outFile, File mutFile) throws FileNotFoundException {
//		PrintWriter out = new PrintWriter(outFile);
//		out.println("<TITLE>FOLDX_commandfile;");
//		out.printf("<BuildModel>%s,%s;\n", outFile, mutFile);
//		out.close();
//	}
	
	private static void writeBuildModelMutScanCommandFiles(File cmdFile, File mutFile, String wtAa, String chain, String pos) throws FileNotFoundException {
		// write command file
		PrintWriter out = new PrintWriter(cmdFile);
		out.println("<TITLE>FOLDX_commandfile;");
		out.printf("<BuildModel>%s,%s;\n", "#", mutFile.getName());
		out.close();
		
		// write individual_list file
		PrintWriter out2 = new PrintWriter(mutFile);
		for(AminoAcid mutAa:AminoAcid.values()) {
			if(mutAa.isStandardAA()) {
				String line = wtAa + chain + pos + mutAa.getOneLetterCode() + ";";
				out2.println(line);
			}
		}
		out2.close();
	}
	
//	private static void writePositionScanCommandFile(File cmdFile, File resultFile, String aa, String chain, String pos) throws FileNotFoundException {
//		PrintWriter out = new PrintWriter(cmdFile);
//		out.println("<TITLE>FOLDX_commandfile;");
//		out.printf("<PositionScan>%s,%s%s%sa;\n", resultFile, aa, chain, pos);
//		out.close();
//	}
	
//	private static void writeMutationsFile(File outFile, String chain, Collection<Mutation> mutations) throws FileNotFoundException {
//		PrintWriter out = new PrintWriter(outFile);
//		for(Mutation m:mutations) {
//			out.printf("%s%s%d%s\n", m.before.getOneLetterCode(),chain,m.position,m.after.getOneLetterCode());
//		}
//		out.close();		
//	}
	
//	private static void runFoldX(File pdbFile, File commandFile, File optionFile) {
//		
//	}
//	
//	private double parseEnergy(File outFile) {
//		return 0;
//	}
//	
//	private double parseEnergyChange(File outFile) {
//		return 0;
//	}
//
//	private Collection<Integer> parseInterfaceResidues(File outFile) {
//		return null;
//	}
	
	/**
	 * Returns the job name which is being used to send FoldX jobs to the cluster. This is used
	 * to match the output files to the submitted jobs.
	 */
	public static String getJobName(String pdbCode, String chain, String position) {
		return String.format("F%s%04d", pdbCode+chain, Integer.parseInt(position));
	}
	
	/*--------------------------------- main --------------------------------*/
	
	/**
	 * Executable for submitting FoldX jobs to the cluster. Results will be
	 * written to stdout.
	 */
	public static void main(String[] args) {
		boolean cleanUp = true;	// if false, temporary files will be kept for debugging
		boolean debug = false; 		// if true, debug output will be written to stdout
		String cmdLogFileName = "cmd.log"; // all executed command will be written to this log
		String foldXLogFileName = "foldx.log"; // all outout of FoldX to stderr or stdout goes to this log
		
		// parse args
		if(args.length < 5) {
			System.err.println("Usage: FoldXRunner pdbFileName position pdbCode chain aaBefore resultDir");
			System.err.println("Example: FoldXRunner RepairPDB_1aie.pdb 1 1aie A E /runs/out");
			System.err.println("Not enough arguments!");
			System.exit(1);
		}
		
		String resultDirName = args[5]; // "/project/StruPPi/henning/projects/mutanom/analysis/foldx_runs/out"
		String beforeAa = args[4]; //"E";
		String chain = args[3]; //"A";
		String pdbCode = args[2]; // "1aie"
		String pos = args[1]; //"326";
		String pdbFileName = args[0]; //"/project/PyMol/pdbs_download/repairPDB/repair/RepairPDB_1aie.pdb";
		
		File resultDir = new File(resultDirName);
		if(!resultDir.isDirectory() || !resultDir.canWrite()) {
			System.err.println("Can not write to results directory: " + resultDir.getAbsolutePath());
			System.exit(1);
		}
		
		File tempDir = createTempDir();
		if(tempDir == null || !tempDir.isDirectory()) {
			System.err.println("Could not create temporary directory: " + tempDir);
			System.exit(1);
		}
		if(debug) System.out.println("Created temp dir " + tempDir);
		
		// Files we will create
		File localPdbFile = new File(tempDir, "mut.pdb"); if(cleanUp) localPdbFile.deleteOnExit();
		File optionFile = new File(tempDir, "options.txt"); if(cleanUp) optionFile.deleteOnExit();
		File commandFile = new File(tempDir, "commands.txt"); if(cleanUp) commandFile.deleteOnExit();
		File mutFile = new File(tempDir, INDIVIDUAL_LIST_FILE); if(cleanUp) mutFile.deleteOnExit();
		File foldXExe = new File(tempDir, FOLDX_EXE); if(cleanUp) foldXExe.deleteOnExit();
		File rotaBase = new File(tempDir, ROTABASE); if(cleanUp) rotaBase.deleteOnExit();
		File cmdLogFile = new File(tempDir, cmdLogFileName); if(cleanUp) cmdLogFile.deleteOnExit();
		File foldXLogFile = new File(tempDir, foldXLogFileName); if(cleanUp) foldXLogFile.deleteOnExit();
		
		String baseName = "mut";
		
		// Files FoldX will create
		File timeFile = new File(tempDir, "Time.txt"); if(cleanUp) timeFile.deleteOnExit();
		File runLogFile = new File(tempDir, "runlog.txt"); if(cleanUp) runLogFile.deleteOnExit();
		File missingFile = new File(tempDir, "missing.txt"); if(cleanUp) missingFile.deleteOnExit();
		File avgFile = new File(tempDir, "Average_BuildModel_" + baseName + ".fxout"); if(cleanUp) avgFile.deleteOnExit();
		File bmFile = new File(tempDir, "BuildModel_" + baseName + ".fxout"); if(cleanUp) bmFile.deleteOnExit();
		File dbmFile = new File(tempDir, "Dif_BuildModel_" + baseName + ".fxout"); if(cleanUp) dbmFile.deleteOnExit();
		File pdblFile = new File(tempDir, "PdbList_BuildModel_" + baseName + ".fxout"); if(cleanUp) pdblFile.deleteOnExit();
		File rawbmFile = new File(tempDir, "Raw_BuildModel_" + baseName + ".fxout"); if(cleanUp) rawbmFile.deleteOnExit();
				
		try {
			writeOptionFile(optionFile, false);
		} catch (FileNotFoundException e) {
			System.err.println("Error writing file: " + optionFile);
			System.exit(1);
		}
		if(debug) System.out.println("Wrote file " + optionFile);
		try {
			writeBuildModelMutScanCommandFiles(commandFile, mutFile, beforeAa, chain, pos);
			//writePositionScanCommandFile(commandFile, resultFile, beforeAa, chain, pos);
		} catch (FileNotFoundException e) {
			System.err.println("Error writing file: " + commandFile);
			System.exit(1);
		}
		if(debug) System.out.println("Wrote file " + commandFile);
		if(debug) System.out.println("Wrote file " + mutFile);

		// copy pdb file to temp dir
		File pdbFile = new File(pdbFileName);
		if(debug) System.out.println("Copying " + pdbFile + " to " + localPdbFile);
		if(!pdbFile.canRead()) {
			System.err.println("Could not find input pdb file " + pdbFileName);
			System.exit(1);
		}
		try {
			BufferedReader in = new BufferedReader(new FileReader(pdbFile));
			PrintWriter out = new PrintWriter(localPdbFile);
			String line;
			while((line = in.readLine()) != null) {
				out.println(line);
			}
			in.close();
			out.close();
		} catch (IOException e1) {
			System.err.println("Error copying file " + pdbFile + " to " + localPdbFile);
		}
		
		// copy executable and rotabase to temp dir (using unix shell)
		if(!new File(FOLDX_PATH_CL + FOLDX_EXE).canRead()) {
			String hostname;
			try {
				hostname=java.net.InetAddress.getLocalHost().getHostName();
			} catch (UnknownHostException e) {
				hostname="unknown";
			}
			System.err.println("Could not find local copy of FoldX on " + hostname + ": " + FOLDX_PATH_CL + FOLDX_EXE);
			System.exit(1);
		}
		if(!new File(FOLDX_PATH_CL + ROTABASE).canRead()) {
			String hostname;
			try {
				hostname=java.net.InetAddress.getLocalHost().getHostName();
			} catch (UnknownHostException e) {
				hostname="unknown";
			}
			System.err.println("Could not find local copy of rotabase on " + hostname + ": " + FOLDX_PATH_CL + ROTABASE);
			System.exit(1);
		}
		
		// This new way to copy files failed because permissions were not handled correctly
//		File rotaSource = new File(FOLDX_PATH_CL, ROTABASE);
//		try {
//			copyFile(rotaSource, rotaBase);
//		} catch(IOException e) {
//			System.err.println("Error copying " + rotaSource.getPath() + " to " + rotaBase.getPath());
//		}
//		File foldxSource = new File(FOLDX_PATH_CL, FOLDX_EXE);
//		try {
//			copyFile(foldxSource, foldXExe);
//		} catch(IOException e) {
//			System.err.println("Error copying " + foldxSource.getPath() + " to " + foldXExe.getPath());
//		}				
		
		String cmd;
		cmd = "cp " + FOLDX_PATH_CL + ROTABASE + " " + rotaBase.getAbsolutePath();
		if(debug) System.out.println("Executing: " + cmd);
		try {
			Runtime.getRuntime().exec(cmd);
		} catch (IOException e) {
			System.err.println("Error executing command line '" + cmd + "':" + e.getMessage());
			System.exit(1);
		}
		cmd = "cp " + FOLDX_PATH_CL + FOLDX_EXE + " " + foldXExe.getAbsolutePath();
		if(debug) System.out.println("Executing: " + cmd);
		try {
			Runtime.getRuntime().exec(cmd);
		} catch (IOException e) {
			System.err.println("Error executing command line '" + cmd + "':" + e.getMessage());
			System.exit(1);
		}
		try {
			Thread.sleep(1000);
		} catch (InterruptedException e) {
			System.err.println("An error occured while waiting for shell: " + e.getMessage());
			System.exit(1);
		} // give some time for the system to perform the copying		
		
		if(!foldXExe.canRead()) {
			System.err.println("Could not find FoldX file: " + foldXExe);
			System.exit(1);
		}
		if(!rotaBase.canRead()) {
			System.err.println("Could not find rotabase file: " + foldXExe);
			System.exit(1);
		}
		
		// run foldx
		PrintStream cmdLog = System.out;
		try {
			cmdLog = new PrintStream(cmdLogFile);
		} catch(IOException e) {
			System.err.println("Error writing to command log " + cmdLogFile + ": " + e.getMessage());
		}
		cmd = String.format("%s -manual %s %s %s", foldXExe.getName(), localPdbFile.getName(), optionFile.getName(), commandFile.getName());
		cmdLog.println(cmd);
		cmdLog.close();
		if(debug) System.out.println("Executing: " + cmd);
		
		try {
			ProcessBuilder pb = new ProcessBuilder(foldXExe.getName(), "-manual", localPdbFile.getName(), optionFile.getName(), commandFile.getName());
			pb.directory(tempDir);
			pb.redirectErrorStream(true);
			Process p = pb.start();
			InputStream processOut = p.getInputStream();
			PrintStream foldXLog = new PrintStream(foldXLogFile);
			int data;
			while((data = processOut.read()) != -1) {
				foldXLog.write(data);
			}
			foldXLog.close();
			processOut.close();
			p.waitFor();
			
			// copying result file to result directory (through NFS)
			File newResultFile = new File(resultDir, getJobName(pdbCode, chain, pos) + ".txt");
			if(debug) System.out.println("Copying " + dbmFile + " (w/o header) to " + newResultFile);
			if(!dbmFile.canRead()) {
				System.err.println("Result file " + dbmFile + " not found.");
				System.exit(1);
			}
			BufferedReader in = new BufferedReader(new FileReader(dbmFile));
			PrintWriter out = new PrintWriter(newResultFile);
			String line;
			while((line = in.readLine()) != null) {
				if(line.startsWith(baseName)) out.println(line);
			}
			in.close();
			out.close();			
		} catch (IOException e) {
			System.err.println("Error executing command line '" + cmd + "':" + e.getMessage());
			System.exit(1);
		} catch (InterruptedException e) {
			System.err.println("An error occured while waiting for FoldX: " + e.getMessage());
			System.exit(1);
		}
		if(debug) System.out.println("done.");
	}
}

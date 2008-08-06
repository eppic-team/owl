import gnu.getopt.Getopt;

import java.io.File;
import java.io.FileNotFoundException;

import tinker.GromacsError;
import tinker.GromacsMDP;
import tinker.GromacsRunner;


public class runMD {

	private static final String PROG_NAME = "runMD";
	
	private static final String GMX_ROOT = "/project/StruPPi/Software/gromacs-3.3.3/";
	private static final String ARCH = GromacsRunner.getArchitecture();
	private static final File GMX_BIN_DIR_32 = new File(GMX_ROOT,"i686/bin");
	private static final File GMX_BIN_DIR_64 = new File(GMX_ROOT,"x86_64/bin");
	
	private static final int DEFAULT_REFINE_TIME = 300;
	private static final int DEFAULT_EQ_TIME = 200;
	private static final String FORCE_FIELD = GromacsRunner.GROMOS96_43a1;;
	
	public static void main(String[] args) {

		String help =
			"Runs a molecular dynamics simulation of a protein structure in explicit solvent using gromacs.\n" +
			"Usage: \n" +
			PROG_NAME+"\n" +
			"   -i <file>   : input PDB file \n" +
			"  [-o <dir>]   : output directory, default: current dir \n"+
			"  [-p <int>]   : number of processes (uses parallel gromacs through MPI)\n" +
			"  [-t <int>]   : md simulation time in picoseconds, default: "+DEFAULT_REFINE_TIME+"\n" +
			"  [-q <int>]   : equilibration time in picoseconds, default: "+DEFAULT_EQ_TIME+"\n" +
			"  [-a]         : the md will be done using a simulated annealing protocol, default: no annealing\n\n"; 
	

		File inputPdb = null;
		File outDir = new File(".");
		int numProc = 1;
		int simulationTime = DEFAULT_REFINE_TIME;
		int equilibrationTime = DEFAULT_EQ_TIME;
		boolean annealing = false;

		Getopt g = new Getopt(PROG_NAME, args, "i:o:p:t:q:ah?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'i':
				inputPdb = new File(g.getOptarg());
				break;
			case 'o':
				outDir = new File(g.getOptarg());
				break;
			case 'p':
				numProc = Integer.parseInt(g.getOptarg());
				break;
			case 't':
				simulationTime = Integer.parseInt(g.getOptarg());
				break;
			case 'q':
				equilibrationTime = Integer.parseInt(g.getOptarg());
				break;
			case 'a':
				annealing = true;
				break;								
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		// input checks
		if (inputPdb == null) {
			System.err.println("An input PDB file (-i) is required.");
			System.err.println(help);
			System.exit(1);
		}
		
		if (!inputPdb.getName().endsWith(".pdb")) {
			System.err.println("Input PDB file must have pdb extension");
			System.exit(1);
		}

		
		// gromacs bin dir
		File gmxBinDir = GMX_BIN_DIR_32;
		if (ARCH.equals("amd64")) {
			gmxBinDir = GMX_BIN_DIR_64;
		} 
		
		// basename of input file
		String basename = inputPdb.getName().substring(0, inputPdb.getName().indexOf(".pdb"));
		
		// log file
		File logFile = new File(outDir,"gromacs.log");
		
		// mdp files
		File emMdp = new File(outDir,"em.mdp");
		File prMdp = new File(outDir,"pr.mdp");
		File mdMdp = new File(outDir,"md.mdp");
		GromacsMDP gmdp = new GromacsMDP();
		try {
			gmdp.setEMValues();
			gmdp.writeToFile(emMdp);
			gmdp.setPRValues(equilibrationTime);
			gmdp.writeToFile(prMdp);
			if (annealing) {
				gmdp.setAnnealValues(simulationTime);
			} else {
				gmdp.setMDValues(simulationTime);
			}
			gmdp.writeToFile(mdMdp);
		} catch (FileNotFoundException e) {
			System.err.println("Couldn't write mdp file, error: "+e.getMessage()+"\nExiting");
			System.exit(1);
		}
		// defining file names for all other files
		File iniGro = new File(outDir,basename+".gro");
		File iniTop = new File(outDir,basename+".top");
		File posreItp = new File(outDir,basename+".posre.itp");
		File watGro = new File(outDir,basename+".wat.gro");
		File ionGro = new File(outDir,basename+".ion.gro");
		File ionLog = new File(outDir,basename+".ion.log");
		File emGro  = new File(outDir,basename+".em.gro");
		File emLog  = new File(outDir,basename+".em.log");
		File emTpr  = new File(outDir,basename+".em.tpr");
		File emTrr  = new File(outDir,basename+".em.trr");
		File emEdr  = new File(outDir,basename+".em.edr");
		File emXtc  = new File(outDir,basename+".em.xtc");
		File prGro  = new File(outDir,basename+".pr.gro");
		File prLog  = new File(outDir,basename+".pr.log");
		File prTpr  = new File(outDir,basename+".pr.tpr");
		File prTrr  = new File(outDir,basename+".pr.trr");
		File prEdr  = new File(outDir,basename+".pr.edr");
		File prXtc  = new File(outDir,basename+".pr.xtc");
		File mdGro  = new File(outDir,basename+".md.gro");
		File mdLog  = new File(outDir,basename+".md.log");
		File mdTpr  = new File(outDir,basename+".md.tpr");
		File mdTrr  = new File(outDir,basename+".md.trr");
		File mdEdr  = new File(outDir,basename+".md.edr");
		File mdXtc  = new File(outDir,basename+".md.xtc");
		
		// running gromacs
		try {
			GromacsRunner gr = new GromacsRunner(gmxBinDir, logFile, numProc);

			// 1. convert pdb 2 gro file
			System.out.println("Converting "+inputPdb+" to gro "+iniGro);
			int charge = gr.convertPdb2Gro(inputPdb, iniGro, iniTop, posreItp, FORCE_FIELD);
			// 2. add solvent
			System.out.println("Adding solvent");
			gr.addSolvent(iniGro, watGro, iniTop);
			// 3. add ions
			System.out.println("Total charge is "+charge+", adding ions");
			gr.addIons(watGro, iniTop, ionGro, emMdp, ionLog, charge);
			// 4. em
			System.out.println("Running Energy Minimization");
			gr.doSimulation(ionGro, iniTop, emMdp, emTpr, emGro, emTrr, emEdr, emXtc, emLog);
			// 5. pr
			System.out.println("Running Position Restrained equilibration");
			gr.doSimulation(emGro, iniTop, prMdp, prTpr, prGro, prTrr, prEdr, prXtc, prLog);
			// 6. md (with annealing)
			System.out.println("Running Molecular Dynamics simulation");
			gr.doSimulation(prGro, iniTop, mdMdp, mdTpr, mdGro, mdTrr, mdEdr, mdXtc, mdLog);
		}
		catch (FileNotFoundException e) {
			System.err.println("Couldn't write to log file "+logFile+", error: "+e.getMessage()+"\nExiting");
			System.exit(1);
		}
		catch (GromacsError e) {
			System.err.println("Gromacs error: "+e.getMessage()+"\nExiting");
			System.exit(1);
		}
	}

}

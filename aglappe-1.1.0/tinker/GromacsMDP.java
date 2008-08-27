package tinker;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Locale;

/**
 * Class representing a set of gromacs simulation parameters, i.e. an mdp file
 * See http://www.gromacs.org/documentation/reference/online/mdp_opt.html
 */
public class GromacsMDP {
	
	/*----------------------------- constants ---------------------------------- */

	// our reference temperature (not a gromacs mdp parameter, but used by other parameters below)
	private static final int REF_TEMPERATURE	= 300;
	
	// our defaults (only the ones differing from gromacs' defaults)
	private static final String TITLE			= "automatically generated MDP file";
	private static final double DT              = 0.002; // 2fs
	private static final int    NSTXOUT         = 0;
	private static final int    NSTVOUT         = 0;
	private static final int    NSTLOG          = 500;
	private static final int    NSTENERGY       = 500;
	private static final int    NSTLIST         = 5;
	private static final String	COULOMBTYPE     = "PME";
	private static final String OPTIMIZE_FFT    = "yes";
	private static final String TC_GRPS         = "protein non-protein";
	private static final String TAU_T           = "0.1 0.1"; 
	private static final String REF_T           = REF_TEMPERATURE+" "+REF_TEMPERATURE; 
	private static final double COMPRESSIBILITY = 4.5e-5;
	private static final String REF_P           = "1.0";
	private static final String GEN_VEL         = "yes";

	// our EM defaults
	private static final String EM_INTEGRATOR 	= "steep";
	private static final int    EM_NSTEPS 		= 5000; 
	
	// our PR defaults
	private static final String PR_INTEGRATOR 	= "md";
	//private static final int 	PR_NSTEPS 		= 100000; // 200ps (at 2fs steps, see default dt)
	private static final int	PR_NSTXTCOUT 	= 500; // 1ps
	private static final String PR_CONSTRAINTS 	= "all-bonds";
	private static final String PR_TCOUPL       = "berendsen";
	private static final String PR_PCOUPL       = "berendsen";
	
	// our MD defaults
	private static final String MD_INTEGRATOR 	= "md";
	//private static final int 	MD_NSTEPS 		= 150000; // 300ps (at 2fs steps)
	private static final int	MD_NSTXTCOUT 	= 500; // 1ps
	private static final String MD_CONSTRAINTS 	= "all-bonds";
	private static final String MD_TCOUPL       = "berendsen";
	
	// our MD annealing defaults (use MD defaults plus these)
	private static final String MD_ANNEALING = "single single";
	private static final String MD_ANNEALING_NPOINTS = "6 6";
	
	
	/*------------------------------- members --------------------------------- */

	// defaults in brackets (), units in square brackets []
	// * means we want to vary the values
	// + means we use a value different from gromacs default but fixed, our default in ()
	// sections we don't use are not filled up
	
	// Preprocessing
	private String title;			// whatever
//	private String cpp;				// (cpp)
//	private String include;			// 
	private String define; 			// (-DFLEXIBLE)				*
	
	// Run control
	private String integrator;		// (md)						*
//	private double tinit;			// (0) [ps]
	private double dt;				// (0.001) [ps] 			+ (0.002)
	private int nsteps; 			// (0) 						*
//	private double init_step;		// (0)
//	private String comm_mode;		// (Linear)
//	private int nstcomm;			// (1) [steps]
//	private String comm_grps;		// ?
	
	// Langevin dynamics
	
	// Energy minimization
//	private double emtol;			// (10.0) [kJ mol-1 nm-1]
//	private double emstep;			// (0.01) [nm]
//	private int nstcgsteep;			// (1000) [steps]
//	private int nbfgscorr;			// (10)
	
	// Shell Molecular Dynamics
	
	// Output control
	private int nstxout;			// (100) [steps]			+ (0)
	private int nstvout;			// (100) [steps]			+ (0)
//	private int nstfout; 			// (0)   [steps]
	private int nstlog;				// (100) [steps]			+ (500)
	private int nstenergy;			// (100) [steps]			+ (500)
	private int nstxtcout; 			// (0) [steps]				* (500)
//	private double xtc_precision;	// (1000) [real] 
//	private String xtc_grps;		// ?
//	private String energygrps;		// ?
	
	// Neighbor searching
	private int nstlist;			// (10) [steps]				+ (5)
//	private String ns_type;			// (grid)
//	private String pbc;				// (xyz)
//	private double rlist;			// (1) [nm]

	// Electrostatics
	private String coulombtype;		// (Cut-off)				+ (PME)
//	private double rcoulomb_switch; // (0) [nm]
//	private double rcoulomb;		// (1) [nm]
//	private double epsilon_r;		// (1)
//	private double epsilon_rf;		// (1) 

	// VdW
//	private String vdwtype;			// (Cut-off)
//	private double rvdw_switch;		// (0) [nm]
//	private double rvdw;			// (1) [nm]
//	private String DispCorr;		// 
	
	// Tables
	
	// Ewald
//	private double fourierspacing;	// (0.12) [nm]
//	private int fourier_nx;			// (0)
//	private int fourier_ny;			// (0)
//	private int fourier_nz;			// (0)
//	private int pme_order;			// (4)
//	private double ewald_rtol;		// (1e-5)
//	private String ewald_geometry;	// (3d)
//	private double epsilon_surface;	// (0)
	private String optimize_fft;	// (no)						+ (yes)
	
	// Temperature coupling
	private String tcoupl;			// (no)						* (berendsen)
	private String tc_grps;			// ()						+ (protein non-protein)
	private String tau_t; 			// () [ps]					+ (0.1 0.1) 
	private String ref_t;			// () [K]					+ (300 300) 
	
	// Pressure coupling
	private String pcoupl;			// (no)						* (berendsen)
//	private String pcoupltype;		// (isotropic)
//	private double 	tau_p;			// (1) [ps]
	private double compressibility; // () [bar^-1]				+ (4.5e-5)
	private String ref_p;			// () [bar]					+ (1.0)

	// Simulated annealing
	private String annealing;		// (no)						* (single single)
	private String annealing_npoints;//()						* (6 6) 
	private String annealing_time;	// ()						* (0 60 120 180 240 300 0 60 120 180 240 300) 
	private String annealing_temp;	// ()						* (300 240 180 120 60 0 300 300 300 300 300 300) 

	// Velocity generation
	private String gen_vel;			// (no)						+ (yes) only meaningful with integrator md
//	private double gen_temp;		// (300) [K]
//	private int gen_seed;			// (173529)

	// Bonds
	private String constraints; 	// (none)					* (all-bonds)
//	private String constraint_algorithm;// (lincs)
//	private String unconstrained_start;	// (no)
//	private double shake_tol;		// (0.0001)
//	private int lincs_order;		// (4)
//	private int lincs_iter;			// (1)
//	private double lincs_warnangle;	// (30) [degrees]
//	private String morse;			// (no)
	
	// Energy group exclusions
	
	// NMR refinement
	
	// Free Energy Perturbation
	
	// Non-equilibrium MD
	
	// Electric fields
	
	// Mixed quantum/classical molecular dynamics
	
	// User defined thingies


	/*------------------------------- constructors --------------------------------- */
	
	public GromacsMDP() {
		setDefaults();
	}
	
	public GromacsMDP(String define, 
			String integrator, 
			int nsteps, int nstxtcout, 
			String annealing, String annealing_npoints, String annealing_time, String annealing_temp, 
			String constraints) {
		
		setDefaults();
		this.define = define;
		this.integrator = integrator;
		this.nsteps = nsteps;
		this.nstxtcout = nstxtcout;
		this.annealing = annealing;
		this.annealing_npoints = annealing_npoints;
		this.annealing_time = annealing_time;
		this.annealing_temp = annealing_temp;
		this.constraints = constraints;
	}
	
	/*------------------------------- private methods --------------------------------- */
	
	private void setDefaults() {
		this.title = TITLE;
		this.dt = DT;
		this.nstxout = NSTXOUT;
		this.nstvout = NSTVOUT;
		this.nstlog = NSTLOG;
		this.nstenergy = NSTENERGY;
		this.nstlist = NSTLIST;
		this.coulombtype = COULOMBTYPE;
		this.optimize_fft = OPTIMIZE_FFT;
		this.tc_grps = TC_GRPS;
		this.tau_t = TAU_T;
		this.ref_t = REF_T;
		this.compressibility = COMPRESSIBILITY;
		this.ref_p = REF_P;
		this.gen_vel = GEN_VEL;		
	}
	
	private void resetValues() {
		this.define = null;
		this.integrator = null;
		this.nsteps = 0;
		this.nstxtcout = 0;
		this.tcoupl = null;
		this.pcoupl = null;
		this.annealing = null;
		this.annealing_npoints = null;
		this.annealing_time = null;
		this.annealing_temp = null;
		this.constraints = null;
	}
	
	/**
	 * 
	 * @param npoints number of annealing steps in simulation
	 * @param totalTime in picoseconds
	 * @return
	 */
	private String getAnnealingTime(int npoints, double totalTime) {
		String annTime = "";
		double step = totalTime/(npoints-1);
		for (double time=0;time<=totalTime;time+=step) {
			annTime += String.format("%.1f", time);
			// we add a space except at the last step
			if (time!=totalTime)  
				annTime+=" ";
		}
		return annTime;
	}
	
	/**
	 * 
	 * @param npoints number of annealing steps in simulation
	 * @param startTemp starting temperature
	 * @return
	 */
	private String getAnnealingTemp(int npoints, int startTemp) {
		String annTemp = "";
		int step = startTemp/(npoints-1);
		for (int temp=startTemp;temp>=0;temp-=step) {
			annTemp += temp;
			// we add a space except at the last step
			if (temp!=0) 
				annTemp+=" ";
		}
		return annTemp;
	}
	
	/*-------------------------------- public methods -------------------------------- */
	
	/**
	 * Sets values for an Energy Minimization simulation
	 * from our set of default values.
	 */
	public void setEMValues() {
		this.resetValues(); // we wipe out existing values just in case
		this.integrator = EM_INTEGRATOR;
		this.nsteps = EM_NSTEPS;
	}
	
	/**
	 * Sets values for a Position Restrained equilibration simulation
	 * from our set of default values.
	 */
	public void setPRValues(int equilibrationTime) {
		this.resetValues(); // we wipe out existing values just in case
		this.tcoupl = PR_TCOUPL;
		this.pcoupl = PR_PCOUPL;
		this.integrator = PR_INTEGRATOR;
		this.nsteps = (int)((double) equilibrationTime/dt);
		this.nstxtcout = PR_NSTXTCOUT;
		this.constraints = PR_CONSTRAINTS;
	}
	
	/**
	 * Sets values for a Molecular Dynamics simulation
	 * from our set of default values.
	 * @param simulationTime total time of the simulation in picoseconds
	 */
	public void setMDValues(int simulationTime) {
		this.resetValues(); // we wipe out existing values just in case
		this.tcoupl = MD_TCOUPL;
		this.integrator = MD_INTEGRATOR;
		this.nsteps = (int)((double) simulationTime/dt);
		this.nstxtcout = MD_NSTXTCOUT;
		this.constraints = MD_CONSTRAINTS;
	}
	
	/**
	 * Sets values for a Molecular Dynamics with annealing simulation
	 * from our set of default values.
	 * @param simulationTime total time of the simulation in picoseconds
	 */
	public void setAnnealValues(int simulationTime) {
		this.resetValues(); // we wipe out existing values just in case
		this.setMDValues(simulationTime);
		this.annealing = MD_ANNEALING;
		this.annealing_npoints = MD_ANNEALING_NPOINTS;
		int npoints1 = Integer.parseInt(annealing_npoints.split(" ")[0]);
		int npoints2 = Integer.parseInt(annealing_npoints.split(" ")[1]);
		this.annealing_time = getAnnealingTime(npoints1, simulationTime)+" "+getAnnealingTime(npoints2, simulationTime);
		this.annealing_temp = getAnnealingTemp(npoints1, REF_TEMPERATURE)+" "+getAnnealingTemp(npoints2, REF_TEMPERATURE);
	}
	
	/**
	 * Writes MDP parameters to given file in gromacs MDP format
	 * @param file
	 * @throws FileNotFoundException
	 */
	public void writeToFile(File file) throws FileNotFoundException {
		PrintWriter out = new PrintWriter(file);
		
		// we only write the values we vary (the ones not in setDefaults) when they have been set (not null)
		out.println("title = "+this.title);
		if (this.define!=null) 
			out.println("define = "+this.define);
		if (this.integrator!=null) 
			out.println("integrator = "+this.integrator);
		out.printf(Locale.US,"dt = %5.3f\n",this.dt);
		if (this.nsteps!=0) 
			out.println("nsteps = "+this.nsteps);
		out.println("nstxout = "+this.nstxout);
		out.println("nstvout = "+ this.nstvout);
		out.println("nstlog = "+this.nstlog);
		out.println("nstenergy = "+this.nstenergy);
		if (this.nstxtcout!=0) 
			out.println("nstxtcout = "+this.nstxtcout);
		out.println("nstlist = "+this.nstlist);
		out.println("coulombtype = "+this.coulombtype);
		out.println("optimize_fft = "+this.optimize_fft);
		if (this.tcoupl!=null)
			out.println("tcoupl = "+this.tcoupl);
		out.println("tc_grps = "+this.tc_grps);
		out.println("tau_t = "+this.tau_t);
		out.println("ref_t = "+this.ref_t);
		if (this.pcoupl!=null)
			out.println("pcoupl = "+this.pcoupl);
		out.printf(Locale.US,"compressibility = %6.1e\n", this.compressibility);
		out.println("ref_p = "+this.ref_p);
		if (this.annealing!=null) 
			out.println("annealing = "+this.annealing);
		if (this.annealing_npoints!=null) 
			out.println("annealing_npoints = "+this.annealing_npoints);
		if (this.annealing_time!=null) 
			out.println("annealing_time = "+this.annealing_time);
		if (this.annealing_temp!=null) 
			out.println("annealing_temp = "+this.annealing_temp);
		out.println("gen_vel = "+this.gen_vel);
		if (this.constraints!=null) 
			out.println("constraints = "+this.constraints);
		out.close();
	}
	
	
	/*----------------------------- main (for testing) ------------------------------ */
	
	public static void main(String[] args) throws Exception {
		File emfile = new File("test_em.mdp");
		File prfile = new File("test_pr.mdp");
		File mdfile = new File("test_md.mdp");
		File mdAnnealFile = new File("test_md_ann.mdp");
		GromacsMDP gmdp = new GromacsMDP();
		gmdp.setEMValues();
		gmdp.writeToFile(emfile);
		gmdp.setPRValues(200);
		gmdp.writeToFile(prfile);
		gmdp.setMDValues(300);
		gmdp.writeToFile(mdfile);
		gmdp.setAnnealValues(300);
		gmdp.writeToFile(mdAnnealFile);
		
	}
}

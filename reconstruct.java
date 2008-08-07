import gnu.getopt.Getopt;
import graphAveraging.ConsensusSquare;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.Formatter;
import java.util.Locale;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.AAinfo;
import proteinstructure.ConformationsNotSameSizeError;
import proteinstructure.PdbLoadError;
import proteinstructure.RIGraph;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbasePdb;

import tinker.TinkerError;
import tinker.TinkerRunner;


public class reconstruct {

	private static final String TINKERBINDIR = "/project/StruPPi/Software/tinker/bin";
	private static final String PRMFILE = "/project/StruPPi/Software/tinker/amber/amber99.prm";

	private static final double DEFAULT_FORCECONSTANT_DISTANCE = TinkerRunner.DEFAULT_FORCECONSTANT_DISTANCE;
	private static final String DEFAULT_CONTACT_TYPE = "Cb";
	private static final double DEFAULT_CUTOFF = 8.0;
	
	private static final TreeMap<Integer, ConsensusSquare> NO_PHIPSI_CONSTRAINTS = null;
	private static final boolean NO_OMEGA_CONSTRAINTS = false;
	
	private static final String PROG_NAME = "reconstruct";
	
	public static void main(String[] args) {
		

		String help = "Reconstructs a protein structure from a contact map using tinker's distgeom.\n" +
			" Two modes of operation:\n" +
			"  a) normal      : NOT IMPLEMENTED YET. specify a sequence and one or more contact map files\n" +
			"  b) benchmarking: specify a pdb code + pdb chain code (-p)\n" +
			"Usage:\n" +
			PROG_NAME+"\n"+
			"  -p <string>   : pdb code + pdb chain code, e.g. 1abcA (benchmarking)\n" +
			"  -b <string>   : base name of output files \n" +
			" [-t <string>]  : contact type. Either 1 contact type: Cb or 2 separated by \n" +
			"                  underscore: Ca_Cg, default: "+DEFAULT_CONTACT_TYPE+" \n" +
			" [-r]           : cross. Use also the cross contact type of the 2 contact types \n" +
			"                  specified. e.g. if Ca_Cg given in -t and -r specified then Ca/Cg \n" +
			"                  contact map will also be used to reconstruct\n" +
			" [-d <floats>]  : distance cutoff(s) comma separated, specify up to 3 distance cutoffs:\n" +
			"                  1st for 1st contact type, 2nd for 2nd contact type 3rd for crossed contact type\n" +
			"                  If 1 specified then it will used for all contact types. Default: "+DEFAULT_CUTOFF+" \n" +
			" [-o <dir>]     : output dir, default: current \n" +
			" [-n <int>]     : number of models to generate, default: 1 \n" +
			" [-m <int>]     : min range, default: 0 \n" +
			" [-M <int>]     : max range, default: 0 \n" +
			" [-f <float>]   : force constant, default: "+DEFAULT_FORCECONSTANT_DISTANCE+" \n" +
			" [-F]           : fast mode: refinement will be done via minimization (faster but \n" +
			"                  worse quality model). Default: slow (refinement via simulated annealing) \n\n";

		String pdbId = null;
		String pdbCode = null;
		String pdbChainCode = null;
		String ct = DEFAULT_CONTACT_TYPE;
		String[] cutoffs = null;
		double cutoff1 = DEFAULT_CUTOFF;
		double cutoff2 = 0.0;
		double cutoff3 = 0.0;
		String outputDir = "."; //default current
		String baseName = null;
		boolean cross = false;
		int n = 1;
		double forceConstant = DEFAULT_FORCECONSTANT_DISTANCE; 
		int minRange = 0;
		int maxRange = 0;
		boolean fast = false;
		
		Getopt g = new Getopt(PROG_NAME, args, "p:t:rd:b:o:n:m:M:f:Fh?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'p':
				pdbId = g.getOptarg();
				break;
			case 't':
				ct = g.getOptarg();
				break;
			case 'r':
				cross = true;
				break;
			case 'd':
				cutoffs = g.getOptarg().split(",");
				break;
			case 'b':
				baseName = g.getOptarg();
				break;
			case 'o':
				outputDir = g.getOptarg();
				break;
			case 'n':
				n = Integer.parseInt(g.getOptarg());
				break;
			case 'm':
				minRange = Integer.parseInt(g.getOptarg());
				break;
			case 'M':
				maxRange = Integer.parseInt(g.getOptarg());
				break;				
			case 'f':
				forceConstant = Double.valueOf(g.getOptarg());
				break;
			case 'F':
				fast = true;
				break;								
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		if (pdbId==null || baseName==null){
			System.err.println("Must specify at least -p and -b");
			System.err.println(help);
			System.exit(1);
		}
		
		if (baseName.contains(".")) {
			System.err.println("Basename can't contain a dot (not allowed by tinker). Exiting");
			System.exit(1);
		}
		
		if (n>999) {
			System.err.println("Maximum number of models is 999. Specify a lower value. Exiting");
			System.exit(1);
		}
		
		// pdb code and chain code
		Pattern p = Pattern.compile("(\\d\\w\\w\\w)(\\w)");
		Matcher m = p.matcher(pdbId);
		if (m.matches()) {
			pdbCode = m.group(1);
			pdbChainCode = m.group(2);
		} else {
			System.err.println("PDB id given not in the correct format, must be in the form 1abcA");
			System.exit(1);
		}
		
		// cutoffs assignment
		if (cutoffs!=null) {
			cutoff1 = Double.parseDouble(cutoffs[0]);
			cutoff2 = cutoff1;
			cutoff3 = cutoff1;
			if (cutoffs.length>1) {
				cutoff2 = Double.parseDouble(cutoffs[1]);
				if (cutoffs.length>2) {
					cutoff3 = Double.parseDouble(cutoffs[2]);
				}
			}
		}

		// contact type assignment
		boolean doublecm = false;
		String ct1 = ct;
		String ct2 = ct;
		String ct3 = null;
		if (ct.contains("_")) {
			ct1 = ct.split("_")[0];
			ct2 = ct.split("_")[1];
			doublecm = true;
		}
		
		if (!AAinfo.isValidContactType(ct1) || !AAinfo.isValidContactType(ct2) || ct.contains("ALL")) {
			System.err.println("Invalid contact type specified. Exiting");
			System.exit(1);
		}
		
		if (cross) {
			ct3 = ct1+"/"+ct2;
		}
		
		
		Pdb pdb = null;
		Pdb mPdb = null;
		try {
			pdb = new PdbasePdb(pdbCode);
			pdb.load(pdbChainCode);
			mPdb = new PdbasePdb(pdbCode);
			mPdb.load(pdbChainCode);
			mPdb.mirror();
		} catch (PdbLoadError e) {
			System.err.println("Error while loading pdb data. Specific error "+e.getMessage());
			System.exit(1);
		} catch (PdbCodeNotFoundError e) {
			System.err.println("Given pdb code "+pdbCode+" couldn't be found in pdbase. Exiting");
			System.exit(1);
		} catch (SQLException e) {
			System.err.println("Problems connecting to database for getting pdb data for "+pdbCode+". Exiting");
			System.exit(1);
		}
		// we also write the file to the out dir so it can be used later for clustering rmsds etc.
		File origPdbFile = new File (outputDir,baseName+".native.pdb");
		try {
			pdb.dump2pdbfile(origPdbFile.getAbsolutePath());
		} catch (IOException e4) {
			System.err.println("Couldn't write original pdb file "+origPdbFile.getAbsolutePath());
			System.err.println("Continuing without it, this is not needed for the rest of the reconstruction process but only for post processing (e.g. comparing rmsds to original)");
		}

		String sequence = pdb.getSequence();

		RIGraph[] graphs =null;
		RIGraph graph1 = pdb.get_graph(ct1, cutoff1);
		if (maxRange>0) graph1.restrictContactsToMaxRange(maxRange);
		if (minRange>0) graph1.restrictContactsToMinRange(minRange);
		RIGraph graph2 = null;
		RIGraph graph3 = null;
		if (doublecm) {
			graph2 = pdb.get_graph(ct2, cutoff2);
			if (maxRange>0) graph2.restrictContactsToMaxRange(maxRange);
			if (minRange>0) graph2.restrictContactsToMinRange(minRange);
			if (cross) {
				graph3 = pdb.get_graph(ct3, cutoff3);
				if (maxRange>0) graph3.restrictContactsToMaxRange(maxRange);
				if (minRange>0) graph3.restrictContactsToMinRange(minRange);
				graphs = new RIGraph[3];
				graphs[0] = graph1;
				graphs[1] = graph2;
				graphs[2] = graph3;
			} else {
				graphs = new RIGraph[2];
				graphs[0] = graph1;
				graphs[1] = graph2;
			}
		} else {
			graphs = new RIGraph[1];
			graphs[0] = graph1;
		}
		
		// defining files
		File reportFile = new File(outputDir,baseName+".report");
		
		// creating TinkerRunner object
		
		TinkerRunner tr = null;
		try {
			tr = new TinkerRunner(TINKERBINDIR, PRMFILE); 
		} catch (FileNotFoundException e3) {
			System.err.println("Couldn't find tinker bin dir "+TINKERBINDIR+". Exiting");
			System.exit(1);
		}
		
		// call reconstruction
		
		try {
			if (fast) {
				// as this is at the moment just a reconstruction benchmarking script it doesn't make sense 
				// at all to use a phi/psi consensus (which would come from templates): we use null and 0 for the 2 phi/psi parameters
				tr.reconstructFast(sequence, graphs, NO_PHIPSI_CONSTRAINTS, NO_OMEGA_CONSTRAINTS, n, forceConstant, 0, outputDir, baseName, false);
			} else {
				tr.reconstruct(sequence, graphs, NO_PHIPSI_CONSTRAINTS, NO_OMEGA_CONSTRAINTS, n, forceConstant, 0, outputDir, baseName, false);
			}
			
		} catch (IOException e) {
			System.err.println("Error while running Tinker reconstruction: " + e.getMessage());
			System.exit(1);
		} catch (TinkerError e) {
			System.err.println("Error while running Tinker reconstruction: " + e.getMessage());
			System.exit(1);
		}
				
		double[] err = tr.getErrorFunctionVal();
		double[] mubv = tr.getMaxUpperBoundViol();
		double[] mlbv = tr.getMaxLowerBoundViol();
		double[] muv = tr.getMaxUpperViol();
		double[] mlv = tr.getMaxLowerViol();
		int[] nubv = tr.getNumUpperBoundViol();
		int[] nlbv = tr.getNumLowerBoundViol();
		int[] nuv = tr.getNumUpperViol();
		int[] nlv = tr.getNumLowerViol();
		double[] rbv = tr.getRmsBoundViol();
		double[] rrv = tr.getRmsRestViol();
		

		// calculate rmsds

		double[] rmsds = new double[n+1];		
		double[] mRmsds = new double[n+1];		
		
		for (int i = 1; i<=n; i++) {
			String ext = new Formatter().format(".%03d",i).toString();
			File outputPdbFile = new File(outputDir, baseName+ext+".pdb");
			try {
				Pdb outputPdb = tr.getStructure(i);
				rmsds[i] = pdb.rmsd(outputPdb, "Ca");
				mRmsds[i] = mPdb.rmsd(outputPdb, "Ca");
			}
			catch (TinkerError e) {
				System.err.println("Error while trying to retrieve results from Tinker: "+ e.getMessage());
			} catch (ConformationsNotSameSizeError e) {
				System.err.println(origPdbFile.getAbsolutePath()+" and "+outputPdbFile.getAbsolutePath()+" don't have the same conformation size, can't calculate rmsd for them.");
			}				
		}					

		// write report file
		
		try {
			PrintWriter reportOut = new PrintWriter(new FileOutputStream(reportFile));
			reportOut.println("#run_id\tcutoff\tcutoff2\tcutoff3\tct\tct2\tct3\tnum_res" +
						"\tresult_id\terror_val" +
						"\tup_bound_viol\tlow_bound_viol\tmax_bound_up\tmax_bound_low\trms_bound" +
						"\tup_viol\tlow_viol\tmax_up\tmax_low\trms_viol\trmsd_to_orig\trmsd_to_mirrored_orig");

			for (int i=1;i<=n;i++){
				String rmsd = String.format(Locale.US,"%6.3f",rmsds[i]);
				String mRmsd = String.format(Locale.US,"%6.3f",mRmsds[i]);
				String errStr = String.format(Locale.US, "%6.3f", err[i]);
				//                     run_id                cutoff      cutoff2      cutoff3      ct1      ct2      ct3         num_res
				reportOut.println(pdbCode+pdbChainCode+"\t"+cutoff1+"\t"+cutoff2+"\t"+cutoff3+"\t"+ct1+"\t"+ct2+"\t"+ct3+"\t"+sequence.length()+"\t"+
				//  result_id     error_val         
						i + "\t" + errStr + "\t" +
				//   up_bound_viol    low_bound_viol   max_bound_up    max_bound_low      rms_bound
						nubv[i] + "\t" + nlbv[i] + "\t" + mubv[i] + "\t" + mlbv[i] + "\t" + rbv[i] + "\t" +
				//      up_viol         low_viol       max_up         max_low         rms_viol
						nuv[i] + "\t" + nlv[i] + "\t" +	muv[i] + "\t" + mlv[i] + "\t"+ rrv[i] + "\t"+
				//      rmsd_to_orig	rmsd_to_mirrored_orig
						rmsd + "\t" + mRmsd);
			}
			reportOut.close();
		} catch (FileNotFoundException e) {
			System.err.println("Couldn't write to report file "+reportFile.getAbsolutePath() +". Error: "+e.getMessage());
		}
		
	}

}

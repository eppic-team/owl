import gnu.getopt.Getopt;
import graphAveraging.ConsensusSquare;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.Locale;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.AAinfo;
import proteinstructure.ConformationsNotSameSizeError;
import proteinstructure.FileRIGraph;
import proteinstructure.GraphFileFormatError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbfilePdb;
import proteinstructure.RIGraph;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbasePdb;

import tinker.TinkerError;
import tinker.TinkerRunner;
import tools.MySQLConnection;


public class reconstruct {

	private static final String TINKERBINDIR = "/project/StruPPi/Software/tinker/bin";
	private static final String PRMFILE = "/project/StruPPi/Software/tinker/amber/amber99.prm";
	
	private static final String PDBASEDB = "pdbase";

	private static final double DEFAULT_FORCECONSTANT_DISTANCE = TinkerRunner.DEFAULT_FORCECONSTANT_DISTANCE;
	private static final String DEFAULT_CONTACT_TYPE = "Cb";
	private static final double DEFAULT_CUTOFF = 8.0;
	private static final String DEFAULT_BASENAME = "rec";
	
	private static final TreeMap<Integer, ConsensusSquare> NO_PHIPSI_CONSTRAINTS = null;
	private static final boolean NO_OMEGA_CONSTRAINTS = false;
	
	private static final String PROG_NAME = "reconstruct";
	
	public static void main(String[] args) {
		

		String help = "\nReconstructs a protein structure from a contact map using tinker's distgeom.\n" +
			" Two modes of operation:\n" +
			"  a) normal      : specify one or more contact map files. The sequence, contact \n" +
			"                   type and cutoff will be taken from the file\n" +
			"  b) benchmarking: specify a pdb code + pdb chain code or a pdb file(-p) and \n" +
			"                   optionally contact type (-t) and cutoff (-d)\n" +
			"Usage:\n" +
			PROG_NAME+" [options] [contact_map_file_1 [contact_map_file_2] [...]] \n"+
			"  -p <string>   : pdb code + pdb chain code, e.g. 1abcA or a pdb file (benchmarking)\n" +
			"                  If in a) i.e. reconstructing from contact map files then the \n" +
			"                  given pdb id/pdb file will be used for rmsd reporting \n" +
			" [-t <string>]  : one or more contact types comma separated (benchmarking). Default: "+DEFAULT_CONTACT_TYPE+" \n" +
			" [-d <floats>]  : one or more distance cutoffs comma separated (benchmarking), \n" +
			"                  matching given contact types. If only one specified then it will be\n" +
			"                  used for all contact types. Default: "+DEFAULT_CUTOFF+" \n" +
			"\n" +
			" [-b <string>]  : base name of output files. Default: "+DEFAULT_BASENAME+" or pdbId given in -p\n" +
			" [-o <dir>]     : output dir. Default: current \n" +
			" [-n <int>]     : number of models to generate. Default: 1 \n" +
			" [-m <int>]     : filter contacts to min range. Default: no filtering \n" +
			" [-M <int>]     : filter contacts to max range. Default: no filtering \n" +
			" [-f <float>]   : force constant. Default: "+DEFAULT_FORCECONSTANT_DISTANCE+" \n" +
			" [-F]           : fast mode: refinement will be done via minimization (faster but \n" +
			"                  worse quality model). Default: slow (refinement via simulate\n" +
			"                  annealing) \n" +
			" [-A]           : if specified reconstruction will be run in parallel (EXPERIMENTAL)\n" +
			" [-g]           : debug mode, prints some debug info\n\n";

		boolean benchmark = true;
		String pdbId = null;
		String pdbCode = null;
		String pdbChainCode = null;
		String[] cts = {DEFAULT_CONTACT_TYPE};
		double[] cutoffs = {DEFAULT_CUTOFF};
		File[] cmFiles = null;
		String outputDir = "."; //default current
		String baseName = null;
		int n = 1;
		double forceConstant = DEFAULT_FORCECONSTANT_DISTANCE; 
		int minRange = 0;
		int maxRange = 0;
		boolean fast = false;
		boolean parallel = false;
		boolean refPdbFromFile = false;
		boolean debug = false;
		
		Getopt g = new Getopt(PROG_NAME, args, "p:b:t:d:o:n:m:M:f:FAgh?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'p':
				pdbId = g.getOptarg();
				break;
			case 'b':
				baseName = g.getOptarg();
				break;				
			case 't':
				cts = g.getOptarg().split(",");
				break;
			case 'd':
				String[] tokens = g.getOptarg().split(",");
				cutoffs=new double[tokens.length];
				for (int i=0;i<tokens.length;i++) 
					cutoffs[i]=Double.parseDouble(tokens[i]);
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
			case 'A':
				parallel = true;
				break;
			case 'g':
				debug = true;
				break;																
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}
		// parsing non-option elements (contact map files)
		if (args.length-g.getOptind()>0) {
			cmFiles = new File[args.length-g.getOptind()];
			for (int index=g.getOptind();index<args.length;index++) {
				cmFiles[index-g.getOptind()]= new File(args[index]);
				benchmark = false;
			}
		}
		
		// input checks
		if (pdbId==null && cmFiles == null) {
			System.err.println("Either a pdb id/file (-p) or at least one contact map file must be given");
			System.err.println(help);
			System.exit(1);
		}
				
		if (n>999) {
			System.err.println("Maximum number of models is 999. Specify a lower value. Exiting");
			System.exit(1);
		}
		
		if (benchmark && cts.length!=cutoffs.length && cutoffs.length!=1) {
			System.err.println("The number of contact types and cutoffs given must be the same. With the exception of only 1 cutoff given that will be taken for all contact types");
			System.exit(1);
		}

		if (benchmark) {
			// assigning default cutoff if only 1 given and multi ct given
			if (cts.length>1 && cutoffs.length==1) {
				double cutoff = cutoffs[0];
				cutoffs = new double[cts.length];
				for (int i=0;i<cts.length;i++) 
					cutoffs[i] = cutoff;
			}
		}
		if (pdbId!=null) {
			// pdb code and chain code
			Pattern p = Pattern.compile("(\\d\\w\\w\\w)(\\w)");
			Matcher m = p.matcher(pdbId);
			if (m.matches()) {
				pdbCode = m.group(1);
				pdbChainCode = m.group(2);
				refPdbFromFile = false;
			} 
			else if (new File(pdbId).exists()) {
				refPdbFromFile = true;
			}
			else {				
				System.err.println("Either PDB id given not in correct format or it is not an existing PDB file: "+pdbId);
				System.exit(1);
			}		
		}		
		
		// setting a default for baseName if it wasn't specified and checking basename
		if (baseName==null) {
			if (benchmark) {
				if (refPdbFromFile) {
					String name = new File(pdbId).getName(); 
					baseName= name.substring(0, name.lastIndexOf("."));
				}
				else baseName=pdbCode+pdbChainCode;
			}
			else baseName = DEFAULT_BASENAME;
		}		
		if (baseName.contains(".")) { 
			System.err.println("Basename ("+baseName+") can't contain a dot (not allowed by tinker). Use a different base name (-b). Exiting");
			System.exit(1);
		}	
		// checking a baseName.pdb file doesn't exist in outdir 
		if (new File(outputDir,baseName+".pdb").exists()) {
			BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
			System.err.print("A file baseName.pdb ("+baseName+".pdb) exists in the output directory. File will be overwritten. Continue (y/n)? ");
			try {
				String answer = br.readLine();
				if (!answer.equals("y")) {
					System.err.println("Use a different basename (-b)");
					System.exit(1);
				}
			} catch (IOException e) {
				System.err.println("Error while trying to read standard input");
				System.exit(1);
			}
		}
				
		// getting pdb data if a pdbId was specified 
		Pdb pdb = null;
		Pdb mPdb = null;
		String sequence = null;
		File origPdbFile = null;
		if (pdbId!=null) {
			try {
				if (refPdbFromFile) {
					origPdbFile = new File(pdbId);
					pdb = new PdbfilePdb(origPdbFile.getAbsolutePath());
					pdb.load(pdb.getChains()[0]);
					mPdb = new PdbfilePdb(origPdbFile.getAbsolutePath());
					mPdb.load(pdb.getChains()[0]);
					mPdb.mirror();
				} else {

					MySQLConnection conn = new MySQLConnection();
					pdb = new PdbasePdb(pdbCode, PDBASEDB, conn);
					pdb.load(pdbChainCode);
					mPdb = new PdbasePdb(pdbCode, PDBASEDB, conn);
					mPdb.load(pdbChainCode);
					mPdb.mirror();
					// we also write the file to the out dir so it can be used later for clustering rmsds etc.
					origPdbFile = new File (outputDir,baseName+".native.pdb");
					try {
						pdb.dump2pdbfile(origPdbFile.getAbsolutePath());
					} catch (IOException e4) {
						System.err.println("Couldn't write original pdb file "+origPdbFile.getAbsolutePath());
						System.err.println("Continuing without it, this is not needed for the rest of the reconstruction process but only for post processing (e.g. comparing rmsds to original)");
					}
				}

			} catch (PdbLoadError e) {
				System.err.println("Error while loading pdb data. Specific error "+e.getMessage());
				System.exit(1);
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Given pdb code "+pdbCode+" couldn't be found in pdbase. Exiting");
				System.exit(1);
			} catch (SQLException e) {
				System.err.println("Problems connecting to database for getting pdb data. Error: "+e.getMessage()+"\nExiting");
				System.exit(1);
			}

			sequence = pdb.getSequence();
		} 

		

		// initialising contact maps
		RIGraph[] graphs = new RIGraph[cts.length];
		
		if (benchmark) {
			for (int i=0;i<cts.length;i++) {
				graphs[i] = pdb.get_graph(cts[i], cutoffs[i]);
			}
		} else {
			graphs = new RIGraph[cmFiles.length];
			cts = new String[cmFiles.length];
			cutoffs = new double[cmFiles.length];
			try {
				for (int i=0;i<cmFiles.length;i++) {
					graphs[i]=new FileRIGraph(cmFiles[i].getAbsolutePath());
					cts[i] = graphs[i].getContactType();
					cutoffs[i] = graphs[i].getCutoff();
					if (cutoffs[i]==0) {
						System.err.println("Contact Map file "+cmFiles[i]+" has no cutoff. Can't continue.");
						System.exit(1);
					}
				}
				if (graphs[0].hasSequence()) {
					sequence = graphs[0].getSequence();
				} else {
					System.err.println("Contact Map file "+cmFiles[0]+" has no sequence. Can't continue.");
					System.exit(1);
				}
				if (graphs.length>1) {
					for (int i=1;i<graphs.length;i++) {
						if (!graphs[i].getSequence().equals(sequence)){
							System.err.println("Sequence of contact map file "+cmFiles[i]+" doesn't coincide with sequence of file "+cmFiles[0]);
							System.exit(1);
						}
					}
				}
			} catch (GraphFileFormatError e) {
				System.err.println("Contact Map file has wrong formatting: "+e.getMessage());
				System.exit(1);
			} catch (IOException e) {
				System.err.println("Couldn't read contact map file: " +e.getMessage());
				System.exit(1);
			}
			
		}
		
		// checking contact types
		for (String ct: cts) {
			if (!AAinfo.isValidContactType(ct) || ct.contains("ALL")) {
				System.err.println("Invalid contact type specified. Exiting");
				System.exit(1);
			}
		}
		
		// restricting to max/min range
		for (RIGraph graph:graphs) {
			if (maxRange>0) graph.restrictContactsToMaxRange(maxRange);
			if (minRange>0) graph.restrictContactsToMinRange(minRange);			
		}

		// defining report file
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
				tr.reconstructFast(sequence, graphs, NO_PHIPSI_CONSTRAINTS, NO_OMEGA_CONSTRAINTS, n, forceConstant, 0, outputDir, baseName, false, parallel);
			} else {
				tr.reconstruct(sequence, graphs, NO_PHIPSI_CONSTRAINTS, NO_OMEGA_CONSTRAINTS, n, forceConstant, 0, outputDir, baseName, false, parallel);
			}
			
		} catch (IOException e) {
			System.err.println("Error while running Tinker reconstruction: " + e.getMessage());
			if (debug) {
				e.printStackTrace();
			}
			System.exit(1);
		} catch (TinkerError e) {
			System.err.println("Error while running Tinker reconstruction: " + e.getMessage());
			if (debug) {
				e.printStackTrace();
			}
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
		
		if (pdbId!=null) {
			for (int i = 1; i<=n; i++) {
				try {
					Pdb outputPdb = tr.getStructure(i);
					rmsds[i] = pdb.rmsd(outputPdb, "Ca");
					mRmsds[i] = mPdb.rmsd(outputPdb, "Ca");
				}
				catch (TinkerError e) {
					System.err.println("Warning: couldn't load tinker pdb output file, error: "+ e.getMessage()+". Can't calculate rmsd for it.");
				} catch (ConformationsNotSameSizeError e) {
					System.err.println(origPdbFile+" and "+tr.getOutPdbFile(i)+" don't have the same conformation size, can't calculate rmsd for them.");
				}				
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
				String runId = "-";
				if (benchmark) {
					if (refPdbFromFile)
						runId = origPdbFile.getName(); 
					else
						runId = pdbCode+pdbChainCode;
				}
				double cutoff1 = cutoffs[0], cutoff2 = 0, cutoff3 = 0;
				String ct1 = cts[0], ct2 = "-", ct3 = "-";
				if (cts.length>1) {
					ct2 = cts[1];
					cutoff2 = cutoffs[1];
					if (cts.length>2) {
						ct3 = cts[2];
						cutoff3 = cutoffs[2];
					}
				}
				//               run_id       cutoff      cutoff2      cutoff3      ct1      ct2      ct3         num_res
				reportOut.println(runId+"\t"+cutoff1+"\t"+cutoff2+"\t"+cutoff3+"\t"+ct1+"\t"+ct2+"\t"+ct3+"\t"+sequence.length()+"\t"+
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
			System.err.println("Couldn't write to report file "+reportFile+". Error: "+e.getMessage());
		}
		
	}

}

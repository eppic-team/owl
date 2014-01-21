import gnu.getopt.Getopt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Locale;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.runners.tinker.TinkerError;
import owl.core.runners.tinker.TinkerRunner;
import owl.core.structure.ContactType;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.features.SecondaryStructure;
import owl.core.structure.graphs.FileRIGraph;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatException;
import owl.core.util.Interval;
import owl.core.util.IntervalSet;
import owl.graphAveraging.ConsensusSquare;
import owl.graphAveraging.PhiPsiAverager;




public class reconstruct {

	private static final String CONFIG_FILE_NAME = "reconstruct.cfg";

	private static final double DEFAULT_FORCECONSTANT_DISTANCE = TinkerRunner.DEFAULT_FORCECONSTANT_DISTANCE;
	private static final int MARGIN_PHIPSI = 3;
	private static final String DEFAULT_CONTACT_TYPE = "Cb";
	private static final double DEFAULT_CUTOFF = 8.0;
	private static final String DEFAULT_BASENAME = "rec";
	
	private static final String PROG_NAME = "reconstruct";
	
	private static String TINKER_BIN_DIR = "/usr/local/bin";
	private static String PRM_FILE = "amber99.prm";
	
	private static String PDB_FTP_URL = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/mmCIF/";
	
	private static final void readTinkerCfgFile(File file) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		Pattern p;
		Matcher m;
		while ((line=br.readLine())!=null) {
			if (line.startsWith("#")) continue;
			p = Pattern.compile("^TINKER_BIN_DIR=(.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				if (new File(m.group(1).trim()).exists()){
					TINKER_BIN_DIR = m.group(1).trim(); 
				} else {
					br.close();
					throw new FileNotFoundException("TINKER_BIN_DIR directory '"+m.group(1).trim()+"' given in config file "+file+" does not exist.");
					//System.err.println("TINKER_BIN_DIR directory '"+m.group(1).trim()+"' given in config file "+file+" does not exist.");
					//System.exit(1);
				}
			}
			p = Pattern.compile("^PRM_FILE=(.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				if (new File(m.group(1).trim()).exists()){
					PRM_FILE = m.group(1).trim(); 
				} else {
					br.close();
					throw new FileNotFoundException("PRM_FILE file '"+m.group(1).trim()+"' given in config file "+file+" does not exist.");
					//System.err.println("PRM_FILE file '"+m.group(1).trim()+"' given in config file "+file+" does not exist.");
					//System.exit(1);
				}
			}
			p = Pattern.compile("^PDB_FTP_URL=(.*)$");
			m = p.matcher(line);
			if (m.matches()) {				
				PDB_FTP_URL = m.group(1).trim(); 
			}			
		}
		br.close();
		
	}
	
	public static void main(String[] args) {
		

		String help = "\nReconstructs a protein structure from a contact map using TINKER's distgeom.\n\n" +
			"This program requires a local installation of the TINKER package. The bin directory and\n" +
			"PRM file used can be specified in the "+CONFIG_FILE_NAME+" file in the current directory\n" +
			"or user's home directory\n\n" +
			" Two modes of operation:\n" +
			"  a) normal      : specify one or more contact map files. The sequence, contact \n" +
			"                   type and cutoff will be taken from the file\n" +
			"  b) benchmarking: specify a pdb code + pdb chain code or a pdb file(-p) and \n" +
			"                   optionally contact type (-t) and cutoff (-d)\n" +
			"Usage:\n\n" +
			" "+PROG_NAME+" [options] [contact_map_file_1 [contact_map_file_2] [...]] \n\n"+
			"  -p <string>     : pdb code + pdb chain code, e.g. 1abcA or a pdb file (benchmarking).\n" +
			"                    The PDB data will be downloaded from the PDB's ftp server.\n" +
			"                    If in a) i.e. reconstructing from contact map files then the \n" +
			"                    given pdb id/pdb file will be used for rmsd reporting \n" +
			" [-t <string>]    : one or more contact types comma separated (benchmarking). Default: "+DEFAULT_CONTACT_TYPE+" \n" +
			" [-d <floats>]    : one or more distance cutoffs comma separated (benchmarking), \n" +
			"                    matching given contact types. If only one specified then it will be\n" +
			"                    used for all contact types. Default: "+DEFAULT_CUTOFF+" \n" +
			"\n" +
			" [-i <intervals>] : use phi/psi restraints from given structure (needs -p). Specify a\n" +
			"                    set of intervals from the given structure from which the phi/psi\n" +
			"                    values will be taken, e.g.: 3-23,30-35,40-50\n" +
			" [-c <file>]      : name of a psipred horizontal file containing secondary structure\n" +
			"                    prediction used to create additional phi/psi constraints. \n" +
			"                    Contraints are NOT applied in fast mode (-F) \n" +
			" [-e]             : restrain omega torsion angles to trans conformation\n" +
			"\n"+
			" [-b <string>]    : base name of output files. Default: "+DEFAULT_BASENAME+" or pdbId given in -p\n" +
			" [-o <dir>]       : output dir. If option -A (parallel) is used then this directory MUST \n" +
			"                    be a globally accessible one (all nodes in cluster must be able to \n" +
			"                    read/write to it).\n" +
			"                    Default: current \n" +
			" [-n <int>]       : number of models to generate. Default: 1 \n" +
			" [-m <int>]       : filter contacts to min range. Default: no filtering \n" +
			" [-M <int>]       : filter contacts to max range. Default: no filtering \n" +
			" [-f <float>]     : force constant. Default: "+DEFAULT_FORCECONSTANT_DISTANCE+" \n" +
			" [-F]             : fast mode: refinement will be done via minimization (faster but \n" +
			"                    worse quality model). Default: slow (refinement via simulate\n" +
			"                    annealing) \n" +
			" [-A]             : if specified reconstruction will be run in parallel using the \n" +
			"                    Sun Grid Engine job scheduler (EXPERIMENTAL)\n" +
			" [-g]             : debug mode, prints some debug info\n\n";

		boolean benchmark = true;
		String pdbId = null;
		String pdbCode = null;
		String pdbChainCode = null;
		String[] cts = {DEFAULT_CONTACT_TYPE};
		double[] cutoffs = {DEFAULT_CUTOFF};
		File[] cmFiles = null;
		IntervalSet intervals = null;
		boolean forceTransOmega = false;
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
		String secondaryStructureFile = null;
		
		Getopt g = new Getopt(PROG_NAME, args, "p:b:t:d:o:i:c:en:m:M:f:FAgh?");
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
			case 'c':
				secondaryStructureFile = g.getOptarg();
				break;
			case 'i':
				tokens = g.getOptarg().split(",");
				intervals = new IntervalSet();
				for (String token:tokens) {
					intervals.add(new Interval(Integer.parseInt(token.split("-")[0]),
											   Integer.parseInt(token.split("-")[1])));
				}
				break;
			case 'e':
				forceTransOmega = true;
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
		
		if(pdbId==null && intervals!=null) {
			System.err.println("Can't specify phi/psi restraints (-i) without specifying a structure (-p).");
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
		
		// reading config file
		boolean cfgFileFound = false;
		File cfgFile = new File(CONFIG_FILE_NAME);
		if (cfgFile.exists()) {
			try {
				readTinkerCfgFile(cfgFile);
				cfgFileFound = true;
			} catch(IOException e) {
				System.err.println("Error while reading config file "+cfgFile+": "+e.getMessage());
			}
		}
		cfgFile = new File(System.getProperty("user.home"),CONFIG_FILE_NAME);
		if (cfgFile.exists()) {
			try {
				readTinkerCfgFile(cfgFile);
				cfgFileFound = true;
			} catch(IOException e) {
				System.err.println("Error while reading config file "+cfgFile+": "+e.getMessage());
			}
		}
		if (!cfgFileFound) {
			System.err.println("Could not find a config file. Using default values:");
			System.err.println("  TINKER_BIN_DIR="+TINKER_BIN_DIR);
			System.err.println("  PRM_FILE="+PRM_FILE);
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
		PdbChain pdb = null;
		PdbChain mPdb = null;
		String sequence = null;
		File origPdbFile = null;
		if (pdbId!=null) {
			try {
				if (refPdbFromFile) {
					origPdbFile = new File(pdbId);
					PdbAsymUnit fullpdb = new PdbAsymUnit(origPdbFile);
					pdb = fullpdb.getFirstChain();
					mPdb = pdb.copy(fullpdb);
					mPdb.mirror();
				} else {

					System.out.println("Downloading PDB entry "+pdbCode+" from PDB's ftp");
					File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".cif");
					cifFile.deleteOnExit();
					PdbAsymUnit.grabCifFile(null, PDB_FTP_URL, pdbCode, cifFile, true);
					PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile);
					pdb = fullpdb.getChain(pdbChainCode);
					mPdb = pdb.copy(fullpdb);
					mPdb.mirror();
					System.out.println("Done");

					// we also write the file to the out dir so it can be used later for clustering rmsds etc.
					origPdbFile = new File (outputDir,baseName+".native.pdb");
					try {
						pdb.writeToPDBFile(origPdbFile);
					} catch (IOException e4) {
						System.err.println("Couldn't write original pdb file "+origPdbFile.getAbsolutePath());
						System.err.println("Continuing without it, this is not needed for the rest of the reconstruction process but only for post processing (e.g. comparing rmsds to original)");
					}
				}

			} catch (PdbLoadException e) {
				System.err.println("Error while loading pdb data. Specific error "+e.getMessage());
				System.exit(1);
			} catch (IOException e) {
				System.err.println("Problems getting pdb data from PDB's ftp. Error: "+e.getMessage()+"\nExiting");
				System.exit(1);				
			} catch (FileFormatException e) {
				System.err.println("Problems loading data from file. Error: "+e.getMessage()+"\nExiting");
				System.exit(1);				
			}

			sequence = pdb.getSequence().getSeq();
			
			// checking that intervals given with -i are valid for this pdb
			if (intervals!=null) {
				int lastResser = pdb.getFullLength();
				for (Interval interv:intervals) {
					if (interv.end>lastResser || interv.beg>lastResser) {
						System.err.println("Specified interval (-i): "+interv+" is not valid for the specified pdb structure "+pdbId);
						System.exit(1);
					}
				}
			}
		} 

		

		// initialising contact maps
		RIGraph[] graphs = new RIGraph[cts.length];
		
		if (benchmark) {
			for (int i=0;i<cts.length;i++) {
				graphs[i] = pdb.getRIGraph(cts[i], cutoffs[i]);
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
			} catch (FileFormatException e) {
				System.err.println("Contact Map file has wrong formatting: "+e.getMessage());
				System.exit(1);
			} catch (IOException e) {
				System.err.println("Couldn't read contact map file: " +e.getMessage());
				System.exit(1);
			}
			
		}
		
		// checking contact types
		for (String ct: cts) {
			if (!ContactType.isValidContactType(ct) || ct.contains("ALL")) {
				System.err.println("Invalid contact type specified. Exiting");
				System.exit(1);
			}
		}
		
		// restricting to max/min range
		for (RIGraph graph:graphs) {
			if (maxRange>0) graph.restrictContactsToMaxRange(maxRange);
			if (minRange>0) graph.restrictContactsToMinRange(minRange);			
		}
		// loading secondary structure prediction
		SecondaryStructure sec = null;
		if (secondaryStructureFile != null) {
			try {
				sec = new SecondaryStructure(new File(secondaryStructureFile));
			} catch (IOException e) {
				System.err.println("Error reading secondary structure prediction. Make sure it's a valid psipred horizontal (.horiz) file. Exiting.");
				System.exit(1);
			}
		}
		// defining report file
		File reportFile = new File(outputDir,baseName+".report");
		
		// creating TinkerRunner object		
		TinkerRunner tr = null;
		try {
			tr = new TinkerRunner(TINKER_BIN_DIR, PRM_FILE); 
		} catch (FileNotFoundException e3) {
			System.err.println(e3.getMessage()+". Are the tinker paths set up in the config file?\nExiting");
			System.exit(1);
		}
		
		// call reconstruction		
		try {
			TreeMap<Integer,ConsensusSquare> phiPsiConstraints = null;
			if (pdbId!=null && intervals!=null) {
				phiPsiConstraints = PhiPsiAverager.getPhiPsiForInterval(pdb, MARGIN_PHIPSI, intervals);
			}
			
			if (sec != null) {
				phiPsiConstraints = sec.getPhiPsiConstraints();
			}
			if (fast) {
				// as this is at the moment just a reconstruction benchmarking script it doesn't make sense 
				// at all to use a phi/psi consensus (which would come from templates): we use null and 0 for the 2 phi/psi parameters
				tr.reconstructFast(sequence, graphs, phiPsiConstraints, forceTransOmega, n, forceConstant, 0, outputDir, baseName, false, parallel);
			} else {
				tr.reconstruct(sequence, graphs, phiPsiConstraints, forceTransOmega, n, forceConstant, 1, outputDir, baseName, false, parallel);
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
		} catch (FileFormatException e) {
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
					PdbChain outputPdb = tr.getStructure(i);
					rmsds[i] = pdb.rmsd(outputPdb, "Ca");
					mRmsds[i] = mPdb.rmsd(outputPdb, "Ca");
				}
				catch (TinkerError e) {
					System.err.println("Warning: couldn't load tinker pdb output file, error: "+ e.getMessage()+". Can't calculate rmsd for it.");
				} catch (IllegalArgumentException e) {
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
				} else {
					// if reconstructing from files (normal mode) we take the name of the first file as id
					runId = cmFiles[0].getName();
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

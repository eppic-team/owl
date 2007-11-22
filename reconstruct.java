import gnu.getopt.Getopt;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.Formatter;
import java.util.Locale;

import proteinstructure.AAinfo;
import proteinstructure.ConformationsNotSameSizeError;
import proteinstructure.RIGraph;
import proteinstructure.Pdb;
import proteinstructure.PdbChainCodeNotFoundError;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbaseInconsistencyError;
import proteinstructure.PdbasePdb;

import tinker.TinkerError;
import tinker.TinkerRunner;


public class reconstruct {

	private static final String TINKERBINDIR = "/project/StruPPi/Software/tinker/bin";
	private static final String PRMFILE = "/project/StruPPi/Software/tinker/amber/amber99.prm";
	
	public static void main(String[] args) {
		
		String programName = reconstruct.class.getName();
		String help = "Usage:\n" +
			programName+" -p <pdb code> -c <pdb chain code> -t <contact_type> [-r] -d <distance cutoff 1> -D <distance cutoff 2> -i <distance cutoff 3> -b <base name> -o <output dir> [-n <number of models>]\n"; 

		String pdbCode = "";
		String pdbChainCode = "";
		String ct = "";
		double cutoff1 = 0.0;
		double cutoff2 = 0.0;
		double cutoff3 = 0.0;
		String outputDir = "";
		String baseName = "";
		boolean cross = false;
		int n = 1;
		//double forceConstant = 100.0;
		
		Getopt g = new Getopt(programName, args, "p:c:d:t:rb:o:d:D:i:n:h?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'p':
				pdbCode = g.getOptarg();
				break;
			case 'c':
				pdbChainCode = g.getOptarg();
				break;
			case 't':
				ct = g.getOptarg();
				break;
			case 'r':
				cross = true;
				break;
			case 'd':
				cutoff1 = Double.valueOf(g.getOptarg());
				break;
			case 'D':
				cutoff2 = Double.valueOf(g.getOptarg());
				break;
			case 'i':
				cutoff3 = Double.valueOf(g.getOptarg());
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
			//case 'f':
			//	forceConstant = Double.valueOf(g.getOptarg());
			//	break;				
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		if (pdbCode.equals("") || pdbChainCode.equals("") || ct.equals("") || cutoff1==0.0 || outputDir.equals("") || baseName.equals("")){
			System.err.println("Must specify at least -p, -c, -t, -d, -o and -b");
			System.err.println(help);
			System.exit(1);
		}
		
		if (baseName.contains(".")) {
			System.err.println("Basename can't contain a dot (not allowed by tinker). Exiting");
			System.exit(1);
		}
		
		if (!AAinfo.isValidContactType(ct) || ct.contains("ALL")) {
			System.err.println("Invalid contact type specified. Exiting");
			System.exit(1);
		}
		
		if (n>999) {
			System.err.println("Maximum number of models is 999. Specify a lower value. Exiting");
			System.exit(1);
		}
		
		boolean doublecm = false;
		String ct1 = ct;
		String ct2 = ct;
		String ct3 = null;
		if (ct.contains("_")) {
			ct1 = ct.split("_")[0];
			ct2 = ct.split("_")[1];
			doublecm = true;
		}
		if (cross) {
			ct3 = ct1+"/"+ct2;
		}
		
		if (cutoff2==0.0) {
			cutoff2 = cutoff1;
		}
		if (cutoff3==0.0) {
			cutoff3 = cutoff1;
		}
		
		Pdb pdb = null;
		try {
			pdb = new PdbasePdb(pdbCode,pdbChainCode);
		} catch (PdbaseInconsistencyError e) {
			System.err.println("Pdbase inconsistency for structure "+pdbCode+". Can't continue, exiting");
			System.exit(1);
		} catch (PdbCodeNotFoundError e) {
			System.err.println("Given pdb code "+pdbCode+" couldn't be found in pdbase. Exiting");
			System.exit(1);
		} catch (SQLException e) {
			System.err.println("Problems connecting to database for getting pdb data for "+pdbCode+". Exiting");
			System.exit(1);
		} catch (PdbChainCodeNotFoundError e) {
			System.err.println("Given pdb chain code "+pdbChainCode+" couldn't be found for pdb code "+pdbCode+". Exiting");
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
		RIGraph graph2 = null;
		RIGraph graph3 = null;
		if (doublecm) {
			graph2 = pdb.get_graph(ct2, cutoff2);		
			if (cross) {
				graph3 = pdb.get_graph(ct3, cutoff3);
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
			tr.reconstruct(sequence, graphs, n, outputDir, baseName, false);
		} catch (IOException e) {
			System.err.println("Error while running Tinker reconstruction: " + e.getMessage());
		} catch (TinkerError e) {
			System.err.println("Error while running Tinker reconstruction: " + e.getMessage());
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
		
		for (int i = 1; i<=n; i++) {
			String ext = new Formatter().format(".%03d",i).toString();
			File outputPdbFile = new File(outputDir, baseName+ext+".pdb");
			try {
				Pdb outputPdb = tr.getStructure(i);
				rmsds[i] = pdb.rmsd(outputPdb, "Ca");
			}
			catch (TinkerError e) {
				System.err.println("Error while trying to retrieve results from Tinker: + e.getMessage()");
			} catch (ConformationsNotSameSizeError e) {
				System.err.println(origPdbFile.getAbsolutePath()+" and "+outputPdbFile.getAbsolutePath()+" don't have the same conformation size, can't calculate rmsd for them.");
			}				
		}					

		// write report file
		
		try {
			PrintWriter reportOut = new PrintWriter(new FileOutputStream(reportFile));
			reportOut.println("run_id\tcutoff\tcutoff2\tcutoff3\tct\tct2\tct3\tnum_res" +
						"\tresult_id\terror_val" +
						"\tup_bound_viol\tlow_bound_viol\tmax_bound_up\tmax_bound_low\trms_bound" +
						"\tup_viol\tlow_viol\tmax_up\tmax_low\trms_viol\trmsd_to_orig");

			for (int i=1;i<=n;i++){
				String rmsd = String.format(Locale.US,"%6.3f",rmsds[i]);
				String errStr = String.format(Locale.US, "%6.3f", err[i]);
				//                         run_id                cutoff      cutoff2      cutoff3      ct1      ct2      ct3         num_res
				reportOut.println(pdbCode+"_"+pdbChainCode+"\t"+cutoff1+"\t"+cutoff2+"\t"+cutoff3+"\t"+ct1+"\t"+ct2+"\t"+ct3+"\t"+sequence.length()+"\t"+
				//  result_id     error_val         
						i + "\t" + errStr + "\t" +
				//   up_bound_viol    low_bound_viol   max_bound_up    max_bound_low      rms_bound
						nubv[i] + "\t" + nlbv[i] + "\t" + mubv[i] + "\t" + mlbv[i] + "\t" + rbv[i] + "\t" +
				//      up_viol         low_viol       max_up         max_low         rms_viol
						nuv[i] + "\t" + nlv[i] + "\t" +	muv[i] + "\t" + mlv[i] + "\t"+ rrv[i] + "\t"+
				//      rmsd_to_orig
						rmsd);
			}
			reportOut.close();
		} catch (FileNotFoundException e) {
			System.err.println("Couldn't write to report file "+reportFile.getAbsolutePath() +". Error: "+e.getMessage());
		}
		
	}

}

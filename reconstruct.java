import gnu.getopt.Getopt;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Formatter;

import proteinstructure.ConformationsNotSameSizeError;
import proteinstructure.Graph;
import proteinstructure.Pdb;
import proteinstructure.PdbChainCodeNotFoundError;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbaseInconsistencyError;
import proteinstructure.PdbasePdb;
import proteinstructure.PdbfileFormatError;
import proteinstructure.PdbfilePdb;

import tinker.ConstraintsMaker;
import tinker.TinkerError;
import tinker.TinkerRunner;


public class reconstruct {

	
	private static final String TINKERBINDIR = "/project/StruPPi/Software/tinker/bin";
	private static final String PRMFILE = "/project/StruPPi/Software/tinker/amber/amber99.prm";
	private static final String PRMTYPE = "amber";
	
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
		String sequence = pdb.getSequence();

		Graph graph1 = pdb.get_graph(ct1, cutoff1);
		Graph graph2 = null;
		Graph graph3 = null;
		if (doublecm) {
			graph2 = pdb.get_graph(ct2, cutoff2);
		}
		if (cross) {
			graph3 = pdb.get_graph(ct3, cutoff3);
		}

		
		// defining files
		File logFile = new File(outputDir,"tinker.log");
		File prmFile = new File(PRMFILE);
		File xyzFile = new File(outputDir,baseName+".xyz");
		File seqFile = new File(outputDir,baseName+".seq");
		File pdbFile = new File(outputDir,baseName+".pdb");
		File keyFile = new File(outputDir,baseName+".key");
		
		// creating TinkerRunner object
		TinkerRunner tr = null;
		try {
			tr = new TinkerRunner(TINKERBINDIR,PRMFILE,logFile);
		} catch (FileNotFoundException e3) {
			System.err.println("Couldn't find tinker bin dir "+TINKERBINDIR+". Exiting");
			System.exit(1);
		}
		
		// 1. run tinker's protein program	
		try {
			tr.runProtein(sequence, outputDir, baseName );
		} catch (IOException e2) {
			System.err.println("Couldn't read file to run 'protein' "+seqFile.getAbsolutePath());
			System.err.println("Exiting");
			System.exit(1);			
		} catch (TinkerError e2) {
			System.err.println("Tinker error while running 'protein', check log file "+logFile.getAbsolutePath()+". Exiting");
			System.exit(1);			
		}

		// 1a. convert xyz file to pdb to be able to map atom serials after
		try {
			tr.runXyzpdb(xyzFile, seqFile, pdbFile);
		} catch (IOException e1) {
			System.err.println("Couldn't read files "+xyzFile.getAbsolutePath()+" or "+seqFile.getAbsolutePath()+" or write to "+pdbFile.getAbsolutePath()+" for running 'xyzpdb'");
			System.err.println("Exiting");
			System.exit(1);
		} catch (TinkerError e1) {
			System.err.println("Tinker error while running xyzpdb, check log file "+logFile.getAbsolutePath()+". Exiting");
			System.exit(1);
		}

		// 2. creating constraints into key file
		ConstraintsMaker cm = null;
		try {
			cm = new ConstraintsMaker(pdbFile,xyzFile,prmFile,PRMTYPE,keyFile);
		} catch (IOException e3) {
			System.err.println("Couldn't read files "+xyzFile.getAbsolutePath()+", "+pdbFile.getAbsolutePath()+" or, "+prmFile.getAbsolutePath()+" write to "+keyFile.getAbsolutePath()+" for creating distance constraints");
			System.err.println("Exiting");
			System.exit(1);
		} catch (PdbfileFormatError e3) {
			System.err.println("pdb file "+pdbFile.getAbsolutePath()+" converted from "+xyzFile.getAbsolutePath()+" doesn't seem to be in the right format. Check log? ("+logFile.getAbsolutePath()+"). Exiting");
			System.exit(1);
		} 
		try {
			cm.createConstraints(graph1);
			if (doublecm) cm.createConstraints(graph2);
			if (cross) cm.createConstraints(graph3);
		} catch (Exception e2) {
			System.err.println("Invalid contact type for writing constraints.");
			System.err.println("Error: "+e2.getMessage());
			System.err.println("Exiting");
			System.exit(1);
		}
		cm.closeKeyFile();

		// 3. run tinker's distgeom
		try {
			System.out.println("Running distgeom...");
			tr.runDistgeom(xyzFile, outputDir, baseName, n);
		} catch (TinkerError e1) {
			System.err.println(e1.getMessage());
			System.err.println("Exiting");
			System.exit(1);			
		} catch (IOException e1) {
			System.err.println("Couldn't read files "+xyzFile.getAbsolutePath()+" or write 'distgeom' output files to output dir "+outputDir);
			System.err.println("Exiting");
			System.exit(1);
		} catch(InterruptedException e) {
			System.err.println("Distgeom was interrupted:" + e.getMessage());
			System.err.println("Exiting.");
			System.exit(1);
		}
		
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
		

		// 4. converting xyz output files to pdb files and calculating rmsds

		double[] rmsds = new double[n+1];		
		
			for (int i = 1; i<=n; i++) {
				String ext = new Formatter().format(".%03d",i).toString();
				File outputXyzFile = new File(outputDir, baseName+ext);
				File outputPdbFile = new File(outputDir, baseName+ext+".pdb");
				try {
					tr.runXyzpdb(outputXyzFile, seqFile, outputPdbFile);
					
					Pdb outputPdb = new PdbfilePdb(outputPdbFile.getAbsolutePath(),"NULL");
					rmsds[i] = pdb.rmsd(outputPdb, "Ca");

				} catch (IOException e) {
					System.err.println("Couldn't read file "+outputXyzFile.getAbsolutePath()+", or "+seqFile.getAbsolutePath()+", or write to "+outputPdbFile.getAbsolutePath()+" while converting with 'xyzpdb'");
					System.err.println("Can't calculate rmsd for it");
				} catch (TinkerError e) {
					System.err.println("Tinker error while running 'xyzpdb' to convert"+outputXyzFile.getAbsolutePath()+", check log file "+logFile.getAbsolutePath());
					System.err.println("Can't calculate rmsd for it");
				}
				catch (PdbfileFormatError e) {
					System.err.println("Output pdb file "+outputPdbFile.getAbsolutePath()+" doesn't seem to be in the correcet format. Can't calculate rmsd for it");
				} catch (PdbChainCodeNotFoundError e) {
					// this shouldn't happen, chain code is hard coded, we throw stack trace and continue if it happens
					e.printStackTrace();
				} catch (ConformationsNotSameSizeError e) {
					System.err.println(pdbFile.getAbsolutePath()+" and "+outputPdbFile.getAbsolutePath()+" don't have the same conformation size, can't calculate rmsd for them.");
				}				
				
				tr.closeLog();
			}					

		
		// 6. report
		System.out.println("accession_code\tchain_pdb_code\tcutoff\tcutoff2\tcutoff3\tct\tct2\tct3" +
				"\tresult_id\tup_bound_viol\tlow_bound_viol\tmax_bound_up\tmax_bound_low\trms_bound" +
				"\tup_viol\tlow_viol\tmax_up\tmax_low\trms_viol\trmsd_to_orig");
		
		for (int i=1;i<=n;i++){
			System.out.println(pdbCode+"\t"+pdbChainCode+"\t"+cutoff1+"\t"+cutoff2+"\t"+cutoff3+"\t"+ct1+"\t"+ct2+"\t"+ct3+"\t"+
					i+"\t"+nubv[i]+"\t"+nlbv[i]+"\t"+mubv[i]+"\t"+mlbv[i]+"\t"+muv[i]+"\t"+mlv[i]+"\t"+rbv[i]+"\t"+
					"\t"+nuv[i]+"\t"+nlv[i]+"\t"+rrv[i]+"\t"+rmsds[i]);
		}
		
	}

}

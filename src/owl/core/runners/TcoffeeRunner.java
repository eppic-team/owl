package owl.core.runners;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;

import owl.core.sequence.Sequence;
import owl.core.sequence.alignment.AlignmentConstructionError;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.util.FileFormatError;


/**
 * A class to run tcoffee.
 * @author duarte
 *
 */
public class TcoffeeRunner {

	// t-coffee seems to have a bug in outputting directly to fasta_aln format: it will add 
	// one extra character to all sequences except for the target. This is why we use clustalw
	private static final String DEFAULT_SEQ2PROF_OUTFORMAT = "clustalw";
	
	private static final boolean DEBUG = false;
	private File tcofProg;
	private File logFile;

	
	public TcoffeeRunner(File tcofProg, File logFile) {
		this.tcofProg = tcofProg;
		this.logFile = logFile;
	}

	/**
	 * Aligns a sequence to a profile (in the form of a multiple sequence alignment)
	 * by using tcoffee's profile comparison
	 * @param seq the sequence object (its tag will be the one used in the output Alignment)
	 * @param profile the multiple sequence alignment representing the profile to align to 
	 * @return
	 * @throws TcoffeeError if tcoffee fails to run
	 * @throws IOException if problem while reading/writing temp files needed to run tcoffee
	 */
	public MultipleSequenceAlignment alignSequence2Profile(Sequence seq, MultipleSequenceAlignment profile) throws TcoffeeError, IOException {
		File inFile = File.createTempFile("tcof.", ".in");
		if (!DEBUG) inFile.deleteOnExit();
		seq.writeToFastaFile(inFile);
		File outFile = File.createTempFile("tcof.", ".out");
		if (!DEBUG) outFile.deleteOnExit();
		File profileFile = File.createTempFile("tcof", "profile");
		if (!DEBUG) profileFile.deleteOnExit();
		PrintStream out = new PrintStream(profileFile);
		profile.writeFasta(out, 80, true);
		out.close();
		runTcoffee(inFile, outFile, DEFAULT_SEQ2PROF_OUTFORMAT, profileFile);
		
		MultipleSequenceAlignment al =  null;
		try {
			al = new MultipleSequenceAlignment(outFile.getAbsolutePath(), MultipleSequenceAlignment.CLUSTALFORMAT);
		} catch (AlignmentConstructionError e) {
			throw new TcoffeeError("Couldn't construct Alignment from Tcoffee output alignment "+outFile+". Error "+e.getMessage());			
		} catch (FileFormatError e) {
			throw new TcoffeeError("Couldn't construct Alignment from Tcoffee output alignment "+outFile+". Error "+e.getMessage());
		}
		
		return al;
	}
	
	private void runTcoffee(File inFile, File outFile, String outFormat, File profileFile) throws TcoffeeError {
		String profStr = "";
		if (profileFile!=null) {
			profStr = "-profile "+profileFile+" -profile_comparison=full50";
		}
		String cmdLine = tcofProg + " "+ inFile + " "+ profStr + " -output=" +outFormat+" -outfile="+outFile;
		
		try {
			PrintWriter tcofLog = new PrintWriter(logFile);
			
			Process proc = Runtime.getRuntime().exec(cmdLine);

			// logging and capturing stdout/stderr
			BufferedReader stdout = new BufferedReader(new InputStreamReader(proc.getInputStream()));
			BufferedReader stderr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

			String line;

			tcofLog.println("#cmd: "+cmdLine);
			tcofLog.println("#################");
			tcofLog.println("# "+tcofProg+" stdout ");
			tcofLog.println("#################");
			while((line = stdout.readLine()) != null) {
				tcofLog.println(line);
			}
			tcofLog.println("#################");
			tcofLog.println("# "+tcofProg+" stderr ");
			tcofLog.println("#################");
			while((line = stderr.readLine()) != null) {
				tcofLog.println(line);
			}


			try {
				int exitValue = proc.waitFor();
				// throwing exception if exit state is not 0 
				if (exitValue!=0) {
					tcofLog.flush();
					throw new TcoffeeError(tcofProg + " exited with value "+exitValue+". Revise log file "+logFile);
				}
			} catch (InterruptedException e) {
				System.err.println("Unexpected error while waiting for "+tcofProg+" to exit. Error: "+e.getMessage());
				System.exit(1);
			}

			tcofLog.close();
			
		} catch (IOException e) {
			throw new TcoffeeError("IO error while trying to run "+tcofProg+": "+e.getMessage());
		}
		
		
	}
	
	public static void main(String[] args) throws Exception {
		File tcofProg = new File("/project/StruPPi/bin/t_coffee");
		File logFile = new File("/tmp/tcoffee.log");
		TcoffeeRunner tcr = new TcoffeeRunner(tcofProg, logFile);
		File seqFile = new File("/project/StruPPi/CASP7/targets/T0290.fa");
		File alFile = new File("/project/StruPPi/CASP8/dryrun/casp7_bla_pb_gtg/bla_max10_e5_s1000_cct4/T0290/T0290.templates_aln.fasta");
		Sequence seq = new Sequence();
		seq.readFromFastaFile(seqFile);
		MultipleSequenceAlignment al = new MultipleSequenceAlignment(alFile.getAbsolutePath(),MultipleSequenceAlignment.FASTAFORMAT);
		MultipleSequenceAlignment outAl = tcr.alignSequence2Profile(seq, al);
		outAl.writeFasta(System.out, 60, true);
		
	}
	
}

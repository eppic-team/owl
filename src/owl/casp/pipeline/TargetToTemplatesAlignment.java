package owl.casp.pipeline;

import gnu.getopt.Getopt;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import owl.core.sequence.Sequence;
import owl.core.sequence.TcoffeeError;
import owl.core.sequence.TcoffeeRunner;
import owl.core.structure.Alignment;
import owl.core.structure.AlignmentConstructionError;
import owl.core.structure.FileFormatError;


public class TargetToTemplatesAlignment {

	protected static final String TARGET2TEMPL_SUFFIX = ".target2templ.fasta";

	protected static final File TCOFPROG = new File("/project/StruPPi/bin/t_coffee");
	
	private static final String PROGRAM_NAME = "target2templatesAln";	
	
	// members
	private Sequence seq;
	private Alignment templatesAln;
	private String baseName;
	private File outDir;
	
	
	public TargetToTemplatesAlignment(Sequence seq, Alignment templatesAln,String baseName, File outDir) {
		this.seq = seq;
		this.templatesAln = templatesAln;
		this.baseName = baseName;
		this.outDir = outDir;
	}

	public Alignment run() throws TcoffeeError, IOException {
		File tcofLogFile = new File(outDir, baseName+".tcof.log");
		
		TcoffeeRunner tcr = new TcoffeeRunner(TCOFPROG,tcofLogFile);
		Alignment al = tcr.alignSequence2Profile(seq, templatesAln);

		return al;
	}
	
	public static void writeTargetToTemplAl(File outDir, String baseName, Alignment al) throws IOException {
		File alOutFile = new File(outDir, baseName+TARGET2TEMPL_SUFFIX);
		PrintStream alOut = new PrintStream(alOutFile);
		al.writeFasta(alOut, 80, true);
	}
	
	public static void main(String[] args) {
		File seqFile = null;
		File tempsAlnFile = null;
		String baseName = null;
		File outDir = new File(".");		

		String help = "Usage: \n" +
		PROGRAM_NAME+"\n" +
		"   -i :  file with input target sequence\n" +
		"   -n :  file with input templates alignment \n"+
		"  [-b]:  basename for output files. Default: basename of input sequence file \n"+
		"  [-o]:  output dir, where output files will be written. Default: current dir \n\n";


		Getopt g = new Getopt(PROGRAM_NAME, args, "i:n:b:o:h?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'i':
				seqFile = new File(g.getOptarg());
				break;
			case 'n':
				tempsAlnFile = new File(g.getOptarg());
				break;				
			case 'b':
				baseName = g.getOptarg();
				break;				
			case 'o':
				outDir = new File(g.getOptarg());
				break;
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		if ((seqFile==null) || (tempsAlnFile==null)) {
			System.err.println("At least a sequence file (-i) and a templates alignment file (-n) must be specified");
			System.err.println(help);
			System.exit(1);
		}
		
		if (baseName==null) {
			baseName = seqFile.getName().substring(0, seqFile.getName().lastIndexOf("."));
		}
		
		
		Sequence seq = new Sequence();
		try {
			seq.readFromFastaFile(seqFile);
		} catch (IOException e) {
			System.err.println("Problem while reading sequence file "+seqFile+". Error: "+e.getMessage());
			System.exit(1);
		} catch (FileFormatError e) {
			System.err.println("Problem while reading sequence file "+seqFile+". Error: "+e.getMessage());
			System.exit(1);
		}
		
		Alignment templatesAln = null;
		try {
			templatesAln = new Alignment(tempsAlnFile.getAbsolutePath(),Alignment.FASTAFORMAT);
		} catch (IOException e) {
			System.err.println("Problem while reading templates alignment file "+tempsAlnFile+". Error: "+e.getMessage());
			System.exit(1);
		} catch (FileFormatError e) {
			System.err.println("Problem while reading templates alignment file "+tempsAlnFile+". Error: "+e.getMessage());
			System.exit(1);			
		} catch (AlignmentConstructionError e) {
			System.err.println("Problem while reading templates alignment file "+tempsAlnFile+". Error: "+e.getMessage());
			System.exit(1);			
		}
		TargetToTemplatesAlignment ttAln = new TargetToTemplatesAlignment(seq,templatesAln,baseName,outDir);
		try {
			Alignment al = ttAln.run();
			writeTargetToTemplAl(outDir, baseName, al);
		} catch (TcoffeeError e) {
			System.err.println("Problem while performing target to templates alignment. Error: "+e.getMessage());
			System.exit(1);
		} catch (IOException e) {
			System.err.println("Problem while performing target to templates alignment. Error: "+e.getMessage());
			System.exit(1);			
		}

	}
}

package owl.casp.pipeline;

import gnu.getopt.Getopt;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;

import owl.core.structure.Alignment;
import owl.core.structure.PaulStructAligner;
import owl.core.structure.PdbLoadError;
import owl.core.structure.StructAlignmentError;
import owl.core.structure.TemplateList;
import owl.core.util.MySQLConnection;


public class TemplatesAlignment {
	
	private static final String PROGRAM_NAME = "templatesAlignment";
	
	protected static final File PAUL_PROG = new File("/project/StruPPi/bin/paul");
	protected static final String DEFAULT_PAULMODE = PaulStructAligner.PAULMODE_SLOW;
	protected static final String ALLPAUL_MODES_STR = PaulStructAligner.PAULMODE_VERYFAST+
													", "+PaulStructAligner.PAULMODE_FAST+
													", "+PaulStructAligner.PAULMODE_SLOW+
													", "+PaulStructAligner.PAULMODE_VERYSLOW;
	
	private static final String MYSQL_SERVER = "talyn";
	private static final String PDBASE_DB = "pdbase";
	private static final String CONTACT_TYPE = "Cb";
	private static final double CUTOFF = 8;
	
	protected static final String USAGE = 
	"  [-p]:  paul mode, one of: "+ALLPAUL_MODES_STR+". Default: "+DEFAULT_PAULMODE+"\n";
	protected static final String GETOPT_OPTIONS = "p:";
	
	// members
	private TemplateList templates;
	private String baseName;
	private File outDir;
	
	private String paulMode;
	
	
	public TemplatesAlignment(TemplateList templates, String paulMode, String baseName, File outDir) {
		this.templates = templates;
		this.baseName = baseName;
		this.outDir = outDir;
		this.paulMode = paulMode;
	}
	
	/**
	 * Runs STEP 2 of our homology modelling pipeline, parameters must be passed 
	 * in the constructor.
	 * Will perform a structural alignment of the templates and return and alignment
	 * writing it also to a file, so that it can be picked up by STEP 3 
	 * 
	 * @throws IOException if I/O problems while running the structural aligner or writing
	 *  final alignment 
	 * @throws StructAlignmentError if the structural aligner fails to run
	 */	
	public Alignment run() throws IOException, StructAlignmentError {
		
		File paulLogFile = new File(outDir,baseName+".paul.log");
		File alOutFile = new File(outDir,baseName+".templates_aln.fasta"); 
		
		PaulStructAligner paulSA = new PaulStructAligner(PAUL_PROG, paulMode, paulLogFile);
		
		Alignment al = paulSA.alignStructures(templates);
		PrintStream alOut = new PrintStream(alOutFile);
		al.writeFasta(alOut, 80, true);
		
		return al; 
	}
	
	public static void main(String[] args) {
		
		File templatesFile = null;
		String baseName = null;
		File outDir = new File(".");
		String paulMode = DEFAULT_PAULMODE;
		

		String help = "Usage: \n" +
		PROGRAM_NAME+"\n" +
		"   -i :  file with input templates list \n"+
		"  [-b]:  basename for output files. Default: basename of input templates file \n"+
		"  [-o]:  output dir, where output files will be written. Default: current dir \n"+
		USAGE+"\n";

		Getopt g = new Getopt(PROGRAM_NAME, args, "i:b:o:"+GETOPT_OPTIONS+"h?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'i':
				templatesFile = new File(g.getOptarg());
				break;
			case 'b':
				baseName = g.getOptarg();
				break;				
			case 'o':
				outDir = new File(g.getOptarg());
				break;
			case 'p':
				paulMode = g.getOptarg();
				break;				
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		if (templatesFile == null) {
			System.err.println("Must specify at least a template file (-i)");
			System.err.println(help);
			System.exit(1);
		}
		if (baseName==null) {
			baseName = templatesFile.getName().substring(0, templatesFile.getName().lastIndexOf("."));
		}
		
		
		MySQLConnection conn = null;
		try {
			conn = new MySQLConnection(MYSQL_SERVER,"");
		} catch(SQLException e) {
			System.err.println("Problems while connecting to MySQL server. Error "+e.getMessage());
			System.exit(1);
		}

		TemplateList templates = null;
		try {
			templates = new TemplateList(templatesFile);
			templates.loadPdbData(conn, PDBASE_DB);
			templates.loadRIGraphs(CONTACT_TYPE, CUTOFF);
		} catch (IOException e) {
			System.err.println("Couldn't read templates file "+templatesFile+". Error: "+e.getMessage());
			System.exit(1);
		} catch (PdbLoadError e) {
			System.err.println("Problems getting PDB data of templates from file "+templatesFile+". Error: "+e.getMessage());
			System.exit(1);
		} catch (SQLException e) {
			System.err.println("Problems getting PDB data of templates from file "+templatesFile+". Error: "+e.getMessage());
			System.exit(1);
		}
	
		
		TemplatesAlignment ta = new TemplatesAlignment(templates, paulMode, baseName, outDir);
		try {
			ta.run();
		} catch (StructAlignmentError e) {
			System.err.println("Problems while performing structural alignment. Error: "+e.getMessage());
			System.exit(1);
		} catch (IOException e) {
			System.err.println("Problems while performing structural alignment. Error: "+e.getMessage());
			System.exit(1);
		}
		
	}
	
}

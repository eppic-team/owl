package owl.casp.pipeline;

import gnu.getopt.Getopt;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.xml.sax.SAXException;

import owl.core.connections.GTGHitList;
import owl.core.connections.GTGParser;
import owl.core.runners.PsipredError;
import owl.core.runners.PsipredRunner;
import owl.core.runners.blast.BlastError;
import owl.core.runners.blast.BlastHitList;
import owl.core.runners.blast.BlastRunner;
import owl.core.runners.blast.BlastUtils;
import owl.core.runners.blast.BlastXMLParser;
import owl.core.sequence.Sequence;
import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadError;
import owl.core.structure.Template;
import owl.core.structure.TemplateList;
import owl.core.util.FileFormatError;
import owl.core.util.MySQLConnection;



/**
 * First step of our homology modelling pipeline
 * 
 * This can be used both as an executable to only run step 1 or as part of an 
 * full automated pipeline, by using the run() method    
 * 
 *
 */
public class TemplateSelection {

	// constants
	private static final String PROGRAM_NAME = "templateSelection";
	
	protected static final String BLASTBIN_DIR = "/project/StruPPi/bin";
	protected static final String PSIPRED_HOMEDIR = "/project/StruPPi/Software/psipred-2.6.1";
	protected static final String BLASTDB_BCKGRD_DIR = "/scratch/local/blast";
	protected static final String BLASTDB_PDB_DIR = "/project/StruPPi/CASP8/blast_dbs";
	
	protected static final String NR_BLASTDB = "uniref90.filtall";
	protected static final String DEFAULT_PDB_BLASTDB = "seqs_pdbase_20080903.fix.reps.fa";
	protected static final int 	DEFAULT_MAXITER = 5;
	protected static final double EVALUE_CUTOFF_PREFILTER = 1e-5;
	protected static final int    TOO_MANY_HITS = 5;
	protected static final double FULLMATCH_IDENTITY_THRESHOLD = 99.9;
	protected static final double FULLMATCH_COVERAGE = 0.95;
	
	protected static final String PSIBLASTOUT_NR_SUFFIX = ".nr-psiblast.out";
	protected static final String PSIBLASTOUT_PDB_SUFFIX = ".pdb-psiblast.out";
	protected static final String PSIBLASTOUT_CLASSIC_PDB_SUFFIX = ".pdb-psiblast.classic.out";
	protected static final String PSIBLAST_PROFILE_SUFFIX = ".profile.chk";
	protected static final String BLASTOUT_PDB_SUFFIX = ".pdb-blast.out"; 
	protected static final String BLASTOUT_CLASSIC_PDB_SUFFIX = ".pdb-blast.classic.out";
	protected static final String PSIPRED_SS2_SUFFIX = ".ss2";
	protected static final String PSIPRED_HORIZ_SUFFIX = ".horiz";
	
	protected static final int 	BLAST_DEFAULT_OUTPUT_TYPE = 7;  // xml output
	protected static final int	BLAST_CLASSIC_OUTPUT = 0;
	protected static final boolean BLAST_NO_FILTERING = true;
	protected static final int	BLAST_PROFILE_OUTPUT = 6;	// needed for conversion to MSA with reformatpsialn
	
	protected static final double DEFAULT_SIMILARITY_GRAPH_RMSD_CUTOFF = 4.0;
	
	protected static final String GTG_RESULTS_DIR = "/project/StruPPi/CASP8/gtg";
	protected static final int    DEFAULT_GTG_SCORE_CUTOFF = 1000;
	
	private static final int	MAX_TITLE_LEN = 60;		// length of title column in report table output on screen
	
	// the 3 possible values for the -l parameter
	protected static final String USE_BLAST_TEMPLATES = "B";
	protected static final String USE_PSIBLAST_TEMPLATES = "P";
	protected static final String USE_GTG_TEMPLATES = "G";
	
	// db constants
	protected static final String PDBASE_DB    = "pdbase";

	// usage string and getopt options
	protected static final String USAGE = 
	"  [-x]:  max number of hits to be taken (after evalue/score cutoff). Default: "+TOO_MANY_HITS+"\n"+
	"  [-m]:  similarity graph rmsd cutoff. Default: "+DEFAULT_SIMILARITY_GRAPH_RMSD_CUTOFF+"\n"+
	"  [-j]:  maximum number of iterations for psi-blast run. Default: "+DEFAULT_MAXITER+"\n"+
	"  [-e]:  e-value cutoff for selecting templates from blast/psiblast. Default: "+EVALUE_CUTOFF_PREFILTER+"\n" +
	"  [-a]:  number of CPUs to be used by blast. Default: 1\n"+
	"  [-B]:  non-redundat sequence blast database in "+BLASTDB_BCKGRD_DIR+". \n" +
	"         Default: "+NR_BLASTDB+"\n"+
	"  [-P]:  PDB blast database in "+BLASTDB_PDB_DIR+".\n" +
	"         Default: "+DEFAULT_PDB_BLASTDB+"\n" +
	"  [-K]:  skip the psi-blast step (psipred requires psiblast, if -K specified no \n" +
	"         psipred will be run). Default: psi-blast runs\n" +
	"  [-S]:  skip the psi-pred secondary structure prediction. Default: psipred runs\n" +
	"  [-g]:  GTG score cutoff. Default: "+ DEFAULT_GTG_SCORE_CUTOFF + "\n" +
	"  [-G]:  GTG directory containing the GTG output files " + GTG_RESULTS_DIR + "\n"+
	"  [-T]:  skip the GTG output parsing step. Default: GTG parsing runs\n" +
	"  [-l]:  specify B for blast, P for psi-blast, G for GTG for the final selected \n" +
	"         templates list. Default: final template list will be the first non-empty \n" +
	"         list of blast, psi-blast or GTG  \n";
	protected static final String GETOPT_OPTIONS = "x:m:j:B:P:e:a:KSg:G:Tl:";
	
	// members
	private File inputSeqFile;
	//private BlastRunner blastRunner;
	private String pdbBlastDb;
	private String nrBlastDb;
	private int blastNumThreads;
	private int maxIter;
	private int maxHits;
	private double eValueCutoff;
	private int gtgScoreCutoff;
	private boolean psiblast;
	private boolean psipred;
	private boolean gtg;
	private File outProfileFile;
	
	private double similarityGraphRmsdCutoff;
	
	private String baseName;
	private File outDir;
	private File gtgDir;
	private String gtgDateFilter;
	
	private String selectTemplates;
	
	private MySQLConnection conn;
	
	private Sequence inputSequence;
	private File outPdbPsiblast;
	private File outPdbBlast;
	private BlastHitList hitsPdbPsiblast;
	private BlastHitList hitsPdbBlast;
	private boolean psipredSucceeded;
	private boolean fullMatch;
	private TemplateList templatesBlast;
	private TemplateList templatesPsiBlast;
	private TemplateList templatesGTG;
	
	/**
	 * Constructs a TemplateSelection by passing all the necessary parameters
	 * @param inputSeqFile input sequence file
	 * @param pdbBlastDb the blast db containing non-redundant PDB sequences
	 * @param nrBlastDb the blast db containing a non-redundant sequence database
	 * @param blastNumThreads the number of CPUs to be used by blast (blast's -a option)
	 * @param maxIter maximum number of iterations for the psiblast step
	 * @param eValueCutoff evalue cutoff for considering a blast hit
	 * @param gtgScoreCutoff
	 * @param maxHits
	 * @param similarityGraphRmsdCutoff
	 * @param psiblast whether to run the psi-blast step or not
	 * @param psipred whether to run the psipred step or not
	 * @param gtg whether to run GTG parsing or not
	 * @param baseName base name of the output files
	 * @param outDir output directory
	 * @param gtgDir the directory with the GTG output files 
	 * @param selectTemplates which TemplateList will be the finally selected, 
	 * valid values: {@link #USE_BLAST_TEMPLATES}, {@link #USE_PSIBLAST_TEMPLATES}, {@link #USE_GTG_TEMPLATES}
	 * @param conn a connection to the database server
	 */
	public TemplateSelection(File inputSeqFile, 
			String pdbBlastDb, String nrBlastDb, int blastNumThreads, int maxIter, double eValueCutoff, int gtgScoreCutoff, 
			int maxHits, double similarityGraphRmsdCutoff, 
			boolean psiblast, boolean psipred, boolean gtg, String baseName, File outDir, File gtgDir, 
			String selectTemplates,
			MySQLConnection conn) {
		
		this.inputSeqFile = inputSeqFile;
		//blastRunner = new BlastRunner(BLASTBIN_DIR, BLASTDB_DIR);
		this.pdbBlastDb = pdbBlastDb;
		Pattern p = Pattern.compile(".*[^0-9](20\\d\\d\\d\\d\\d\\d)[^0-9].*");
		Matcher m = p.matcher(pdbBlastDb);
		this.gtgDateFilter = new SimpleDateFormat("yyyymmdd").format(new Date());
		if (m.matches()) {
			this.gtgDateFilter = m.group(1);
		}
		this.nrBlastDb = nrBlastDb;
		this.blastNumThreads = blastNumThreads;
		this.baseName = baseName;
		this.outDir = outDir;
		this.maxIter = maxIter;
		this.eValueCutoff = eValueCutoff;
		this.gtgScoreCutoff = gtgScoreCutoff;
		this.maxHits = maxHits;
		this.similarityGraphRmsdCutoff = similarityGraphRmsdCutoff;
		this.psiblast = psiblast;
		this.psipred = psipred;
		this.gtg = gtg;
		this.gtgDir = gtgDir;
		this.selectTemplates = selectTemplates;
		this.conn = conn;
		
		this.psipredSucceeded = true;
		this.fullMatch = false;
		
		this.templatesBlast = null;
		this.templatesPsiBlast = null;
		this.templatesGTG = null;
	}
	
	/**
	 * Blast against PDB 
	 * @return
	 * @param queryLength
	 * @throws BlastError
	 * @throws IOException
	 */
	private BlastHitList doBlastAgainstPdb() throws BlastError, IOException {
		System.out.println("### BLAST");
		System.out.println("Blasting against PDB sequence db "+pdbBlastDb + "...");

		outPdbBlast = new File(outDir,baseName+BLASTOUT_PDB_SUFFIX);
		File classicOutFile = new File(outDir,baseName+BLASTOUT_CLASSIC_PDB_SUFFIX);
		BlastRunner blastRunner = new BlastRunner(BLASTBIN_DIR, BLASTDB_PDB_DIR);
		blastRunner.runBlastp(inputSeqFile, pdbBlastDb, outPdbBlast, BLAST_DEFAULT_OUTPUT_TYPE, BLAST_NO_FILTERING, blastNumThreads);
		blastRunner.runBlastp(inputSeqFile, pdbBlastDb, classicOutFile, BLAST_CLASSIC_OUTPUT, BLAST_NO_FILTERING, blastNumThreads);
		
		try {
			BlastXMLParser blastParser = new BlastXMLParser(outPdbBlast);
			hitsPdbBlast = blastParser.getHits();
		} catch (SAXException e) {
			// if this happens it means that blast doesn't format correctly its XML, i.e. has a bug
			System.err.println("Unexpected error: "+e.getMessage());
		}
		
		if (hitsPdbBlast.size()>0) System.out.println("Blast done, best e-value was "+hitsPdbBlast.getBestHit().getEvalueMaxScoringHsp());
		System.out.println("Number of hits: "+hitsPdbBlast.size());
		hitsPdbBlast.applyCutoff(eValueCutoff);
		if (hitsPdbBlast.size()>=maxHits) {
			System.out.println("More than "+maxHits+" hits. Taking only "+maxHits);
			hitsPdbBlast.filterByMaxRank(maxHits);
//			System.out.print("Trying lower cutoffs: ");
//			while (hits.size()>=maxHits) {
//				eValueCutoff *= 1e-1;
//				hits.applyCutoff(eValueCutoff);
//				System.out.printf("cutoff: %2.1e, hits: %s; ",eValueCutoff,hits.size());
//			}
//			System.out.println();
		}
		System.out.println("Number of hits after cutoff: "+hitsPdbBlast.size());
		hitsPdbBlast.print();

		return hitsPdbBlast;
	}
	
	/**
	 * Psiblast against nr sequence database and then with the profile against PDB
	 * @param queryLength
	 * @return
	 * @throws BlastError
	 * @throws IOException
	 */
	private BlastHitList doPsiblastAgainstNrAndPdb() throws BlastError, IOException {
		System.out.println("### PSI-BLAST");
		System.out.println("Psi-blasting with maximum of "+maxIter+" iterations...");
		
		File outNrPsiblast = new File(outDir,baseName+PSIBLASTOUT_NR_SUFFIX);
		outPdbPsiblast = new File(outDir,baseName+PSIBLASTOUT_PDB_SUFFIX);
		File classicOutPdbPsiblast = new File(outDir, baseName+PSIBLASTOUT_CLASSIC_PDB_SUFFIX);
		outProfileFile = new File(outDir,baseName+PSIBLAST_PROFILE_SUFFIX);

		BlastRunner blastRunner = new BlastRunner(BLASTBIN_DIR, BLASTDB_BCKGRD_DIR);
		blastRunner.runPsiBlast(inputSeqFile, nrBlastDb, outNrPsiblast, maxIter, outProfileFile, null, BLAST_PROFILE_OUTPUT, BLAST_NO_FILTERING, blastNumThreads);
		blastRunner = new BlastRunner(BLASTBIN_DIR, BLASTDB_PDB_DIR);
		blastRunner.runPsiBlast(inputSeqFile, pdbBlastDb, outPdbPsiblast, 1, null, outProfileFile, BLAST_DEFAULT_OUTPUT_TYPE, BLAST_NO_FILTERING, blastNumThreads);
		blastRunner.runPsiBlast(inputSeqFile, pdbBlastDb, classicOutPdbPsiblast, 1, null, outProfileFile, BLAST_CLASSIC_OUTPUT, BLAST_NO_FILTERING, blastNumThreads);

		try {
			BlastXMLParser blastParser = new BlastXMLParser(outPdbPsiblast);
			hitsPdbPsiblast = blastParser.getHits();
		} catch (SAXException e) {
			// if this happens it means that blast doesn't format correctly its XML, i.e. has a bug
			System.err.println("Unexpected error: "+e.getMessage());
		}
		
		if (hitsPdbPsiblast.size()>0) System.out.println("Psi-blast done, best e-value was "+hitsPdbPsiblast.getBestHit().getEvalueMaxScoringHsp());
		System.out.println("Number of hits: "+hitsPdbPsiblast.size());
		hitsPdbPsiblast.applyCutoff(eValueCutoff);
		if (hitsPdbPsiblast.size()>=maxHits) {
			System.out.println("More than "+maxHits+" hits. Taking only "+maxHits);
			hitsPdbPsiblast.filterByMaxRank(maxHits);
//			System.out.print("Trying lower cutoffs: ");
//			while (hits.size()>=maxHits) {
//				eValueCutoff *= 1e-1;
//				hits.applyCutoff(eValueCutoff);
//				System.out.printf("cutoff: %2.1e, hits: %s; ",eValueCutoff,hits.size());
//			}
//			System.out.println();
		}
		System.out.println("Number of hits after cutoff: "+hitsPdbPsiblast.size());
		hitsPdbPsiblast.print();
		
		
		return hitsPdbPsiblast;

	}
	
	/**
	 * Runs psipred
	 * @throws PsipredError
	 * @throws IOException
	 */
	private void doPsipred() throws PsipredError, IOException {
		System.out.println("Running psipred for secondary structure");
		PsipredRunner psipred = new PsipredRunner(PSIPRED_HOMEDIR);
		File outSs2File   = new File(outDir, baseName+PSIPRED_SS2_SUFFIX);
		File outHorizFile = new File(outDir, baseName+PSIPRED_HORIZ_SUFFIX);
		psipred.run(inputSeqFile, outSs2File, outHorizFile, outProfileFile, BLASTBIN_DIR);
		System.out.println("Psipred output in files "+outSs2File+" and "+outHorizFile);
	}

	/**
	 * Writes cluster graph files (a gdl file and a matrix file) for the given templates
	 * @param templates
	 * @param suffix the name appended as an extension to the output file names 
	 */
	private void doWriteClusterGraphFile(TemplateList templates, String suffix) {
		// writing cluster graph file
		try {
			BlastUtils.writeClusterGraph(templates, conn, PDBASE_DB, new File(outDir,baseName+"."+suffix+".gdl"), new File(outDir,baseName+"."+suffix+".matrix"), similarityGraphRmsdCutoff);
		} catch (IOException e){
			System.err.println("Error while getting rmsds  for cluster graph file: "+e.getMessage());
		} catch (PdbLoadError e) {
			System.err.println("Error while getting PDB data for cluster graph file: "+e.getMessage());
		} catch (SQLException e) {
			System.err.println("Error while getting PDB data for cluster graph file: "+e.getMessage());	
		}
	}
	
	/**
	 * Gets the GTG hits by reading the precomputed output file
	 * @param queryLength
	 * @return the GTGHitList or an empty list if problems occur while parsing
	 * @throws IOException
	 */
	private GTGHitList doGetGTGHits(int queryLength) throws IOException {
		Pattern p = Pattern.compile("T\\d\\d\\d\\d");
		Matcher m = p.matcher(baseName);
		if (!m.matches()) {
			System.err.println("Specified basename "+baseName+" doesn't look like a CASP target. Skipping GTG.");
			return new GTGHitList(); // basename doesn't look like a CASP target we return empty list 
		}

		File gtgOut = new File(gtgDir, baseName+".gtg");

		System.out.println("### GTG");
		System.out.println("Getting GTG hits from file "+gtgOut);
		
		GTGHitList hits = null;
		try {
			GTGParser gtgParser = new GTGParser(gtgOut);
			hits = gtgParser.getHits();		
		} catch (IOException e) {
			System.err.println("Error while parsing GTG file: "+e.getMessage());
			return new GTGHitList();
		} catch (FileFormatError e) {
			System.err.println("Error while parsing GTG file: "+e.getMessage());
			return new GTGHitList();			
		}
		
		int hitsBeforeFilter = hits.size();
		System.out.println("Filtering by date "+gtgDateFilter);
		try {
			hits.filterByMaxReleaseDate(conn, PDBASE_DB, gtgDateFilter);
		} catch (SQLException e) {
			System.err.println("Error while filtering out the GTG hits by date: "+gtgDateFilter+". Error: "+e.getMessage());
		}
		System.out.println((hitsBeforeFilter-hits.size())+" hits filtered out.");
		
		if (hits.size()>0) System.out.println("Best GTG score: "+hits.getBestHit().getTotalScore());
		System.out.println("Number of hits: "+hits.size());

		hits.setQueryLength(queryLength);
		hits.applyCutoff(gtgScoreCutoff);

		System.out.println("Number of hits after cutoff "+gtgScoreCutoff+": "+hits.size());
		
		if (hits.size()>=maxHits) {
			System.out.println("More than "+maxHits+" hits. Taking only "+maxHits);
			hits.filterByMaxRank(maxHits);
//			System.out.print("Trying higher cutoffs: ");
//			while (hits.size()>=maxHits) {
//				gtgScoreCutoff += 1000;
//				hits.applyCutoff(gtgScoreCutoff);
//				System.out.printf("cutoff: %s, hits: %s ; ",gtgScoreCutoff,hits.size());
//			}
//			System.out.println();
		}

		// before printing we reassign the residue numbers of the subject, 
		// otherwise the graphical alignment print out would be refering to the observed sequence
		try {
			hits.reassignSubjectSerials(conn, PDBASE_DB);
		} catch (SQLException e) {
			System.err.println("Error while reassigning serials for GTG matches: "+e.getMessage());
		} catch (PdbCodeNotFoundException e) {
			System.err.println("Error while reassigning serials for GTG matches: "+e.getMessage());
		} catch (PdbLoadError e) {
			System.err.println("Error while reassigning serials for GTG matches: "+e.getMessage());
		}
		
		hits.print();
		
		return hits;

	}
	
	/**
	 * Inner class to store a record of the summary table
	 */
	private class SummaryTableRecord {
		String id;
		int rankBlast;
		double eValBlast;
		int rankPsiBlast;
		double eValPsiBlast;
		int rankGTG;
		int scoreGTG;
		String scopSccs;
		String title;
		public SummaryTableRecord(String id, int rankBlast, double eValBlast, int rankPsiBlast, double eValPsiBlast, int rankGTG, int scoreGTG, String scopSccs, String title) {
			this.id = id;
			this.rankBlast = rankBlast;
			this.eValBlast = eValBlast;
			this.rankPsiBlast = rankPsiBlast;
			this.eValPsiBlast = eValPsiBlast;
			this.rankGTG = rankGTG;
			this.scoreGTG = scoreGTG;
			this.scopSccs = scopSccs;
			this.title = title;
		}		
	}

	/**
	 * Inner class to store the summary table 
	 */
	private class SummaryTable extends ArrayList<SummaryTableRecord> {
		private static final long serialVersionUID = 1L;

		public void sortOnBlast() {
			Collections.sort(this, new Comparator<SummaryTableRecord>() {
				public int compare(SummaryTableRecord o1, SummaryTableRecord o2) {
					return new Integer(o1.rankBlast).compareTo(o2.rankBlast);
				}
			});
		}
		
		public void sortOnPsiBlast() {
			Collections.sort(this, new Comparator<SummaryTableRecord>() {
				public int compare(SummaryTableRecord o1, SummaryTableRecord o2) {
					return new Integer(o1.rankPsiBlast).compareTo(o2.rankPsiBlast);
				}
			});
		}

		public void sortOnGTG() {
			Collections.sort(this, new Comparator<SummaryTableRecord>() {
				public int compare(SummaryTableRecord o1, SummaryTableRecord o2) {
					return new Integer(o1.rankGTG).compareTo(o2.rankGTG);
				}
			});
		}
		
		/**
		 * Writes the table to an output stream. The last column 'title' will be restricted to the given length.
		 * Negative maxLen means no length restriction.
		 * @param Out
		 * @param maxLen
		 */
		public void writeTable(PrintStream Out, int maxLen) {
			String titleStr;
			Out.printf("%5s\t%9s\t%8s\t%9s\t%8s\t%9s\t%8s\t%30s\t%s\n","id","blast", "eVal", "psiblast", "eVal", "GTG", "score", "scop id", "title");
			// setting the title
			for (SummaryTableRecord record:this) {
				if(record.title == null) {
					titleStr = "";
				} else {
					if(maxLen < 0) {
						titleStr = record.title;
					} else {
						titleStr = record.title.length() <=MAX_TITLE_LEN?record.title:record.title.substring(0, MAX_TITLE_LEN-3) + "...";
					}
				}
				Out.printf("%5s\t%9s\t%8s\t%9s\t%8s\t%9s\t%8s\t%30s\t%s\n",
						record.id,
						record.rankBlast==Integer.MAX_VALUE?"-":record.rankBlast, 
						record.eValBlast==Double.MAX_VALUE?"-":String.format("%8.1e",record.eValBlast),
						record.rankPsiBlast==Integer.MAX_VALUE?"-":record.rankPsiBlast,
						record.eValPsiBlast==Double.MAX_VALUE?"-":String.format("%8.1e",record.eValPsiBlast),
						record.rankGTG==Integer.MAX_VALUE?"-":record.rankGTG,
						record.scoreGTG==0?"-":record.scoreGTG,
						record.scopSccs,
						titleStr);
			}
		}
		
		public void print() {
			writeTable(System.out, MAX_TITLE_LEN);
		}
	}
	
	/**
	 * Prints a tabular summary of blast, psiblast and GTG hits, including scop 
	 * assignments and writes the table also to a basename.report file
	 * The table is sorted on the ranks of the first non-empty list (blast, psiblast or gtg)
	 * @param templatesBlast
	 * @param templatesPsiBlast
	 * @param templatesGTG
	 * @throws FileNotFoundException if can't find report file to write to it
	 */
	private void doTabularSummary(TemplateList templatesBlast, TemplateList templatesPsiBlast, TemplateList templatesGTG) throws FileNotFoundException {
		TreeSet<String> ids = getTemplateListUnion(templatesBlast, templatesPsiBlast, templatesGTG);

		boolean pdbDataAvailable = true;
		try {
			if (templatesBlast!=null) templatesBlast.loadPdbData(conn, PDBASE_DB);
			if (templatesPsiBlast!= null) templatesPsiBlast.loadPdbData(conn, PDBASE_DB);
			if (templatesGTG!=null) templatesGTG.loadPdbData(conn, PDBASE_DB);
		} catch (SQLException e) {
			System.err.println("Couldn't get SCOP identifiers. Error "+e.getMessage());
			pdbDataAvailable = false;
		} catch (PdbLoadError e) {
			System.err.println("Couldn't get SCOP identifiers. Error "+e.getMessage());
			pdbDataAvailable = false;
		}
		
		SummaryTable table = new SummaryTable();
		for (String id:ids) {
			String scop = "";
			String title = "";
			if (pdbDataAvailable) {
				if (templatesBlast.contains(id)) {
					Template template = templatesBlast.getTemplate(id); 
					scop = template.getScopSccsString();
					title = template.getTitle();
				} else if (templatesPsiBlast!=null && templatesPsiBlast.contains(id)) {
					Template template = templatesPsiBlast.getTemplate(id); 
					scop = template.getScopSccsString();
					title = template.getTitle();
				} else if (templatesGTG!=null && templatesGTG.contains(id)) {
					Template template = templatesGTG.getTemplate(id); 
					scop = template.getScopSccsString();
					title = template.getTitle();
				}
			}
			
			double eValBlast = Double.MAX_VALUE;
			double eValPsiBlast = Double.MAX_VALUE;
			int scoreGTG = 0;
			if (templatesBlast.contains(id)) {
				Template template = templatesBlast.getTemplate(id); 
				eValBlast = template.getBlastHit().getEvalueMaxScoringHsp();
			} 
			if (templatesPsiBlast!=null && templatesPsiBlast.contains(id)) {
				Template template = templatesPsiBlast.getTemplate(id); 
				eValPsiBlast = template.getBlastHit().getEvalueMaxScoringHsp();
			} 
			if (templatesGTG!=null && templatesGTG.contains(id)) {
				Template template = templatesGTG.getTemplate(id); 
				scoreGTG = template.getGTGHit().getTotalScore();	
			}			
			
			int rankBlast = templatesBlast.getRank(id);
			int rankPsiBlast = Integer.MAX_VALUE;
			if (templatesPsiBlast!=null) rankPsiBlast = templatesPsiBlast.getRank(id);
			int rankGTG = Integer.MAX_VALUE;
			if (templatesGTG!=null) rankGTG = templatesGTG.getRank(id);
			 
			table.add(new SummaryTableRecord(id,rankBlast, eValBlast, rankPsiBlast, eValPsiBlast, rankGTG, scoreGTG, scop, title));
		}
		if (templatesBlast!=null && templatesBlast.size()>0) {
			table.sortOnBlast();
		} else if (templatesPsiBlast!=null && templatesPsiBlast.size()>0) {
			table.sortOnPsiBlast();
		} else if (templatesGTG!=null && templatesGTG.size()>0) {
			table.sortOnGTG();
		}
		table.print();
		
		File reportFile = new File(outDir, baseName+".report");
		PrintStream reportPS = new PrintStream(reportFile);
		reportPS.printf("#maxHits: %2d, eValCutoff: %6.0e, gtgCutoff: %5d\n", maxHits, eValueCutoff, gtgScoreCutoff);
		table.writeTable(reportPS, -1);	// -1 meaning unrestricted output length
		reportPS.close();
		System.out.println("Summary table wrote to "+reportFile);
	}
	
	/**
	 * Gets the union of all templateIds of the 3 given TemplateList
	 * @param templatesBlast
	 * @param templatesPsiBlast
	 * @param templatesGTG
	 * @return
	 */
	private TreeSet<String> getTemplateListUnion(TemplateList templatesBlast, TemplateList templatesPsiBlast, TemplateList templatesGTG) {
		TreeSet<String> ids = new TreeSet<String>(); 
		for (Template template: templatesBlast) {
			ids.add(template.getId());
		}
		if (templatesPsiBlast!=null) {
			for (Template template: templatesPsiBlast) {
				ids.add(template.getId());
			}			
		}
		if (templatesGTG!=null) {
			for (Template template: templatesGTG) {
				ids.add(template.getId());
			}			
		}		
		return ids;
	}
	
	private TemplateList doTemplateSelection(TemplateList templatesBlast, TemplateList templatesPsiBlast, TemplateList templatesGTG) {
		// TODO this is just a manual selection process, ideally we want a full automated (and clever!) selection
		
		// if not specific selection criterium specified, then try first blast, then psiblast and finally GTG
		if (selectTemplates.equals("")) {
			if (templatesBlast!=null && templatesBlast.size()>0) {
				// we want to catch the special case of a full match, we don't want modelling for that case!
				Template firstTemp = templatesBlast.iterator().next();
				if (firstTemp.getBlastHit().getTotalPercentIdentity()>FULLMATCH_IDENTITY_THRESHOLD &&
						(firstTemp.getBlastHit().getQueryCoverage()>FULLMATCH_COVERAGE)) {
					
					fullMatch = true;
					System.out.println("###########################################################");
					System.out.println("An IDENTICAL MATCH TEMPLATE (100% identity, above "+FULLMATCH_COVERAGE+" coverage) was found: "+firstTemp.getId());
					System.out.println("###########################################################");
				}
				return templatesBlast;
			} else if (templatesPsiBlast!=null && templatesPsiBlast.size()>0) {
				return templatesPsiBlast;
			} else if (templatesGTG!=null && templatesGTG.size()>0) {
				return templatesGTG;
			}
		}
		
		// a selection criterium was specified, we return the corresponding templates
		if (selectTemplates.equals(USE_BLAST_TEMPLATES)) {
			return templatesBlast;
		}
		if (selectTemplates.equals(USE_PSIBLAST_TEMPLATES)) {
			return templatesPsiBlast;
		}
		if (selectTemplates.equals(USE_GTG_TEMPLATES)) {
			return templatesGTG;
		}
		return new TemplateList(); // if all other failed, we return an empty list
	}
	
	/**
	 * Runs STEP 1 of our homology modelling pipeline, parameters must be passed 
	 * in the constructor.
	 * As well as returning a list of selected templates it writes out to a file a 
	 * list with the selected template's ids, so that the list can be used to be 
	 * picked up by STEP 2. 
	 * A psipred secondary structure prediction is also performed
	 * A few other files are also written: a gdl file with a cluster graph of the 
	 * selected templates and several blast output files
	 * @return the selected templates
	 * @throws BlastError if blast fails to run
	 * @throws IOException if I/O problems while running blast or psipred or writing results
	 * @throws FileFormatError if input sequence file not in FASTA format
	 */
	public TemplateList run() throws BlastError, IOException, FileFormatError {
		// getting query sequence for the length
		inputSequence = new Sequence();
		inputSequence.readFromFastaFile(inputSeqFile);

		// blasting
		BlastHitList hitsBlast = this.doBlastAgainstPdb();
		templatesBlast = new TemplateList(hitsBlast);
		templatesBlast.setSource(TemplateList.SRC_BLAST);
		File idsFile = new File(outDir,baseName+".blast.templates");
		if (hitsBlast.size()>0) {
			System.out.println("Writing blast template ids to file "+idsFile);
			templatesBlast.writeIdsToFile(idsFile);
			doWriteClusterGraphFile(templatesBlast, "blast");
		}
		
		// psiblasting
		BlastHitList hitsPsiBlast = null;
		templatesPsiBlast = null;
		if (psiblast) {			
			hitsPsiBlast = this.doPsiblastAgainstNrAndPdb();
			templatesPsiBlast = new TemplateList(hitsPsiBlast);
			templatesPsiBlast.setSource(TemplateList.SRC_PSIBLAST);
			idsFile = new File(outDir,baseName+".psiblast.templates");
			if (hitsPsiBlast.size()>0) {
				System.out.println("Writing psiblast template ids to file "+idsFile);
				templatesPsiBlast.writeIdsToFile(idsFile);
				doWriteClusterGraphFile(templatesPsiBlast, "psiblast");
			}
		}
		
		// gtg
		GTGHitList hitsGTG = null;
		templatesGTG = null;
		if (gtg) {
			hitsGTG = doGetGTGHits(inputSequence.getLength());
			templatesGTG = new TemplateList(hitsGTG);
			templatesGTG.setSource(TemplateList.SRC_GTG);
			if (hitsGTG.size()>0) {
				idsFile = new File(outDir,baseName+".GTG.templates");
				System.out.println("Writing GTG template ids to file "+idsFile);
				templatesGTG.writeIdsToFile(idsFile);
				doWriteClusterGraphFile(templatesGTG, "GTG");
			}
		}
		
		// running secondary structure prediction with psipred from psi-blast profile
		if (psiblast && psipred) {
			try {
				doPsipred();
			} catch (PsipredError e) {
				System.err.println("Psipred failed to run. Error: "+e.getMessage());
				psipredSucceeded = false;
			}
		}
		
		// selection of TemplateList and writing final selection to file
		TemplateList templates = doTemplateSelection(templatesBlast, templatesPsiBlast, templatesGTG);
		
		idsFile = new File(outDir,baseName+".templates");
		System.out.println("Writing final selected template ids to file "+idsFile);
		templates.writeIdsToFile(idsFile);
		
		// print summary in tabular form
		doTabularSummary(templatesBlast, templatesPsiBlast, templatesGTG);
		
		return templates;
	}
	
	/**
	 * Gets the input sequence read from the input file 
	 * @return
	 */
	public Sequence getInputSequence() {
		return inputSequence;
	}
	
	/**
	 * Gets the output pdb psiblast file (in default format xml)
	 * @return
	 */
	public File getOutPdbPsiblast() {
		return outPdbPsiblast;
	}

	/**
	 * Gets the output pdb blast file (in default format xml)
	 * @return
	 */
	public File getOutPdbBlast() {
		return outPdbBlast;
	}
	
	/**
	 * Returns the blast hit list of the psiblast against pdb
	 * @return
	 */
	public BlastHitList getHitsPdbPsiblast() {
		return hitsPdbPsiblast;
	}
	
	/**
	 * Returns the blast hit list of the blast against pdb
	 * @return
	 */
	public BlastHitList getHitsPdbBlast() {
		return hitsPdbBlast;
	}
	
	/**
	 * Returns true if psipred ran correctly, false otherwise
	 * @return
	 */
	public boolean hasPsipredSucceeded() {
		return psipredSucceeded;
	}

	/**
	 * Returns true if full match template (100% identity and above 95% coverage)
	 * was found.
	 * @return
	 */
	public boolean foundFullMatch() {
		return fullMatch;
	}
	
	public TemplateList getTemplatesBlast() {
		return templatesBlast;
	}
	
	public TemplateList getTemplatesPsiBlast() {
		return templatesPsiBlast;
	}
	
	public TemplateList getTemplatesGTG() {
		return templatesGTG;
	}
	
	/*-------------------------------- main -----------------------------------*/
	
	public static void main(String[] args) throws Exception {

		File inputSeqFile = null;
		String pdbBlastDb = DEFAULT_PDB_BLASTDB;
		String nrBlastDb = NR_BLASTDB;
		int blastNumThreads = 1;
		String baseName = null;
		File outDir = new File(".");
		int maxHits = TOO_MANY_HITS;
		double similarityGraphRmsdCutoff = DEFAULT_SIMILARITY_GRAPH_RMSD_CUTOFF;
		int maxIter = DEFAULT_MAXITER;
		double eValueCutoff = EVALUE_CUTOFF_PREFILTER;
		int gtgScoreCutoff = DEFAULT_GTG_SCORE_CUTOFF;
		boolean psiblast = true;
		boolean psipred = true;
		boolean gtg = true;
		File gtgDir = new File(GTG_RESULTS_DIR);
		String selectTemplates = ""; 

		String help = "Usage: \n" +
		PROGRAM_NAME+"\n" +
		"   -i :  file with input target sequence in FASTA format \n"+
		"  [-b]:  basename for output files. Default: basename of input sequence file \n"+
		"  [-o]:  output dir, where output files will be written. Default: current dir \n"+
		USAGE + "\n";
		
		Getopt g = new Getopt(PROGRAM_NAME, args, "i:b:o:"+GETOPT_OPTIONS+"h?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'i':
				inputSeqFile = new File(g.getOptarg());
				break;
			case 'b':
				baseName = g.getOptarg();
				break;				
			case 'o':
				outDir = new File(g.getOptarg());
				break;
			case 'x':
				maxHits = Integer.parseInt(g.getOptarg());
				break;				
			case 'm':
				 similarityGraphRmsdCutoff = Double.parseDouble(g.getOptarg());
				break;								
			case 'j':
				maxIter = Integer.parseInt(g.getOptarg());
				break;
			case 'B':
				nrBlastDb = g.getOptarg();
				break;				
			case 'P':
				pdbBlastDb = g.getOptarg();
				break;
			case 'e':
				eValueCutoff = Double.parseDouble(g.getOptarg());
				break;
			case 'a':
				blastNumThreads = Integer.parseInt(g.getOptarg());
				break;				
			case 'K':
				psiblast = false;
				break;
			case 'S':
				psipred = false;
				break;
			case 'g':
				gtgScoreCutoff = Integer.parseInt(g.getOptarg());
				break;																
			case 'G':
				gtgDir = new File(g.getOptarg());
				break;
			case 'T':
				gtg = false;
				break;
			case 'l':
				selectTemplates = g.getOptarg();
				if (!selectTemplates.equals(USE_BLAST_TEMPLATES) && 
						!selectTemplates.equals(USE_PSIBLAST_TEMPLATES) && 
						!selectTemplates.equals(USE_GTG_TEMPLATES)) {
					System.err.println("Invalid value specified for -l option, allowed values are: "+USE_BLAST_TEMPLATES+", "+USE_PSIBLAST_TEMPLATES+" or "+USE_GTG_TEMPLATES);
					System.exit(1);
				}
				break;																				
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		if (inputSeqFile==null) {
			System.err.println("Must specify input sequence file (-i)");
			System.err.println(help);
			System.exit(1);
		}
		
		if (baseName==null) {
			baseName = inputSeqFile.getName().substring(0, inputSeqFile.getName().lastIndexOf("."));
		}
		
		if (!psiblast && selectTemplates.equals(USE_PSIBLAST_TEMPLATES)) {
			System.err.println("Can't specify option use psi-blast templates (-l P) if skip psi-blast specified (-K)");
			System.exit(1);
		}
		if (!gtg && selectTemplates.equals(USE_GTG_TEMPLATES)){
			System.err.println("Can't specify option use GTG templates (-l G) if skip GTG specified (-T)");
			System.exit(1);			
		}
		
		MySQLConnection conn = new MySQLConnection();
		
		TemplateSelection ts = 
			new TemplateSelection(inputSeqFile, pdbBlastDb, nrBlastDb, blastNumThreads, maxIter, eValueCutoff, 
					gtgScoreCutoff, maxHits, 
					similarityGraphRmsdCutoff, 
					psiblast, psipred, gtg, 
					baseName, outDir, gtgDir, 
					selectTemplates, 
					conn);
		
		ts.run();

	}


}

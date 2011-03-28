package owl.casp.pipeline;

import gnu.getopt.Getopt;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.TreeMap;

import org.xml.sax.SAXException;

import owl.core.runners.TcoffeeException;
import owl.core.runners.blast.BlastException;
import owl.core.runners.blast.BlastHit;
import owl.core.runners.blast.BlastHitList;
import owl.core.runners.blast.BlastXMLParser;
import owl.core.runners.tinker.TinkerError;
import owl.core.runners.tinker.TinkerRunner;
import owl.core.sequence.Sequence;
import owl.core.sequence.alignment.AlignmentConstructionException;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.structure.Pdb;
import owl.core.structure.PdbLoadException;
import owl.core.structure.TemplateList;
import owl.core.structure.alignment.StructAlignmentException;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatException;
import owl.core.util.MySQLConnection;
import owl.graphAveraging.ConsensusSquare;
import owl.graphAveraging.GraphAverager;
import owl.graphAveraging.GraphAveragerException;
import owl.graphAveraging.PhiPsiAverager;



public class RunAll {
	
	private static final String PROGRAM_NAME = "model_it";
	
	// graph averaging constants
	private static final String[] DEFAULT_CONTACT_TYPES = {"Ca", "Cg"};
	private static final double[] DEFAULT_CUTOFFS = {8.0, 8.0};
	private static final double DEFAULT_CCT = 0.5;

	// phi/psi averaging constants
	private static final double PHIPSI_CONSENSUS_THRESHOLD = 0.5;
	private static final int    PHIPSI_CONSENSUS_INTERVAL = 20;

	// tinker reconstruction constants
	private static final String FORCEFIELD_FILE = 		"/project/StruPPi/Software/tinker/amber/amber99.prm";
	private static final String TINKER_BIN_DIR = 		"/project/StruPPi/Software/tinker/bin";
	private static final double FORCE_CONSTANT_DIST =   100.0;
	private static final double FORCE_CONSTANT_TORSION = 1.0;
	private static final int DEFAULT_NUM_TINKER_MODELS = 40;
	
	// dssp constants
	private static final String DSSP_EXE = "/project/StruPPi/bin/dssp";
	
	// template alignment constants
	private static final String CONTACT_TYPE_PAUL ="Cb";
	private static final double CUTOFF_PAUL = 8;
	
	// sec struct comparison file
	private static final String SSCOMPARE_FILE_SUFFIX = ".ss.compare";
	
	// usage string and getopt options string
	private static final String USAGE = 
	"  [-D]:  don't perform tinker reconstruction. Default: does perform reconstruction \n" +
	"  [-F]:  don't use phi/psi constraints for reconstruction. Default: phi/psi \n" +
	"         constraints used\n" +
	"  [-O]:  don't enforce trans conformation for omega torsion angle (peptide bond). \n" +
	"         Default: enforcing\n" +
	"  [-s]:  contact conservation threshold (CCT). Default: "+DEFAULT_CCT +"\n" +
	"  [-n]:  number of tinker models. To get accurate results at least 20 must be \n" +
	"         specified. Default: "+DEFAULT_NUM_TINKER_MODELS+"\n"+
	"  [-A]:  if specified reconstruction will be run in parallel (EXPERIMENTAL)\n";
	private static final String GETOPT_OPTIONS = "DFOs:n:A";
	
	public static void main(String[] args) {
		
		File inputSeqFile = null;
		String pdbBlastDb = TemplateSelection.DEFAULT_PDB_BLASTDB;
		String nrBlastDb = TemplateSelection.NR_BLASTDB;
		int blastNumThreads = 1;
		String baseName = null;
		File outDir = new File(".");
		int maxHits = TemplateSelection.TOO_MANY_HITS;
		double similarityGraphRmsdCutoff = TemplateSelection.DEFAULT_SIMILARITY_GRAPH_RMSD_CUTOFF;
		int maxIter = TemplateSelection.DEFAULT_MAXITER;
		double eValueCutoff = TemplateSelection.EVALUE_CUTOFF_PREFILTER;
		int gtgScoreCutoff = TemplateSelection.DEFAULT_GTG_SCORE_CUTOFF;
		boolean psiblast = true;
		boolean psipred = true;
		boolean gtg = true;
		File gtgDir = new File(TemplateSelection.GTG_RESULTS_DIR);
		String selectTemplates = ""; 
		
		File templatesFile = null;
		File target2templatesFile = null;
		
		String paulMode = TemplatesAlignment.DEFAULT_PAULMODE;
		
		boolean reconstruct = true;
		boolean usePhiPsiConstraints = true;
		boolean forceTransOmega = true;
		double consensusThreshold = DEFAULT_CCT;
		int numberTinkerModels = DEFAULT_NUM_TINKER_MODELS;
		boolean parallel = false;

		String help = "Usage: \n" +
		PROGRAM_NAME+"\n" +
		"   -i :  file with input target sequence in FASTA format \n"+
		"  [-b]:  basename for output files. Default: basename of input sequence file \n"+
		"  [-o]:  output dir, where output files will be written. Default: current dir \n"+
		TemplateSelection.USAGE+
		"\n"+
		"  [-t]:  specify a templates list file: no template selection, use given list of \n" +
		"         templates\n" +
		"  [-L]:  specify a target to template(s) alignment fasta file. Requires -t. Will do graph\n" +
		"         averaging and reconstruction with given alignment. No template selection and alignment\n" +
		"         steps\n" +
		"\n"+
		TemplatesAlignment.USAGE+
		"\n" +
		USAGE+
		"\n";
		
		Getopt g = new Getopt(PROGRAM_NAME, args, "i:b:o:"+TemplateSelection.GETOPT_OPTIONS+"t:L:"+TemplatesAlignment.GETOPT_OPTIONS+GETOPT_OPTIONS+"h?");
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
				if (!selectTemplates.equals(TemplateSelection.USE_BLAST_TEMPLATES) && 
						!selectTemplates.equals(TemplateSelection.USE_PSIBLAST_TEMPLATES) && 
						!selectTemplates.equals(TemplateSelection.USE_GTG_TEMPLATES)) {
					System.err.println("Invalid value specified for -l option, allowed values are: "
							+TemplateSelection.USE_BLAST_TEMPLATES+", "+TemplateSelection.USE_PSIBLAST_TEMPLATES+" or "+TemplateSelection.USE_GTG_TEMPLATES);
					System.exit(1);
				}
				break;
			case 't':
				templatesFile = new File(g.getOptarg());
				break;
			case 'L':
				target2templatesFile = new File(g.getOptarg());
				break;								
			case 'p':
				paulMode = g.getOptarg();
				break;
			case 'D':
				reconstruct = false;
				break;
			case 'F':
				usePhiPsiConstraints = false;
				break;
			case 'O':
				forceTransOmega = false;
				break;
			case 's':
				consensusThreshold = Double.parseDouble(g.getOptarg());
				break;
			case 'n':
				numberTinkerModels = Integer.parseInt(g.getOptarg());
				break;
			case 'A':
				parallel = true;
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
		
		if (!psiblast && selectTemplates.equals(TemplateSelection.USE_PSIBLAST_TEMPLATES)) {
			System.err.println("Can't specify option use psi-blast templates (-l P) if skip psi-blast specified (-K)");
			System.exit(1);
		}
		if (!gtg && selectTemplates.equals(TemplateSelection.USE_GTG_TEMPLATES)){
			System.err.println("Can't specify option use GTG templates (-l G) if skip GTG specified (-T)");
			System.exit(1);			
		}
		if (target2templatesFile!=null && templatesFile==null) {
			System.err.println("Option -L requires a templates file (-t)");
			System.exit(1);
		}
		
		
		MySQLConnection conn = connectToDb();

		TemplateList templates = null;	
		MultipleSequenceAlignment target2templatesAln = null;
		MultipleSequenceAlignment templatesAln = null;
		Sequence inputSeq = null;
		
		// STEP 1: template selection

		if (templatesFile==null) {
		    System.out.println("################################################");
		    System.out.println("# STEP 1:  TEMPLATE SELECTION");
		    System.out.println("################################################");
		
			TemplateSelection ts = new TemplateSelection(inputSeqFile, 
					pdbBlastDb, nrBlastDb, blastNumThreads, maxIter, eValueCutoff, 
					gtgScoreCutoff, 
					maxHits, similarityGraphRmsdCutoff, psiblast, psipred, gtg, baseName, outDir, gtgDir, 
					selectTemplates, conn);

			try {
				templates = ts.run();

			} catch (FileFormatException e) {
				System.err.println("Input sequence file doesn't conform to FASTA format. Error: "+e.getMessage()+"\nExiting");
				System.exit(1);
			} catch (BlastException e) {
				System.err.println("Blast failed to run. Error: "+e.getMessage()+"\nExiting");
				System.exit(1);
			} catch (IOException e) {
				System.err.println("Problem while reading/writing files for blast/psipred. Error: "+e.getMessage()+"\nExiting");
				System.exit(1);			
			} catch (InterruptedException e) {
				System.err.println("Thread interrupted: "+e.getMessage());
				System.exit(1);
			}

			if (ts.foundFullMatch()) {
				// TemplateSelection already printed out info about this
				System.out.println("Full match was found. No need to model!");
				System.exit(0);
			}

			inputSeq = ts.getInputSequence();
		
		} else { // if -t or -L option, skip STEP 1 and load templates from file
			try {
				templates = new TemplateList(templatesFile);
				inputSeq = new Sequence();
				inputSeq.readFromFastaFile(inputSeqFile);
			} catch (IOException e) {
				System.err.println("Couldn't read templates or sequence file. Error: "+e.getMessage());
				System.exit(1);
			} catch (FileFormatException e) {
				System.err.println("Couldn't read input sequence file. Error: "+e.getMessage());
				System.exit(1);
			}
			// load pdb data 
			try {
				templates.loadPdbData(conn, TemplateSelection.PDBASE_DB);
			} catch (SQLException e) {
				System.err.println("Problems getting PDB data for templates alignment. Error: "+e.getMessage()+"\nCan't continue");
				System.exit(1);			
			} catch (PdbLoadException e) {
				System.err.println("Problems getting PDB data for templates alignment. Error: "+e.getMessage()+"\nCan't continue");
				System.exit(1);			
			}
		}
		
		if (target2templatesFile==null) {
			// STEP 2: alignment
			System.out.println("################################################");
			System.out.println("# STEP 2: ALIGNMENT");
			System.out.println("################################################");

			// check templates
			if (templates.size()==0) {
				System.out.println("No templates found, can't model!");
				System.exit(0);
			} else if (templates.size()==1) {
				target2templatesAln =getAlignment(templates, outDir, baseName, inputSeq.getSeq()); // gets the blast/psiblast alignment (depending of where templates come from)

				if (target2templatesAln == null){ // non of the above returned an alignment
					// We don't have blast or psiblast alignments: 
					//    we will perform tcoffee alignment against the single template sequence instead 
					//    of a profile alignment of the multi-templates
					// We create a trivial "alignment" of the 1 sequence to later use it with tcoffee
					templatesAln = getTrivialAlignment(templates);
				}
			} 		

			// align templates (if more than 1)
			if (templates.size()>1) {
				templates.loadRIGraphs(CONTACT_TYPE_PAUL, CUTOFF_PAUL);

				TemplatesAlignment ta = new TemplatesAlignment(templates, paulMode, baseName, outDir);

				try {
					System.out.println("Aligning templates structurally");
					templatesAln = ta.run();
				} catch (IOException e) {
					System.err.println("Problems while running the structural aligner program. Error "+e.getMessage()+"\nExiting");
					System.exit(1);
				} catch (StructAlignmentException e) {
					System.err.println("Problems while running the structural aligner program. Error "+e.getMessage()+"\nExiting");
					System.exit(1);
				}
			}

			// align target to templates
			if (target2templatesAln==null) { // if !=null then we are using a blas/psiblast alignment so target2templatesAln is set and we don't want to run tcoffee
				TargetToTemplatesAlignment tta = new TargetToTemplatesAlignment(inputSeq, templatesAln, baseName, outDir);
				try {
					System.out.println("Aligning target to templates");
					target2templatesAln = tta.run();
				} catch (TcoffeeException e) {
					System.err.println("Problems while running tcoffee for target to template alignment. Error "+e.getMessage()+"\nExiting");
					System.exit(1);
				} catch (IOException e) {
					System.err.println("Problems while running tcoffee for target to template alignment. Error "+e.getMessage()+"\nExiting");
					System.exit(1);
				} catch (InterruptedException e) {
					System.err.println("Thread was interrupted: "+e.getMessage());
					System.exit(1);
				}
			}
			// writing alignment
			try {
				TargetToTemplatesAlignment.writeTargetToTemplAl(outDir, baseName, target2templatesAln);
			} catch (IOException e) {
				System.err.println("Couldn't write the target to templates alignment file. Error: "+e.getMessage());
			}
		} else {
			try {
				// we load from file and don't check its contents yet, the GraphAverager will do that later
				target2templatesAln = new MultipleSequenceAlignment(target2templatesFile.getAbsolutePath(),"FASTA");
			} catch (AlignmentConstructionException e) {
				System.err.println("Couldn't construct the alignment from given file "+target2templatesFile+". Error: "+e.getMessage());
				System.exit(1);
			} catch (IOException e) {
				System.err.println("Couldn't read the alignment from given file "+target2templatesFile+". Error: "+e.getMessage());
				System.exit(1);
			} catch (FileFormatException e) {
				System.err.println("Wrong format of target to template alignment file "+target2templatesFile+". Error: "+e.getMessage());
				System.exit(1);
			}
		}

		// STEP 3: average graph and reconstruction
	    System.out.println("################################################");
	    System.out.println("# STEP 3: GRAPH AVERAGING AND RECONSTRUCTION");
	    System.out.println("################################################");

		// average graph
		// array to store the consensus graphs (one per contact type/cutoff) for the reconstruction
		RIGraph[] consensusGraphs = new RIGraph[DEFAULT_CONTACT_TYPES.length];
		System.out.println("Performing graph averaging");
		try {
			for (int i=0; i<DEFAULT_CONTACT_TYPES.length;i++) {
				templates.loadRIGraphs(DEFAULT_CONTACT_TYPES[i], DEFAULT_CUTOFFS[i]);
				GraphAverager ga = new GraphAverager(target2templatesAln, templates, inputSeq.getName(), inputSeq.getSeq());
				consensusGraphs[i] = ga.getConsensusGraph(consensusThreshold);
				String ctStr = DEFAULT_CONTACT_TYPES[i].replace("/", ":");
				ga.getAverageGraph().writeToFile(new File(outDir,baseName+"."+ctStr+"_"+DEFAULT_CUTOFFS[i]+".avrgd.cm").getAbsolutePath());
			}
		} catch (GraphAveragerException e) {
			System.err.println("Problems performing graph averaging. Error: "+e.getMessage()+"\nCan't continue with reconstruction");
			System.exit(1);
		} catch (IOException e) {
			// this is not fatal, we don't exit
			System.err.println("Couldn't write averaged graph file. Error "+e.getMessage());
		}

		// writing secondary structure alignment file
		writeSecStructMatching(templates, target2templatesAln, outDir, baseName, inputSeq.getName());

		// reconstruct
		if (reconstruct) {
			System.out.println("Reconstructing "+numberTinkerModels+" models");			
			// average phi/psi restraints
			TreeMap<Integer, ConsensusSquare> phiPsiConsensus = null;
			if (usePhiPsiConstraints) {
				System.out.println("Getting phi/psi consensus from templates");
				PhiPsiAverager phiPsiAvrger = new PhiPsiAverager(templates,target2templatesAln); // here we get only PDB data out of templates, so doesn't matter which graphs templates contains by now
				phiPsiConsensus = phiPsiAvrger.getConsensusPhiPsiOnTarget(PHIPSI_CONSENSUS_THRESHOLD, PHIPSI_CONSENSUS_INTERVAL, inputSeq.getName());
			}

			if (forceTransOmega) {
				System.out.println("Omega torsion angles will be restrained to trans conformation");
			}

			
			Pdb pdb = null;
			try {
				TinkerRunner tr = new TinkerRunner(TINKER_BIN_DIR,FORCEFIELD_FILE);
				if (parallel) System.out.println("Running distgeom in parallel through SGE");
				tr.reconstruct(inputSeq.getSeq(), consensusGraphs, phiPsiConsensus, forceTransOmega,
						numberTinkerModels, FORCE_CONSTANT_DIST, FORCE_CONSTANT_TORSION, 
						outDir.getPath(), baseName, true, parallel);
				pdb = tr.getStructure(tr.pickByLeastBoundViols());
			} catch (TinkerError e1) {
				System.err.println("Problem while running tinker for reconstruction. Error: "+e1.getMessage());
				System.exit(1);
			} catch (IOException e1) {
				System.err.println("Problem while running tinker for reconstruction. Error: "+e1.getMessage());
				System.exit(1);
			}
			
			File outpdbfile = new File(outDir,baseName+".reconstructed.pdb");
			
			try {
				pdb.writeToPDBFile(outpdbfile.getAbsolutePath());
				System.out.println("Done reconstruction. Final selected model written to " + outpdbfile);
			} catch (IOException e) {
				System.err.println("Done reconstruction. Couldn't write final selected model to pdb file! "+outpdbfile);
				System.exit(1);
			}

		} else {
			System.out.println("Skipping reconstruction (option -D)");
		}
		
	}
	
	/*-------------------------------- private methods ----------------------------------------*/
 	
	private static MySQLConnection connectToDb() {
		MySQLConnection conn;
		try {		
			conn = new MySQLConnection();

		} catch (SQLException e) {
			conn = null;
			System.err.println("Problems connecting to database. Error: "+e.getMessage()+"\nExiting");
			System.exit(1);
		}
		return conn;
	}
	
	private static MultipleSequenceAlignment getAlignment(TemplateList templates, File outDir, String baseName, String fullQuerySeq) {
		BlastHit hit = null;
		if (!templates.getSource().equals(TemplateList.SRC_OTHER)) {
			System.out.println("Target to template alignment from "+templates.getSource());
			hit = templates.getTemplate(0).getBlastHit();
			return hit.getMaxScoringHsp().getAlignmentFullSeqsWithPDBTag(fullQuerySeq, templates.getTemplate(0).getPdb().getSequence());
		} else { // we've read templates from a file (-t option)
			// we try to find the blast or psiblast alignments in the outdir, if one of them exists we use it
			File blastFile = new File(outDir,baseName+TemplateSelection.PSIBLASTOUT_PDB_SUFFIX);
			File psiblastFile = new File(outDir,baseName+TemplateSelection.BLASTOUT_PDB_SUFFIX);
			if (blastFile.exists()) {
				BlastXMLParser bxp;
				try {
					bxp = new BlastXMLParser(blastFile);
					BlastHitList hits = bxp.getHits();
					if (hits.contains(BlastHit.templateIdToSubjectId(templates.getTemplate(0).getId()))) {
						System.out.println("Target to template alignment from blast file "+blastFile);
						hit = hits.getHit(BlastHit.templateIdToSubjectId(templates.getTemplate(0).getId()));
						return hit.getMaxScoringHsp().getAlignmentFullSeqsWithPDBTag(fullQuerySeq, templates.getTemplate(0).getPdb().getSequence());
					} 
				} catch (SAXException e) {
					System.err.println("Problem parsing blast xml file "+blastFile+". Can't use blast alignment. Error: "+e.getMessage());
				} catch (IOException e) {
					System.err.println("Problem parsing blast xml file "+blastFile+". Can't use blast alignment. Error: "+e.getMessage());
				}

			} 
			if (psiblastFile.exists()) {
				BlastXMLParser bxp;
				try {
					bxp = new BlastXMLParser(psiblastFile);
					BlastHitList hits = bxp.getHits();
					if (hits.contains(BlastHit.templateIdToSubjectId(templates.getTemplate(0).getId()))) {
						System.out.println("Target to template alignment from psi-blast file "+psiblastFile);
						hit = hits.getHit(BlastHit.templateIdToSubjectId(templates.getTemplate(0).getId()));
						return hit.getMaxScoringHsp().getAlignmentFullSeqsWithPDBTag(fullQuerySeq, templates.getTemplate(0).getPdb().getSequence());
					}
				} catch (SAXException e) {
					System.err.println("Problem parsing psi-blast xml file "+psiblastFile+". Can't use psi-blast alignment. Error: "+e.getMessage());
				} catch (IOException e) {
					System.err.println("Problem parsing psi-blast xml file "+psiblastFile+". Can't use psi-blast alignment. Error: "+e.getMessage());
				}				
			}
		}
		// everything else failed, we return null: couldn't find this subject in the blast/psiblast files 
		return null;
	}
		
	private static MultipleSequenceAlignment getTrivialAlignment(TemplateList templates) {
		String[] tags = {templates.getTemplate(0).getId()};
		String[] seqs = {templates.getTemplate(0).getPdb().getSequence()};
		MultipleSequenceAlignment al;
		try {
			al = new MultipleSequenceAlignment(tags, seqs);
		} catch (AlignmentConstructionException e) {
			// this shouldn't happen
			al = null;
			System.err.println("Unexpected error: "+e.getMessage());
			System.exit(1);
		}
		return al;
		
	}
	
	private static void	writeSecStructMatching(TemplateList templates, MultipleSequenceAlignment target2templatesAln, File outDir, String baseName, String targetTag) {;
		target2templatesAln.addSecStructAnnotation(templates, DSSP_EXE);
		PrintStream out;
		try {
			File psipredFile = new File(outDir, baseName+TemplateSelection.PSIPRED_HORIZ_SUFFIX);
			if (!psipredFile.canRead()) throw new FileNotFoundException("Can't read psipred file "+psipredFile);
			out = new PrintStream(new File(outDir, baseName+SSCOMPARE_FILE_SUFFIX));
			target2templatesAln.writeWithSecStruct(out, targetTag, psipredFile, false);
			out.close();
		} catch (IOException e) {
			System.err.println("Couldn't write secondary structure comparison file. Error "+e.getMessage());
		}

	}

}

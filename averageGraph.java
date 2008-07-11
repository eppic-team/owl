import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.Alignment;
import proteinstructure.PairwiseSequenceAlignment;
import proteinstructure.Pdb;
import proteinstructure.PdbasePdb;
import proteinstructure.PredEval;
import proteinstructure.RIGraph;
import proteinstructure.Template;
import proteinstructure.TemplateList;
import proteinstructure.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import sequence.Sequence;
import tinker.TinkerRunner;
import tools.MySQLConnection;
//import tinker.TinkerRunner;
import gnu.getopt.Getopt;
import graphAveraging.ConsensusSquare;
import graphAveraging.GraphAverager;
import graphAveraging.PhiPsiAverager;


public class averageGraph {
	
	/*------------------------------ constants ------------------------------*/
	private static final double DEFAULT_THRESHOLD =		0.5;
	
	private static final String MUSCLE_BIN = 			"muscle";
	
	private static final String FORCEFIELD_FILE = 		"/project/StruPPi/Software/tinker/amber/amber99.prm";
	private static final String TINKER_BIN_DIR = 		"/project/StruPPi/Software/tinker/bin";
	private static final double FORCE_CONSTANT_DIST =   100.0;
	private static final double FORCE_CONSTANT_TORSION = 1.0;
	
	private static final double PHIPSI_CONSENSUS_THRESHOLD = 0.5;
	private static final int    PHIPSI_CONSENSUS_INTERVAL = 20;
	
	private static final String	PDB_DB = 				"pdbase";
	private static final String	DB_HOST = 				"white";								
	private static final String	DB_USER = 				MySQLConnection.getUserName();
	private static final String	DB_PWD = 				"nieve";
	
	private static final String CASP_CONTACT_TYPE =		"Cb";
	private static final double CASP_CUTOFF = 			8.0;
	
	private static enum Modes {	BENCHMARK, PREDICT, ALIGN	};
	
	private static final String CASP_METHOD_STR = "1. Template detection by Blast, Psiblast and Global Trace Graph (Heger et al. 2007) " +
			"\n2. Graph based sequence to multiple structure alignment " +
			"\n3. Model building by distance geometry using Tinker (Ponder et al. 1987) " +
			"\n4. Physics-based simulated annealing refinement with explicit solvent using Gromacs (Berendsen et al. 1995)";
	
	/*------------------------- private methods ------------------------------*/
	
	/**
	 * Creates a temporary sequence file with both target and template sequences and runs muscle
	 * @param aliFile
	 * @param outDir
	 * @param basename
	 * @param targetSeq
	 * @param targetTag
	 * @param templates
	 */
	private static void runMuscle(File aliFile, String outDir, String basename, String targetSeq, String targetTag, TemplateList templates) {
		
		// we have to create a temporary seq file for muscle's input
		File tmpSeqFile = new File(outDir,basename+".tmp.fasta");
		tmpSeqFile.deleteOnExit();
		String[] seqs = new String[templates.size()+1];
		String[] tags = new String[templates.size()+1];
		seqs[0]=targetSeq;
		tags[0]=targetTag;
		int i = 1;
		for (Template template:templates) {
			seqs[i] = template.getPdb().getSequence();
			tags[i] = template.getId();;
			i++;
		}	
		try {
			Sequence.writeSeqs(tmpSeqFile, seqs, tags);
		} catch (IOException e) {
			System.err.println("Couldn't write file "+tmpSeqFile.getAbsolutePath()+" with sequences for input of muscle. Exiting. ");
			System.exit(1);
		}

		try {
			Process muscleProc = Runtime.getRuntime().exec(MUSCLE_BIN+" -in "+tmpSeqFile.getCanonicalPath()+" -out "+aliFile.getCanonicalPath());
			if (muscleProc.waitFor()!=0) {
				System.err.println("muscle finished with an error (exit status "+muscleProc.exitValue()+"). Couldn't calculate alignment. Exiting");
				System.exit(1);
			}
		} catch (IOException e) {
			System.err.println("Couldn't run muscle. Error: "+e.getMessage()+". Exiting");
			System.exit(1);
		} catch (InterruptedException e) {
			System.err.println("Unexptected error: "+e.getMessage()+". Exiting.");
			System.exit(1);
		}
		
	}
	
	@SuppressWarnings("unused")
	private static String readCaspMethodFromFile(File caspMethodFile) throws IOException {
		String methodStr = "";
		BufferedReader br = new BufferedReader(new FileReader(caspMethodFile));
		String line;
		while ((line=br.readLine())!=null){
			methodStr += line+"\n";
		}
		return methodStr;
	}

	private static void checkSequences(Alignment ali, String targetTag, String targetSeq, Modes mode, File seqFile, String pdbCodeTarget, String pdbChainCodeTarget) {
		if(!ali.getSequenceNoGaps(targetTag).equals(targetSeq)) {
			if (mode==Modes.PREDICT) {
				System.err.println("Target sequence in alignment does not match target sequence from file "+seqFile);
			} else if (mode==Modes.BENCHMARK) {
				System.err.println("Target sequence in alignment does not match target sequence taken from db with pdb id "+pdbCodeTarget+pdbChainCodeTarget);
			} else if (mode==Modes.ALIGN) {
				// in ALIGN mode the 2 sequences shouldn't disagree
				System.err.println("Unexpected error, target dummy sequence and alignment dummy sequence disagree. This shouldn't be happening. Please report the bug!");
			}
			System.err.println("Trying to align sequences: ");
			try {
				PairwiseSequenceAlignment alCheck = new PairwiseSequenceAlignment(targetSeq,ali.getSequenceNoGaps(targetTag),"graph","alignment");
				alCheck.printAlignment();
			} catch (PairwiseSequenceAlignmentException e) {
				System.err.println("Error while creating alignment check, can't display an alignment, error: "+e.getMessage()+". The 2 sequences are: ");
				System.err.println("file:     "+targetSeq);
				System.err.println("alignment: "+ali.getSequenceNoGaps(targetTag));
			}
			System.exit(1);
		}
	}
	
	/*----------------------------- main --------------------------------*/
	public static void main(String[] args) throws Exception {
				
		String help =
				"Performs graph averaging. Three modes of operation: \n" +
				"a) alignment only:   specify only a list of templates (-P) \n" +
				"b) prediction:       specify a sequence file (-f) \n" +
				"c) benchmarking:     specify a pdb code+pdb chain code (-p) \n\n" + 
				"Usage: \n" +
				averageGraph.class.getName()+"\n" +
				"   -p <string> : target pdb code+target chain code (benchmarking), e.g. 1bxyA \n\n" +
				"   -f <file>   : file with target sequence to be predicted in FASTA format (prediction) \n\n"+
				"   -P <file>   : file with list of PDB ids (pdb codes+pdb chain codes in 1 column)\n" +
				"                 or a list of PDB files to be used as templates (or a mix of ids \n" +
				"                 and files)\n" +
				"   -t <string> : comma separated list of contact types \n" +
				"   -d <floats> : comma separated list of distance cutoffs (one per contact type) \n" +
				"   -b <string> : basename for output files (averaged graph, averaged graph with voters \n" +
				"                 and consensus graphs) \n"+
				"  [-a] <file>  : input alignment file. Default: perform multiple sequence alignment of \n" +
				"                 target and templates with muscle \n" +				
				"  [-s] <floats>: comma separated list of contact conservation thresholds (CCT) e.g. 0.5 \n" +
				"                 will predict an edge in target when present in half of the templates. \n" +
				"                 Default: "+DEFAULT_THRESHOLD+" is used\n"+
				"  [-o] <dir>   : output dir, where output files will be written. Default: current dir \n"+
				"  [-r] <int>   : if specified tinker's distgeom will be run to reconstruct the consensus \n" +
				"                 graph creating the specified number of models and finally outputting one \n" +
				"                 pdb file with the chosen model. If more than 1 CCT were specified, then \n" +
				"                 the first one is taken. This can take very long!\n" +
				"  [-R]         : used with -r above will run the reconstruction with given number of models \n" +
				"                 keeping the whole tinker output in the output dir\n" +
				"  [-c] <string>: write final reconstructed model also in CASP TS format using as AUTHOR the \n" +
				"                 specified string (with underscores instead of hyphens!). The target tag in the \n" +
				"                 target sequence file must comply with the CASP target naming convention, \n" +
				"                 e.g. T0100 \n" +
				"  [-m] <file>  : CURRENTLY DISABLED: take METHOD text from given file as the CASP TS METHOD text. If text in file \n" +
				"                 contains new lines it will be splitted over several METHOD lines \n" +
				"  [-F] <number>: use phi/psi consensus values as torsion angle constraints for the reconstruction.\n" +
				"                 The default consensus threshold is fixed at 50%.\n" +
				"                 The specified number will be taken as the angle interval for defining the consensus.\n" +
				"                 Default: "+PHIPSI_CONSENSUS_INTERVAL+"\n" +
				"  [-O]         : enforce trans conformation for the omega torsion angle. Default: not enforcing\n\n"+
				"A set of templates must always be specified (-P). A multiple sequence alignment of \n" +
				"target and templates should be specified as well (-a). If one is not given, then an \n" +
				"alignment is calculated with muscle. If no target sequence is given, a dummy sequence\n" + 
				"with the length of alignment is created. In this case reconstruction is not possible.\n" +
				"For reconstruction, please note that at least 20 models should be specified to \n" +
				"get a reasonable final selected model. \n\n";
		
		String[] cts = null;
		double[] cutoffs = null;
		
		File aliFile = null;
		String outDir = "."; // default we set to current
		String basename = "";
		
		String pdbCodeTarget =  "";
		String pdbChainCodeTarget = "";
		File templatesFile = null;
		
		File seqFile = null;
		
		Modes mode = Modes.ALIGN;
		
		double[] consensusThresholds = {DEFAULT_THRESHOLD};
		
		boolean reconstruct = false;
		int numberTinkerModels = 0;
		boolean keepModels = false;
		
		boolean usePhiPsiConstraints = false;
		int phiPsiConsensusInterval = PHIPSI_CONSENSUS_INTERVAL;
		
		boolean forceTransOmega = false;
		
		boolean casp = false;
		String caspAuthorStr = null;
		//File caspMethodFile = null;
		
		Getopt g = new Getopt(averageGraph.class.getName(), args, "p:P:d:t:a:s:o:b:f:r:Rc:m:F:Oh?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'p':
				pdbCodeTarget = g.getOptarg().substring(0, 4);
				pdbChainCodeTarget = g.getOptarg().substring(4);
				mode = Modes.BENCHMARK;
				break;
			case 'P':
				templatesFile = new File(g.getOptarg());
				break;
			case 'd':
				String[] tokens  = g.getOptarg().split(",");
				cutoffs = new double[tokens.length];
				for (int i=0;i<cutoffs.length;i++) { 
					cutoffs[i] = Double.parseDouble(tokens[i]);
				}
				break;
			case 't':
				cts = g.getOptarg().split(",");
				break;				
			case 'a':
				aliFile = new File(g.getOptarg());
				break;
			case 's':
				tokens  = g.getOptarg().split(",");
				consensusThresholds = new double[tokens.length];
				for (int i=0;i<consensusThresholds.length;i++) { 
					consensusThresholds[i] = Double.parseDouble(tokens[i]);
				}
				break;				
			case 'o':
				outDir = g.getOptarg();
				break;
			case 'b':
				basename = g.getOptarg();
				break;
			case 'f':
				seqFile = new File(g.getOptarg());
				mode = Modes.PREDICT;
				break;
			case 'r':
				numberTinkerModels = Integer.parseInt(g.getOptarg());
				reconstruct = true;
				break;
			case 'R':
				keepModels = true;
				break;				
			case 'c':
				caspAuthorStr = g.getOptarg().replace('_', '-');
				casp = true;
				break;
			case 'm':
				//caspMethodFile = new File(g.getOptarg());
				break;				
			case 'F':
				usePhiPsiConstraints = true;
				phiPsiConsensusInterval = Integer.parseInt(g.getOptarg());
				break;
			case 'O':
				forceTransOmega = true;
				break;																
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}
		
		// input checks
		if (templatesFile==null || cts==null || cutoffs==null || basename.equals("")) {
			System.err.println("Some missing option");
			System.err.println(help);
			System.exit(1);
		}
		if (mode==Modes.ALIGN && reconstruct) {
			System.err.println("Cannot reconstruct in align mode. Either provide a target sequence (-f) or pdb and chain code (-p) for benchmarking");
			System.err.println(help);
			System.exit(1);			
		}
		if (mode==Modes.ALIGN && casp) {
			System.err.println("The casp option (-c) is incompatible with align mode.");
			System.err.println(help);
			System.exit(1);						
		}
		if (mode==Modes.ALIGN && aliFile==null) {
			System.err.println("Alignment must be provided (-a) in align mode.");
			System.err.println(help);
			System.exit(1);			
		}
		if (seqFile!=null && !pdbCodeTarget.equals("")){
			// in this case mode has been assigned first to BENCHMARK ant then to PREDICT. It's better to check directly the variables seqFile and pdbCodeTarget
			System.err.println("Options -f (prediction), and -p (benchmark) are exclusive");
			System.err.println(help);
			System.exit(1);
		}
		if (mode==Modes.BENCHMARK && casp) {
			System.err.println("Options -p (benchmark) and -c (casp output) are incompatible");
			System.err.println(help);
			System.exit(1);			
		}
		// check that we specified same number of contact types and cutoffs
		if (cts.length!=cutoffs.length) {
			System.err.println("Specified list of contact types differs in length from list of cutoffs. Exiting");
			System.exit(1);			
		}
				
		MySQLConnection conn = new MySQLConnection(DB_HOST,DB_USER,DB_PWD);
		
		String targetSeq = null;
		String targetTag = null;
		Pdb targetPdb = null;
		
		if (mode==Modes.BENCHMARK) {
			// 1) benchmark: from a known structure sequence we repredict it based on template structures
			targetPdb = new PdbasePdb(pdbCodeTarget);
			targetPdb.load(pdbChainCodeTarget);
			targetSeq = targetPdb.getSequence();
			targetTag = pdbCodeTarget+pdbChainCodeTarget;			
		} else {
			// 2) prediction: from a sequence with unknown structure, we predict the structure based on template structures
			// or
			// 3) inspection: just create average and consensus graph for the set of templates using a dummy sequence (so nothing to do at this point)
			if (mode==Modes.PREDICT) {
				Sequence seq = new Sequence();
				try {
					seq.readFromFastaFile(seqFile);
				} catch (IOException e) {
					System.err.println("Couldn't read sequence file "+seqFile.getAbsolutePath()+": "+e.getMessage()+". Exiting");
					System.exit(1);
				}
				targetTag = seq.getName();
				// we take the sequence that we will use later to create the predicted graph
				targetSeq = seq.getSeq();
			}
			
			if (casp) { // in ALIGN mode targetTag is not defined, but we already checked that both options casp && ALIGN mode were not specified together
				Pattern p = Pattern.compile("T\\d\\d\\d\\d");
				Matcher m = p.matcher(targetTag);
				if (!m.matches()) {
					System.err.println("Target tag '"+targetTag+"' found in sequence file "+seqFile+" does not look like a CASP target name. If this is not a CASP prediction don't use the -c switch. Exiting.");
					System.exit(1);
				}
			}

		}
		
		// getting template structures

		TemplateList templates = null;
		try {
			templates = new TemplateList(templatesFile);
			templates.loadPdbData(conn, PDB_DB);
		} catch (IOException e) {
			System.err.println("Error while reading templates file "+templatesFile+": "+ e.getMessage()+"\nExiting");
			System.exit(1);
		}
		if (templates.size()==0) {
			System.err.println("Couldn't find any pdb code+pdb chain code in templates file\nExiting");
			System.exit(1);
		}
				
		// if an alignment file was not specified, perform alignment
		if (aliFile == null){  
			aliFile = new File(outDir,basename+".muscle_ali.fasta");
			// do alignment with muscle
			System.out.println("Performing alignment with muscle");
			runMuscle(aliFile, outDir, basename, targetSeq, targetTag, templates);
		}

		// read the alignment from file
		System.out.println("Reading alignment from "+aliFile);
		Alignment ali = new Alignment(aliFile.getCanonicalPath(), "FASTA");
		
		// checking that targetSeq and sequence from alignment match, we exit if not 
		if (mode==Modes.BENCHMARK || mode==Modes.PREDICT) { // in ALIGN mode targetTag and targetSeq are not defined, so we can't check
			checkSequences(ali, targetTag, targetSeq, mode, seqFile, pdbCodeTarget, pdbChainCodeTarget);
		}

		
		System.out.println("Averaging...");
		
		if (mode==Modes.BENCHMARK) {
			// printing headers for table of statistics
			System.out.printf("%10s\t","ct_cutoff"); 
			PredEval.printHeaders();
		}
		
		// array to store the consensus graphs (one per contact type/cutoff) for later use them in the reconstruction section
		RIGraph[] consensusGraphs = new RIGraph[cts.length];
		
		for (int ctIdx=0;ctIdx<cts.length;ctIdx++) {
			// if in benchmark we get the original graph to later calculate accuracy/coverage
			RIGraph originalGraph = null;
			if (mode==Modes.BENCHMARK) {
				originalGraph = targetPdb.get_graph(cts[ctIdx], cutoffs[ctIdx]);
			}
			
			// we get graphs for our templates
			templates.loadRIGraphs(cts[ctIdx], cutoffs[ctIdx]);
		
			String ctStr = cts[ctIdx].replace("/", ":");
			File avrgdGraphFile = new File(outDir,basename+"."+ctStr+"_"+cutoffs[ctIdx]+".avrgd.cm");
			File avrgdVotersGraphFile = new File(outDir,basename+"."+ctStr+"_"+cutoffs[ctIdx]+".avrgd.voters.cm");

			GraphAverager ga = null;
			if (mode==Modes.ALIGN) {
				ga = new GraphAverager(ali, templates);
			} else {
				ga = new GraphAverager(ali, templates, targetTag);
			}
			
			RIGraph averagedGraph = ga.getAverageGraph();
			//System.out.println("Writing average graph to " + avrgdGraphFile + " and average graph with voters to " + avrgdVotersGraphFile);
			averagedGraph.write_graph_to_file(avrgdGraphFile.getAbsolutePath());
			ga.writeAverageGraphWithVoters(avrgdVotersGraphFile.getAbsolutePath());

			for (double consensusThreshold: consensusThresholds) {
				File consGraphFile = new File(outDir,basename+"."+ctStr+"_"+cutoffs[ctIdx]+".CCT"+(String.format("%2.0f",consensusThreshold*100))+".cm");
				RIGraph consensusGraph = ga.getConsensusGraph(consensusThreshold);
				//System.out.printf("Writing consensus graph at CCT %2.0f to %s \n",consensusThreshold*100,consGraphFile);
				consensusGraph.write_graph_to_file(consGraphFile.getAbsolutePath());

				if (mode==Modes.BENCHMARK) {
					PredEval eval = consensusGraph.evaluatePrediction(originalGraph);
					System.out.printf("%6s_%3.1f\t%3.1f",cts[ctIdx],cutoffs[ctIdx],consensusThreshold);
					eval.printRow();
					//eval.printSummary();
				}

			}

			// printing pairwise overlap statistics (only if not in benchmark mode, otherwise it overlaps from the output of the benchmarking)
			if (mode!=Modes.BENCHMARK) {
				System.out.println("Contact type: "+cts[ctIdx]+", cutoff: "+cutoffs[ctIdx]);
				ga.printPairwiseOverlaps();
				int overlap = ga.getSumOfPairsOverlap();
				System.out.println("Sum of pairs contact overlap: " + overlap);
				// double consensus = ga.getEnsembleConsensusScore();
				// System.out.println("Ensemble consensus score:" + consensus);
			}
			
			// for reconstruction we take the first given consensus threshold value
			consensusGraphs[ctIdx] = ga.getConsensusGraph(consensusThresholds[0]);
		}
		
		// writing Cb 8 graph in CASP RR format if casp output was specified
		if (casp) {
			// we get the Cb graphs
			templates.loadRIGraphs(CASP_CONTACT_TYPE, CASP_CUTOFF);
			
			GraphAverager ga = new GraphAverager(ali, templates, targetTag);
			RIGraph averagedGraph = ga.getAverageGraph();
			File caspRRFile = new File(outDir,basename+".CASP.RR");
			int targetNum = Integer.parseInt(targetTag.substring(1)); // note: target tag must be like T0100, otherwise this fails!
			averagedGraph.setTargetNum(targetNum);
			averagedGraph.setCaspModelNum(1);
			averagedGraph.setAuthorStr(caspAuthorStr);
			averagedGraph.writeToCaspRRFile(caspRRFile.getAbsolutePath());
			System.out.println("CASP RR file with Cb 8 contact prediction written to " + caspRRFile);
		}
		
		
		// reconstruct
		if (reconstruct) {
			TreeMap<Integer, ConsensusSquare> phiPsiConsensus = null;
			if (usePhiPsiConstraints) {
				System.out.println("Getting phi/psi consensus from templates for reconstruction");
				PhiPsiAverager phiPsiAvrger = new PhiPsiAverager(templates,ali); // here we get only PDB data out of templates, so doesn't matter which graphs templates contains by now
				phiPsiConsensus = phiPsiAvrger.getConsensusPhiPsiOnTarget(PHIPSI_CONSENSUS_THRESHOLD, phiPsiConsensusInterval, targetTag);
			}
			
			System.out.println("Reconstructing");
			
			TinkerRunner tr = new TinkerRunner(TINKER_BIN_DIR,FORCEFIELD_FILE);
			
			Pdb pdb = null;
			if (keepModels) {
				tr.reconstruct(targetSeq, consensusGraphs, phiPsiConsensus, forceTransOmega, numberTinkerModels, FORCE_CONSTANT_DIST, FORCE_CONSTANT_TORSION, outDir, basename, false);
				pdb = tr.getStructure(tr.pickByLeastBoundViols());
			} else {
				pdb = tr.reconstruct(targetSeq, consensusGraphs, phiPsiConsensus, forceTransOmega, numberTinkerModels, FORCE_CONSTANT_DIST, FORCE_CONSTANT_TORSION);
			}
			
			File outpdbfile = new File(outDir,basename+".reconstructed.pdb");
			pdb.dump2pdbfile(outpdbfile.getAbsolutePath());
			System.out.println("Done reconstruction. Final selected model written to " + outpdbfile);
			if (casp) {
				File outcasptsfile = new File(outDir,basename+".reconstructed.casp");
				int targetNum = Integer.parseInt(targetTag.substring(1)); // note: target tag must be like T0100, otherwise this fails!
				pdb.setTargetNum(targetNum);
				pdb.setCaspModelNum(1);
				pdb.setCaspAuthorStr(caspAuthorStr);
				pdb.setCaspMethodStr(CASP_METHOD_STR);
				//pdb.setCaspMethodStr(readCaspMethodFromFile(caspMethodFile));
				pdb.setParents(templates.getIds());
				pdb.writeToCaspTSFile(outcasptsfile);
				System.out.println("Model written also to CASP TS file " + outcasptsfile);
			}

			if (mode==Modes.BENCHMARK) {
				System.out.printf("rmsd to native: %5.2f\n",pdb.rmsd(targetPdb, "Ca"));
			}
		}
		
		conn.close();
	}

}

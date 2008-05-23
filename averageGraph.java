import java.io.File;
import java.io.IOException;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.Alignment;
import proteinstructure.Pdb;
import proteinstructure.PdbasePdb;
import proteinstructure.PredEval;
import proteinstructure.RIGraph;
import proteinstructure.TemplateList;
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
	private static final double PHIPSI_CONSENSUS_THRESHOLD = 0.5;
	private static final int    PHIPSI_CONSENSUS_INTERVAL = 20;
	
	private static final String	PDB_DB = 				"pdbase";
	private static final String	DB_HOST = 				"white";								
	private static final String	DB_USER = 				MySQLConnection.getUserName();
	private static final String	DB_PWD = 				"nieve";
	
	private static final String CASP_CONTACT_TYPE =		"Cb";
	private static final double CASP_CUTOFF = 			8.0;
	
	/*------------------------- private methods ------------------------------*/
	
	/**
	 * Creates a temporary sequence file with both target and template sequences and runs muscle
	 * @param aliFile
	 * @param outDir
	 * @param basename
	 * @param targetSeq
	 * @param targetTag
	 * @param templateGraphs
	 */
	private static void runMuscle(File aliFile, String outDir, String basename, String targetSeq, String targetTag, Pdb[] templatePdbs) {
		
		// we have to create a temporary seq file for muscle's input
		File tmpSeqFile = new File(outDir,basename+".tmp.fasta");
		tmpSeqFile.deleteOnExit();
		String[] seqs = new String[templatePdbs.length+1];
		String[] tags = new String[templatePdbs.length+1];
		seqs[0]=targetSeq;
		tags[0]=targetTag;
		for (int i=0; i<templatePdbs.length;i++) {
			seqs[i+1] = templatePdbs[i].getSequence();
			tags[i+1] = templatePdbs[i].getPdbCode()+templatePdbs[i].getPdbChainCode();
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
	
	/*----------------------------- main --------------------------------*/
	public static void main(String[] args) throws Exception {
				
		String help =
				"Performs graph averaging. Three modes of operation: \n" +
				"a) benchmarking: specify a pdb code+pdb chain code (-p) \n" +
				"b) prediction:   specify a sequence file (-f) \n" +
				"c) inspection:   specify only a list of templates (-P) \n\n" + 
				"Usage: \n" +
				averageGraph.class.getName()+"\n" +
				"   -p <string> : target pdb code+target chain code (benchmarking), e.g. 1bxyA \n\n" +
				"   -f <file>   : file with target sequence to be predicted in FASTA format (prediction) \n\n"+
				"   -P <file>   : file with list of templates' pdb codes+pdb chain codes in 1 column\n" +
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
				"  [-c] <string>: write final reconstructed model also in CASP TS format using as AUTHOR the \n" +
				"                 specified string (with underscores instead of hyphens!). The target tag in the \n" +
				"                 target sequence file must comply with the CASP target naming convention, \n" +
				"                 e.g. T0100 \n" +
				"  [-m] <string>: use given string as the CASP TS METHOD text. If string contains new lines it \n" +
				"                 will be splitted over several METHOD lines \n" +
				"  [-F] <number>: use phi/psi consensus values as torsion angle constraints for the reconstruction.\n" +
				"                 The default consensus threshold is fixed at 50%.\n" +
				"                 The specified number will be taken as the angle interval for defining the consensus.\n" +
				"                 Default: "+PHIPSI_CONSENSUS_INTERVAL+"\n\n"+
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
		
		boolean benchmark = false;
		
		double[] consensusThresholds = {DEFAULT_THRESHOLD};
		
		boolean reconstruct = false;
		int numberTinkerModels = 0;
		
		boolean usePhiPsiConstraints = false;
		int phiPsiConsensusInterval = PHIPSI_CONSENSUS_INTERVAL;
		
		boolean casp = false;
		String caspAuthorStr = null;
		String caspMethodStr = null;
		
		Getopt g = new Getopt(averageGraph.class.getName(), args, "p:P:d:t:a:s:o:b:f:r:c:m:F:h?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'p':
				pdbCodeTarget = g.getOptarg().substring(0, 4);
				pdbChainCodeTarget = g.getOptarg().substring(4);
				benchmark = true;
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
				benchmark = false;
				break;
			case 'r':
				numberTinkerModels = Integer.parseInt(g.getOptarg());
				reconstruct = true;
				break;
			case 'c':
				caspAuthorStr = g.getOptarg().replace('_', '-');
				casp = true;
				break;
			case 'm':
				caspMethodStr = g.getOptarg();
				break;				
			case 'F':
				usePhiPsiConstraints = true;
				phiPsiConsensusInterval = Integer.parseInt(g.getOptarg());
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
		if (seqFile==null && !benchmark && reconstruct) {
			System.err.println("Cannot reconstruct in inspect mode. Either provide a target sequence (-f) or pdb- and chain code (-p/-c) for benchmarking");
			System.err.println(help);
			System.exit(1);			
		}
		if (seqFile==null && !benchmark && reconstruct && aliFile == null) {
			System.err.println("Alignment has to be provided (-a) in inspect mode.");
			System.err.println(help);
			System.exit(1);			
		}
		if (seqFile!=null && !pdbCodeTarget.equals("")){
			System.err.println("Options -f (prediction), and -p/-c (benchmark) are exclusive");
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
		
		if (benchmark) {
			// 1) benchmark: from a known structure sequence we repredict it based on template structures
			targetPdb = new PdbasePdb(pdbCodeTarget);
			targetPdb.load(pdbChainCodeTarget);
			targetSeq = targetPdb.getSequence();
			targetTag = pdbCodeTarget+pdbChainCodeTarget;			
		} else {
			// 2) prediction: from a sequence with unknown structure, we predict the structure based on template structures
			// or
			// 3) inspection (seqFile=null): just create average and consensus graph for the set of templates using a dummy sequence
			if(seqFile != null) {
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
			
			if (casp) {
				Pattern p = Pattern.compile("T\\d\\d\\d\\d");
				Matcher m = p.matcher(targetTag);
				if (!m.matches()) {
					System.err.println("Target tag '"+targetTag+"' found in sequence file "+seqFile+" does not look like a CASP target name. If this is not a CASP prediction don't use the -c switch. Exiting.");
					System.exit(1);
				}
			}

		}
		
		// getting template structures
		//TODO eventually we should read using RIGEnsemble.loadFromListFile adding this case: list file is a list of pdbCodes+pdbChainCodes
		//     Then RIGEnsemble should also have an alignment to be able to use it for ensembles not sharing same sequence
		//     By using loadFromListFile we would get the benefit of reading templatesFile that contain a list of pdb files/casp RR files/cm files/cif files
		String[] codesTemplates = null;
		try {
			codesTemplates = TemplateList.readIdsListFile(templatesFile); 
		} catch (IOException e) {
			System.err.println("Error while reading templates file "+templatesFile+": "+ e.getMessage()+"\nExiting");
			System.exit(1);
		}
		if (codesTemplates.length==0) {
			System.err.println("Couldn't find any pdb code+pdb chain code in templates file\nExiting");
			System.exit(1);
		}
		
		Pdb[] templatePdbs = new Pdb[codesTemplates.length];
		for (int i=0;i<codesTemplates.length;i++) {
			Pdb pdb = new PdbasePdb(codesTemplates[i].substring(0, 4), PDB_DB, conn);
			pdb.load(codesTemplates[i].substring(4));
			templatePdbs[i]=pdb;
		}
		
		// if an alignment file was not specified, perform alignment
		if (aliFile == null){  
			aliFile = new File(outDir,basename+".muscle_ali.fasta");
			// do alignment with muscle
			System.out.println("Performing alignment with muscle");
			runMuscle(aliFile, outDir, basename, targetSeq, targetTag, templatePdbs);
		}

		// read the alignment from file
		System.out.println("Reading alignment from "+aliFile);
		Alignment ali = new Alignment(aliFile.getCanonicalPath(), "FASTA");
		if(seqFile == null) {
			targetSeq = GraphAverager.makeDummySequence(ali.getAlignmentLength());
			targetTag = GraphAverager.makeDummyTag();
			ali.addSequence(targetTag, targetSeq);
		}
		
		System.out.println("Averaging...");
		
		if (benchmark) {
			// printing headers for table of statistics
			System.out.printf("%10s\t","ct_cutoff"); 
			PredEval.printHeaders();
		}
		
		// array to store one graph per contact type for later use them in the reconstruction section
		RIGraph[] graphsForReconstruction = new RIGraph[cts.length];
		
		for (int ctIdx=0;ctIdx<cts.length;ctIdx++) {
			// if in benchmark we get the original graph to later calculate accuracy/coverage
			RIGraph originalGraph = null;
			if (benchmark) {
				originalGraph = targetPdb.get_graph(cts[ctIdx], cutoffs[ctIdx]);
			}

			TreeMap<String, RIGraph> templateGraphs = new TreeMap<String, RIGraph>();

			for (int i=0;i<codesTemplates.length;i++) {
				RIGraph graph = templatePdbs[i].get_graph(cts[ctIdx], cutoffs[ctIdx]);
				templateGraphs.put(graph.getPdbCode()+graph.getPdbChainCode(),graph);
			}

			//System.out.println("Contact type: "+cts[ctIdx]+", cutoff: "+cutoffs[ctIdx]);

			String ctStr = cts[ctIdx].replace("/", ":");
			File avrgdGraphFile = new File(outDir,basename+"."+ctStr+"_"+cutoffs[ctIdx]+".avrgd.cm");
			File avrgdVotersGraphFile = new File(outDir,basename+"."+ctStr+"_"+cutoffs[ctIdx]+".avrgd.voters.cm");

			GraphAverager ga = new GraphAverager(targetSeq, ali, templateGraphs, targetTag);
			ga.printPairwiseOverlaps();
			int overlap = ga.getSumOfPairsOverlap();
			System.out.println("Sum of pairs contact overlap: " + overlap);
//			double consensus = ga.getEnsembleConsensusScore();
//			System.out.println("Ensemble consensus score:" + consensus);
			
			RIGraph averagedGraph = ga.getAverageGraph();
			//System.out.println("Writing average graph to " + avrgdGraphFile + " and average graph with voters to " + avrgdVotersGraphFile);
			averagedGraph.write_graph_to_file(avrgdGraphFile.getAbsolutePath());
			ga.writeAverageGraphWithVoters(avrgdVotersGraphFile.getAbsolutePath());

			for (double consensusThreshold: consensusThresholds) {
				File consGraphFile = new File(outDir,basename+"."+ctStr+"_"+cutoffs[ctIdx]+".CCT"+(String.format("%2.0f",consensusThreshold*100))+".cm");
				RIGraph consensusGraph = ga.getConsensusGraph(consensusThreshold);
				//System.out.printf("Writing consensus graph at CCT %2.0f to %s \n",consensusThreshold*100,consGraphFile);
				consensusGraph.write_graph_to_file(consGraphFile.getAbsolutePath());

				if (benchmark) {
					PredEval eval = consensusGraph.evaluatePrediction(originalGraph);
					System.out.printf("%6s_%3.1f\t%3.1f",cts[ctIdx],cutoffs[ctIdx],consensusThreshold);
					eval.printRow();
					//eval.printSummary();
				}

			}
			
			// for reconstruction we take the first given consensus threshold value
			graphsForReconstruction[ctIdx] = ga.getConsensusGraph(consensusThresholds[0]);
		}
		
		// writing Cb 8 graph in CASP RR format if casp output was specified
		if (casp) {
			TreeMap<String, RIGraph> templateGraphs = new TreeMap<String, RIGraph>();

			for (int i=0;i<codesTemplates.length;i++) {
				RIGraph graph = templatePdbs[i].get_graph(CASP_CONTACT_TYPE, CASP_CUTOFF);
				templateGraphs.put(graph.getPdbCode()+graph.getPdbChainCode(),graph);
			}
			
			GraphAverager ga = new GraphAverager(targetSeq, ali, templateGraphs, targetTag);
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
			System.out.println("Getting phi/psi consensus from templates for reconstruction");
			TreeMap<Integer, ConsensusSquare> phiPsiConsensus = null;
			if (usePhiPsiConstraints) {
				// we are re-reading from db the PDB data, this is really inefficient (we already have them in templatePdbs)
				// TODO make templatePdbs a TemplateList 
				TemplateList templates = new TemplateList(codesTemplates);
				templates.loadPdbData(conn, PDB_DB);
				PhiPsiAverager phiPsiAvrger = new PhiPsiAverager(templates,ali);
				phiPsiConsensus = phiPsiAvrger.getConsensusPhiPsiOnTarget(PHIPSI_CONSENSUS_THRESHOLD, phiPsiConsensusInterval, targetTag);
			}
			
			System.out.println("Reconstructing");
			
			TinkerRunner tr = new TinkerRunner(TINKER_BIN_DIR,FORCEFIELD_FILE);
			Pdb pdb = tr.reconstruct(targetSeq, graphsForReconstruction, phiPsiConsensus, numberTinkerModels);
			File outpdbfile = new File(outDir,basename+".reconstructed.pdb");
			pdb.dump2pdbfile(outpdbfile.getAbsolutePath());
			System.out.println("Done reconstruction. Final selected model written to " + outpdbfile);
			if (casp) {
				File outcasptsfile = new File(outDir,basename+".reconstructed.casp");
				int targetNum = Integer.parseInt(targetTag.substring(1)); // note: target tag must be like T0100, otherwise this fails!
				pdb.setTargetNum(targetNum);
				pdb.setCaspModelNum(1);
				pdb.setCaspAuthorStr(caspAuthorStr);
				pdb.setCaspMethodStr(caspMethodStr);
				pdb.setParents(codesTemplates);
				pdb.writeToCaspTSFile(outcasptsfile);
				System.out.println("Model written also to CASP TS file " + outcasptsfile);
			}

			if (benchmark) {
				System.out.printf("rmsd to native: %5.2f\n",pdb.rmsd(targetPdb, "Ca"));
			}
		}
		
		conn.close();
	}

}

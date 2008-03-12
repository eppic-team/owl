import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.Alignment;
import proteinstructure.Pdb;
import proteinstructure.PdbasePdb;
import proteinstructure.PredEval;
import proteinstructure.RIGraph;
import tinker.TinkerRunner;
//import tinker.TinkerRunner;
import gnu.getopt.Getopt;
import graphAveraging.GraphAverager;


public class averageGraph {
	
	/*------------------------------ constants ------------------------------*/
	private static final double DEFAULT_THRESHOLD =		0.5;
	
	private static final String MUSCLE_BIN = 			"muscle";
	
	private static final String FORCEFIELD_FILE = 		"/project/StruPPi/Software/tinker/amber/amber99.prm";
	private static final String TINKER_BIN_DIR = 		"/project/StruPPi/Software/tinker/bin";
	
	/*------------------------- private methods ------------------------------*/
	/**
	 * Reads a sequence file in FASTA format
	 * @param seqFile
	 * @return an array with 2 members: FASTA tag and sequence
	 */
	private static String[] readSeq(File seqFile) {
		String tag = "";
		String seq = "";
		try {
			BufferedReader fileIn = new BufferedReader(new FileReader(seqFile));
			String nextLine;
			// read sequences
			while((nextLine = fileIn.readLine()) != null) {
				Pattern p = Pattern.compile("^>(.*)$");
				Matcher m = p.matcher(nextLine);
				if (m.find()) {
					tag = m.group(1).trim();
				} else {
					seq += nextLine.trim();
				}
			}
		} catch (IOException e) {
			System.err.println("Couldn't read sequence file "+seqFile.getAbsolutePath()+": "+e.getMessage()+". Exiting");
			System.exit(1);
		}
		String[] tagAndseq = {tag,seq};
		return tagAndseq;
	}
	
	/**
	 * Writes given sequences and tags to given sequence file in FASTA format
	 * @param seqFile
	 * @param seqs
	 * @param tags
	 */
	private static void writeSeqs(File seqFile, String[] seqs, String[] tags) {
		try {
			PrintStream Out = new PrintStream(new FileOutputStream(seqFile));
			int len = 80;
			for (int seqIdx=0;seqIdx<seqs.length;seqIdx++) { 
				Out.println(">"+tags[seqIdx]);
				for(int i=0; i<seqs[seqIdx].length(); i+=len) {
					Out.println(seqs[seqIdx].substring(i, Math.min(i+len,seqs[seqIdx].length())));
				}		
			}
			Out.close();
		} catch (IOException e) {
			System.err.println("Couldn't write file "+seqFile.getAbsolutePath()+" with sequences for input of muscle. Exiting. ");
			System.exit(1);
		}
	}
	
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
		writeSeqs(tmpSeqFile, seqs, tags);		

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
	
	/**
	 * Reads a list file containing 1 column of pdbCodes+pdbChainCodes, e.g. 1bxyA
	 * @param templatesFile
	 * @return
	 */
	private static String[] readTemplatesFile(File templatesFile) {
		ArrayList<String> codesAL = new ArrayList<String>(); 
		try {
			BufferedReader fileIn = new BufferedReader(new FileReader(templatesFile));
			String line;
			int lineCount=0;
			while((line = fileIn.readLine()) != null) {
				lineCount++;
				if (line.length()!=0 && !line.startsWith("#")) {
					Pattern p = Pattern.compile("^\\d\\w\\w\\w\\w");
					Matcher m = p.matcher(line);
					if (m.matches()) {
						codesAL.add(line);
					} else {
						System.err.println("Line "+lineCount+" in templates file doesn't look like a pdb code+pdb chain code");
					}
				}
			}
		} catch (FileNotFoundException e) {
			System.err.println("Couldn't find templates file "+templatesFile+" . Exiting.");
			System.exit(1);
		} catch (IOException e) {
			System.err.println("Error while reading templates file "+templatesFile+". Exiting");
			System.exit(1);
		}
		if (codesAL.isEmpty()) {
			System.err.println("Couldn't find any pdb code+pdb chain code in templates file. Exiting");
			System.exit(1);
		}
		String[] codes = new String[codesAL.size()];
		codesAL.toArray(codes);
		return codes;
	}
	
	/*----------------------------- main --------------------------------*/
	public static void main(String[] args) throws Exception {
				
		String help = "Usage: \n" +
				averageGraph.class.getName()+"\n" +
				"  -p: target pdb code+target chain code (benchmarking), e.g. 1bxyA \n" +
				"\n"+
				"  -f: file with target sequence to be predicted in FASTA format (prediction) \n"+
				"\n"+
				"  -P: file with list of templates' pdb codes+pdb chain codes in 1 column\n" +
				"  -t: comma separated list of contact types \n" +
				"  -d: comma separated list of distance cutoffs (one per contact type) \n" +
				"  -b: basename for output files (averaged graph, averaged graph with voters and consensus graphs) \n"+
				"  [-a]: input alignment file, if not specified, a multiple sequence alignment of target and templates will be calculated with muscle \n" +				
				"  [-s]: comma separated list of contact conservation thresholds (CCT) e.g. 0.5 will predict an edge in target when present in half of the templates. If not specified "+DEFAULT_THRESHOLD+" is used\n"+
				"  [-o]: output dir, where output files will be written. If not specified current dir will be used \n"+
				"  [-r]: if specified tinker's distgeom will be run to reconstruct the consensus graph creating the specified number of models and finally outputting 1 pdb file with the chosen model. If more than 1 CCT were specified, then the first one is taken. This can take very long!\n"+
				"Performs graph averaging. Two modes of operation: \n" +
				"a) benchmarking: specify a pdb code/pdb chain code (-p/-c) \n" +
				"b) prediction:   specify a sequence file (-f) \n" +
				"A set of templates must always be specified (-P/-C). Also as an input a multiple sequence alignment of target and templates should be specified (-a). If one is not given, then a an alignment is calculated with muscle. \n" +
				"For reconstruction, please not that at least 20 models should be specified to get a reasonable final selected model. \n";
		
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
		
		Getopt g = new Getopt(averageGraph.class.getName(), args, "p:P:d:t:a:s:o:b:f:r:h?");
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
		if (seqFile==null && (pdbCodeTarget.equals("") && pdbChainCodeTarget.equals(""))) {
			System.err.println("Either a sequence file or a target pdb code and chain code (for benchmarking) must be specified");
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
		
		
		
		String targetSeq = null;
		String targetTag = null;
		Pdb targetPdb = null;
		
		if (benchmark) {
			// 1) benchmark: from a known structure sequence we repredict it based on a msa with known structures
			targetPdb = new PdbasePdb(pdbCodeTarget);
			targetPdb.load(pdbChainCodeTarget);
			targetSeq = targetPdb.getSequence();
			targetTag = pdbCodeTarget+pdbChainCodeTarget;			
		} else {
			// 2) prediction: from a sequence with unknown structure, we predict the structure based on a msa with known structures
			String[] tagAndSeq = readSeq(seqFile);
			targetTag = tagAndSeq[0];
			// we take the sequence that we will use later to create the predicted graph
			targetSeq = tagAndSeq[1];
		}
		
		// getting template structures
		//TODO eventually we should read using RIGEnsemble.loadFromListFile adding this case: list file is a list of pdbCodes+pdbChainCodes
		//     Then RIGEnsemble should also have an alignment to be able to use it for ensembles not sharing same sequence
		//     By using loadFromListFile we would get the benefit of reading templatesFile that contain a list of pdb files/casp RR files/cm files/cif files
		String[] codesTemplates = readTemplatesFile(templatesFile); 
		Pdb[] templatePdbs = new Pdb[codesTemplates.length];
		for (int i=0;i<codesTemplates.length;i++) {
			Pdb pdb = new PdbasePdb(codesTemplates[i].substring(0, 4));
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
		
		
		// reconstruct
		if (reconstruct) {
			System.out.println("Reconstructing");
			TinkerRunner tr = new TinkerRunner(TINKER_BIN_DIR,FORCEFIELD_FILE);
			Pdb pdb = tr.reconstruct(targetSeq, graphsForReconstruction, numberTinkerModels);
			File outpdbfile = new File(outDir,basename+".reconstructed.pdb");
			pdb.dump2pdbfile(outpdbfile.getAbsolutePath());
			System.out.println("Done reconstruction. Final selected model written to " + outpdbfile);
			if (benchmark) {
				System.out.printf("rmsd to native: %5.2f\n",pdb.rmsd(targetPdb, "Ca"));
			}
		}
	}

}

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
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
	
	private static final double DEFAULT_THRESHOLD =		0.5;
	
	private static final String MUSCLE_BIN = 			"muscle";
	
	private static final String FORCEFIELD_FILE = 		"/project/StruPPi/Software/tinker/amber/amber99.prm";
	private static final String TINKER_BIN_DIR = 		"/project/StruPPi/Software/tinker/bin";
	
	private static final int 	NUMBER_TINKER_MODELS = 	40;
	
	private static String[] readSeq(String seqFile) {
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
			System.err.println("Couldn't read sequence file "+seqFile+": "+e.getMessage()+". Exiting");
			System.exit(1);
		}
		String[] tagAndseq = {tag,seq};
		return tagAndseq;
	}
	
	private static void writeSeqs(String seqFile, String[] seqs, String[] tags) {
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
			System.err.println("Couldn't write file "+seqFile+" with sequences for input of muscle. Exiting. ");
			System.exit(1);
		}
	}
	
	public static void main(String[] args) throws Exception {
				
		String help = "Usage: \n" +
				averageGraph.class.getName()+"\n" +
				"  -p: target pdb code (benchmarking) \n" +
				"  -c: target chain code (benchmarking) \n" +
				"\n"+
				"  -f: file with target sequence to be predicted in FASTA format (prediction) \n"+
				"\n"+
				"  -P: comma separated list of templates' pdb codes \n" +
				"  -C: comma separated list of templates' pdb chain codes \n" +
				"  -t: contact type \n" +
				"  -d: distance cutoff \n" +
				"  -b: basename for output files (averaged graph, averaged graph with voters and consensus graphs) \n"+
				"  [-a]: input alignment file, if not specified, a multiple sequence alignment of target and templates will be calculated with muscle \n" +				
				"  [-s]: comma separated list of contact conservation thresholds (CCT) e.g. 0.5 will predict an edge in target when present in half of the templates. If not specified "+DEFAULT_THRESHOLD+" is used\n"+
				"  [-o]: output dir, where output files will be written. If not specified current dir will be used \n"+
				"  [-r]: if specified tinker's distgeom will be run to reconstruct the consensus graph, outputting 1 pdb file with the chosen model among "+NUMBER_TINKER_MODELS+" models. If more than 1 CCT were specified, then the first one is taken. This can take very long!\n"+
				"Performs graph averaging. Two modes of operation: \n" +
				"a) benchmarking: specify a pdb code/pdb chain code (-p/-c) \n" +
				"b) prediction:   specify a sequence file (-f) \n" +
				"A set of templates must always be specified (-P/-C). Also as an input a multiple sequence alignment of target and templates should be specified (-a). If one is not given, then a an alignment is calculated with muscle. \n";
		
		String ct = "";
		double cutoff = 0.0;
		
		String aliFile = "";
		String outDir = "."; // default we set to current
		String basename = "";
		
		String pdbCodeTarget =  "";
		String pdbChainCodeTarget = "";
		String[] pdbCodesTemplates = null; 
		String[] pdbChainCodesTemplates = null;
		
		String seqFile = "";
		
		boolean benchmark = false;
		boolean reconstruct = false;
		
		double[] consensusThresholds = {DEFAULT_THRESHOLD};
		
		
		Getopt g = new Getopt(averageGraph.class.getName(), args, "p:c:P:C:d:t:a:s:o:b:f:rh?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'p':
				pdbCodeTarget = g.getOptarg();
				benchmark = true;
				break;
			case 'c':
				pdbChainCodeTarget = g.getOptarg();
				break;				
			case 'P':
				pdbCodesTemplates = g.getOptarg().split(",");
				break;
			case 'C':
				pdbChainCodesTemplates = g.getOptarg().split(",");
				break;
			case 'd':
				cutoff = Double.parseDouble(g.getOptarg());
				break;
			case 't':
				ct = g.getOptarg();
				break;				
			case 'a':
				aliFile = g.getOptarg();
				break;
			case 's':
				String[] tokens  = g.getOptarg().split(",");
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
				seqFile = g.getOptarg();
				benchmark = false;
				break;
			case 'r':
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
		if (pdbCodesTemplates==null || pdbChainCodesTemplates==null || ct.equals("") || cutoff==0.0 || basename.equals("")) {
			System.err.println("Some missing option");
			System.err.println(help);
			System.exit(1);
		}
		if (seqFile.equals("") && (pdbCodeTarget.equals("") && pdbChainCodeTarget.equals(""))) {
			System.err.println("Either a sequence file or a target pdb code and chain code (for benchmarking) must be specified");
			System.err.println(help);
			System.exit(1);			
		}
		if (!seqFile.equals("") && !pdbCodeTarget.equals("")){
			System.err.println("Options -f (prediction), and -p/-c (benchmark) are exclusive");
			System.err.println(help);
			System.exit(1);
		}
		
		// check that we specified same number of template pdb codes and template pdb chain codes
		if (pdbCodesTemplates.length!=pdbChainCodesTemplates.length) {
			System.err.println("Specified list of pdb codes differs in length from list of pdb chain codes. Exiting");
			System.exit(1);
		}
		
		String targetSeq = null;
		String targetTag = null;
		RIGraph originalGraph = null;
		Pdb targetPdb = null;
		
		if (benchmark) {
			// 1) benchmark: from a known structure sequence we repredict it based on a msa with known structures
			targetPdb = new PdbasePdb(pdbCodeTarget);
			targetPdb.load(pdbChainCodeTarget);
			targetSeq = targetPdb.getSequence();
			targetTag = pdbCodeTarget+pdbChainCodeTarget;
			originalGraph = targetPdb.get_graph(ct, cutoff);
		}
		else {
			// 2) prediction: from a sequence with unknown structure, we predict the structure based on a msa with known structures
			String[] tagAndSeq = readSeq(seqFile);
			targetTag = tagAndSeq[0];
			// we take the sequence that we will use later to create the predicted graph
			targetSeq = tagAndSeq[1];
		}
		
				
		TreeMap<String, RIGraph> templateGraphs = new TreeMap<String, RIGraph>();
		
		for (int i=0;i<pdbCodesTemplates.length;i++) {
			Pdb pdb = new PdbasePdb(pdbCodesTemplates[i]);
			pdb.load(pdbChainCodesTemplates[i]);
			RIGraph graph = pdb.get_graph(ct, cutoff);
			templateGraphs.put(graph.getPdbCode()+graph.getPdbChainCode(),graph);
		}
		
		// if an alignment file was not specified, perform alignment
		if (aliFile.equals("")){  
			aliFile = new File(outDir,basename+".muscle_ali.fasta").getAbsolutePath();
			// do alignment with muscle
			System.out.println("Performing alignment with muscle");
			if (seqFile.equals("")) { // if seqFile doesn't exist (not given as input)
				// we have to create a temporary seq file for muscle's input
				File seqF = new File(outDir,basename+".tmp.fasta");
				seqFile = seqF.getAbsolutePath();
				seqF.deleteOnExit();
				String[] seqs = new String[templateGraphs.size()+1];
				String[] tags = new String[templateGraphs.size()+1];
				int i=0;
				seqs[i]=targetSeq;
				tags[i]=targetTag;
				for (String tag:templateGraphs.keySet()) {
					i++;
					seqs[i] = templateGraphs.get(tag).getSequence();
					tags[i] = tag;
				}	
				writeSeqs(seqFile, seqs, tags);
			}
			System.out.println(MUSCLE_BIN+" -in "+seqFile+" -out "+aliFile);
			Process muscleProc = Runtime.getRuntime().exec(MUSCLE_BIN+" -in "+seqFile+" -out "+aliFile);
			if (muscleProc.waitFor()!=0) {
				System.err.println("muscle finished with an error (exit status "+muscleProc.exitValue()+"). Couldn't calculate alignment. Exiting");
				System.exit(1);
			}
		}
		
		// read the alignment from file
		System.out.println("Reading alignment from "+aliFile);
		Alignment ali = new Alignment(aliFile, "FASTA");
		
		
		// averaging
		String templatesStr="";
		for (String tag:templateGraphs.keySet()){
			templatesStr+= tag+" ";
		}
		System.out.println("Calculating average graph of "+targetTag+ " based on "+templatesStr);
		System.out.println("Contact type for graphs is "+ct+", cutoff "+cutoff);
		
		String avrgdGraphFile = new File(outDir,basename+".avrgd.cm").getAbsolutePath();
		String avrgdVotersGraphFile = new File(outDir,basename+".avrgd.voters.cm").getAbsolutePath();;

		GraphAverager ga = new GraphAverager(targetSeq, ali, templateGraphs, targetTag);
		RIGraph averagedGraph = ga.getAverageGraph();
		System.out.println("Writing average graph to file " + avrgdGraphFile + " and average graph with voters to " + avrgdVotersGraphFile);
		averagedGraph.write_graph_to_file(avrgdGraphFile);
		ga.writeAverageGraphWithVoters(avrgdVotersGraphFile);
		
		for (double consensusThreshold: consensusThresholds) {
			String consGraphFile = new File(outDir,basename+".CCT"+(String.format("%2.0f",consensusThreshold*100))+".cm").getAbsolutePath();;
			RIGraph consensusGraph = ga.getConsensusGraph(consensusThreshold);
			System.out.printf("Writing consensus graph at CCT %2.0f to file %s \n",consensusThreshold*100,consGraphFile);
			consensusGraph.write_graph_to_file(consGraphFile);
			
			if (benchmark) {
				PredEval eval = consensusGraph.evaluatePrediction(originalGraph);
				System.out.println("## Prediction with CCT: "+consensusThreshold*100+"%");
				//eval.print();
				System.out.println("Number of native contacts:    "+eval.original);
				System.out.println("Number of predicted contacts: "+eval.predicted);
				System.out.printf("Accuracy: %4.3f\n", eval.accuracy);
				System.out.printf("Coverage: %4.3f\n", eval.coverage);
			}

		}
		
		// reconstruct
		//TODO ideally we would like to be able to reconstruct using several consensus graphs: Ca+Cg+C/Cg, at the moment only implemented reconstruction with one graph
		if (reconstruct) {
			RIGraph[] reconstGraphs = {ga.getConsensusGraph(consensusThresholds[0])};
			TinkerRunner tr = new TinkerRunner(TINKER_BIN_DIR,FORCEFIELD_FILE);
			Pdb pdb = tr.reconstruct(targetSeq, reconstGraphs, NUMBER_TINKER_MODELS);
			String outpdbfile = new File(outDir,basename+".reconstructed.pdb").getAbsolutePath();
			pdb.dump2pdbfile(outpdbfile);
		}
	}

}

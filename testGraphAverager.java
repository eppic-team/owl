import graphAveraging.GraphAverager;

import java.io.*;
import java.util.*;
import proteinstructure.*;
import tinker.TinkerError;
import tinker.TinkerRunner;

public class testGraphAverager {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// Plan:
		// use getopt to parse the following parameters
		// - original structure (optional, for benchmarking)
		// - file with list of structures for averaging (required)
		// - contact type & distance cutoff (both required)
		// - consensus edge threshold (defaults to 0.5?)
		// - run tinker or not
		// - parameters for tinker model picking (optional)
		// - number of tinker models to output (only if run tinker is set, defaults to one)
		// - output benchmark result (graph comparison or structure comparison, only if original structure is specified)
		// - output file for resulting graph (optional)
		
		// take pdb files from command line (later: use getopt)
		// create trivial alignment (check alignment)
		// create graphAverager
		// create consensus graph
		// compare with graph for native structure (first argument)
		
		// do tinker reconstruction
		// compare with native structure
		
		// benchmarking:
		// - use score of best model vs. native (either graph based or structure based)
		// - use matlab to optimize input parameters
		// - create matlab executable for optimizing command line parameters:
		//   * specify command line to be run and constant parameters
		//   * specify parameter ranges to be optimized (discrete numerical, continuous numerical with steps, discrete set)
		//   * specify regexp to parse result from output (stdout or file)
		//   * specify optimization procedure (use default matlab ones, e.g. each parameter separate, genetic algorithm, simul. annealing, ...)
		//   * specify convergence threshold or timeout
		//   * run on cluster?
		
		// Input options:
		// - edgeType
		// - distCutoff
		// - targetPdb
		// - targetChainCode
		// - list of predictions files
		// - edge threshold
		// - min seq separation (for benchmarking and reconstruction)
		
		// Output options:
		// -c	write consensus graph
		// -w	write weighted average graph
		// -t	write target graph
		// -p	write result pdb file
		// -a	write accuracy of result vs. native
		// -o	write coverage of results vs. native
		// -r	write rmsd of result vs. native
		// -b	write table of benchmarking results
		// -v	verbose output
		
		// read command line parameters
		if(args.length < 5) {
			System.out.println("Usage: testGraphAverager <targetPdbFile> <ChainCode> <ListOfPredictionFiles> <edgeThreshold> <minSeqSep>");
			System.exit(1);
		}
		
		// constants
		String maxClusterExecutable = "/project/StruPPi/bin/maxcluster";
		String TINKERBINDIR = "/project/StruPPi/Software/tinker/bin";
		String PRMFILE = "/project/StruPPi/Software/tinker/amber/amber99.prm";
		
		// input variables
		File targetFile = new File(args[0]);
		File listFile = new File(args[2]);
		String chainCode = args[1];
		String thresholdStr = args[3];
		String minSeqSepStr = args[4];
		
		// parameters (should all become command line parameters eventually)
		String contactType = "Cb";												// edge type for graph generation (and reconstruction)
		double distCutoff = 8.0;												// distance cutuff for graph generation (and reconstruction)
		double graphAveragingThreshold = Double.parseDouble(thresholdStr); 		// take edges which occur in at least this fraction of models
		int minSeqSep = Integer.parseInt(minSeqSepStr);							// consider only edges with at least this sequence separation
		String outConsensusGraphFile = "consensusGraph.cm";						// name of consensus graph file
		String outAverageGraphFile = "averageGraph.cm";							// name of average graph file
		String outTargetGraphFile = "targetGraph.cm";							// name of target graph file
		String outStructFile = "predictedStructure.pdb";						// name of structure file
		int numberOfTinkerModels = 100;											// number of models to be reconstructed by Tinker
		
		boolean verbose = false;												// if false, only results will be written to stdout
		boolean outputConsensusGraph = false;									// whether to write the resulting consensus graph to a file
		boolean outputAverageGraph = false;										// whether to write the weighted average graph to a file
		boolean outputTargetGraph = false;										// whether to write the target graph to a file
		boolean outputSingleLineResult = false;									// whether to output graph measures of consensus graph
		boolean outputFullResultTable = false;									// whether to output a table with graph measures for all inputs
		boolean doReconstruct = true;											// whether to use tinker to create a 3D prediction
		boolean outputSingleLine3dResult = true;
		boolean outputPredictedStructure = false;								// whether to write reconstructed structure to a file
		
		// local variables
		String targetFileName = targetFile.getAbsolutePath();					// target file name
		Pdb target = null;														// target structure
		RIGraph targetGraph = null;												// target graph
		Vector<String> templateFileNames = new Vector<String>();				// vector of template file names
		Vector<RIGraph> models = new Vector<RIGraph>();								// vector of template graphs
		RIGraph averageGraph;														// weighted graph with fraction of occurrance for each edge
		RIGraph consensusGraph;													// consensus graph after applying threshold
		
		if(!targetFile.canRead()) {
			System.err.println("Can not read from file " + targetFileName);
			System.exit(1);
		}
		if(!listFile.canRead()) {
			System.err.println("Can not read from file " + listFile.getAbsolutePath());
			System.exit(1);
		}
		
		// read target graph
		if(verbose) System.out.println("Reading input files...");
		try {
			target = new PdbfilePdb(targetFile.getAbsolutePath());
			target.load(chainCode);
			targetGraph = target.get_graph(contactType, distCutoff);
		} catch(PdbLoadError e) {
			System.err.println("Error while trying to load pdb data from file " + targetFile.getAbsolutePath()+", specific error: "+e.getMessage());
			System.exit(1);
		}
		
		// read list of predictions
		try {
			BufferedReader in = new BufferedReader(new FileReader(listFile));
			String line;
			File file;
			Pdb pdb;
			RIGraph graph;
			while ((line =  in.readLine()  ) != null) {
				file = new File(line);
				if(!file.canRead()) {
					System.err.println("File " + line + " not found.");
					System.exit(1);
				} else {
					templateFileNames.add(line);
					try {
						pdb = new PdbfilePdb(file.getAbsolutePath());
						pdb.load(Pdb.NULL_CHAIN_CODE);
						graph = pdb.get_graph(contactType, distCutoff);
						models.add(graph);
						if(verbose) System.out.print(".");
					} catch(PdbLoadError e) {
						System.err.println("Error while trying to load pdb data from file " + file.getAbsolutePath()+", specific error: "+e.getMessage());
						System.exit(1);
					}
				}
			}
			in.close();

		} catch(FileNotFoundException e) {
			System.err.println("File " + listFile.getAbsolutePath() + " not found.");
			System.exit(1);
		} catch(IOException e) {
			System.err.println("Error reading from file " + listFile.getAbsolutePath());
			System.exit(1);
		}
		if(verbose) System.out.println();
		
		// create trivial alignment and template graphs
		TreeMap<String,String> sequences = new TreeMap<String, String>();
		TreeMap<String, RIGraph> templateGraphs = new TreeMap<String, RIGraph>();
		sequences.put(targetFile.getAbsolutePath(), targetGraph.getSequence());
		for(int i=0; i < models.size(); i++) {
			sequences.put(templateFileNames.get(i), models.get(i).getSequence());
			templateGraphs.put(templateFileNames.get(i), models.get(i));
		}
		Alignment al = null;
		try {
			al = new Alignment(sequences);
		} catch (AlignmentConstructionError e) {
			System.err.println("Could not create alignment: " + e.getMessage());
		}
		
		// create GraphAverager
		GraphAverager grav = new GraphAverager(targetGraph.getSequence(), al, templateGraphs, targetFile.getAbsolutePath());
		if(verbose) System.out.println("Calculating average...");
		consensusGraph = grav.doAveraging(graphAveragingThreshold);
		averageGraph = grav.getAverageGraph();

		// compare consensus graph with target
		PredEval eval = consensusGraph.evaluatePrediction(targetGraph,minSeqSep);	
		
		// generate result table (consensus vs. individual models)
		ArrayList<PredEval> evals = new ArrayList<PredEval>();
		String targetTitle = "Average";
		eval.title = targetTitle;
		evals.add(eval);
		for (int i = 0; i < models.size(); i++) {
			RIGraph g = models.get(i);
			String title = templateFileNames.get(i).split("/")[templateFileNames.get(i).split("/").length-1];
			eval = g.evaluatePrediction(targetGraph,minSeqSep);
			eval.title = title;
			evals.add(eval);
		}
		// sort results
		Collections.sort(evals, new Comparator<PredEval>() {
			public int compare(PredEval e, PredEval e2) {
				return -((new Double(e.accuracy+e.coverage)).compareTo(new Double(e2.accuracy+e2.coverage)));
			}
		});
		// find best model accuracy and coverage (excluding average graph)
		double bestAcc;
		double bestCov;
		if(evals.get(1).title!=targetTitle) {
			bestAcc = evals.get(1).accuracy;
			bestCov = evals.get(1).coverage;
		} else {
			bestAcc = evals.get(2).accuracy;
			bestCov = evals.get(2).coverage;			
		}
		
		// print full result table
		if(outputFullResultTable) {
			// TODO: add gdt to native order by gdt (and number of models in header)
			System.out.printf("Results (ContactType=%s DistCutoff=%.1f AveragingThreshold=%s MinSeqSep=%d):\n", contactType, distCutoff, thresholdStr,minSeqSep);
			System.out.println();
			System.out.print(" # "); PredEval.printHeaders();
			int c = 1;
			for(PredEval e:evals) {
				if(e.title == targetTitle) System.out.println();
				System.out.printf("%2d ",c++);
				e.printRow();
				if(e.title == targetTitle) System.out.println();
			}
		}
		
		// print accuracy and coverage of consensus graph and best model vs. target graph
		if(outputSingleLineResult) {
			System.out.printf("%.2f\t%.2f\t%.2f\t%.2f\n",eval.accuracy, eval.coverage, bestAcc, bestCov); // picked_gdt, best_pred_gdt, best_templ_gdt
		}
		
		// write resulting consensus graph to file
		if(outputConsensusGraph) {
			try {
				consensusGraph.write_graph_to_file(outConsensusGraphFile);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		// write weighted average graph to file
		if(outputAverageGraph) {
			try {
				averageGraph.write_graph_to_file(outAverageGraphFile);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}			
		}
		
		// write target graph to file
		if(outputTargetGraph) {
			try{
				targetGraph.write_graph_to_file(outTargetGraphFile);

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		// run reconstruction for consensus graph
		if(doReconstruct) {
						
			RIGraph[] graphs = {consensusGraph};
			String sequence = targetGraph.getSequence();
			
			TinkerRunner tr = null;
			try {
				tr = new TinkerRunner(TINKERBINDIR, PRMFILE);
			} catch (FileNotFoundException e3) {
				System.err.println("Couldn't find tinker bin dir "+TINKERBINDIR+". Exiting");
				System.exit(1);
			}
			
			try {
				Pdb resultPdb = tr.reconstruct(sequence, graphs, numberOfTinkerModels);
				if(outputPredictedStructure) {
					resultPdb.dump2pdbfile(outStructFile);
					System.err.println("Output of reconstruction written to " + outStructFile);
				}
			} catch (IOException e) {
				System.err.println("Error while running Tinker reconstruction: " + e.getMessage());
			} catch (TinkerError e) {
				System.err.println("Error while running Tinker reconstruction: " + e.getMessage());
			}	
			
			// benchmark 3D prediction

			// compare all tinker models with real
			double[] tinkerModelGdts = null;
			try {
				tinkerModelGdts = tr.getGdtsToNative(targetFileName, maxClusterExecutable);
			} catch(IOException e) {
				System.err.println("Error while calculating GDT scores: " + e.getMessage());
			}
			// find best gdt from all tinker models
			double bestGdt = 0;
			//int bestModel = 0;
			for (int i = 0; i < tinkerModelGdts.length; i++) {
				if(tinkerModelGdts[i] > bestGdt) {
					//bestModel = i;
					bestGdt = tinkerModelGdts[i];
				}
			}
			// get gdt of picked tinker model (picked by least bound violations)
			int pickedIdx = tr.pickByLeastBoundViols();
			double pickedGdt = tinkerModelGdts[pickedIdx-1];
			
			// get template gdt scores
			//double[] templateGdtScores = new double[templateFileNames.size()];
			double bestTemplateGdt = 0;
			MaxClusterRunner maxCluster = null;		
			try {
				maxCluster = new MaxClusterRunner(maxClusterExecutable);
				for(String tpl:templateFileNames) {
					double gdtScore = maxCluster.calculateGdt(tpl, targetFileName);
					if(gdtScore > bestTemplateGdt) {
						bestTemplateGdt = gdtScore;
					}
				}
			} catch(IOException e) {
				System.err.println("Error running maxcluster:" + e.getMessage());
			}
			if(outputSingleLine3dResult) {
				// row output: best template gdt, best predicted gdt, picked predicted gdt
				System.out.printf("%6.3f\t%6.3f\t%6.3f\n", bestGdt, pickedGdt, bestTemplateGdt);
			}
			
		}
	}
}

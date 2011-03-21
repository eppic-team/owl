package owl.cccp;

import owl.graphAveraging.GraphAverager;
import owl.graphAveraging.GraphAveragerException;
import gnu.getopt.Getopt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;

import owl.core.util.FileFormatException;
import owl.core.structure.Pdb;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbfilePdb;
import owl.core.structure.graphs.GraphComparisonResult;
import owl.core.structure.graphs.RIGEnsemble;
import owl.core.structure.graphs.RIGraph;

import owl.core.sequence.Sequence;
import owl.core.util.MySQLConnection;

public class benchmarkCCCP {

	/*------------------------------ constants ------------------------------*/
	
	// file paths
	public static final String modelDir = "/project/StruPPi/CASP7/decoys/all_full/";
	public static final String seqDir = "/project/StruPPi/CASP7/targets/";
	public static final String nativeDir = "/project/StruPPi/CASP7/answers/baker/CASP7_natives";
	
	// casp constants
	public static final String DEFAULT_CONTACT_TYPE = "Cb";
	public static final double DEFAULT_DIST_CUTOFF = 8.0;
	
	// database constants
	public static final String DB_HOST = "talyn";
	public static final String DB_NAME = "casp7";	
	public static final String TARGET_TABLE = "casp7.targets";
	public static final String table = "casp7.mod_qual_est";
	
	// benchmarking constants
	public static final int	LR_CONTACTS = 24;	// contacts with at least this sequence separation are considered long range
	public static final int ALL_CONTACTS = 1;	// sequence separation without filtering
	
	/*--------------------------- type definitions --------------------------*/
	public class ContGraphComparisonResult {
		
		// comparison against native
		public double acc;
		public double cov;
		public double accLR;
		public double covLR;
		
		// comparison against input templates
		public int rank;
		public double bestAcc;
		public double bestCov;
		public String bestModel;
		
		public int rankLR;
		public double bestAccLR;
		public double bestCovLR;
		public String bestModelLR;
	}
	
	
	
	// benchmark target difficulty prediction
	// - Input: List of targets
	// - Interm: score / classification per target
	// - Compare to average/median GDT / blast score / Casp classification
	// - Result: Correlation, Scatterplot
	public static void benchmarkDiffPred(File targetListFile) {
		// The first version of this function simply returns a vector of consensus scores
		// which can be evaluated using matlab or R. The second version could compute and
		// return the correlation coefficient directly given the list of targets and the
		// vector of GDT scores. This function could then be run in a loop to optimize
		// certain parameters.
	
		try {
			BufferedReader in = new BufferedReader(new FileReader(targetListFile));
			String line;
			while((line=in.readLine()) != null) {
				// find target directory
				File dir = new File(modelDir,line);
				if(!dir.isDirectory()) {
					System.err.println(dir.getPath() + " does not exist. Skipping.");
					continue;
				}
				// load sequence
				String parent = line.substring(0, 5);
				File seqFile = new File(seqDir, parent + ".fa");
				if(!seqFile.canRead()) {
					System.err.println(seqFile.getPath() + " can not be read. Skipping.");
					continue;
				}
				Sequence seq = new Sequence();
				try {
					seq.readFromFastaFile(seqFile);
				} catch (FileFormatException e) {
					System.err.println("Failed to read Fasta file " + seqFile.getPath() + ":" + e.getMessage() + ". Skipping.");
					continue;
				}

				System.out.printf("%s\t%s\t%s\t%s\t", line, parent, dir.getPath(), seqFile.getPath());
				
				// load RIGEnsemble
				RIGEnsemble rigs = new RIGEnsemble();
				int filesRead = rigs.loadFromFileList(dir, seq);
				GraphAverager ga = null;
				try {
					ga = new GraphAverager(rigs);
				} catch (GraphAveragerException e) {
					System.err.println("Could not create graphAverager: " + e.getMessage());
				}
				double consensus = ga.getEnsembleConsensusScore();
				
				System.out.printf("%d\t%.2f\n", filesRead, consensus);
			}		
		} catch (FileNotFoundException e) {
			System.err.println("Could not find file with targets " + targetListFile.getPath());
		} catch (IOException e) {
			System.err.println("Error reading from target file:" + e.getMessage());
		}
	}
	
	// benchmark model quality prediction (whole structure)
	// - Input: List of targets
	// - Interm: score per model
	// - Compare to GDT / RMSD
	// - Result: Correlation per target, Correlation over all models, Scatterplots
	public static void benchmarkModQualEstimation(File targetListFile, MySQLConnection conn) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(targetListFile));
			String line;
			while((line=in.readLine()) != null) {
				estimateQuality(line,conn);
			}		
		} catch (FileNotFoundException e) {
			System.err.println("Could not find file with targets " + targetListFile.getPath());
		} catch (IOException e) {
			System.err.println("Error reading from target file:" + e.getMessage());
		}
	}
	
	// Estimates the quality of each model for a given target and return results as an array.
	// The correlation of this vector to the real model quality (e.g. GDT score) then gives
	// the quality of the estimation. This can be investigated per target or on average.
	// Returns a map (model_id -> quality_estimate) or null if an error occured.
	public static HashMap<String, Double> estimateQuality(String target, MySQLConnection conn) {
		HashMap<String, Double> result = new HashMap<String, Double>();
		// find target directory
		File dir = new File(modelDir,target);
		if(!dir.isDirectory()) {
			System.err.println("Target directory " + dir.getPath() + " does not exist. Skipping.");
			return null;
		}
		// load sequence
		String parent = target.substring(0, 5);
		File seqFile = new File(seqDir, parent + ".fa");
		if(!seqFile.canRead()) {
			System.err.println(seqFile.getPath() + " can not be read. Skipping.");
			return null;
		}
		Sequence seq = new Sequence();
		try {
			seq.readFromFastaFile(seqFile);
		} catch (IOException e) {
			System.err.println(seqFile.getPath() + " can not be read. Skipping.");
			return null;
		} catch (FileFormatException e) {
			System.err.println("Format error in " + seqFile.getPath() + ":" + e.getMessage() + ". Skipping.");
			return null;
		}
		
		System.out.println(target);
		
		// load RIGEnsemble
		RIGEnsemble rigs = new RIGEnsemble();
		GraphAverager ga = null;
		int filesRead = 0;
		try {
			filesRead = rigs.loadFromFileList(dir, seq);
			ga = new GraphAverager(rigs);
		} catch (IOException e) {
			System.err.println("Target directory " + dir.getPath() + " does not exist. Skipping.");
			return null;
		} catch (GraphAveragerException e) {
			System.err.println("Could not create graphAverager: " + e.getMessage());
			return null;
		}
		
		// for each model, estimate quality
		for (int i = 0; i < rigs.getEnsembleSize(); i++) {
			RIGraph rig = rigs.getRIG(i);
			String filename = rigs.getFileName(i);
			String gaTag = String.format("%03d", i); // this is how the graphs are identified in the graphAverager TODO: use filename?
			double cons = ga.getConsensusScore(gaTag, true, true);
			result.put(filename, cons);
			//System.out.printf("%s\t%d\t%d\t%d\t%.3f\n", new File(filename).getName(), rig.getObsLength(), rig.getFullLength(), rig.getEdgeCount(), cons);
			writeQualityEstimateToDatabase(conn, new File(filename).getName(), rig.getObsLength(), rig.getFullLength(), rig.getEdgeCount(), cons);
		}
		System.out.println(filesRead);
		
		return result;
	}
	
	/*---------------------------- private methods --------------------------*/
		
	public static double getCPscore(RIGraph pred, RIGraph real) {
		GraphComparisonResult pe = pred.evaluatePrediction(real);		
		return 0.5*(pe.accuracy + pe.coverage);
	}
	
	public double getCPLRscore(RIGraph pred, RIGraph real) {
		return 0;
	}
	
	public double getQPscore(double[] pred, double[] real) {
		return 0;
	}
		
	public double get3Dscore(Pdb pred, Pdb real) {
		return 0;
	}
	
	public ContGraphComparisonResult evalContactPrediction(RIGraph consensusGraph, RIGraph targetGraph, RIGEnsemble e) {
		ContGraphComparisonResult result = new ContGraphComparisonResult();
		
		// -------------------- All contacts -----------------------------
		int minSeqSep = ALL_CONTACTS;
		
		// compare consensus graph with target
		GraphComparisonResult eval = consensusGraph.evaluatePrediction(targetGraph,minSeqSep);
		result.acc = eval.accuracy;
		result.cov = eval.coverage;
		
		// generate result table (consensus vs. individual models)
		ArrayList<GraphComparisonResult> evals = new ArrayList<GraphComparisonResult>();
		String targetTitle = "Prediction";
		eval.title = targetTitle;
		evals.add(eval);
		RIGraph[] rigs = e.getRIGs();
		String[] filenames = e.getFilenames();
		for (int i = 0; i < rigs.length; i++) {
			RIGraph g = rigs[i];
			String title = filenames[i].split("/")[filenames[i].split("/").length-1];	// basename
			eval = g.evaluatePrediction(targetGraph,minSeqSep);
			eval.title = title;
			evals.add(eval);
		}
		// sort results
		Collections.sort(evals, new Comparator<GraphComparisonResult>() {
			public int compare(GraphComparisonResult e, GraphComparisonResult e2) {
				return -((new Double(e.accuracy+e.coverage)).compareTo(new Double(e2.accuracy+e2.coverage)));
			}
		});
		
		// find ranking
		int rank = 0;
		for (int i = 0; i < evals.size(); i++) {
			if(evals.get(i).title.equals(targetTitle)) {
				rank = i+1;
				break;
			}
		}
		result.rank = rank;
		
		// find best model accuracy and coverage (excluding average graph)
		double bestAcc;
		double bestCov;
		String bestModel;
		if(!evals.get(1).title.equals(targetTitle)) {
			bestAcc = evals.get(1).accuracy;
			bestCov = evals.get(1).coverage;
			bestModel = evals.get(1).title;
		} else {
			bestAcc = evals.get(2).accuracy;
			bestCov = evals.get(2).coverage;
			bestModel = evals.get(2).title;
		}
		result.bestAcc = bestAcc;
		result.bestCov = bestCov;
		result.bestModel = bestModel;
		

		// -------------------- Only long range -----------------------------
		
		minSeqSep = LR_CONTACTS;
		
		// compare consensus graph with target
		eval = consensusGraph.evaluatePrediction(targetGraph,minSeqSep);
		result.accLR = eval.accuracy;
		result.covLR = eval.coverage;
		
		// generate result table (consensus vs. individual models)
		evals = new ArrayList<GraphComparisonResult>();
		targetTitle = "Prediction";
		eval.title = targetTitle;
		evals.add(eval);
		rigs = e.getRIGs();
		filenames = e.getFilenames();
		for (int i = 0; i < rigs.length; i++) {
			RIGraph g = rigs[i];
			String title = filenames[i].split("/")[filenames[i].split("/").length-1];	// basename
			eval = g.evaluatePrediction(targetGraph,minSeqSep);
			eval.title = title;
			evals.add(eval);
		}
		// sort results
		Collections.sort(evals, new Comparator<GraphComparisonResult>() {
			public int compare(GraphComparisonResult e, GraphComparisonResult e2) {
				return -((new Double(e.accuracy+e.coverage)).compareTo(new Double(e2.accuracy+e2.coverage)));
			}
		});
		
		// find ranking
		rank = 0;
		for (int i = 0; i < evals.size(); i++) {
			if(evals.get(i).title.equals(targetTitle)) {
				rank = i+1;
				break;
			}
		}
		result.rankLR = rank;
		
		// find best model accuracy and coverage (excluding average graph)
		if(!evals.get(1).title.equals(targetTitle)) {
			bestAcc = evals.get(1).accuracy;
			bestCov = evals.get(1).coverage;
			bestModel = evals.get(1).title;
		} else {
			bestAcc = evals.get(2).accuracy;
			bestCov = evals.get(2).coverage;
			bestModel = evals.get(2).title;
		}
		result.bestAccLR = bestAcc;
		result.bestCovLR = bestCov;
		result.bestModelLR = bestModel;
		
		return result;
	}
	
	private static void writeQualityEstimateToDatabase(MySQLConnection conn, String model, int obsLen, int fullLen, int edgeCount, double qualityEstimate) {
		String sql = String.format("INSERT INTO %s VALUES(\"%s\",%d,%d,%d,%.3f)", table, model, obsLen, fullLen, edgeCount, qualityEstimate);
		try {
			conn.executeSql(sql);
		} catch (SQLException e) {
			System.err.println("Error executing SQL query: " + e.getMessage());
		}
	}
	
	// benchmark model quality prediction (per residue)
	// - Input: List of targets
	// - Interm: score per residue
	// - Compare to distance from native?
	// - Result: Correlation per model, per target and overall, Scatterplot
	
	// benchmark contact prediction
	// - Input: List of targets
	// - Interm: Acc/Cov per target and SR/LR
	// - Result: Avg acc+cov, Ranking (our prediction vs. inputs)
	
	// benchmark 3D prediction
	// - Input: List of targets
	// - Interm: RMSD/GDT to native for each target
	// - Result: Avg GDT/RMSD (our prediction vs. inputs)
	
	public void runBenchmarks(String[] targetList, boolean diff, boolean mq, boolean mqpr, boolean cp, boolean ts, CccpParams params) {
		String[] targets = getAllTargets();
		
		System.out.printf("Target\tnModels\tLength\tnCont\tnPred\t\tacc\tcov\tacc+cov\trank\tbest ac\t\tLR:acc\tLR:cov\tLR:a+c\tLR:rank\tbest ac\t\tbest model\tbest model LR\n");
		int numTargets = 0;
		int sumRanks = 0;
		int sumRanksLR = 0;
		double sumAccCov = 0;
		double sumAccCovLR = 0;
		for(String target:targets) {
			File nativeFile = new File(nativeDir, target + ".pdb");
			File seqFile = new File(seqDir, target + ".fa");
			File targetModDir = new File(modelDir, target);

			// get target sequence
			Sequence seq = new Sequence();
			try {
				seq.readFromFastaFile(seqFile);
			} catch(IOException e) {
				System.err.println("An IOException occured while reading file " + seqFile.getAbsolutePath() + ": " + e.getMessage());
			} catch(FileFormatException e) {
				System.err.println("A FastFileFormatError occured while reading file " + seqFile.getAbsolutePath() + ": " + e.getMessage());
			}
			int seqLen = seq.getLength();
			
			// get native structure and graph
			Pdb nativePdb = new PdbfilePdb(nativeFile.getAbsolutePath());
			try {
				nativePdb.load("A");
			} catch (PdbLoadException e1) {
				System.err.println("Failed to load native structure " + nativeFile.getAbsolutePath() + ". Skipping target.");
				continue;
			}
			RIGraph nativeGraph = nativePdb.getRIGraph(DEFAULT_CONTACT_TYPE, DEFAULT_DIST_CUTOFF);
			int numNativeCont = nativeGraph.getEdgeCount();
			
			// do contact prediction
			int numPredCont = 0;
			int nModels = 0;
			double acc,cov,accLR,covLR;
			double accCov = 0, accCovLR = 0, bestAccCov, bestAccCovLR;
			String bestModel, bestModelLR;
			int rank = 0, rankLR = 0;
			ContGraphComparisonResult eval;
			try {
				ConsensusPredictor consPred = new ConsensusPredictor(targetModDir, seqFile, params);
				RIGraph predGraph = consPred.getGraphWithPredictedNumContacts();
				eval = evalContactPrediction(predGraph, nativeGraph, consPred.getTemplateEnsemble());
				
				// evaluate
				nModels = consPred.getNumModels();
				numPredCont = predGraph.getEdgeCount();
				
				rank = eval.rank;
				acc = eval.acc;
				cov = eval.cov;
				accCov = 0.5 * (eval.acc+eval.cov);
				bestAccCov = 0.5 * (eval.bestAcc + eval.bestCov);
				bestModel = eval.bestModel;
				
				rankLR = eval.rankLR;
				accLR = eval.accLR;
				covLR = eval.covLR;
				accCovLR = 0.5 * (eval.accLR + eval.covLR);
				bestAccCovLR = 0.5 * (eval.bestAccLR + eval.bestCovLR);
				bestModelLR = eval.bestModelLR;
				
				// sum scores
				numTargets++;
				sumRanks += rank;
				sumAccCov += accCov;
				sumRanksLR += rankLR;
				sumAccCovLR += accCovLR;
			} catch (IOException e) {
				System.out.println("An IOException occured for target " + target + ": " + e.getMessage() + ". Skipping target.");
				continue;
			}			
			System.out.printf("%s\t%4d\t%4d\t%5d\t%5d\t\t%5.3f\t%5.3f\t%5.3f\t%3d\t%5.3f\t\t%5.3f\t%5.3f\t%5.3f\t%3d\t%5.3f\t\t%s\t%s\n", target, nModels, seqLen, numNativeCont, numPredCont, acc, cov, accCov, rank, bestAccCov, accLR, covLR, accCovLR, rankLR, bestAccCovLR, bestModel, bestModelLR);
		}
		System.out.printf("Number of targets:    %4d\n",targets.length);
		System.out.printf("Targets evaluated:    %4d\n",numTargets);
		System.out.printf("Mean acc+cov:         %5.3f\n",sumAccCov / numTargets);
		System.out.printf("Mean rank:            %4.0f\n",1.0 * sumRanks / numTargets);
		System.out.printf("Mean acc+cov LR:      %5.3f\n",sumAccCovLR / numTargets);
		System.out.printf("Mean rank LR:         %4.0f\n",1.0 * sumRanksLR / numTargets);
		
		
	}
	
	/**
	 * Returns an array of all Casp7 targets for which we have native structures and which were not cancelled.
	 */
	private static String[] getAllTargets() {
		LinkedList<String> targets = new LinkedList<String>();
		String[] ta = new String[0];
		
		try {
		MySQLConnection conn = new MySQLConnection(DB_HOST, DB_NAME);
		// load targets from database
		String query = "SELECT target_id FROM " + TARGET_TABLE + " WHERE has_native=true";
		Statement st = conn.createStatement();
		ResultSet rs = st.executeQuery(query);
		while(rs.next()) {
			targets.add(rs.getString(1));
		}
		} catch (SQLException e) {
			System.err.println("An SQL error occured: " + e.getMessage());
		}
		return targets.toArray(ta);
	}
	
	public enum BenchmarkMode {diff, mq, mqpr, cp, ts};
	
	public static void main(String[] args) {
//		if(args.length < 1) {
//			System.out.println("Usage: benchmarkCCCP <option>");
//			System.out.println("-d  target difficulty estimation");
//			System.out.println("-w  whole-structure model quality prediction");
//			System.out.println("-r  per-residue model quality prediction");
//			System.out.println("-c  contact prediction");
//			System.out.println("-t  tertiary structure prediction");			
//			System.exit(1);
//		}
//		String opt = args[0];
//		BenchmarkMode mode = null;
//		if(opt.equals("-d")) {
//			mode = BenchmarkMode.diff;
//		} else
//		if(opt.equals("-w")) {
//			mode = BenchmarkMode.mq;
//		} else 
//		if(opt.equals("-r")) {
//			mode = BenchmarkMode.mqpr;
//		} else
//		if(opt.equals("-c")) {
//			mode = BenchmarkMode.cp;
//		} else 
//		if(opt.equals("-t")) {
//			mode = BenchmarkMode.ts;		
//		} else {
//			System.err.println("Unknown option " + opt + ".");
//			System.exit(1);
//		}
		
		if(args.length < 1) {
			System.out.println("Usage: benchmarkCCCP <option(s)>");
			System.out.println("-v  		plain vanilla, use defaults");
			System.out.println("-f  		use first models only");
			System.out.println("-p <val> 	filter out poor predictions");
			System.out.println("-m <val> 	predict more edges");
			System.out.println("-l <val> 	predict less edges");			
			System.exit(1);
		}
		
		String programName = benchmarkCCCP.class.getName();
		Getopt g = new Getopt(programName, args, "vfp:m:l:");
		int c;
		double q = 0;
		double f = 0;
				
		CccpParams params = new CccpParams();
		
		while((c = g.getopt()) != -1) {
			switch(c) {
			case 'v': break;
			case 'f': params.useOnlyFirstModels = true; break;
			case 'l':
			case 'm': params.useQuantileEdges = true; q = Double.parseDouble(g.getOptarg()); break;
			case 'p': params.filterByBestConsensus = true; f = Double.parseDouble(g.getOptarg()); break;
			default:
				System.exit(1);
			}
		}
		if(q < 0 || q > 1) {
			System.err.println("Quantile needs to be between 0 and 1. Exiting.");
			System.exit(1);
		}
		if(f < 0 || f > 1) {
			System.err.println("Filtering precentage needs to be between 0 and 1. Exiting.");
			System.exit(1);
		}
		
		for(String opt:args) {
			if(opt.equals("-v")) {
				break;
			} else
			if(opt.equals("-f")) {
				params.useOnlyFirstModels = true;
			} else 
			if(opt.equals("-p")) {
				params.filterByBestConsensus = true;
				params.keepBestConsensusPercentage = 0.75;
			} else
			if(opt.equals("-m")) {
				params.useQuantileEdges = true;
				params.numContactsQuantile = 0.75;
			} else 
			if(opt.equals("-l")) {
				params.useQuantileEdges = true;
				params.numContactsQuantile = 0.25;
			} else {
				System.err.println("Unknown option " + opt + ".");
				System.exit(1);
			}			
		}
		
		// running benchmark
		benchmarkCCCP bm = new benchmarkCCCP();
		bm.runBenchmarks(getAllTargets(), true, true, true, true, true, params);
		
//		switch(mode) {
//		case diff: 
//			System.err.println("#Running benchmark...");
//			String targetList = "/project/StruPPi/CASP7/targets/single_dom_targets.txt";
//			benchmarkDiffPred(new File(targetList));
//			System.err.println("#done.");
//			break;
//			
//		case mq:
//			targetList = "/project/StruPPi/CASP7/targets/single_dom_targets.txt";
//			
//			// connect to database
//			MySQLConnection conn;
//			try {
//				conn = new MySQLConnection("white","casp7");
//				//HashMap<String, Double> ret = estimateQuality("T0283", conn);	// single target, results -> DB
//				benchmarkModQualEstimation(new File(targetList), conn);				// multiple targets, resutls -> DB
//				System.out.println("done.");
//			} catch (SQLException e) {
//				System.err.println("Could not connect to database: " + e.getMessage());
//				e.printStackTrace();
//			}
//			break;
//			
//		case mqpr:
//			break;
//			
//		case cp:
//			break;
//			
//		case ts:
//			break;
//		}

	}
	
}

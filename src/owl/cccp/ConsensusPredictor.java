package owl.cccp;

import edu.uci.ics.jung.graph.util.Pair;
import owl.graphAveraging.ConsensusSquare;
import owl.graphAveraging.GraphAverager;
import owl.graphAveraging.GraphAveragerException;
//import graphAveraging.PhiPsiAverager;
//import graphAveraging.ConsensusSquare;

import java.io.*;
import java.util.*;

import org.ggf.drmaa.DrmaaException;

import owl.core.util.FileFormatException;
import owl.core.structure.PdbChain;
import owl.core.structure.graphs.RIGEnsemble;
import owl.core.structure.graphs.RIGraph;
import owl.core.structure.features.SecondaryStructure;

import owl.core.sequence.Sequence;
import owl.core.runners.tinker.TinkerError;
import owl.core.runners.tinker.TinkerRunner;
import owl.core.util.Goodies;

/**
 * This class provides methods for the different predictions that can be made based on our
 * consensus contact strategy given a set of server predictions for a particular Casp target.
 * The main method then allows to run the predictions from command line and produce the
 * output which will then be submitted to Casp.
 * @author stehr
 *
 */
@SuppressWarnings("unused")
public class ConsensusPredictor {

	/*------------------------------ constants ------------------------------*/

	// tinker constants
	private static final String FORCEFIELD_FILE = 		"/project/StruPPi/Software/tinker/amber/amber99.prm";
	private static final String TINKER_BIN_DIR = 		"/project/StruPPi/Software/tinker/bin";
	//private static final String TINKER_BIN_DIR = 		"/project/StruPPi/Software/tinker/bin-modified/large";
	private static final String DSSP_EXECUTABLE = "/project/StruPPi/Software/dssp/dsspcmbi";
	private static final String DSSP_PARAMETERS = "--";
	
	// casp constants
	//private static final String CASP_GROUP_CODE = "3715-2099-9446";	// CASP8
	private static final String CASP_GROUP_CODE = "6915-8127-6629"; // CASP9
	private static final String CASP_METHOD_QA = "Contact based consensus score compared to server models";
	private static final String CASP_METHOD_RR = "Consensus contact prediction based on server models";
	private static final String CASP_METHOD_TS = "Distance geometry structure based on consensus contact prediction";
	private static final int CASP_MAX_CONTACT_LINES = 32000; // casp format does not allow more lines in RR files
	
	// benchmark constants
	public static final String CASP7_modelDir = "/project/StruPPi/CASP7/decoys/all_full/";
	public static final String CASP7_seqDir = "/project/StruPPi/CASP7/targets/";
	public static final String CASP7_nativeDir = "/project/StruPPi/CASP7/answers/baker/CASP7_natives";
	
	// other constants
	public static final int 	NUM_REPORTED_PARENTS = 10;
	
	/*--------------------------- type definitions --------------------------*/
	public enum CaspTargetClass{HA_TBM, TBM, TBM_FM, FM, Other};
	public class GlobalQualities extends TreeMap<String, Double> {private static final long serialVersionUID = 1L;};
	public class LocalQualities extends TreeMap<String, Double[]> {private static final long serialVersionUID = 1L;};
	
	/*--------------------------- member variables --------------------------*/
	
	File modelDir;
	File seqFile;
	Sequence seq;
	RIGEnsemble ensemble;
	public GraphAverager ga;
	CccpParams params;							// selected parameters which should be accessible to automatic optimization
	Map<String,Double> consensusScores = null;	// consensus scores for each template, calculated on demand
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Create a new consensusPredictor for a given target sequence and directory with predictions.
	 * @param models
	 * @param seq
	 */
	public ConsensusPredictor(File modelDir, File seqFile, CccpParams params) throws IOException {
		
		this.modelDir = modelDir;
		this.seqFile = seqFile;
		this.params = params;
		
		// read sequence
		seq = new Sequence();
		try {
			seq.readFromFastaFile(seqFile);
		} catch (FileFormatException e) {
			throw new IOException(e);
		}
		System.err.println(seq.getName());		
		
		// create RIGEnsemble
		ensemble = new RIGEnsemble();
		if(params.useDsspConsensusSecondaryStructure) ensemble.useDssp(DSSP_EXECUTABLE, DSSP_PARAMETERS);
		if(params.useOnlyFirstModels) ensemble.loadOnlyFirstModels();
		ensemble.loadFromFileList(modelDir, seq);
		
		// create graphAverager		
		try {
		this.ga = new GraphAverager(ensemble);
		} catch (GraphAveragerException e) {
			System.err.println("Could not create graphAverager: " + e.getMessage());
		}
		
		// filter out bad models
		// TODO: Move to graph averager?
		if(params.filterByBestConsensus) {

			consensusScores = getConsensusScores();
			int toKeep = (int) Math.ceil(params.keepBestConsensusPercentage * ensemble.getEnsembleSize());
			int toTrash = ensemble.getEnsembleSize() - toKeep;
			LinkedHashMap<String,RIGraph> goodGraphs = new LinkedHashMap<String, RIGraph>();
			for (int i = 0; i < ensemble.getEnsembleSize(); i++) {
				String filename = ensemble.getFileName(i);
				RIGraph rig = ensemble.getRIG(i);
				goodGraphs.put(filename, rig);
			}
			LinkedHashMap<String, Double> sortedScores;
			int trashed = 0;
			
			sortedScores = Goodies.sortMapByValue(consensusScores, Goodies.ASCENDING);
			
			for(String filename:sortedScores.keySet()) {
				goodGraphs.remove(filename);
				trashed++;
				if(trashed >= toTrash) break;
			}
			
			// now create new ensemble and new graph averager for filtered set
			ensemble = new RIGEnsemble();
			ensemble.loadFromGraphMap(goodGraphs);
			try {
				this.ga = new GraphAverager(ensemble);
			} catch (GraphAveragerException e) {
				System.err.println("Could not create graphAverager for filtered set: " + e.getMessage());
			}
		}		
	}

	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * @return the number of files in the model directory (disregarding file type)
	 */
	public int getNumFilesInModelDir() {
		return modelDir.list().length;
	}
	
	/**
	 * @return the number of input models
	 */
	public int getNumModels() {
		return ensemble.getEnsembleSize();
	}
	
	/**
	 * @return the name of the input sequence
	 */
	public String getTargetName() {
		return seq.getName();
	}
	
	/**
	 * Returns the ensemble of templates used for this ConsensusPredictor
	 * @return
	 */
	public RIGEnsemble getTemplateEnsemble() {
		return this.ensemble;
	}
	
	/* --- methods to create Casp relevant results ---*/
	/**
	 * @return the consensus graph (= contact prediction)
	 */
	public RIGraph getConsensusGraph() {
		RIGraph cons = null;
		if(ga != null) {
			cons = ga.getConsensusGraph(params.consensusThreshold);
		} else {
			System.err.println("Error. GraphAverager is null!");
		}
		return cons;
	}
	
	/**
	 * Returns a graph with the predicted number of contacts, all set to weight 1.
	 * The strategy for predictin the number of contacts (e.g. median) depends on the settings in the params object.
	 * @return a graph with the top n contacts with highest consensus, where n is the predicted number of contacts 
	 */
	public RIGraph getGraphWithPredictedNumContacts() {
		RIGraph cons = null;
		if(ga != null) {
			int n = 0;
			if(params.useQuantileEdges) {
				n = ga.getQuantNumContacts(params.numContactsQuantile);
			} else {
				n = ga.getMedNumContacts();
			}
			cons = ga.getGraphWithTopContacts(n, true);	// set edge weights to one
		} else {
			System.err.println("Error. GraphAverager is null!");
		}
		return cons;
	}
	
	/**
	 * @return a graph with the top 32000 (or less) contacts with highest consensus (limit for Casp RR format) 
	 */
	public RIGraph getAvGraphWithTop32000Contacts() {
		RIGraph av = null;
		if(ga != null) {
			av = ga.getGraphWithTopContacts(CASP_MAX_CONTACT_LINES, false);	// leave edge weights as they are
		} else {
			System.err.println("Error. GraphAverager is null!");
		}
		return av;
	}	
	
	/**
	 * @return an estimation of the quality of the predicted residues in Angstrom deviation
	 */
	public double[] getPredQualEstimPerResidue() {
		return null;
	}
	
	/**
	 * @return an estimation of the probability that the predicted contact is a true contact
	 */
	public Map<Pair<Integer>, Double> getPredQualEstimPerEdge() {
		return null;
	}
	
	/**
	 * Calculates a 3D model based on the given contact graph.
	 * @param predGraph the contact prediction based on which the structure will be reconstructed
	 * @return a prediction of the 3D structure
	 */	
	public PdbChain get3DPrediction(RIGraph predGraph) throws IOException {
		return get3DPrediction(predGraph, null);
	}
	
	/**
	 * Calculates a 3D model based on the given contact graph using the given secondary structure
	 * annotation as constraints for the reconstruction.
	 * @param predGraph the contact prediction based on which the structure will be reconstructed
	 * @param ss the secondary structure object used for the reconstruction, needs to match predGraph's sequence
	 * @return a prediction of the 3D structure
	 */
	public PdbChain get3DPrediction(RIGraph predGraph, SecondaryStructure ss) throws IOException {
		RIGraph[] rigs = {predGraph};
		TinkerRunner tr = new TinkerRunner(TINKER_BIN_DIR,FORCEFIELD_FILE);
		TreeMap<Integer,ConsensusSquare> angleConstraints = null;
		if(ss != null) {
			if(!predGraph.getSequence().equals(ss.getSequence())) {
				throw new IOException("Sequences of contact graph and secondary structure annotation do not match");
			}
			angleConstraints = ss.getPhiPsiConstraints();
		}
		//PhiPsiAverager ppa = new PhiPsiAverager(null,null);	// TODO: Create alignment, pass PDBs
		//TreeMap<Integer,ConsensusSquare> angleConstraints = ppa.getConsensusPhiPsi(0.5, 20);
		PdbChain pdb = null;
		try {
			if(params.tinkerFastMode) {
				pdb = tr.reconstructFast(seq.getSeq(), rigs, angleConstraints, params.forceTransOmega, params.numTinkerModels, params.useParallelTinker);
			} else {
				pdb = tr.reconstruct(seq.getSeq(), rigs, angleConstraints, params.forceTransOmega, params.numTinkerModels, 100, 1, params.useParallelTinker);
			}
		} catch (TinkerError e) {
			throw new IOException(e.getMessage());
		} catch (FileFormatException e) {
			throw new IOException(e.getMessage());
		}
		return pdb;
	}
	
	/**
	 * @return estimations for the overall model quality of the input models
	 */
	public GlobalQualities getModQualEstims() {
		
		GlobalQualities result = new GlobalQualities();
		
		if(consensusScores == null) {					// use cached scores
			consensusScores = getConsensusScores();		// or generate from scratch
		}
		// for each model, map consensus score to quality estimation
		for(String filename:consensusScores.keySet()) {
			double cons = consensusScores.get(filename);

			// normalize to [0;1]
			double predGdt = (params.cons2qualSlope * cons + params.cons2qualOffset) / 100;
			if(predGdt < 0) {
				predGdt = 0;
			}
			if(predGdt > 1) {
				predGdt = 1;
			}
			result.put(filename, predGdt);		
		}
		return result;
	}
	
	/** 
	 * Returns the consensus scores for all models.
	 */
	private LinkedHashMap<String,Double> getConsensusScores() {
		LinkedHashMap<String,Double> scores = new LinkedHashMap<String,Double>();
		// for each model, estimate quality
		for (int i = 0; i < ensemble.getEnsembleSize(); i++) {
			String gaTag = String.format("%03d", i); // this is how the graphs identified in the graphAverager TODO: use filename?
			double cons = ga.getConsensusScore(gaTag, true, true);
			String filename = ensemble.getFileName(i);
			scores.put(filename, cons);
		}
		return scores;
	}
	
	/**
	 * @return estimations for the per-residue quality of the input models
	 */
	public LocalQualities getModQualEstimsPerResidue() {
		return null;
	}
	
	/* --- methods to create other informative results ---*/

	/**
	 * Run all predictions and output an informative summary report of the results.
	 */
	public void predictionReport() {
		// Report for target ...
		// Sequence length
		// Total models in directory
		// Models which could be read
		// Minimum contacts
		// Maximum contacts
		// Average contacts
		// Median contacts
		// The overall consensus is, predicted average GDT is, I think this target is a ... target
		// Contact prediction written to ...
		// 3D prediction written to ...
		// Quality prediction written to ...
	}
	
	/**
	 * @return the predicted casp classification of the target
	 */
	public CaspTargetClass getTargetClassification() {
		return CaspTargetClass.Other;
	}
	
	/**
	 * @return the estimated difficulty of the target in average GDT_TS score
	 */
	public double getTargetDiffEstimation() {
		double consensus = ga.getEnsembleConsensusScore();
		double predGdt = 0.064*consensus - 25;
		if(predGdt < 0 || predGdt > 100) {
			System.err.println("This target has an unusually high or low ensemble consensus score of " + consensus + " and should be checked manually.");
		}
		return predGdt;
	}
	
	/**
	 * @return the average graph for visual inspection of the consensus
	 */
	public RIGraph getAverageGraph() {
		return ga.getAverageGraph();
	}
	
	/**
	 * @return an estimation of the quality of the predicted graph/structure between 0 and 1
	 */
	public double getPredictionQualityEstim() {
		return 0.0;
	}
	
	/**
	 * @return the consensus secondary structure prediction
	 */
	public SecondaryStructure getConsensusSecStrucPrediction() {
		return ga.getConsensusSecondaryStructure(params.consensusSecondaryStructureThreshold);
	}
	
	/* --- methods to write results to Casp and other formats --- */

	public void writeSecStrucPred(SecondaryStructure pred, File outFile) {
	}
	
	public void writeAverageGraph(RIGraph avGraph, File outFile) throws IOException {
		avGraph.writeToFile(outFile.getAbsolutePath());
	}
	
	public void writeCaspRRFile(RIGraph consensusGraph, File outFile) throws IOException {
		int targetNum = Integer.parseInt(this.getTargetName().substring(1));
		consensusGraph.setTargetNum(targetNum);
		consensusGraph.setCaspModelNum(1);
		consensusGraph.setAuthorStr(CASP_GROUP_CODE);
		consensusGraph.setMethodStr(CASP_METHOD_RR);
		// TODO: add AUTHOR and METHOD tags
		consensusGraph.writeToCaspRRFile(outFile.getAbsolutePath());
	}
	
	public void writeCaspQA1File(GlobalQualities glQuals, File outFile) throws IOException {
		
		// write a casp quality prediction file without per-residue scores (QA1)
		PrintWriter out = new PrintWriter(outFile);
		out.printf("PFRMAT QA\n");
		out.printf("TARGET %s\n", this.getTargetName());
		out.printf("AUTHOR %s\n", CASP_GROUP_CODE);
		out.printf("METHOD %s\n", CASP_METHOD_QA);
		out.printf("MODEL 1\n");
		out.printf("QMODE 1\n");
		for(String model:glQuals.keySet()) {
			double qualPred = glQuals.get(model);
			out.printf("%s %.3f\n", new File(model).getName(), qualPred);
		}
		out.printf("END");
		out.close();
	}
	
	public void writeCaspQA2File(GlobalQualities glQuals, LocalQualities locQuals, File outFile) {
	}
	
	public void writeCaspTSFile(PdbChain pred, File outFile) throws IOException {
		int targetNum = Integer.parseInt(this.getTargetName().substring(1)); // note: target tag must be like T0100, otherwise this fails!
		pred.setTargetNum(targetNum);
		pred.setCaspModelNum(1);
		pred.setCaspAuthorStr(CASP_GROUP_CODE);
		pred.setCaspMethodStr(CASP_METHOD_TS);
		pred.writeToCaspTSFile(outFile);	// TODO: add AUTHOR and METHOD tags
	}
	
	/**
	 * Writes a list of frequently used parents to stdout.
	 */
	public void writeParents() {
		System.out.println("Frequently used parents:");
		Map<String,Integer> parents = this.ga.getParentFrequencies();
		int maxNum = NUM_REPORTED_PARENTS;
		int c = 0;
		for(String p:parents.keySet()) {
			System.out.printf("%s\t%d\n", p, parents.get(p));
			if(++c >= maxNum) break;
		}
	}
	
	/**
	 * Writes all results to files using default files names.
	 */
	public void writeAllCaspFiles(File outDir) throws IOException {
		
		String target = this.getTargetName();
		System.out.println("============");
		System.out.println("Target " + target);
		System.out.println("============");
		System.out.println("Target length: " + this.seq.getLength());
		System.out.println("Number of files in model dir: " + this.getNumFilesInModelDir());		
		System.out.println("Number of models read: " + ensemble.getEnsembleSize());
		System.out.println("Min num contacts: " + ga.getMinNumContacts());
		System.out.println("Max num contacts: " + ga.getMaxNumContacts());
		System.out.println("Med num contacts: " + ga.getMedNumContacts());
		System.out.println("Avg num contacts: " + ga.getAvgNumContacts());
		System.out.printf("Estimated average GDT for this target: %.2f\n", this.getTargetDiffEstimation());
		System.out.println("");
		writeParents();
		System.out.println("");
		
		// Consensus secondary structure
		System.out.println("Consensus secondary structure:");
		SecondaryStructure consSS = getConsensusSecStrucPrediction();
		consSS.print();
		System.out.println("");
		
		// Quality prediction
		File qsFile = new File(outDir, target + ".qs");
		GlobalQualities glQual = this.getModQualEstims();
		this.writeCaspQA1File(glQual, qsFile);
		System.out.println("Model quality predictions written to " + qsFile.getPath());
		
		// Average graph (for visual inspection of ensemble)
		File avFile = new File(outDir, target + ".cm");
		RIGraph avGraph = this.getAverageGraph();
		this.writeAverageGraph(avGraph, avFile);
		System.out.println("Average graph written to " + avFile.getPath());

		// Contact prediction with average graph limited to 32000 edges (limit for Casp RR format)
		File rrFile = new File(outDir, target + ".rr");
		RIGraph contPred = this.getAvGraphWithTop32000Contacts();
		this.writeCaspRRFile(contPred, rrFile);
		System.out.println("Contact prediction written to " + rrFile.getPath());
		
//		// Contact prediction with consensus threshold
//		File rrFile = new File(outDir, target + ".rr");
//		RIGraph contPred = this.getConsensusGraph();
//		this.writeCaspRRFile(contPred, rrFile);
//		System.out.println("Contact prediction written to " + rrFile.getPath());
		
		// Contact prediction with given number of contacts
		File cmFile2 = new File(outDir, target + ".2.cm");
		RIGraph contPred2 = this.getGraphWithPredictedNumContacts();
		this.writeAverageGraph(contPred2, cmFile2);
		System.out.println("Edges predicted: " + contPred2.getEdgeCount());

		// 3D prediction, optionally w/ secondary structure
		if(params.tinkerFastMode) {
			System.out.println("Reconstructing " + params.numTinkerModels + " models in fast mode...");
		} else {
			System.out.println("Reconstructing " + params.numTinkerModels + " models...");			
		}
		File tsFile = new File(outDir, target + ".ts");
		PdbChain struPred = null;
		if(params.useDsspConsensusSecondaryStructure) {
			struPred = this.get3DPrediction(contPred2, consSS);	// use file with median number of contacts and consensus secondary structure
		} else {
			struPred = this.get3DPrediction(contPred2);	// without secondary structure
		}
		this.writeCaspTSFile(struPred, tsFile);
		System.out.println("3D prediction written to " + tsFile.getPath());	
		
//		// write additional 3D prediction w/o secondary structure
//		if(params.useDsspConsensusSecondaryStructure) {
//			struPred = null;	// allow garbage collection
//			if(params.tinkerFastMode) {
//				System.out.println("Reconstructing " + params.numTinkerModels + " models in fast mode...");
//			} else {
//				System.out.println("Reconstructing " + params.numTinkerModels + " models...");			
//			}
//			File tsFile2 = new File(outDir, target + ".2.ts");
//			PdbChain struPred2 = this.get3DPrediction(contPred2);	// use file with median number of contacts and consensus secondary structure
//			this.writeCaspTSFile(struPred2, tsFile2);
//			System.out.println("3D prediction written to " + tsFile2.getPath());
//		}
		
	}
	
	/**
	 * Writes all results (where it makes sense) to database using default tables and columns.
	 */
	public void writeAllToDb() {
	}
		
	/*--------------------------------- main --------------------------------*/
	
	/**
	 * Input:   A directory with server predictions, a target sequence file, an output dir
	 * Output:  1. A target difficulty estimation (in expected average GDT_TS and target classification) (for our information and for internal use)
	 * 			2. An estimated quality (GDT and RMSD) for each of the models
	 * 			3. One Casp quality prediction file per model containing the above estimation plus estimations per residue
	 * 			4. A contact prediction (Casp RR format, including probabilities for each contact)
	 * 			5. A 3D prediction (Casp TS format, with confidence values for each residue?)
	 * 			6. An estimation of the quality of the 3D prediction
	 * 			7. A consensus secondary structure prediction (just for our information and for manual modelling)
	 * 			8. An average graph as .cm or .png (just for our information and for manual modelling)
	 * 			Bonus: A prediction of who will win Casp (based on our quality estimation)
	 * @param args
	 */
	public static void main(String[] args) {
		// Input: seqFile, modelDir, outDir
		// create ConsensusPredictor instance
		// call all evaluation functions
		// write all results to Casp files
		// output some guesses of how difficult this target is and how good the prediction is likely to be
		// with optional parameters, perform only some of the predictions (default is all)
		
		if(args.length < 3) {
			System.out.println("Usage: ConsensusPredictor <seqFile> <modelDir> <outDir>");
			System.exit(1);
		}
		
		File seqFile = new File(args[0]);
		File modelDir = new File(args[1]);
		File outDir = new File(args[2]);
		
//		String target = "T0354";
//		File seqFile = new File(CASP7_seqDir, target + ".fa");
//		File modelDir = new File(CASP7_modelDir);
//		File outDir = new File(".");
		
		// this will cause the program to fail fast if SGE is n/a
		try {
			Object dummy = new org.ggf.drmaa.SimpleJobTemplate().getJobName();
		} catch (DrmaaException e1) {
			// do nothing
		}

		
		try {
			ConsensusPredictor cp = new ConsensusPredictor(modelDir, seqFile, new CccpParams());
			cp.writeAllCaspFiles(outDir);
		} catch(IOException e) {
			System.out.println("Error processing target: " + e.getMessage());
		}
		
	}

}

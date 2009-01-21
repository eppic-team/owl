package embed;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;
import java.util.TreeSet;

import edu.uci.ics.jung.graph.util.Pair;
import graphAveraging.GraphAverager;
import graphAveraging.GraphAveragerError;

import proteinstructure.IntPairSet;
import proteinstructure.Pdb;
import proteinstructure.PdbasePdb;
import proteinstructure.RIGEdge;
import proteinstructure.RIGEnsemble;
import proteinstructure.RIGraph;

public class Distiller {
	
	private static final int DIAGONALS_TO_SKIP = 3;
	
	
	/**
	 * Class to store a IntPairSet (set of contacts) together with its score.
	 * Can be sorted based on the score order. 
	 * @author duarte
	 *
	 */
	private class SetScore implements Comparable<SetScore>{
		public double score;
		public IntPairSet set;
		public SetScore(double score, IntPairSet set) {
			this.score = score;
			this.set = set;
		}
		public int compareTo(SetScore o) {
			if (this.score<o.score) return -1;
			if (this.score>o.score) return 1;
			return 0;
		}
	}
	
	
	/*-------------------------- members  -------------------------------*/
	
	private Bound[][] cmapBounds;
	private int totalNumberContacts;
	private double finalContacts;
	private int cmapSize;
	private RIGraph rig;
	
	private ArrayList<SetScore> allSampledSets;
	
	/*----------------------- constructors  -----------------------------*/
	
	/**
	 * 
	 * @param rig
	 * @param finalContacts
	 */
	public Distiller(RIGraph rig, double finalContacts) {
		this.rig = rig;
		this.cmapBounds = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		this.totalNumberContacts = rig.getEdgeCount();
		this.finalContacts = finalContacts;
		this.cmapSize = cmapBounds.length;
	}
	
	/*------------------------- privates  -------------------------------*/
	
	/**
	 * Calculates the sum square deviation of the given all pairs bounds matrix 
	 * against the contact map bounds measuring the deviation of the upper bounds 
	 * to the ones in the contact map: sum((u'i-ui)2)
	 * @param bounds
	 * @return
	 */
	private double computeSumSquareDeviation(Bound[][] bounds) {
		double sumSquareDev = 0;
		for (int i=0;i<cmapSize;i++) {
			for (int j=i+1;j<cmapSize;j++) {
				if (cmapBounds[i][j]!=null) {
					sumSquareDev += ((bounds[i][j].upper-cmapBounds[i][j].upper)*(bounds[i][j].upper-cmapBounds[i][j].upper));
				}
			}
		}
		return sumSquareDev;

	}
	
	/**
	 * Measures the error function for the given sparseBounds matrix (a subset of the contact map) by first 
	 * inferring bounds for all pairs through the triangle inequality and then 
	 * measuring how well the all pairs matrix fits the contact map.  
	 * Thus a high error value means the given sparseBounds matrix has low information content
	 * about the rest of the contact map and low error value that the matrix has high 
	 * information content  
	 * @param sparseBounds
	 * @return
	 */
	private double getError(Bound[][] sparseBounds) {
		// infer bounds for all pairs through triangle inequality
		Bound[][] bounds = inferAllBounds(sparseBounds);
		return (computeSumSquareDeviation(bounds))/(double) cmapSize; 
	}
	
	/**
	 * Infer bounds for all pairs from a spare bounds matrix via the triangle inequality
	 * @param sparseBounds
	 * @return
	 */
	private Bound[][] inferAllBounds(Bound[][] sparseBounds) {
		BoundsSmoother bs = new BoundsSmoother(sparseBounds);
		return bs.getInitialBoundsAllPairs();
	}

	/**
	 * Gets all pairs of the given bounds matrix except for the first 
	 * minSeqSeparation diagonals, the indices of the pairs are the matrix 
	 * indices (0 to cmapSize-1)
	 * @param bounds the bounds matrix
	 * @param minSeqSeparation only contacts above this sequence separation will be returned 
	 * @return
	 */
	private IntPairSet getPairs(Bound[][] bounds, int minSeqSeparation) {
		IntPairSet allCMPairs = new IntPairSet();
		for (int i=0;i<cmapSize;i++) {
			for (int j=i+1+minSeqSeparation;j<cmapSize;j++) {
				if (bounds[i][j]!=null)  {
					allCMPairs.add(new Pair<Integer>(i,j));
				}
			}
		}
		return allCMPairs;
	}
	
	/**
	 * Randomly samples one subset of numSampledContacts contacts from contacts with 
	 * sequence separation above {@link #DIAGONALS_TO_SKIP} 
	 * @param numSampledContacts
	 * @return
	 */
	private Bound[][] sampleSubset(int numSampledContacts) {
		Bound[][] subset = new Bound[cmapSize][cmapSize];
		
		Random rand = new Random();
		IntPairSet allcmpairs = getPairs(cmapBounds, DIAGONALS_TO_SKIP);
		TreeSet<Integer> sampledIndices = new TreeSet<Integer>();
		for (int i=0;i<numSampledContacts;i++) {
			// the TreeSet takes care of not having duplicates indices
			// add() returns false if the value is a duplicate, so if it returns false 
			// we loop until the return value is true, i.e. we have a new index
			while (!sampledIndices.add(rand.nextInt(allcmpairs.size())));
		}
		// the iterator returns the indices in order since it is from a TreeSet
		Iterator<Integer> indexIterator = sampledIndices.iterator();
		int index = 0;
		int sampledIndex = indexIterator.next();
		for (Pair<Integer> pair:allcmpairs) {
			if (index==sampledIndex) {
				Bound b = cmapBounds[pair.getFirst()][pair.getSecond()];
				subset[pair.getFirst()][pair.getSecond()] = new Bound(b.lower, b.upper);
				if (indexIterator.hasNext()) {
					sampledIndex = indexIterator.next();
				} else {
					break;
				}
			}
			index++;
		}
		Reconstructer.addBackboneRestraints(subset);
		return subset;
	}
	
	/*-------------------------- publics  -------------------------------*/
	
	/**
	 * Samples random subsets of the contact map scoring them and storing all results in the 
	 * allSampledSets ArrayList (sorted by scores). Get the max and min scoring subsets with {@link #getMaxErrorSetScore()} and 
	 * {@link #getMinErrorSetScore()}. Write all scores out with {@link #writeScores(PrintStream)}
	 * @param numSamples
	 */
	public void distillRandomSampling(int numSamples) {
		 
		allSampledSets = new ArrayList<SetScore>();		

		int numSampledContacts = (int)(totalNumberContacts*finalContacts);
		System.out.println("Sampling "+numSamples+" subsets of "+numSampledContacts+" contacts ");
		
		for (int i=0;i<numSamples;i++) {
			Bound[][] boundsSubset = sampleSubset(numSampledContacts);
			double error = getError(boundsSubset);
		
			allSampledSets.add(new SetScore(error,getPairs(boundsSubset, DIAGONALS_TO_SKIP)));			
			//System.out.printf("sample "+i+": "+getPairs(boundsSubset,1).size()+" contacts. Error value: %4.2e\n",error);
			if (i!=0 && i%10==0) System.out.print(".");
			if (i!=0 && i%1000==0) System.out.println(" "+i);
		}
				
		Collections.sort(allSampledSets);
	}
	
	/**
	 * Returns the SetScore with the minimum score of all samples
	 * @return
	 */
	public SetScore getMinErrorSetScore() {
		return allSampledSets.get(0);
	}

	/**
	 * Returns the SetScore with the maximum score of all samples
	 * @return
	 */
	public SetScore getMaxErrorSetScore() {
		return allSampledSets.get(allSampledSets.size()-1);
	}
	
	/**
	 * Writes to the given PrintStream all the scores in a 1 column text file
	 * @param out
	 */
	public void writeScores(PrintStream out) {
		for (SetScore ss:allSampledSets) {
			out.printf("%8.3f\n",ss.score);
		}
	}
	
	/**
	 * Returns a RIGraph containing the union of edges of the best percentileForBestSet percent 
	 * subsets out of all the sampled subsets. The edges are weighted by the fraction of occurrence 
	 * in the subsets.
	 * @param percentileForBestSet
	 * @return
	 */
	public RIGraph getRIGEnsembleForBestSets(double percentileForBestSet) {
		int numberBestSets = (int)(allSampledSets.size()*percentileForBestSet);
		System.out.println("Getting the 1st percentile of the best scores: "+numberBestSets);
		double sumScore = 0;
		RIGEnsemble rigs = new RIGEnsemble(rig.getContactType(), rig.getCutoff());
		for (int i=0;i<numberBestSets;i++) {
			SetScore setScore = allSampledSets.get(i);
			sumScore+=setScore.score;
			rigs.addRIG(createRIGraphFromIntPairSet(rig.getSequence(), setScore.set, rig.getContactType(), rig.getCutoff()));
		}
		System.out.printf("Average error value for best sets: %8.3f\n",(sumScore/(double)numberBestSets));
		GraphAverager ga = null;
		try { 
			ga = new GraphAverager(rigs);
		} catch (GraphAveragerError e) {
			// this shouldn't happen
			System.err.println("Unexpected error while creating the average RIG: "+e.getMessage());
		}
		return ga.getAverageGraph();
	}
		
	/*------------------------- statics  -----------------------------*/
	
	private static RIGraph createRIGraphFromIntPairSet(String sequence, IntPairSet set, String contactType, double distCutoff) {
		RIGraph graph = new RIGraph(sequence);
		for (Pair<Integer> pair:set) {
			graph.addEdge(new RIGEdge(1.0), graph.getNodeFromSerial(pair.getFirst()+1), graph.getNodeFromSerial(pair.getSecond()+1));
		}
		graph.setContactType(contactType);
		graph.setCutoff(distCutoff);
		return graph;
	}
	
	/*-------------------------- main  -------------------------------*/
	
	public static void main(String[] args) throws Exception {
		
		if (args.length<2) {
			System.err.println("Usage: Distiller <num_samples> <out_dir>");
		}
		int numSamples = Integer.parseInt(args[0]);
		String outDir = args[1];

		double finalContactPercent = 0.05;
		double percentileForBestSet = 0.01;
		
		String pdbCode = "1sha";
		String pdbChainCode = "A";
		String ct = "Ca";
		double cutoff = 8.0;

		File essenceFile = new File(outDir,pdbCode+pdbChainCode+"_essence.cm");
		File maxErrorFile = new File(outDir,pdbCode+pdbChainCode+"_maxError.cm");
		File scoresFile = new File(outDir,pdbCode+pdbChainCode+".scores");
		File ensembleFile = new File(outDir,pdbCode+pdbChainCode+"_ensemble_best_"+String.format("%4.2f",percentileForBestSet)+".cm");
		
		Pdb pdb = new PdbasePdb(pdbCode);
		pdb.load(pdbChainCode);
		RIGraph graph = pdb.get_graph(ct, cutoff);
		System.out.println("Total contacts: "+graph.getEdgeCount());
		
		Distiller dist = new Distiller(graph, finalContactPercent);
		
		dist.distillRandomSampling(numSamples);

		SetScore maxErrorSetScore = dist.getMaxErrorSetScore();
		SetScore minErrorSetScore = dist.getMinErrorSetScore();

		System.out.println();
		
		System.out.printf("Max error: %8.3f \n", maxErrorSetScore.score);
		createRIGraphFromIntPairSet(pdb.getSequence(), maxErrorSetScore.set, ct, cutoff).write_graph_to_file(maxErrorFile.getAbsolutePath());

		System.out.printf("Min error: %8.3f \n", minErrorSetScore.score);
		createRIGraphFromIntPairSet(pdb.getSequence(), minErrorSetScore.set, ct, cutoff).write_graph_to_file(essenceFile.getAbsolutePath());
				
		dist.writeScores(new PrintStream(scoresFile));
		
		RIGraph ensembleOfBests = dist.getRIGEnsembleForBestSets(percentileForBestSet);
		ensembleOfBests.write_graph_to_file(ensembleFile.getAbsolutePath());
	}
}

package embed;

//import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;
import java.util.TreeSet;

import edu.uci.ics.jung.graph.util.Pair;

import proteinstructure.IntPairSet;
import proteinstructure.Pdb;
import proteinstructure.PdbasePdb;
import proteinstructure.RIGEdge;
import proteinstructure.RIGraph;

public class Distiller {

	private static final int DIAGONALS_TO_SKIP = 3;
	
	//private RIGraph rig;
	private Bound[][] cmapBounds;
	private int totalNumberContacts;
	private double finalContacts;
	private int cmapSize;
	
	public Distiller(RIGraph rig, double finalContacts) {
		//this.rig = rig;
		this.cmapBounds = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		this.totalNumberContacts = rig.getEdgeCount();
		this.finalContacts = finalContacts;
		this.cmapSize = cmapBounds.length;
	}
	
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
	 * Scores the given sparseBounds matrix (a subset of the contact map) by first 
	 * inferring bounds for all pairs through the triangle inequality and then 
	 * measuring how well the all pairs matrix fits the contact map.  
	 * Thus a high score means the given sparseBounds matrix has high information content
	 * about the rest of the contact map and low score that the matrix has low 
	 * information content  
	 * @param sparseBounds
	 * @return
	 */
	private double score(Bound[][] sparseBounds) {
		// infer bounds for all pairs through triangle inequality
		Bound[][] bounds = inferAllBounds(sparseBounds);
		return 1/computeSumSquareDeviation(bounds); 
	}
	
	private Bound[][] inferAllBounds(Bound[][] sparseBounds) {
		BoundsSmoother bs = new BoundsSmoother(sparseBounds);
		return bs.getInitialBoundsAllPairs();
	}

//	public IntPairSet distillGreedy() {
//		HashMap<Pair<Integer>, Double> scores = new HashMap<Pair<Integer>, Double>();
//		
//		IntPairSet allCMPairs = getPairs(cmapBounds, DIAGONALS_TO_SKIP);
//
//		Pair<Integer> pairToRemove = findHighestScoringEdge(allCMPairs, BoundsSmoother.copyBounds(cmapBounds));
//		//removePair(BoundsSmoother.copyBounds(cmapBounds));
//				
//		
//		return null;
//	}
//	
//	private Pair<Integer> findHighestScoringEdge(IntPairSet pairs, Bound[][] sparseBounds) {
//		// we proceed removing one edge at a time and scoring the resulting all pairs matrix
//		// (inferred via triangle inequality) to see how well it fits the contact map. 
//		// A high score for a matrix (with a removed edge) means good fit to the contact map
//		// Thus removing the edge that gives a matrix with the highest score means we remove
//		// the edge which was least important in getting a good score.
//		Pair<Integer> highestScEdge = null;
//		double maxScore = 0;
//		for (Pair<Integer> pair:pairs) {
//			// we remove the pair
//			removePair(sparseBounds,pair);
//			// and score
//			double score = score(sparseBounds);
//			if (score>maxScore) {
//				maxScore = score;
//				highestScEdge = pair;
//			}
//		}
//		return highestScEdge;
//	}
//	
//	private void removePair(Bound[][] bounds, Pair<Integer> pair) {
//		bounds[pair.getFirst()][pair.getSecond()] = null;
//	}
	
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
	
	public IntPairSet distillRandomSampling(int numSamples) {
		double maxScore = 0;
		int indexMaxScoringSample = -1;
		Bound[][] bestScoringBounds = null;
		for (int i=0;i<numSamples;i++) {
			Bound[][] boundsSubset = sampleSubset();
			double score = score(boundsSubset);
			System.out.printf("sample "+i+": "+getPairs(boundsSubset,1).size()+" contacts. Score: %4.2e\n",score);
			if (score>maxScore) {
				bestScoringBounds = boundsSubset;
				maxScore = score;
				indexMaxScoringSample = i;
			}
		}
		System.out.printf("Best score: %4.2e for sample %d\n", maxScore, indexMaxScoringSample);
		return getPairs(bestScoringBounds, DIAGONALS_TO_SKIP);
	}
	
	private Bound[][] sampleSubset() {
		Bound[][] subset = new Bound[cmapSize][cmapSize];

		
		int numSampledContacts = (int)(totalNumberContacts*finalContacts);
		System.out.println("Sampling "+numSampledContacts+" contacts");

		Random rand = new Random();
		IntPairSet allcmpairs = getPairs(cmapBounds, DIAGONALS_TO_SKIP);
		TreeSet<Integer> sampledIndices = new TreeSet<Integer>();
		for (int i=0;i<numSampledContacts;i++) {
			// the TreeSet takes care of not having duplicates indices
			// add() returns false if the value is a duplicate, so if it returns false 
			// we loop until the return value is true, i.e. we have a new index
			while (!sampledIndices.add(rand.nextInt(allcmpairs.size())));
		}
		// the iterator returns the indices in order as it is from a TreeSet
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
	
	/*-------------------------- main  -------------------------------*/
	
	public static void main(String[] args) throws Exception {
		String pdbCode = "1sha";
		String pdbChainCode = "A";
		String ct = "Ca";
		double cutoff = 8.0;

		Pdb pdb = new PdbasePdb(pdbCode);
		pdb.load(pdbChainCode);
		RIGraph graph = pdb.get_graph(ct, cutoff);
		System.out.println("Total contacts: "+graph.getEdgeCount());
		
		Distiller dist = new Distiller(graph, 0.05);
		//dist.distillGreedy();
		IntPairSet bestSet = dist.distillRandomSampling(1000);
		RIGraph essence = new RIGraph(pdb.getSequence());
		for (Pair<Integer> pair:bestSet) {
			essence.addEdge(new RIGEdge(1.0), essence.getNodeFromSerial(pair.getFirst()+1), essence.getNodeFromSerial(pair.getSecond()+1));
			//System.out.println(pair.getFirst()+"\t"+pair.getSecond());
		}
		essence.write_graph_to_file("/project/StruPPi/jose/embed/"+pdbCode+pdbChainCode+"_essence.cm");
	}
}

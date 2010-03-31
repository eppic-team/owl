package embed;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.Random;
import java.util.TreeSet;

import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbasePdb;
import owl.core.structure.RIGEdge;
import owl.core.structure.RIGEnsemble;
import owl.core.structure.RIGraph;
import owl.core.util.IntPairSet;

import edu.uci.ics.jung.graph.util.Pair;
import graphAveraging.GraphAverager;
import graphAveraging.GraphAveragerError;


public class Distiller {
	
	private static final int DIAGONALS_TO_SKIP = 4;
	protected static final String SCORE_PRINT_FORMAT = "%9.4f";
	
	/*-------------------------- members  -------------------------------*/
	
	private Bound[][] cmapBounds;
	private int totalNumberContacts;
	private int cmapSize;
	private RIGraph rig;
	
	private ArrayList<SetScore> allSampledSets;
	
	/*----------------------- constructors  -----------------------------*/
	
	/**
	 * Constructs a Distiller 
	 * @param rig the residue interaction graph
	 */
	public Distiller(RIGraph rig) {
		this.rig = rig;
		this.cmapBounds = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		this.totalNumberContacts = rig.getEdgeCount();
		this.cmapSize = cmapBounds.length;
	}
	

	/*------------------------- privates  -------------------------------*/

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
	public Bound[][] sampleSubset(int numSampledContacts) {
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
	 * @param contactsToSample fraction of total contacts we want to be selected for the subset 
	 * @throws PdbLoadError 
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundError 
	 */
	public void distillRandomSampling(int numSamples, double contactsToSample) throws PdbCodeNotFoundError, SQLException, PdbLoadError {
		 
		allSampledSets = new ArrayList<SetScore>();		

		int numSampledContacts = (int)(totalNumberContacts*contactsToSample);
		System.out.println("Sampling "+numSamples+" subsets of "+numSampledContacts+" contacts, with sequence separation above "+DIAGONALS_TO_SKIP);
		
		for (int i=0;i<numSamples;i++) {
			Bound[][] boundsSubset = sampleSubset(numSampledContacts);
			double error = Scorer.getCMError(boundsSubset, rig);
		
			// we have to skip the first diagonal that we've added to the random subset
			allSampledSets.add(new SetScore(error,getPairs(boundsSubset, 1)));			
			//System.out.printf("sample "+i+": "+getPairs(boundsSubset,1).size()+" contacts. Error value: %4.2e\n",error);
			if (i!=0 && i%100==0) System.out.print(".");
			if (i!=0 && i%10000==0) System.out.println(" "+i);
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
			out.printf(""+SCORE_PRINT_FORMAT+"\n",ss.score);
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
		System.out.printf("Average error value for best sets: "+SCORE_PRINT_FORMAT+"\n",(sumScore/(double)numberBestSets));
		GraphAverager ga = null;
		try { 
			ga = new GraphAverager(rigs);
		} catch (GraphAveragerError e) {
			// this shouldn't happen
			System.err.println("Unexpected error while creating the average RIG: "+e.getMessage());
		}
		return ga.getAverageGraph();
	}
	
	/**
	 * Writes all sampled subsets and scores to two files in table format with a subsetId column
	 * linking the two tables. Useful to import into a relational database.
	 * @param edgesTableFile
	 * @param scoresTableFile
	 * @param startingID the starting subset id to use
	 * @throws IOException
	 */
	public void writeSubsetsToTables(File edgesTableFile, File scoresTableFile, int startingID) throws IOException {
		PrintWriter pwE = new PrintWriter(edgesTableFile);
		PrintWriter pwS = new PrintWriter(scoresTableFile);
		int i = startingID;
		for (SetScore ss: allSampledSets) {
			pwS.printf("%d\t"+SCORE_PRINT_FORMAT+"\n",i,ss.score);
			for (Pair<Integer> pair:ss.set) {
				pwE.println(i+"\t"+(pair.getFirst()+1)+"\t"+(pair.getSecond()+1));
			}
			i++;
		}
		pwE.close();
		pwS.close();
	}
		
	/*------------------------- statics  -----------------------------*/
	
	public static RIGraph createRIGraphFromIntPairSet(String sequence, IntPairSet set, String contactType, double distCutoff) {
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
		
		if (args.length<5) {
			System.err.println("Usage: Distiller <pdbCode+pdbChainCode> <num_samples> <out_dir> <starting_subset_id> <subset_percent_size>");
			System.exit(1);
		}
		String pdbId = args[0]; 
		int numSamples = Integer.parseInt(args[1]);
		String outDir = args[2];
		int startingID = Integer.parseInt(args[3]);
		double contactsToSample = Double.parseDouble(args[4])/(double)100;
		
		String pdbCode = pdbId.substring(0, 4);
		String pdbChainCode = pdbId.substring(4,5);
		String ct = "Ca";
		double cutoff = 8.0;

		File edgesTableFile = new File(outDir,pdbCode+pdbChainCode+"_"+startingID+".subsets");
		File scoresTableFile = new File(outDir,pdbCode+pdbChainCode+"_"+startingID+".scores");
		
		Pdb pdb = new PdbasePdb(pdbCode);
		pdb.load(pdbChainCode);
		RIGraph graph = pdb.getRIGraph(ct, cutoff);
		System.out.println("Total contacts: "+graph.getEdgeCount());
		
		Distiller dist = new Distiller(graph);
		
		dist.distillRandomSampling(numSamples, contactsToSample);
		
		dist.writeSubsetsToTables(edgesTableFile, scoresTableFile, startingID);
		
		SetScore maxErrorSetScore = dist.getMaxErrorSetScore();
		SetScore minErrorSetScore = dist.getMinErrorSetScore();

		System.out.println();
		
		System.out.printf("Max error: "+SCORE_PRINT_FORMAT+" \n", maxErrorSetScore.score);
		System.out.printf("Min error: "+SCORE_PRINT_FORMAT+" \n", minErrorSetScore.score);	

	}
}

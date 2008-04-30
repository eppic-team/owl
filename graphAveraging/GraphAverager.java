package graphAveraging;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Locale;
import java.util.TreeMap;
import java.util.TreeSet;

import proteinstructure.Alignment;
import proteinstructure.AlignmentConstructionError;
import proteinstructure.IntPairComparator;
import proteinstructure.PairwiseSequenceAlignment;
import proteinstructure.RIGEdge;
import proteinstructure.RIGEnsemble;
import proteinstructure.RIGNode;
import proteinstructure.RIGraph;
import proteinstructure.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;

import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;


/**
 * Provides methods for an ensemble of RIGs: Calculating a consensus graph,
 * and consensus scores for individual members or the whole ensemble. 
 * @author stehr
 * @date 2007-12-19
 */
public class GraphAverager {

	/*--------------------------- member variables --------------------------*/
	
	private Alignment al;
	private TreeMap<String,RIGraph> templateGraphs;	// identified by tag-string
	private String targetTag;		// id of target sequence in alignment
	private String sequence;		// sequence of the final consensus graph
	private String contactType;		// contact type of the final consensus graph
	private double distCutoff;		// cutoff of the final consensus graph
	
	private HashMap<Pair<Integer>,Vote> contactVotes;
	
	/*---------------------------- inner classes ----------------------------*/
	/**
	 * Inner class to store vote counts together with voters.
	 * Used in contactVotes Map
	 */
	private class Vote {
		private int voteCount;
		private TreeSet<String> voters;
		
		public Vote(int voteCount, TreeSet<String> voters) {
			this.voteCount = voteCount;
			this.voters = voters;
		}
		
		public int getVoteCount() {
			return voteCount;
		}
		
		public TreeSet<String> getVoters() {
			return voters;
		}
	}
	
	/**
	 * Extension of RIGEdge for an average graph. So that we can add 
	 * to the RIGEdges the voters information
	 *
	 */
	private class AveragedRIGEdge extends RIGEdge {
		private TreeSet<String> voters;
		
		public AveragedRIGEdge(double weight, TreeSet<String> voters) {
			this.setWeight(weight);
			this.voters = voters;
		}
		
		public TreeSet<String> getVoters() {
			return this.voters;
		}
	}
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Creates a GraphAverager given a RIGEnsemble
	 * @param rigs the RIGEnsemble for which graph will be averaged 
	 */
	public GraphAverager(RIGEnsemble rigs) {
		this.targetTag = Long.toString(new Date().getTime()); // a unique identifier
		
		this.templateGraphs = new TreeMap<String,RIGraph>();
		for (int i=0;i<rigs.getEnsembleSize();i++){
			this.templateGraphs.put(String.format("%03d", i),rigs.getRIG(i));
		}
		this.sequence = rigs.getRIG(0).getSequence();
		RIGraph firstGraph = templateGraphs.get(templateGraphs.firstKey());
		this.contactType = firstGraph.getContactType();
		this.distCutoff = firstGraph.getCutoff();		
		
		TreeMap<String,String> sequences = new TreeMap<String, String>();
		sequences.put(targetTag, sequence);	
		for(String id:templateGraphs.keySet()) {
			sequences.put(id, templateGraphs.get(id).getSequence());
		}
		// create trivial alignment
		try {
			this.al = new Alignment(sequences);
		} catch (AlignmentConstructionError e) {
			System.err.println("Could not create alignment: " + e.getMessage());
		}
		
		checkSequences();	
		countVotes(); // does the averaging by counting the votes and putting them into contactVotes		
	}
	
	/**
	 * Create a graph averager given an alignment of the target sequence to the template graphs.
	 * @param sequence the target sequence
	 * @param al a multiple alignment of the target sequence and the template graphs
	 * @param templateGraphs a collection of template graphs to be averaged
	 * @param targetTag the identifier of the target sequence in the alignment
	 */
	public GraphAverager(String sequence, Alignment al, TreeMap<String,RIGraph> templateGraphs, String targetTag) {
		this.al = al;
		this.templateGraphs = templateGraphs;
		this.targetTag = targetTag;
		this.sequence = sequence;
		RIGraph firstGraph = templateGraphs.get(templateGraphs.firstKey());
		this.contactType = firstGraph.getContactType();
		this.distCutoff = firstGraph.getCutoff();
		
		checkSequences();	
		countVotes(); // does the averaging by counting the votes and putting them into contactVotes
	}
	
	/**
	 * Create a graph averager assuming that the target and template structures all have the same size.
	 * A trivial alignment will be internally created. Fails if sizes do not match.
	 * @param sequence the target sequence
	 * @param templateGraphs a collection of template graphs to be averaged
	 */
	public GraphAverager(String sequence, TreeMap<String,RIGraph> templateGraphs) {
		this.templateGraphs = templateGraphs;
		this.targetTag = Long.toString(new Date().getTime()); // a unique identifier
		this.sequence = sequence;
		RIGraph firstGraph = templateGraphs.get(templateGraphs.firstKey());
		this.contactType = firstGraph.getContactType();
		this.distCutoff = firstGraph.getCutoff();		
		
		TreeMap<String,String> sequences = new TreeMap<String, String>();
		sequences.put(targetTag, sequence);	
		for(String id:templateGraphs.keySet()) {
			sequences.put(id, templateGraphs.get(id).getSequence());
		}
		// create trivial alignment
		try {
			this.al = new Alignment(sequences);
		} catch (AlignmentConstructionError e) {
			System.err.println("Could not create alignment: " + e.getMessage());
		}

		
		checkSequences();	
		countVotes(); // does the averaging by counting the votes and putting them into contactVotes		
	}
	
	/*---------------------------- private methods --------------------------*/
	
	/**
	 * Checks that tags and sequences are consistent between this.al and this.templateGraphs and between this.al  and this.graph/this.targetTag 
	 *
	 */
	private void checkSequences(){
		if (!al.hasTag(targetTag)){
			System.err.println("Alignment doesn't seem to contain the target sequence, check the FASTA tags");
			//TODO throw exception
		}
		for (String tag:templateGraphs.keySet()){
			if (!al.hasTag(tag)){
				System.err.println("Alignment is missing template sequence "+tag+", check the FASTA tags");
				// TODO throw exception
			}
		}
		if (templateGraphs.size()!=al.getNumberOfSequences()-1){
			System.err.println("Number of sequences in alignment is different from number of templates +1 ");
			// TODO throw exception
		}
		if(!al.getSequenceNoGaps(targetTag).equals(this.sequence)) {
			System.err.println("Target sequence in alignment does not match sequence in target graph");
			System.err.println("Trying to align sequences of alignment vs graph: ");
			try {
				PairwiseSequenceAlignment alCheck = new PairwiseSequenceAlignment(this.sequence,al.getSequenceNoGaps(targetTag),"graph","alignment");
				alCheck.printAlignment();
			} catch (PairwiseSequenceAlignmentException e) {
				System.err.println("Error while creating alignment check, can't display an alignment, error: "+e.getMessage()+". The 2 sequences are: ");
				System.err.println("graph:     "+sequence);
				System.err.println("alignment: "+al.getSequenceNoGaps(targetTag));
			}
			// TODO throw exception
		}
		for (String tag:templateGraphs.keySet()){
			if(!al.getSequenceNoGaps(tag).equals(templateGraphs.get(tag).getSequence())) {
				System.err.println("Sequence of template graph "+tag+" does not match sequence in alignment");
				System.err.println("Trying to align sequences of alignment vs graph: ");
				try {
					PairwiseSequenceAlignment alCheck = new PairwiseSequenceAlignment(templateGraphs.get(tag).getSequence(),al.getSequenceNoGaps(tag),"graph","alignment");
					alCheck.printAlignment();
				} catch (PairwiseSequenceAlignmentException e) {
					System.err.println("Error while creating alignment check, can't display an alignment. The 2 sequences are: ");
					System.err.println("graph:     "+templateGraphs.get(tag).getSequence());
					System.err.println("alignment: "+al.getSequenceNoGaps(tag));
				}
				// TODO throw exception
			}			
		}
	}
	
	/**
	 * Counts the votes for each possible alignment edge and puts all the votes in contactVotes Map
	 *
	 */
	private void countVotes() {
		
		contactVotes = new HashMap<Pair<Integer>, Vote>();

		// we get the first graph in templates to see if they are directed or undirected
		boolean directed = templateGraphs.get(templateGraphs.firstKey()).isDirected();
		
		// we go through all positions in the alignment
		for (int i=1; i<=al.getAlignmentLength(); i++){
			for (int j=1; j<=al.getAlignmentLength(); j++) {
 
				if (directed) {
					// in directed case we only want to skip loop edges 
					if (i==j) continue;
				} else {
					// in undirected case we want to skip half of the matrix
					// this way the final contactVotes Map will contain Pair<Integer> always with j>i
					if (i>=j) continue;
				}
				
				int vote = 0;
				TreeSet<String> voters = new TreeSet<String>(); 
				// scanning all templates to see if they have this contact
				for (String tag:templateGraphs.keySet()){			
					RIGraph thisGraph = templateGraphs.get(tag);
					int iSeqIdx = al.al2seq(tag, i);
					int jSeqIdx = al.al2seq(tag, j);
					
					// if each of the ends map to a gap in this sequence we skip it
					if ((iSeqIdx!=-1) && (jSeqIdx!=-1) && thisGraph.containsEdgeIJ(iSeqIdx, jSeqIdx)) {  
						RIGEdge thisGraphCont = thisGraph.getEdgeFromSerials(iSeqIdx, jSeqIdx);
						Pair<RIGNode> pair = thisGraph.getEndpoints(thisGraphCont);
						if (thisGraph.containsEdgeIJ(pair.getFirst().getResidueSerial(), pair.getSecond().getResidueSerial())) {
							vote++;
							voters.add(tag);
						}
					}
				}
				// putting vote in contactVotes Map
				if (vote>0){
					contactVotes.put(new Pair<Integer>(i,j), new Vote(vote,voters));
				}				
			}
		}		
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Returns the number of templates.
	 * @return
	 */
	public int getNumberOfTemplates() {
		return this.templateGraphs.size();
	}
	
	/**
	 * Calculates the consensus graph from the set of template graphs. An edge is contained
	 * in the consensus graph if the fractions of template graphs it is contained in is above
	 * the given threshold.
	 * The output is a new RIGraph object created from the given sequence in the constructor 
	 * to which we add the averaged edges.  
	 * @param threshold the threshold above which an edge is taken to be a consensus edge
	 * @return the new graph
	 */
	public RIGraph getConsensusGraph(double threshold) {
		
		RIGraph graph = new RIGraph(this.sequence);
		graph.setContactType(this.contactType);
		graph.setCutoff(this.distCutoff);
		int numTemplates = templateGraphs.size();
		
		// if vote above threshold we take the contact for our target
		int voteThreshold = (int) Math.ceil((double)numTemplates*threshold); // i.e. round up of 50%, 40% or 30% (depends on threshold given)
		for (Pair<Integer> alignCont:contactVotes.keySet()){
			if (contactVotes.get(alignCont).getVoteCount()>=voteThreshold) {
				int target_i_res = al.al2seq(targetTag,alignCont.getFirst());
				int target_j_res = al.al2seq(targetTag,alignCont.getSecond());
				if (target_i_res!=-1 && target_j_res!=-1) { // we can't add contacts that map to gaps!!
					graph.addEdge(new RIGEdge(), graph.getNodeFromSerial(target_i_res), graph.getNodeFromSerial(target_j_res), EdgeType.UNDIRECTED);
				}
			}
		}
		return graph;
	}
	
	/**
	 * Returns a RIGraph containing the union of edges of the template graphs 
	 * weighted by the fraction of occurrence in the templates.
	 * The sequence of the graph is initialized to the sequence of the template.
	 * @return
	 */
	public RIGraph getAverageGraph() {
		RIGraph graph = new RIGraph(this.sequence);
		graph.setContactType(this.contactType);
		graph.setCutoff(this.distCutoff);
		int numTemplates = templateGraphs.size();
		
		for (Pair<Integer> alignCont:contactVotes.keySet()){
			double weight = 1.0 * contactVotes.get(alignCont).getVoteCount() / numTemplates;
			int target_i_res = al.al2seq(targetTag,alignCont.getFirst());
			int target_j_res = al.al2seq(targetTag,alignCont.getSecond());
			if (target_i_res!=-1 && target_j_res!=-1) { // we can't add contacts that map to gaps!!
				graph.addEdge(new AveragedRIGEdge(weight,contactVotes.get(alignCont).getVoters()), graph.getNodeFromSerial(target_i_res), graph.getNodeFromSerial(target_j_res), EdgeType.UNDIRECTED);				
			}
		}		
		return graph;
	}

	/**
	 * Writes the averaged graph (see {@link #getAverageGraph()} to outfile  
	 * with edge weights and voters: i.e. templates that voted for the edge 
	 * 
	 * @param outfile
	 */
	public void writeAverageGraphWithVoters(String outfile) throws IOException {
		//TODO  We might want in the future to allow this format (with list of voters for each edge) as our aglappe format. 
		// 		At the moment we write it to a file without headers because they are not compatible (although would be 
		// 		quite simple to make them compatible)
		
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		RIGraph graph = getAverageGraph();
		
		// printing column names
		Out.print("#i\tj\tweight");
		for (String tag:templateGraphs.keySet()) {
			Out.print("\t"+tag);
		}
		Out.println();
		
		// we use temp TreeMaps to be able to order the output
		TreeMap<Pair<Integer>,Double> pairs2weights = new TreeMap<Pair<Integer>,Double>(new IntPairComparator());
		TreeMap<Pair<Integer>,TreeSet<String>> pairs2voters = new TreeMap<Pair<Integer>,TreeSet<String>>(new IntPairComparator()); 
		for (RIGEdge cont:graph.getEdges()){
			AveragedRIGEdge avgCont = (AveragedRIGEdge) cont; 
			Pair<RIGNode> pair = graph.getEndpoints(cont);
			int i_resser=pair.getFirst().getResidueSerial();
			int j_resser=pair.getSecond().getResidueSerial();
			double weight=cont.getWeight();
			pairs2weights.put(new Pair<Integer>(i_resser,j_resser),weight);
			pairs2voters.put(new Pair<Integer>(i_resser,j_resser), avgCont.getVoters());
			
		}
		for (Pair<Integer> pair:pairs2weights.keySet()) { 
			Out.printf(Locale.US,pair.getFirst()+"\t"+pair.getSecond()+"\t%6.3f",pairs2weights.get(pair));
			for (String tag:templateGraphs.keySet()) {
				if (pairs2voters.get(pair).contains(tag)) {
					Out.print("\t1");
				} else {
					Out.print("\t0");
				}
			}
			Out.println();
		}
	}
	
	/**
	 * Calculates the consensus score for the graph g with tag t with respect to the current ensemble of graphs.
	 * g has to have the same size (i.e. number of nodes) as the graphs in the ensemble.
	 * This is done by summing over each edge e in g, the number of graphs in the ensemble which also
	 * contain edge e. This value can be normalized by the number of nodes and the size of the ensemble.
	 * The normalization makes the values comparable across different targets.
	 * This score is a rough estimate of model quality. 
	 * @param g2
	 * @return the consensus score or -1 if something went wrong
	 */
	public double getConsensusScore(String t, boolean normalizeByNumNodes, boolean normalizeByNumTemplates) {
		double score = 0;
		RIGraph g = templateGraphs.get(t);
		if(g != null) {
			for(RIGEdge e:g.getEdges()) {
				// TODO: Use index mapping from alignment
				Pair<RIGNode> n = g.getEndpoints(e);
				Pair<Integer> i = new Pair<Integer>(n.getFirst().getResidueSerial(), n.getSecond().getResidueSerial());
				int count = contactVotes.get(i).getVoteCount();
				score = score + count;
			}
		} else return -1;
		if(normalizeByNumNodes) {
			score = score / al.getAlignmentLength();
		}
		if(normalizeByNumTemplates) {
			score = score / templateGraphs.size();
		}
		return score;
	}
	
	/**
	 * Return the consensus score for the whole ensemble (summed over all members).
	 * This is a rough estimate of target difficulty.
	 * @return
	 */
	public double getEnsembleConsensusScore() {
		double score = 0;
		for(String tag:templateGraphs.keySet()) {
			double s = getConsensusScore(tag, true, true);
			if(s > 0) score = score + s;
		}
		return score;
	}
	
	/**
	 * Filters the graphs in this ensemble by the consensus score. All graphs
	 * whith a consensus score equal to or greater than the given threshold are kept,
	 * others are thrown away.
	 * @param minScore the filtering threshold
	 */
	public void filterByConsensusScore(double minScore) {
		// filter
		Iterator<String> it = templateGraphs.keySet().iterator();
		while(it.hasNext()) {
			String tag = it.next();
			double cs = getConsensusScore(tag, true, true);
			if(cs < minScore) {
				// delete from templateGraphs, keep in alignment (to avoid side effects)
				it.remove();
			}
		}
		// update contact votes
		countVotes();
	}
	
	// Run graph averaging from the command line
	public static void main(String[] args) {
		// Take a list of pdb/cm files, calculate the average & evtl. reconstruct using Tinker
		// 0. Create graphs (if necessary)
		// 1. select models (using overall consensus score) [or do this outside?]
		// 2. create consensus graph
		// 3. reconstruct (if necessary)
		
		// define getopt options
	}
}

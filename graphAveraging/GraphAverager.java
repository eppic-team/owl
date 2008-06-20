package graphAveraging;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

import proteinstructure.Alignment;
import proteinstructure.AlignmentConstructionError;
import proteinstructure.IntPairComparator;
import proteinstructure.PairwiseSequenceAlignment;
import proteinstructure.RIGEdge;
import proteinstructure.RIGEnsemble;
import proteinstructure.RIGNode;
import proteinstructure.RIGraph;
import proteinstructure.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import tools.Goodies;

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
	
	private Alignment al;			// alignment containing the graphs below
	private TreeMap<String,RIGraph> templateGraphs;	// template graphs identified by tag-string
	private String targetTag;		// id of target sequence in alignment, or null if no target specified
	private String sequence;		// sequence of the final consensus graph, dummy sequence if targetTag=null
	private String contactType;		// contact type of the final consensus graph
	private double distCutoff;		// cutoff of the final consensus graph
	
	private int[] numContacts;		// sorted array with number of contacts for each template graph (see countVotes)
	private HashMap<Pair<Integer>,Vote> contactVotes; // number of votes for each contact (see countVotes)
	
	/*---------------------------- inner classes ----------------------------*/
	/**
	 * Inner class to store vote counts together with voters.
	 * Used in contactVotes Map
	 */
	private class Vote {
		private int voteCount;			// how many graphs actually contain this contact
		private int potentialVoteCount; // how many graphs could make this contact (because both ends map to non-gaps)
		private TreeSet<String> voters;
		
		public Vote(int voteCount, int potentialVoteCount, TreeSet<String> voters) {
			this.voteCount = voteCount;
			this.potentialVoteCount = potentialVoteCount;
			this.voters = voters;
		}
		
		public int getVoteCount() {
			return voteCount;
		}
		
		public int getPotentialVoteCount() {
			return this.potentialVoteCount;
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
	 * Creates a GraphAverager given a RIGEnsemble. In this case the alignment is
	 * trivial and the target sequence is the same as all the input sequences.
	 * @param rigs the RIGEnsemble for which graph will be averaged 
	 * @throws GraphAveragerError
	 */
	public GraphAverager(RIGEnsemble rigs) throws GraphAveragerError {
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
	 * Create a graph averager given a multiple alignment of the target sequence to the template graphs.
	 * @param al a multiple alignment of the target sequence and the template graphs
	 * @param templateGraphs a collection of template graphs to be averaged
	 * @param targetTag the identifier of the target sequence in the alignment
	 * @throws GraphAveragerError
	 */
	public GraphAverager(Alignment al, TreeMap<String,RIGraph> templateGraphs, String targetTag) throws GraphAveragerError {
		this.al = al;
		this.templateGraphs = templateGraphs;
		this.targetTag = targetTag;
		this.sequence = al.getSequenceNoGaps(targetTag);
		RIGraph firstGraph = templateGraphs.get(templateGraphs.firstKey());
		this.contactType = firstGraph.getContactType();
		this.distCutoff = firstGraph.getCutoff();
		
		checkSequences();	
		countVotes(); // does the averaging by counting the votes and putting them into contactVotes
	}
	
	/**
	 * Create a graph averager for a set of templates with the given alignment. A dummy target sequence will
	 * be created internally which has the length of the alignment. The resulting average or consensus graphs
	 * will then have coordinates corresponding to the columns in the alignment.
	 * @param al
	 * @param templateGraphs
	 * @throws GraphAveragerError
	 */
	public GraphAverager(Alignment al, TreeMap<String,RIGraph> templateGraphs) throws GraphAveragerError {
		this.sequence = makeDummySequence(al.getAlignmentLength());
		this.targetTag = makeDummyTag();
		try {
			this.al = al.copyAndAdd(this.targetTag, this.sequence);
		} catch (AlignmentConstructionError e) {
			throw new GraphAveragerError(e);
		}
		this.templateGraphs = templateGraphs;
		RIGraph firstGraph = templateGraphs.get(templateGraphs.firstKey());
		this.contactType = firstGraph.getContactType();
		this.distCutoff = firstGraph.getCutoff();
		
		checkSequences();	
		countVotes(); // does the averaging by counting the votes and putting them into contactVotes
	}

	
	/*---------------------------- private methods --------------------------*/
	
	/**
	 * @return a dummy sequence of the given length.
	 */
	private String makeDummySequence(int length) {
		char dummyChar = 'X';
		StringBuilder buf = new StringBuilder(length);
		for (int i = 0; i < length; i++) {
			buf.append(dummyChar);
		}
		return buf.toString();
	}
	
	/**
	 * @return a randomly generated sequence tag
	 */
	private String makeDummyTag() {
		return "dummyTag";
	}
	
	/**
	 * Checks that tags and sequences are consistent between this.al and this.templateGraphs and between this.al  and this.graph/this.targetTag 
	 * @throws GraphAveragerError
	 */
	private void checkSequences() throws GraphAveragerError {
		if (!al.hasTag(targetTag)){
			throw new GraphAveragerError("Alignment doesn't seem to contain the target sequence, check the FASTA tags");
		}
		for (String tag:templateGraphs.keySet()){
			if (!al.hasTag(tag)){
				throw new GraphAveragerError("Alignment is missing template sequence "+tag+", check the FASTA tags");
			}
		}
		// we check that at the number of graphs is not bigger than sequences in alignment -1 
		// that means we do allow alignments that contain more sequences 
		if (templateGraphs.size()>al.getNumberOfSequences()-1){
			throw new GraphAveragerError("Number of sequences in alignment is different from number of templates +1 ");
		}
		// now we check if every id from the graphs is really present in the alignment
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
				throw new GraphAveragerError("Sequence of template graph "+tag+" does not match sequence in alignment");
			}			
		}
	}
	
	/**
	 * Counts the votes for each possible alignment edge and puts all the votes in contactVotes Map.
	 * Also initializes the numContacts array with the number of contacts for each template.
	 *
	 */
	private void countVotes() {
		
		// count num contacts for each template
		numContacts = new int[this.getNumberOfTemplates()];
		int c = 0;
		for(String tag:templateGraphs.keySet()) {
			RIGraph g = templateGraphs.get(tag);
			numContacts[c] = g.getEdgeCount();
			c++;
		}
		Arrays.sort(numContacts);
		
		// count number of votes for each contact
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
				
				int votes = 0;
				int potentialVotes = 0;
				TreeSet<String> voters = new TreeSet<String>(); 
				// scanning all templates to see if they have this contact
				for (String tag:templateGraphs.keySet()){			
					RIGraph thisGraph = templateGraphs.get(tag);
					int iSeqIdx = al.al2seq(tag, i);
					int jSeqIdx = al.al2seq(tag, j);
					
					// if either of the ends maps to a gap in this sequence we skip it
					if ((iSeqIdx!=-1) && (jSeqIdx!=-1)) {
						potentialVotes++;
						if(thisGraph.containsEdgeIJ(iSeqIdx, jSeqIdx)) {
							RIGEdge thisGraphCont = thisGraph.getEdgeFromSerials(iSeqIdx, jSeqIdx);
							Pair<RIGNode> pair = thisGraph.getEndpoints(thisGraphCont);
							// TODO: Why is this check necessary? Isn't this always true?
							if (thisGraph.containsEdgeIJ(pair.getFirst().getResidueSerial(), pair.getSecond().getResidueSerial())) {
								votes++;
								voters.add(tag);
							}
						}
					}
				}
				// putting vote in contactVotes Map
				if (votes>0){
					contactVotes.put(new Pair<Integer>(i,j), new Vote(votes, potentialVotes, voters));
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
	 * @return The average number of contacts over all templates.
	 */
	public double getAvgNumContacts() {
		double sum = 0.0;
		for (int i = 0; i < numContacts.length; i++) {
			sum += numContacts[i];
		}
		return sum / getNumberOfTemplates();
	}
	
	/**
	 * @return The median number of contacts over all templates
	 */
	public int getMedNumContacts() {
		int medIdx = getNumberOfTemplates() / 2;
		return numContacts[medIdx];
	}
	
	/**
	 * Returns the number of contacts such that a fraction t of templates has less contacts.
	 * and 1-t has less. The index is rounded to the nearest integer and t is forced into [0;1].
	 * @param t the quantile at which the number of contacts is returned
	 * @return the number of contacts at the t quantile (rounded).
	 */
	public int getQuantNumContacts(double t) {
		if(t > 1) t = 1;
		if(t < 0) t = 0;
		int qntIdx = (int) Math.round(t * getNumberOfTemplates()); 
		return numContacts[qntIdx];
	}
	
	/**
	 * @return the minimum number of contacts over all templates
	 */
	public int getMinNumContacts() {
		return numContacts[0];
	}
	
	/**
	 * @return the maximum number of contacts over all templates
	 */
	public int getMaxNumContacts() {
		return numContacts[getNumberOfTemplates()-1];
	}

	/**
	 * Returns a map containing all parents used in the templates for this graph averager mapped
	 * to their frequencies.
	 * @return a map from parent strings to frequencies.
	 */
	public Map<String, Integer> getParentFrequencies() {
		HashMap<String, Integer> parents = new HashMap<String, Integer>();
		for(String tag:templateGraphs.keySet()) {
			RIGraph g = templateGraphs.get(tag);
			String[] ps = g.getParents();
			if(ps != null) {
				for(String p:ps) {
					if(parents.containsKey(p)) {
						int count = parents.get(p);
						parents.put(p,count+1);
					} else {
						parents.put(p,1);
					}
				}
			}
		}
		Map<String, Integer> orderedParents = Goodies.sortMapByValue(parents, Goodies.DESCENDING);
		return orderedParents;	
	}
		
	/**
	 * @return the overlap (=number of shared contacts) between two templates under the given alignment
	 */
	public int getPairwiseOverlap(String tag1, String tag2) {
		RIGraph rig1 = templateGraphs.get(tag1);
		RIGraph rig2 = templateGraphs.get(tag2);
		int sharedEdges = 0;
		for(RIGEdge e:rig1.getEdges()) {
			Pair<RIGNode> eps = rig1.getEndpoints(e);
			int i = eps.getFirst().getResidueSerial();
			int j = eps.getSecond().getResidueSerial();
			// map i,j to alignment
			int ali = this.al.seq2al(tag1, i);
			int alj = this.al.seq2al(tag1, j);
			// map to second sequence
			int i2 = al.al2seq(tag2, ali);
			int j2 = al.al2seq(tag2, alj);
			// check whether rig2 contains this edge
			if(rig2.containsEdgeIJ(i2, j2)) sharedEdges++;
		}		
		return sharedEdges;
	}
	
	/**
	 * @return the sum of pairwise overlaps between all templates as a measure of alignment quality
	 */
	public int getSumOfPairsOverlap() {
		int sum = 0;
		for(String tag1:templateGraphs.keySet()) {
			for(String tag2:templateGraphs.keySet()) {
				if(tag1.compareTo(tag2) < 0) {
					int ol1 = getPairwiseOverlap(tag1, tag2);
					sum += ol1;
					// for debugging only:
					int ol2 = getPairwiseOverlap(tag2, tag1);
					if(ol1 != ol2) {
						System.err.printf("Error in GraphAverager.getSumOfPairsOverlap(): overlap(%s,%s) = %d != %d = overlap(%s,%s)",
								tag1, tag2, ol1, ol2, tag2, tag1);
					}
				}
			}
		}
		return sum;
	}
	
	/**
	 * Prints the overlap values for all pairs of templates to stdout.
	 */
	public void printPairwiseOverlaps() {
		System.out.println("Pairwise contact overlaps:");
		for(String tag1:templateGraphs.keySet()) {
			for(String tag2:templateGraphs.keySet()) {
				if(tag1.compareTo(tag2) <= 0) {
					int ol1 = getPairwiseOverlap(tag1, tag2);
					System.out.printf("%s\t%s\t%d\n", tag1, tag2, ol1);
				}
			}
		}
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
		
		boolean legacyMode = false;	// take threshold as fraction of total number of templates
		RIGraph graph = new RIGraph(this.sequence);
		graph.setContactType(this.contactType);
		graph.setCutoff(this.distCutoff);
		int numTemplates = templateGraphs.size();
		
		// if vote above threshold we take the contact for our target
		// old: int voteThreshold = (int) Math.ceil((double)numTemplates*threshold); // i.e. round up of 50%, 40% or 30% (depends on threshold given)
		for (Pair<Integer> alignCont:contactVotes.keySet()){
			int potentialContacts = contactVotes.get(alignCont).getPotentialVoteCount();
			int actualContacts = contactVotes.get(alignCont).getVoteCount();
			if(legacyMode) potentialContacts = numTemplates;
			double weight = 1.0 * actualContacts / potentialContacts;
			if (weight > threshold) {
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
	 * The sequence of the graph is the target sequence.
	 * @return the average graph
	 */
	public RIGraph getAverageGraph() {

		boolean legacyMode = false;	// take threshold as fraction of total number of templates		
		RIGraph graph = new RIGraph(this.sequence);
		graph.setContactType(this.contactType);
		graph.setCutoff(this.distCutoff);
		int numTemplates = templateGraphs.size();
		
		for (Pair<Integer> alignCont:contactVotes.keySet()){
			int potentialContacts = contactVotes.get(alignCont).getPotentialVoteCount();
			int actualContacts = contactVotes.get(alignCont).getVoteCount();
			if(legacyMode) potentialContacts = numTemplates;
			double weight = 1.0 * actualContacts / potentialContacts;
			int target_i_res = al.al2seq(targetTag,alignCont.getFirst());
			int target_j_res = al.al2seq(targetTag,alignCont.getSecond());
			if (target_i_res!=-1 && target_j_res!=-1) { // we can't add contacts that map to gaps!!
				graph.addEdge(new AveragedRIGEdge(weight,contactVotes.get(alignCont).getVoters()), graph.getNodeFromSerial(target_i_res), graph.getNodeFromSerial(target_j_res), EdgeType.UNDIRECTED);				
			}
		}		
		return graph;
	}

	/**
	 * Returns a RIGraph with numContacts contacts where the contacts are picked from
	 * the union of all templates in order of consensus score (i.e. fraction of templates
	 * confirming the contact). 
	 * @param numContacts the number of contacts picked
	 * @return the graph with top contacts
	 */
	public RIGraph getGraphWithTopContacts(int numContacts) {
		
		// order edges by weight
		RIGraph av = getAverageGraph();
		Collection<RIGEdge> edges = av.getEdges();
		ArrayList<RIGEdge> edgeList = new ArrayList<RIGEdge>(edges);
		Collections.sort(edgeList, new Comparator<RIGEdge>() {
			public int compare(RIGEdge e1, RIGEdge e2) {
			return -1 * Double.compare(e1.getWeight(), e2.getWeight());
		}
		});
		
		// create new graph
		RIGraph graph = new RIGraph(this.sequence);
		graph.setContactType(this.contactType);
		graph.setCutoff(this.distCutoff);
		
		// add top edges to graph
		for (int i = 0; i < Math.min(edgeList.size(), numContacts); i++) {
			RIGEdge e = edgeList.get(i);
			e.setWeight(1.0);
			graph.addEdge(e, av.getEndpoints(e));
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
				int count = 0;
				if(contactVotes.containsKey(i)) {
					count = contactVotes.get(i).getVoteCount();
				} else {
					System.err.println("Severe error in GraphAverager.getConsensusScore(): contactVotes does not contain key " + i + "(this may be a serious bug!)");
				}
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
	 * TODO: How does this compare to SumOfPairsOverlap?
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

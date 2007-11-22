package proteinstructure;

import java.util.TreeMap;

import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

//TODO test the class. After rewriting with JUNG2 hasn't been tested
public class GraphAverager {

	private Alignment al;
	private TreeMap<String,RIGraph> templateGraphs;
	private String targetTag;
	private int numTemplates;
	private String sequence;		// sequence of the final consensus graph
	private String contactType;		// contact type of the final consensus graph
	private double distCutoff;		// cutoff of the final consensus graph
	
	private TreeMap<Pair<Integer>,Integer> contactVotes;
	
	public GraphAverager(String sequence, Alignment al, TreeMap<String,RIGraph> templateGraphs, String targetTag) {
		this.al = al;
		this.templateGraphs = templateGraphs;
		this.targetTag = targetTag;
		this.sequence = sequence;
		RIGraph firstGraph = templateGraphs.get(templateGraphs.firstKey());
		this.contactType = firstGraph.getContactType();
		this.distCutoff = firstGraph.getCutoff();
		
		this.numTemplates = templateGraphs.size();
		checkSequences();
		
		countVotes(); // does the averaging by counting the votes and putting them into contactVotes
		
	}
	
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
			// TODO throw exception
		}
		for (String tag:templateGraphs.keySet()){
			if(!al.getSequenceNoGaps(tag).equals(templateGraphs.get(tag).getSequence())) {
				System.err.println("Sequence of template graph "+tag+" does not match sequence in alignment");
				// TODO throw exception
			}			
		}
	}
	
	/**
	 * Counts the votes for each possible alignment edge and puts all the votes in contactVotes TreeMap
	 *
	 */
	private void countVotes() {
		
		contactVotes = new TreeMap<Pair<Integer>, Integer>();
		
		// we go through all positions in the alignment
		for (int i=0; i<al.getAlignmentLength(); i++){
			for (int j=0; j<al.getAlignmentLength(); j++) {
 
				int vote = 0; 
				// scanning all templates to see if they have this contact
				for (String tag:templateGraphs.keySet()){			
					RIGraph thisGraph = templateGraphs.get(tag);
					//NOTE that order in which we give the nodes in findEdge doesn't matter (but ONLY if graph is undirected!)
					RIGEdge thisGraphCont = thisGraph.findEdge(thisGraph.getNodeFromSerial(al.al2seq(tag, i)), thisGraph.getNodeFromSerial(al.al2seq(tag, j)));
					Pair<RIGNode> pair = thisGraph.getEndpoints(thisGraphCont);
					if (thisGraph.findEdge(thisGraph.getNodeFromSerial(pair.getFirst().getResidueSerial()),thisGraph.getNodeFromSerial(pair.getSecond().getResidueSerial()))!=null) {
						vote++;
					}
				}
				// putting vote in contactVotes TreeMap
				if (vote>0){
					contactVotes.put(new Pair<Integer>(i,j), vote);
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
	public RIGraph doAveraging(double threshold) {
		
		RIGraph graph = new RIGraph(this.sequence);
		graph.setContactType(this.contactType);
		graph.setCutoff(this.distCutoff);
		
		// if vote above threshold we take the contact for our target
		int voteThreshold = (int) Math.ceil((double)numTemplates*threshold); // i.e. round up of 50%, 40% or 30% (depends on threshold given)
		for (Pair<Integer> alignCont:contactVotes.keySet()){
			if (contactVotes.get(alignCont)>=voteThreshold) {
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
	 * Returns a RIGraph containing the union of edges of the template graphs weighted by the fraction of occurance in the templates.
	 * The sequence of the graph is initalized to the sequence of the template.
	 * @return
	 */
	public RIGraph getAverageGraph() {
		RIGraph graph = new RIGraph(this.sequence);
		graph.setContactType(this.contactType);
		graph.setCutoff(this.distCutoff);
		
		for (Pair<Integer> alignCont:contactVotes.keySet()){
			double weight = 1.0 * contactVotes.get(alignCont) / numTemplates;
			int target_i_res = al.al2seq(targetTag,alignCont.getFirst());
			int target_j_res = al.al2seq(targetTag,alignCont.getSecond());
			if (target_i_res!=-1 && target_j_res!=-1) { // we can't add contacts that map to gaps!!
				graph.addEdge(new RIGEdge(weight), graph.getNodeFromSerial(target_i_res), graph.getNodeFromSerial(target_j_res), EdgeType.UNDIRECTED);				
			}
		}		
		return graph;
	}
}

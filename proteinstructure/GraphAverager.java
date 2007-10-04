package proteinstructure;

import java.util.TreeMap;

public class GraphAverager {

	private Alignment al;
	private TreeMap<String,Graph> templateGraphs;
	private String targetTag;
	private int numTemplates;
	private String sequence;
	
	private TreeMap<Edge,Integer> contactVotes;
	
	public GraphAverager(String sequence, Alignment al, TreeMap<String,Graph> templateGraphs, String targetTag) {
		this.al = al;
		this.templateGraphs = templateGraphs;
		this.targetTag = targetTag;
		this.sequence = sequence;
		
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
		
		contactVotes = new TreeMap<Edge, Integer>();
		
		// we go through all positions in the alignment
		for (int i=0; i<al.getAlignmentLength(); i++){
			for (int j=0; j<al.getAlignmentLength(); j++) {
 
				int vote = 0; 
				// scanning all templates to see if they have this contact
				for (String tag:templateGraphs.keySet()){					
					Edge thisGraphCont = new Edge(al.al2seq(tag, i),al.al2seq(tag, j));
					if (templateGraphs.get(tag).containsContact(thisGraphCont)) {
						vote++;
					}
				}
				// putting vote in contactVotes TreeMap
				if (vote>0){
					contactVotes.put(new Edge(i,j), vote);
				}				
			}
		}		
	}
	
	/**
	 * Calculates the consensus graph from the set of template graphs. An edge is contained
	 * in the consensus graph if the fractions of template graphs it is contained in is above
	 * the given threshold. The resulting consensus edges are added to the output graph passed
	 * to the construtor on creation. Note that access to this graph is by reference, so the
	 * original graph is modified.
	 * @param threshold the threshold above which an edge is taken to be a consensus edge
	 */
	public Graph doAveraging(double threshold) {
		
		Graph graph = new Graph(this.sequence);
		
		// if vote above threshold we take the contact for our target
		int voteThreshold = (int) Math.ceil((double)numTemplates*threshold); // i.e. round up of 50%, 40% or 30% (depends on threshold given)
		for (Edge alignCont:contactVotes.keySet()){
			if (contactVotes.get(alignCont)>=voteThreshold) {
				Edge targetGraphCont = new Edge(al.al2seq(targetTag,alignCont.i),al.al2seq(targetTag,alignCont.j));
				if (targetGraphCont.i!=-1 && targetGraphCont.j!=-1){ // we can't add contacts that map to gaps!!
					graph.addEdge(targetGraphCont);
				}
			}
		}
		return graph;
		
	}
	
}

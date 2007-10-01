package proteinstructure;

import java.util.TreeMap;

public class GraphAverager {

	private static final double DEFAULT_THRESHOLD = 0.5;
	
	private Alignment al;
	private TreeMap<String,Graph> templateGraphs;
	private String targetTag;
	private double threshold;
	
	private Graph graph;
	
	public GraphAverager(Graph graph, Alignment al, TreeMap<String,Graph> templateGraphs, String targetTag) {
		this.graph = graph;
		this.al = al;
		this.templateGraphs = templateGraphs;
		this.targetTag = targetTag;
		this.threshold = DEFAULT_THRESHOLD;
		checkSequences();
		
	}
	
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
		if(!al.getSequenceNoGaps(targetTag).equals(graph.getSequence())) {
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
	
	public void setThreshold(double threshold){
		this.threshold = threshold;
	}

	public Graph predict() {

		int numTemplates = templateGraphs.size(); 
		
		TreeMap<Edge,Integer> contactVotes = new TreeMap<Edge, Integer>();
		
		// we take all contacts of first template and assign votes based on their presence in the other templates

		for (int i=0; i<al.getAlignmentLength(); i++){
			for (int j=0; j<al.getAlignmentLength(); j++) {
 
				int vote = 0; 
				// scanning other templates to see if they have this contact
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
		
		return this.graph;
	}
	
}

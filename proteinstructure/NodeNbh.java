package proteinstructure;

import java.util.TreeMap;

public class NodeNbh extends TreeMap<Integer,String> {
	

	private static final long serialVersionUID = 1L;
	
	public static final String centralLetter="x";

	// central residue
	public int central_resser;
	public String central_resType;
	
	/**
	 * Specific NodeNbh: is a neighbourhood in a specific structure e.g. ABCxEFG with x=D in position 25
	 * @param resser
	 * @param resType
	 */
	public NodeNbh(int resser, String resType){
		super();
		this.central_resser=resser;
		this.central_resType=resType;
	}

	public String getMotifFullGaps(){
		String motif="";
		if(!this.isEmpty()) {
			int min=Math.min(central_resser, this.firstKey());
			int max=Math.max(central_resser, this.lastKey());
			for (int i=min;i<=max;i++){
				if (this.containsKey(i)){
					motif+=AAinfo.threeletter2oneletter(this.get(i));
				} else if (i==central_resser){
					motif+=centralLetter;
				} else {
					motif+="_";
				}
			}
		}
		return motif;
	}
	
	public String getMotif(){
		String motif="";
		if(!this.isEmpty()) {
			int min=Math.min(central_resser, this.firstKey());
			int max=Math.max(central_resser, this.lastKey());
			int gapSize=0;
			String gap="";
			for (int i=min;i<=max;i++){
				if (this.containsKey(i)){
					motif+=gap;
					motif+=AAinfo.threeletter2oneletter(this.get(i));
					gapSize=0;
					gap="";
				} else if (i==central_resser){
					motif+=gap;
					motif+=centralLetter;
					gapSize=0;
					gap="";
				} else {
					gapSize++;
					gap="_{"+gapSize+"}";
				}
			}
		}
		return motif;
	}

	public String getMotifNoGaps(){
		String motif="";
		if(!this.isEmpty()) {
			int min=Math.min(central_resser, this.firstKey());
			int max=Math.max(central_resser, this.lastKey());
			for (int i=min;i<=max;i++){
				if (this.containsKey(i)){
					motif+=AAinfo.threeletter2oneletter(this.get(i));
				} else if (i==central_resser){
					motif+=centralLetter;
				}
			}
		}
		return motif;
	}

	public String getMotifReducedAlphabet(Graph graph) {
		String motif="";
		if(!this.isEmpty()) {
			int min=Math.min(this.central_resser, this.firstKey());
			int max=Math.max(this.central_resser, this.lastKey());
			for (int i=min;i<=max;i++){
				if (this.containsKey(i)){
					Edge edge = graph.contacts.getEdge(Math.min(i, this.central_resser),Math.max(i, this.central_resser));
					double weight = edge.weight;
					if (weight>0) { // SC dominated
						motif+=AAinfo.threeletter2oneletter(this.get(i));
					} else if (weight<0) { //BB dominated
						// for the sec structure we take the "seen" node (so not the central but the other), following Michael's convention
						SecStrucElement elem = graph.getSecondaryStructure().getSecStrucElement(i);
						char ssType = elem==null?'O':elem.getType();
						char ssLetter = 'o';
						if (ssType=='S') ssLetter = 'b';
						else if (ssType=='H') ssLetter= 'z';
						else if (ssType=='O') ssLetter = 'o';
						motif+=ssLetter;
					}
				} else if (i==this.central_resser){
					motif+=NodeNbh.centralLetter;
				}
			}
		}
		return motif;		
	}
	
	public String getCommaSeparatedResSerials(){
		String ressers="";
		for (int resser:this.keySet()){
			if (this.lastKey()!=resser) {
				ressers += resser+","; // we are not in last element, we put a comma
			} else {
				ressers += resser; // for last element no comma
			}
		}
		return ressers;
	}

	public String toString(){
		if (this.isEmpty()) return "";
		else return this.getMotif();
	}
	
	/**
	 * Returns a copy (deep) of this NodeNbh as a new NodeNbh object
	 * @return
	 */
	public NodeNbh getCopy(){
		NodeNbh copy = new NodeNbh(this.central_resser,this.central_resType);
		copy.putAll(this);
		return copy;
	}
}

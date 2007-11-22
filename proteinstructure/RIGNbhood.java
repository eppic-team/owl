package proteinstructure;

import java.util.TreeMap;



public class RIGNbhood extends TreeMap<Integer,RIGNode> {
	

	private static final long serialVersionUID = 1L;
	
	public static final String centralLetter="x";

	// central residue
	private RIGNode centralResidue;
	
	/**
	 * Construct a RIGNbhood by passing the central residue as a RIGNode object
	 * @param centralResidue
	 */
	public RIGNbhood(RIGNode centralResidue){
		super();
		this.centralResidue = centralResidue;
	}

	public RIGNode getCentralResidue() {
		return centralResidue;
	}
	
	public String getMotifFullGaps(){
		String motif="";
		if(!this.isEmpty()) {
			int min=Math.min(centralResidue.getResidueSerial(), this.firstKey());
			int max=Math.max(centralResidue.getResidueSerial(), this.lastKey());
			for (int i=min;i<=max;i++){
				if (this.containsKey(i)){
					motif+=AAinfo.threeletter2oneletter(this.get(i).getResidueType());
				} else if (i==centralResidue.getResidueSerial()){
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
			int min=Math.min(centralResidue.getResidueSerial(), this.firstKey());
			int max=Math.max(centralResidue.getResidueSerial(), this.lastKey());
			int gapSize=0;
			String gap="";
			for (int i=min;i<=max;i++){
				if (this.containsKey(i)){
					motif+=gap;
					motif+=AAinfo.threeletter2oneletter(this.get(i).getResidueType());
					gapSize=0;
					gap="";
				} else if (i==centralResidue.getResidueSerial()){
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
			int min=Math.min(centralResidue.getResidueSerial(), this.firstKey());
			int max=Math.max(centralResidue.getResidueSerial(), this.lastKey());
			for (int i=min;i<=max;i++){
				if (this.containsKey(i)){
					motif+=AAinfo.threeletter2oneletter(this.get(i).getResidueType());
				} else if (i==centralResidue.getResidueSerial()){
					motif+=centralLetter;
				}
			}
		}
		return motif;
	}

	public String getMotifReducedAlphabet(RIGraph graph) {
		String motif="";
		if(!this.isEmpty()) {
			int min=Math.min(centralResidue.getResidueSerial(), this.firstKey());
			int max=Math.max(centralResidue.getResidueSerial(), this.lastKey());
			for (int i=min;i<=max;i++){
				if (this.containsKey(i)){
					//TODO check if the following is correct. If JUNG's SparseGraph is correctly implemented then order of vertices given below shouldn't matter
					RIGEdge edge = graph.findEdge(centralResidue,this.get(i));
					double weight = edge.getWeight();
					if (weight>0) { // SC dominated
						motif+=AAinfo.threeletter2oneletter(this.get(i).getResidueType());
					} else if (weight<0) { //BB dominated
						// for the sec structure we take the "seen" node (so not the central but the other), following Michael's convention
						char ssType = this.get(i).getSecStrucElement().getType();
						char ssLetter = 'o';
						if (ssType == 0 || ssType == 'O') ssLetter='o';
						if (ssType=='S') ssLetter = 'b';
						else if (ssType=='H') ssLetter= 'z';
						motif+=ssLetter;
					}
				} else if (i==this.centralResidue.getResidueSerial()){
					motif+=RIGNbhood.centralLetter;
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
	 * Returns a copy (deep) of this NodeNbh 
	 * @return
	 */
	public RIGNbhood copy(){
		RIGNbhood copy = new RIGNbhood(this.centralResidue.copy());
		for (int resser:this.keySet()) {
			copy.put(resser,this.get(resser).copy());
		}
		return copy;
	}
}

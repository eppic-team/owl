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
					motif+=AA.threeletter2oneletter(this.get(i));
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
					motif+=AA.threeletter2oneletter(this.get(i));
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
					motif+=AA.threeletter2oneletter(this.get(i));
				} else if (i==central_resser){
					motif+=centralLetter;
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

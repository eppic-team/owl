package proteinstructure;

import java.util.TreeSet;


public class EdgeSet extends TreeSet<Edge> {

	private static final long serialVersionUID = 1L;

	public EdgeSet() {
		super();
	}
	
	/**
	 * Gets the maximum range of this EdgeSet
	 * i.e. the sequence separation for the pair with maximum sequence separation
	 * @return
	 */
	public int getMaxRange() {
		int max=0;
		for (Edge cont:this){
			max = Math.max(max, cont.getRange());
		}
		return max;
	}

	/**
	 * Gets the maximum node serial in this EdgeSet
	 * @return
	 */
	public int getMaxNode(){
		int max=0;
		for (Edge cont:this){
			int contactMax=Math.max(cont.i, cont.j);
			max = Math.max(max,contactMax);
		}
		return max;
	}
}

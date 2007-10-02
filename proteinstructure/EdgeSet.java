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
	
	/**
	 * Returns a deep copy of this edgeSet
	 */
	public EdgeSet copy() {
		EdgeSet copy = new EdgeSet();
		for(Edge e:this) {
			copy.add(e.copy());
		}
		return copy;
	}
	
	/**
	 * Filters the edgeset to remove everything above/below/equal to the given cutoff.
	 * If sign < 0 remove edges with weight below the cutoff,
	 * if sign > 0 remove edges with weight above the cutoff,
	 * if sign==0 remove edges with weight equal to the cutoff.
	 * @param cutoff the filter cutoff
	 * @param sign whether to remove edges below/above or equal to the threshold
	 */
	public void filterEdges(double cutoff, int sign) {
		EdgeSet tempEdges = this.copy();
		for(Edge e:tempEdges) {
			if(sign > 0) {
				if(e.weight > cutoff) this.remove(e);
			} else
			if(sign < 0) {
				if(e.weight < cutoff) this.remove(e);
			} else {
				if(e.weight == cutoff) this.remove(e);
			}
		}
	}
}

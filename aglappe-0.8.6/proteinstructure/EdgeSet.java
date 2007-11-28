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
	 * Returns true if this EdgeSet equals the given other EdgeSet. Edges are considered
	 * equal if they connect the same node nums.
	 * @param other the EdgeSet to compare with
	 * @return true if the two sets are equal, false otherwise
	 */
	public boolean equals(EdgeSet other) {
		if(this.size() != other.size()) return false;
		for(Edge e:this) {
			if(!other.contains(e)) {
				return false;
			}
		}
		// just to be safe:
		for(Edge e:other) {
			if(!this.contains(e)) {
				return false;
			}
		}
		return true;
	}
	
	/**
	 * Filters the edgeset to remove everything above/below/equal to the given weight cutoff.
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
	
	/**
	 * Gets set of all incident nodes to edges of this edge set.
	 * @return node set of incident nodes 
	 */
	public NodeSet getIncidentNodes() {
	    NodeSet nodes = new NodeSet();
	    for( Edge e : this ) {
		nodes.add(new Node(e.i));
		nodes.add(new Node(e.j));
	    }
	    return nodes;
	}
	
	/**
	 * Returns the Edge object given the 2 residue serials i_num and j_num.
	 * If no Edge exists for the given i_num, j_num then returns null 
	 * NOTE: in non-directed case, edges are stored only in 1 direction (j>i): thus i,j must be passed correctly (j>i)
	 * @param i_num
	 * @param j_num
	 * @return
	 */
	//TODO this is a temporary fix, we need to redesign the EdgeSet to be a Map which has a get method
	public Edge getEdge(int i_num, int j_num){
		for (Edge cont:this) {
			if (cont.equals(new Edge(i_num,j_num))) {
				return cont;
			}
		}
		return null;
	}
}

package proteinstructure;
import java.lang.Comparable;

public class Edge implements Comparable {
	
	public static final double DEFAULT_WEIGHT = 1;
	
	public int i;
	public int j;
	public double weight;	// Note: two edges are considered the same if i and j are the same, disregarding the weight

	/**
	 * Creates a new edge from i to j with default edge weight.
	 * @param i
	 * @param j
	 */
	public Edge(int i,int j){
		this.i=i;
		this.j=j;
		this.weight=DEFAULT_WEIGHT;
	}
	
	/**
	 * Create a new edge from i to j with the given edge weight.
	 * @param i
	 * @param j
	 * @param weight
	 */
	public Edge(int i,int j, double weight){
		this.i=i;
		this.j=j;
		this.weight=weight;
	}
	
	/**
	 * Create an edge from node i to node j. Note that this is equivalent to
	 * Edge((i.getNum(), j.getNum()) and is added for convenience only.
	 * @param i node i
	 * @param j node j
	 */
	public Edge(Node i,Node j) {
		this(i.getNum(), j.getNum());
	}

	/**
	 * Create an edge from node i to node j with the given edge weight. Note that this is equivalent to
	 * Edge((i.getNum(), j.getNum(), weight) and is added for convenience only.
	 * @param i node i
	 * @param j node j
	 */
	public Edge(Node i,Node j, double weight) {
		this(i.getNum(), j.getNum());
	}
	
	/**
	 * Return a deep copy if this edge
	 * @return
	 */
	public Edge copy() {
		return new Edge(i, j, weight);
	}
	
	public int compareTo(Object o) {
		Edge other = (Edge) o;
		if (this.i>other.i){
			return 1;
		} 
		else if (this.i<other.i){
			return -1;
		}
		else { // only remains case this.i==other.i
			if (this.j>other.j){
				return 1;
			}
			else if (this.j<other.j){
				return -1;
			}
		}				
		return 0; // if none of the conditions before returned, then both i and j are equal 
	}

	public boolean equals(Object o){
		Edge other = (Edge) o;
		if (this.i==other.i && this.j==other.j){
			return true;
		}
		return false;
	}

	public String toString() {
		return this.i+" "+this.j;
	}
	
	public int hashCode() {
		return i*100000+j; // hash function found after a lot of experimenting!! do not touch!
	}

	/**
	 * Gets range (i.e. sequence separation) of contact
	 * @return
	 */
	public int getRange(){
		return Math.abs(this.i-this.j);
	}
	
}


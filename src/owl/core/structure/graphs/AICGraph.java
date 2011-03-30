package owl.core.structure.graphs;

import owl.core.structure.Atom;
import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * An atom inter-chain interaction graph.
 * Note that for the clash methods, a clash distance must be passed. A reasonable value for 
 * it is anything under the disulfide bond length (2.05A).
 * 
 * @author duarte_j
 *
 */
public class AICGraph  extends SparseGraph<Atom,AICGEdge> {

	private static final long serialVersionUID = 1L;
	
	private double distCutoff;
	
	public double getDistCutoff() {
		return distCutoff;
	}
	
	public boolean hasClashes(double clashDistance) {
		for (AICGEdge edge:this.getEdges()) {
			if (edge.getDistance()<clashDistance) {
				return true;
			}
		}
		return false;
	}
	
	public int getNumClashes(double clashDistance) {
		int count = 0;
		for (AICGEdge edge:this.getEdges()) {
			if (edge.getDistance()<clashDistance) {
				count++;
			}
		}
		return count;		
	}
	
	/**
	 * Equality based on having exact same set of edges between atoms of same serial and same atom code.
	 */
	public boolean equals(Object other) {
		if (! (other instanceof AICGraph)) return false;
		
		AICGraph o = (AICGraph) other;
		
		if (this.getEdgeCount()!=o.getEdgeCount()) {
			return false;
		}
		
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);
			if (o.findEdge(pair.getFirst(), pair.getSecond())==null) {
				return false;
			}
		}
		return true;
	}

}

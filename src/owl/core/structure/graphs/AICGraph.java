package owl.core.structure.graphs;

import owl.core.structure.Atom;
import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * An atom inter-chain interaction graph.
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
	
	/**
	 * Equality based on having exact same set of edges between atoms of same serial and same atom code.
	 */
	public boolean equals(Object other) {
		if (! (other instanceof AICGraph)) return false;
		
		AICGraph o = (AICGraph) other;
		
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);
			if (!o.containsEdge(o.findEdge(pair.getFirst(), pair.getSecond()))) {
				return false;
			}
		}
		for (AICGEdge edge:o.getEdges()) {
			Pair<Atom> pair = o.getEndpoints(edge);
			if (!this.containsEdge(this.findEdge(pair.getFirst(), pair.getSecond()))) {
				return false;
			}
		}
		return true;
	}
	
}

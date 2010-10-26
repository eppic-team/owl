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
	
	// in theory the shortest distance between any 2 non-H atoms would be that of a disulfide bond (2.05A)
	public static final double CLASH_DISTANCE_CUTOFF = 1.5; 

	private double distCutoff;
	
	public double getDistCutoff() {
		return distCutoff;
	}
	
	public boolean hasClashes() {
		for (AICGEdge edge:this.getEdges()) {
			if (edge.getDistance()<CLASH_DISTANCE_CUTOFF) {
				return true;
			}
		}
		return false;
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

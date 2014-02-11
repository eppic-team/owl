package owl.core.util;

import java.util.HashSet;
import java.util.TreeSet;

import edu.uci.ics.jung.graph.util.Pair;

public class IntPairSet extends HashSet<Pair<Integer>> {

	private static final long serialVersionUID = 1L;

	/**
	 * Gets set of all incident nodes to edges of this edge set.
	 * 
	 * TODO this is a temporary fix to replace getIncidentNodes function in former EdgeSet class
	 * 		where should this go? Should IntPairSet be actually a simple graph (where nodes are ints)?  
	 * @return set of incident nodes 
	 */
	public TreeSet<Integer> getIncidentNodes() {
		TreeSet<Integer> nodes = new TreeSet<Integer>();
		for( Pair<Integer> e : this ) {
			nodes.add(e.getFirst());
			nodes.add(e.getSecond());
		}
		return nodes;
	}
}

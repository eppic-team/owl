package owl.core.util;

import java.util.HashSet;
import java.util.TreeSet;

import edu.uci.ics.jung.graph.util.Pair;

public class FloatPairSet extends HashSet<Pair<Float>> {

	private static final long serialVersionUID = 1L;

	public TreeSet<Float> getIncidentNodes() {
		TreeSet<Float> nodes = new TreeSet<Float>();
		for( Pair<Float> e : this ) {
			nodes.add(e.getFirst());
			nodes.add(e.getSecond());
		}
		return nodes;
	}
}

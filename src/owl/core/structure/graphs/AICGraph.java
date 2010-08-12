package owl.core.structure.graphs;

import owl.core.structure.Atom;
import edu.uci.ics.jung.graph.SparseGraph;

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
}

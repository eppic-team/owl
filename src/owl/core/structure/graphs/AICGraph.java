package owl.core.structure.graphs;

import java.util.ArrayList;
import java.util.List;

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
	
	public static final double DISULFIDE_BRIDGE_DIST = 2.05;
	public static final double DISULFIDE_BRIDGE_DIST_SIGMA = 0.1;
	
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
	 * Returns true if at least one pair of atoms in this atom interchain graph forms 
	 * disulfide bridge. Conditions for it are: they both belong to CYS residues, they are
	 * both SG atoms and the distance between them is 2.05A+-0.1A
	 * @return
	 */
	public boolean hasDisulfideBridges() {
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);
			Atom atomi = pair.getFirst();
			Atom atomj = pair.getSecond();
			if (    atomi.getParentResidue().getShortCode()=='C' &&
					atomj.getParentResidue().getShortCode()=='C' &&
					atomi.getCode().equals("SG") &&
					atomj.getCode().equals("SG") &&
					edge.getDistance()<(DISULFIDE_BRIDGE_DIST+DISULFIDE_BRIDGE_DIST_SIGMA) && 
					edge.getDistance()>(DISULFIDE_BRIDGE_DIST-DISULFIDE_BRIDGE_DIST_SIGMA)) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Returns a list of all pairs of atoms that form a disulfide bond in this interface.
	 * @return
	 */
	public List<Pair<Atom>> getDisulfidePairs() {
		List<Pair<Atom>> disPairs = new ArrayList<Pair<Atom>>();
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);
			Atom atomi = pair.getFirst();
			Atom atomj = pair.getSecond();
			if (    atomi.getParentResidue().getShortCode()=='C' &&
					atomj.getParentResidue().getShortCode()=='C' &&
					atomi.getCode().equals("SG") &&
					atomj.getCode().equals("SG") &&
					edge.getDistance()<(DISULFIDE_BRIDGE_DIST+DISULFIDE_BRIDGE_DIST_SIGMA) && 
					edge.getDistance()>(DISULFIDE_BRIDGE_DIST-DISULFIDE_BRIDGE_DIST_SIGMA)) {
				disPairs.add(pair);
			}
		}
		return disPairs;
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

package owl.core.structure.graphs;

import owl.core.structure.Atom;
import owl.core.structure.Residue;
import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

public class RICGraph extends SparseGraph<RICGNode,RICGEdge> {

	private static final long serialVersionUID = 1L;

	public RICGraph(AICGraph atomGraph) {
		
		populateGraph(atomGraph);

	}
	
	private void populateGraph(AICGraph atomGraph) {
		
		for (AICGEdge aicEdge:atomGraph.getEdges()) {

			Pair<Atom> pair = atomGraph.getEndpoints(aicEdge);

			Residue iResidue = pair.getFirst().getParentResidue();
			Residue jResidue = pair.getSecond().getParentResidue();
						
			RICGNode iNode = new RICGNode(iResidue.getSerial(),iResidue.getLongCode(),iResidue.getParent().getPdbChainCode());
			RICGNode jNode = new RICGNode(jResidue.getSerial(),jResidue.getLongCode(),jResidue.getParent().getPdbChainCode());
			
			boolean isClash = aicEdge.getDistance()<AICGraph.CLASH_DISTANCE;
			boolean isDisulfide = atomGraph.isDisulfideInteraction(pair, aicEdge.getDistance());
			boolean isHBond = atomGraph.isHbondInteraction(pair, aicEdge.getDistance());
			
			if (this.findEdge(iNode, jNode)==null) {

				RICGEdge ricEdge = new RICGEdge();
				ricEdge.setClash(isClash);
				ricEdge.setDisulfide(isDisulfide);
				ricEdge.setnAtoms(0);
				ricEdge.setnHBonds((isHBond?1:0));
				
				this.addEdge(ricEdge, 
					iNode, 
					jNode, 
					EdgeType.UNDIRECTED);
				
			} else {
				
				RICGEdge ricEdge = this.findEdge(iNode, jNode);
				
				// only 1 atom with clash or disulfide will set it for the whole residue,
				// we don't want to reset the state once it is already true
				if (!ricEdge.isClash()) ricEdge.setClash(isClash);
				if (!ricEdge.isDisulfide()) ricEdge.setDisulfide(isDisulfide);
				ricEdge.setnAtoms(ricEdge.getnAtoms()+1);
				ricEdge.setnHBonds(ricEdge.getnHBonds()+(isHBond?1:0));
				
			}

		}
	}
}

package proteinstructure;

import java.util.HashMap;
import java.util.TreeMap;

import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * An Atom Interaction Graph
 *
 */
public class AIGraph extends ProtStructGraph<AIGNode,AIGEdge> {
	
	private static final long serialVersionUID = 1L;

	private boolean crossed; // true if this AIGraph has been obtained from a cross contact type (the ones with "/")
	private double distCutoff;
	
	public AIGraph() {
		super();
		this.crossed = false;
		this.distCutoff = 0;
	}
	
	public double getCutoff() {
		return distCutoff;
	}
	
	public void setCutoff(double distCutoff) {
		this.distCutoff = distCutoff;
	}
	
	public boolean isCrossed() {
		return crossed;
	}
	
	public void setCrossed(boolean crossed) {
		this.crossed = crossed;
	}
	
	/**
	 * Returns a RIGraph by collapsing atom contacts into residue contacts,
	 * using the number of atom edges per residue as the atom weights for the RIGEdges
	 * TODO eventually we can pass a parameter for other ways of assigning atom contact weights
	 * @return
	 */
	public RIGraph getRIGraph() {
		EdgeType et = EdgeType.UNDIRECTED;
		
		// TODO we still use here DIRECTED as default for crossed, eventually this should change by taking another parameter "boolean directed", so crossed could have DIRECTED/UNDIRECTED versions 
		if (this.isCrossed()) {
			et = EdgeType.DIRECTED;
		}
		RIGraph resGraph = new RIGraph();
		resGraph.setPdbCode(this.pdbCode);
		resGraph.setChainCode(this.chainCode);
		resGraph.setPdbChainCode(this.pdbChainCode);
		resGraph.setModel(this.model);
		resGraph.setSequence(this.sequence);
	
		TreeMap<Integer,RIGNode> rignodes = new TreeMap<Integer,RIGNode>();
		for (AIGNode atomNode:this.getVertices()) {
			RIGNode resNode = atomNode.getParent();
			int resser = resNode.getResidueSerial();
			rignodes.put(resser, resNode); // we put in the map each RIGNode several times, that should be fine
		}
		
		// putting the RIGnodes into the RIGraph 
		for (int resser:rignodes.keySet()){
			resGraph.addVertex(rignodes.get(resser));
		}
		// now also adding unobserved residues as RIGNodes (tagged with observed=false) from the sequence
		for (int resser=1; resser<=sequence.length(); resser++) {
			if (!rignodes.containsKey(resser)) {
				String resType = AAinfo.oneletter2threeletter(Character.toString(sequence.charAt(resser-1)));
				RIGNode resNode = new RIGNode(resser,resType);
				resNode.setObserved(false);
				rignodes.put(resser,resNode);
				resGraph.addVertex(resNode);
			}
		}
		
		resGraph.setSerials2NodesMap(rignodes);
		
		// collapsing atomPairs into resPairs and counting atom contacts to assign atom weights
		HashMap<Pair<RIGNode>,Integer> pairs2weights = new HashMap<Pair<RIGNode>, Integer>();  
		for (AIGEdge atomEdge: this.getEdges()){
			Pair<AIGNode> atomPair = this.getEndpoints(atomEdge);
			RIGNode v1 = atomPair.getFirst().getParent();
			RIGNode v2 = atomPair.getSecond().getParent();
			Pair<RIGNode> resPair = new Pair<RIGNode>(v1,v2);
			if (v1!=v2) {
				if (!pairs2weights.containsKey(resPair)) {
					//NOTE the pairs2weights map takes care of eliminating duplicate residue pairs (Maps don't accept duplicate as keys)
					pairs2weights.put(resPair, 1);
				} else {
					pairs2weights.put(resPair,pairs2weights.get(resPair)+1);
				}
			}
		}

		// putting the RIGEdges in the resGraph
		for (Pair<RIGNode> resPair:pairs2weights.keySet()) {
			//TODO put distance in the RIGEdge objects (only makes sense in single atom contact types)
			RIGEdge e = new RIGEdge(pairs2weights.get(resPair));
			resGraph.addEdge(e, resPair, et);//(e, pair, et);
		}
		
		return resGraph;
	}
	
	public int getContactRange(AIGEdge edge) {
		Pair<AIGNode> pair = this.getEndpoints(edge);
		return Math.abs(pair.getFirst().getParent().getResidueSerial()-pair.getSecond().getParent().getResidueSerial());
	}

}

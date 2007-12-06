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
	public RIGraph getRIGraph(boolean directed) {
		EdgeType et = EdgeType.UNDIRECTED;
		
		if (directed) {
			et = EdgeType.DIRECTED;
		}
		RIGraph resGraph = new RIGraph();
		resGraph.setPdbCode(this.pdbCode);
		resGraph.setChainCode(this.chainCode);
		resGraph.setPdbChainCode(this.pdbChainCode);
		resGraph.setModel(this.model);
		resGraph.setSequence(this.sequence);
		resGraph.setSecondaryStructure(this.secondaryStructure);
	
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
		
		resGraph.setSerials2NodesMap(rignodes);
		
		// collapsing atomPairs into resPairs and counting atom contacts to assign atom weights
		HashMap<Pair<RIGNode>,Integer> pairs2weights = new HashMap<Pair<RIGNode>, Integer>();
		HashMap<Pair<RIGNode>,Double> pairs2distances = new HashMap<Pair<RIGNode>,Double>();
		for (AIGEdge atomEdge: this.getEdges()){
			Pair<AIGNode> atomPair = this.getEndpoints(atomEdge);
			RIGNode v1 = atomPair.getFirst().getParent();
			RIGNode v2 = atomPair.getSecond().getParent();
			Pair<RIGNode> resPair = new Pair<RIGNode>(v1,v2);
			if (v1!=v2) {
				if (!pairs2weights.containsKey(resPair)) {
					//NOTE the pairs2weights map takes care of eliminating duplicate residue pairs (Maps don't accept duplicate as keys)
					pairs2weights.put(resPair, 1);
					pairs2distances.put(resPair, atomEdge.getDistance());
				} else {
					pairs2weights.put(resPair,pairs2weights.get(resPair)+1);
					pairs2distances.put(resPair, Math.min(pairs2distances.get(resPair), atomEdge.getDistance()));
				}
			}
		}

		// putting the RIGEdges in the resGraph
		for (Pair<RIGNode> resPair:pairs2weights.keySet()) {
			// if undirected and edge already exists
			if (!directed && (resGraph.findEdge(resPair.getFirst(), resPair.getSecond())!=null)) {
				//increase weight
				RIGEdge e = resGraph.findEdge(resPair.getFirst(), resPair.getSecond());
				e.setAtomWeight(e.getAtomWeight()+pairs2weights.get(resPair));
				e.setDistance(Math.min(e.getDistance(), pairs2distances.get(resPair)));
			} else {
				//add edge
				RIGEdge e = new RIGEdge(pairs2weights.get(resPair));
				e.setDistance(pairs2distances.get(resPair));
				resGraph.addEdge(e, resPair, et);//(e, pair, et);				
			}
		}
		
		return resGraph;
	}
	
	public int getContactRange(AIGEdge edge) {
		Pair<AIGNode> pair = this.getEndpoints(edge);
		return Math.abs(pair.getFirst().getParent().getResidueSerial()-pair.getSecond().getParent().getResidueSerial());
	}
	
	public boolean addGraph(AIGraph graph) {
		//NOTE:The checks below would make sense only for adding RIGraphs
		//In AIGraphs we have as nodes only the selected atoms and not all atoms
		/*
		if (this.getVertexCount()!=graph.getVertexCount()) {
			return false;
		}
		Iterator<Integer> it = graph.getSerials().iterator();
		for (int serial:this.getSerials()) {			
			AIGNode node = this.getNodeFromSerial(serial);
			AIGNode node2 = graph.getNodeFromSerial(it.next());
			if (!node.equals(node2)){
				return false;
			}
		}*/

		boolean change = false;
		
		TreeMap<Integer,RIGNode> rignodes = new TreeMap<Integer,RIGNode>();
		for (AIGNode atomNode : this.getVertices()) {
			rignodes.put(atomNode.getParent().getResidueSerial(), atomNode.getParent());
		}
		
		for (AIGNode atomNode : graph.getVertices()) {
			RIGNode v = null;
			if (!rignodes.containsKey(atomNode.getParent().getResidueSerial())) {
				v = atomNode.getParent();
				rignodes.put(v.getResidueSerial(), v);
				change = true;
			} else {
				v = rignodes.get(atomNode.getParent().getResidueSerial()); 
			}
			if (!this.serials2nodes.containsKey(atomNode.getAtomSerial())) {
				change = true;
				AIGNode v1 = new AIGNode (atomNode.getAtomSerial(), atomNode.getAtomName(), v);
				this.addVertex(v1);
				this.serials2nodes.put(atomNode.getAtomSerial(), v1);
			}
		}
		
		for (AIGEdge atomEdge: graph.getEdges()){
			Pair<AIGNode> atomPair = graph.getEndpoints(atomEdge);
			AIGNode v1 = this.getNodeFromSerial(atomPair.getFirst().getAtomSerial());
			AIGNode v2 = this.getNodeFromSerial(atomPair.getSecond().getAtomSerial());
			// This condition is to take care of not adding multiple instances of the same atomic edge
			if (this.findEdge(v1, v2)==null) {				
				change = true;
				this.addEdge(atomEdge.copy(), v1, v2, graph.getEdgeType(atomEdge));
			}
		}
		this.setCrossed((this.crossed || graph.crossed));
		return change;
	}
	
	/**
	 * Removes a vertex from this graph.
	 * This overridden function also updates the serials2nodes map.
	 * @return true if vertex is present in this graph and thus can be removed, false if vertex is not present in this graph 
	 */
	@Override
	public boolean removeVertex(AIGNode vertex) {
		serials2nodes.remove(vertex.getAtomSerial());
		return super.removeVertex(vertex);
	}
	
}

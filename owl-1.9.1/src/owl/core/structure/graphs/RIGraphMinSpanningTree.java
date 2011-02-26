package owl.core.structure.graphs;

import org.apache.commons.collections15.Factory;
import org.apache.commons.collections15.Transformer;

import owl.core.structure.Pdb;
import owl.core.structure.PdbasePdb;

import edu.uci.ics.jung.algorithms.shortestpath.PrimMinimumSpanningTree;

/**
 * Wrapper for JUNG's PrimMinimumSpanningTree class to be used to find the 
 * Minimum Spanning Tree of a RIGraph. Weights will be the ones set in the RIGEdges.
 * 
 *
 */
public class RIGraphMinSpanningTree {
	
		
	Factory<RIGraph> RIGraphFactory = new Factory<RIGraph>() {
		public RIGraph create() {
			return new RIGraph();
		}
	};

	Transformer<RIGEdge, Double> RIGEdgeWeightTransformer = new Transformer<RIGEdge, Double>() {
		public Double transform(RIGEdge input) {
			return input.getWeight();
		}
	};

	PrimMinimumSpanningTree<RIGNode, RIGEdge> mst;

	public RIGraphMinSpanningTree() {
		this.mst = new PrimMinimumSpanningTree<RIGNode, RIGEdge>(RIGraphFactory, RIGEdgeWeightTransformer);
	}

	/**
	 * Given a RIGraph returns a new RIGraph that is the Minimum Spanning Tree
	 * of the input graph.
	 * Weights are taken from the weights set in the RIGEdges.
	 * @param graph
	 * @return
	 */
	public RIGraph getMST(RIGraph graph) { 
		RIGraph mstGraph = (RIGraph) this.mst.transform(graph);
		mstGraph.setPdbCode(graph.getPdbCode());
		mstGraph.setPdbChainCode(graph.getPdbChainCode());
		mstGraph.setChainCode(graph.getChainCode());
		mstGraph.setModel(graph.getModel());
		mstGraph.setTargetNum(graph.getTargetNum());
		mstGraph.setGroupNum(graph.getGroupNum());
		mstGraph.setCaspModelNum(graph.getCaspModelNum());
		mstGraph.setContactType(graph.getContactType());
		mstGraph.setCutoff(graph.getCutoff());
		mstGraph.setSequence(graph.getSequence());
		
		mstGraph.setSerials2NodesMap();
		
		return mstGraph;
	}

	
	/**
	 * For testing the class
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		
		Pdb pdb = new PdbasePdb("1bxy");
		pdb.load("A");
		RIGraph graph = pdb.getRIGraph("Cb", 8);
		
		System.out.println("Vertex count: "+graph.getVertexCount());
		System.out.println("Edge count: "+graph.getEdgeCount());

		RIGraphMinSpanningTree pmst = new RIGraphMinSpanningTree();
		RIGraph mstGraph = pmst.getMST(graph);
		
		System.out.println("Vertex count: "+mstGraph.getVertexCount());
		System.out.println("Edge count: "+mstGraph.getEdgeCount());

		// code to write a gdl file of the final MST graph
//		GraphIOGDLFile<RIGNode, RIGEdge> graphIO = new GraphIOGDLFile<RIGNode, RIGEdge>();
//		graphIO.writeGdlFile(mstGraph, "test.gdl", 
//				new Transformer<RIGNode,Integer>() {
//					public Integer transform(RIGNode node) {
//						return node.getResidueSerial();
//					}
//				}, 
//				new Transformer<RIGNode,String>() {
//					public String transform(RIGNode node) {
//						return node.getResidueSerial()+node.getResidueType();
//					}
//				},
//				new Transformer<RIGNode,String>() {
//					public String transform(RIGNode node) {
//						return null;
//					}
//				});
		

	}

}

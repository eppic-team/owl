package owl.core.structure;

//import java.util.ArrayList;
//import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import edu.uci.ics.jung.graph.SparseMultigraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * The graph of the interfaces connections, in principle an infinite periodic graph, 
 * but here we store only the "unit cell" of it. If an edge goes through the unit cell
 * boundary then it is marked specially.
 *  
 * We don't store parallel interfaces in it.
 * 
 * @see InterfaceGraphEdge
 * @author duarte_j
 *
 */
public class InterfaceGraph {

	private SparseMultigraph<SubunitId, InterfaceGraphEdge> graph;
	
	public InterfaceGraph(ChainInterfaceList cil) {
		graph  = new SparseMultigraph<SubunitId, InterfaceGraphEdge>();
		for (int id:cil.getNonParallelInterfacesIds()) {
			graph.addEdge(new InterfaceGraphEdge(id,cil.get(id).isWithinUnitCell()), 
					new Pair<SubunitId>(cil.get(id).getFirstSubunitId(),cil.get(id).getSecondSubunitId()),
					EdgeType.UNDIRECTED);
		}
	}
	
	public void printInfo() {
		System.out.println("Distinct nodes: "+graph.getVertexCount());
		for (SubunitId node:graph.getVertices()) {
			System.out.println(node);
		}
		System.out.println("Distinct edges: "+graph.getEdgeCount());
		for (InterfaceGraphEdge edge:graph.getEdges()) {
			System.out.println(edge+": "+graph.getEndpoints(edge).getFirst()+"-"+graph.getEndpoints(edge).getSecond());
		}
	}
	
	public Set<InterfaceGraphEdge> getEdgesNotWithinUnitCell() {
		Set<InterfaceGraphEdge> edges = new TreeSet<InterfaceGraphEdge>();
		for (InterfaceGraphEdge edge:graph.getEdges()) {
			if (!edge.isWithinUnitCell()) edges.add(edge);
		}
		return edges;
	}
	
//	public List<Integer> getInducedInterfaces(Assembly ass) {
//		List<Integer> list = new ArrayList<Integer>();
//
//		Set<SubunitId> nodes = ass.getSubunits();
//		
//		for (InterfaceGraphEdge edge:graph.getEdges()) {
//			Pair<SubunitId> pair = graph.getEndpoints(edge);
//			if (!ass.contains(edge.getId()) && 
//				nodes.contains(pair.getFirst()) && nodes.contains(pair.getSecond()) &&
//				edge.isWithinUnitCell()) {
//				// TODO finish! this is not finished at all!
//				list.add(edge.getId());
//			}
//		}
//		
//		return list;
//	}
}

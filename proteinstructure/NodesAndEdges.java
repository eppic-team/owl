package proteinstructure;

/**
 * A light weight graph class consisting of an EdgeSet and a NodeSet.
 * Name chosen because Graph is already taken for a ResidueInteractionGraph.
 * @author Henning Stehr
 */
public class NodesAndEdges {

	/*--------------------------- member variables --------------------------*/
	protected NodeSet nodes;
	protected EdgeSet edges;
	
	/*----------------------------- constructors ----------------------------*/
	/**
	 * Create a new graph given a node set and an edge set
	 */
	public NodesAndEdges(NodeSet n, EdgeSet e) {
		nodes = n;
		edges = e;
	}
	
	/**
	 * Create a new graph with empty node- and edge sets
	 */
	public NodesAndEdges() {
		nodes = new NodeSet();
		edges = new EdgeSet();
	}
	
	/*---------------------------- public methods ---------------------------*/
	/**
	 * Returns a deep copy of this graph
	 */
	public NodesAndEdges copy() {
		return new NodesAndEdges(nodes.copy(), edges.copy());
	}
	
	/**
	 * Returns the node set of the graph
	 * @return the node set
	 */
	public NodeSet getNodes() { return nodes; }
	
	/**
	 * Returns the edge set of the graph
	 * @return the edge set
	 */
	public EdgeSet getEdges() { return edges; }

	/**
	 * Returns the number of nodes
	 * @return the number of nodes
	 */
	public int getNumNodes() { return nodes.size(); }
	
	/**
	 * Returns the number of edges
	 * @return the number of edges
	 */
	public int getNumEdges() { return edges.size(); }
	
	/**
	 * Check whether all edges map to nodes from the node set.
	 * @return true if graph is consistent, false otherwise
	 */
	public boolean isConsistent() {
		boolean consistent = true;
		for(Edge e:edges) {
			if(!nodes.contains(e.i) || !nodes.contains(e.j)) {
				consistent = false;
			}
		}
		return consistent;
	}
	
	/**
	 * Output information about the graph to stdout
	 * @param printNodes whether to output all nodes
	 * @param printEdges whether to output all edges
	 */
	public void printGraph(boolean printNodes, boolean printEdges) {
		System.out.println("Nodes: " + nodes.size() + " Edges: " + edges.size());
		if(printNodes) {
			System.out.println("Nodes:");
			for(int node:nodes) {
				System.out.print(node + " ");
			}
			System.out.println();
		}
		if(printEdges) {
			System.out.println("Edges:");
			for(Edge e:edges) {
				System.out.print("(" + e.i + "," + e.j + "," + e.weight + ") ");
			}
			System.out.println();
		}
	}
	
}

package owl.core.structure;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.vecmath.Point3i;

import owl.core.util.CombinationsGenerator;

import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * The graph of the interfaces connection topology of a given PDB crystal structure,
 * in principle an infinite periodic graph. Here we store only the "unit cell" of it. 
 * <p> 
 * We don't store parallel interfaces in graph (in terms of the graph those are loop 
 * edges).
 * </p>
 * <p>
 * Nodes are SubnitIds: chain identifier (identifier within assymetric unit) + operator 
 * identifier (identifier within unit cell) + translation identifier. 
 * Two unit-cell-translation-related chains will have the same chain identifier and operator
 * but different translation identifiers
 * </p>
 * <p>
 * Edges are InterfaceGraphEdges: interface id (1 to number of unique interfaces in crystal)   
 * </p>
 * 
 * @see InterfaceGraphEdge
 * @see SubunitId
 * 
 * @author duarte_j
 *
 */
public class InterfaceGraph {

	private SparseGraph<SubunitId, InterfaceGraphEdge> graph;
	private SpaceGroup sg;
	
	/**
	 * Constructs the complete InterfaceGraph given a full list of all unique 
	 * interfaces in crystal.
	 * @param cil
	 */
	public InterfaceGraph(ChainInterfaceList cil) {
		graph  = new SparseGraph<SubunitId, InterfaceGraphEdge>();
		for (int id:cil.getNonParallelInterfacesIds()) {
			graph.addEdge(new InterfaceGraphEdge(cil.get(id)), 
					new Pair<SubunitId>(cil.get(id).getFirstSubunitId(),cil.get(id).getSecondSubunitId()),
					EdgeType.UNDIRECTED);
		}
		sg = cil.get(0).getFirstMolecule().getParent().getSpaceGroup();
	}
	
	/**
	 * Constructs an empty InterfaceGraph.
	 * Use {@link #addEdge(InterfaceGraphEdge, SubunitId, SubunitId)} to add edges to it.
	 */
	public InterfaceGraph(SpaceGroup sg) {
		graph = new SparseGraph<SubunitId, InterfaceGraphEdge>();
		this.sg = sg;
	}
	
	public boolean addEdge(InterfaceGraphEdge edge, SubunitId v1, SubunitId v2) {
		return graph.addEdge(edge, v1, v2, EdgeType.UNDIRECTED);
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
	
	/**
	 * Tells whether this graph contains any node in cells different
	 * from the original unit cell
	 * @return
	 */
	public boolean hasNodesInAdjacentCells() {
		for (SubunitId node:graph.getVertices()){
			if (!node.getTransl().equals(new Point3i(0,0,0))) {
				return true;
			}
			
		}
		return false;
	}
	
	/**
	 * Returns the set of all nodes that are in cells different from the
	 * original unit cell
	 * @return
	 */
	public Set<SubunitId> getAdjacentCellsNodes() {
		Set<SubunitId> nodes = new HashSet<SubunitId>();
		for (SubunitId node:graph.getVertices()){
			if (!node.getTransl().equals(new Point3i(0,0,0))) {
				nodes.add(node);
			}
		}
		return nodes;
	}
	
	public Set<InterfaceGraph> getAllAssemblies() {
		Set<InterfaceGraph> assemblies = new HashSet<InterfaceGraph>(); //TODO must define equals() and hashCode() !
		
		// we then enumerate all assemblies with 1 interface, 2 interfaces, 3 interfaces .... up to n
		System.out.println("Theoretical total assemblies: "+(((int)Math.pow(2, graph.getEdgeCount()))-1));
		for (int n=1;n<=graph.getEdgeCount();n++) {
			Set<InterfaceGraph> sizenassemblies = getAssembliesSizeN(n);
			
			assemblies.addAll(sizenassemblies);
		}
		return assemblies;
	}
	
	private Set<InterfaceGraph> getAssembliesSizeN(int n) {
		Set<InterfaceGraph> assemblies = new HashSet<InterfaceGraph>();
		Collection<InterfaceGraphEdge> edges = graph.getEdges();
		InterfaceGraphEdge[] alledges = new InterfaceGraphEdge[edges.size()];
		edges.toArray(alledges);
		CombinationsGenerator cg = new CombinationsGenerator (alledges.length, n);
		// 1 we enumerate all n-combinations of edges
		while (cg.hasMore()) {
			int[] indices = cg.getNext();
			InterfaceGraphEdge[] seledges = new InterfaceGraphEdge[n];
			for (int i = 0; i < indices.length; i++) {
				seledges[i] = seledges[indices[i]];
			}
			
			InterfaceGraph subgraph = new InterfaceGraph(this.sg);

			HashSet<SubunitId> selnodes = new HashSet<SubunitId>();
			for (InterfaceGraphEdge edge:seledges) {
				Pair<SubunitId> pair = graph.getEndpoints(edge);
				selnodes.add(pair.getFirst());
				selnodes.add(pair.getSecond());
			}
			// 2 for each combination of size n we find any other induced edges (all edges that connect any node connected by the engaged edges)
			for (SubunitId v1:selnodes) {
				for (SubunitId v2:selnodes) {
					if (v1.compareTo(v2)<=0) continue; // to skip the repeated and self pairs
					for (InterfaceGraphEdge edge:graph.findEdgeSet(v1, v2)) {
						subgraph.addEdge(edge, v1, v2);
					}
				}
			}
					
			// if the assembly is infinite then it's not a valid assembly, we skip it
			if (subgraph.isInfinite()) continue;

			// now we need to find whether the corresponding symmetry edges to the engaged ones add any more nodes
			subgraph.addSymmetryEdges();

			// once the symmetry nodes are also added then we can add the subgraph to the assemblies set, 
			// taking care of duplicates (because of that assemblies needs to be a set) 
			assemblies.add(subgraph);

		}
		System.out.println("size "+n+": "+assemblies.size());
		return assemblies;
	}
	
	/**
	 * Tells whether this InterfaceGraph corresponds to that of an infinite assembly.
	 * The condition is: there exists a path (length>=1) between a subunit and its symmetry 
	 * related equivalent, e.g. A0 in original cell and A0 in cell (-1,-1,0)
	 * As the subgraphs are always connected this will be true whenever there are two symmetry 
	 * related equivalent subunits in the same graph.
	 * @return
	 */
	private boolean isInfinite() {
		Set<SubunitId> nodes = getAdjacentCellsNodes();
		for (SubunitId extnode:nodes) {
			for (SubunitId node:graph.getVertices()) {
				if (extnode.isSymRelatedEquivalent(node)) return true;
			}
		}
		return false;
	}

	/**
	 * Based on the current connectivity of the graph, inferred symmetry edges are found
	 * and added to this graph (including their nodes)
	 * Rules that apply here:
	 *   In-Jm => Im-Jn
	 *   In-Jn => Im-Jm
	 *   TODO this is not true in general, it really depends on the space group, plus there are possibly many more rules
	 *   TODO must found a formal way to find these rules 
	 */
	private void addSymmetryEdges() {
		// nothing to do in this case
		if (graph.getEdgeCount()<2) return;
		
		for (InterfaceGraphEdge edge:graph.getEdges()) {
			Pair<SubunitId> pair = graph.getEndpoints(edge);
			SubunitId inode = pair.getFirst();
			SubunitId jnode = pair.getSecond();
			List<Pair<SubunitId>> relatedPairs = getRelatedEdges(inode,jnode);
			for (Pair<SubunitId> relatedPair:relatedPairs) {
				// TODO check: we modify the graph here adding more node/edges, that in turns will alter the output of getRelatedEdges... doesn't this go into an infinite loop?
				if (graph.containsVertex(relatedPair.getFirst()) || graph.containsVertex(relatedPair.getSecond())) {
					// note that if the edge was already there, then nothing happens as our graph does not admit parallel edges
					graph.addEdge(edge, relatedPair.getFirst(), relatedPair.getSecond(), EdgeType.UNDIRECTED);
				}
			}
		}
	}
	
	/**
	 * Rules that apply here:
	 *   In-Jm => Im-Jn
	 *   In-Jn => Im-Jm
	 *   TODO this is not true in general, it really depends on the space group, plus there are possibly many more rules 
	 *   TODO must found a formal way to find these rules 
 
	 * @param inode
	 * @param jnode
	 * @return
	 */
	private List<Pair<SubunitId>> getRelatedEdges(SubunitId inode, SubunitId jnode) {
		List<Pair<SubunitId>> pairs = new ArrayList<Pair<SubunitId>>();
		
		// 1 if both are same chain then there's no related edges 
		//   e.g. A0-A1 
		if (inode.getPdbChainCode()==jnode.getPdbChainCode()) return pairs;
		
		// 2 if they are different chains and same operator then all other operators between those 2 chains are related
		//   e.g. A0-B0 => A1-B1, A2-B2, ...
		if (inode.getTransformId()==jnode.getTransformId()) {
			for (int n=0;n<sg.getNumOperators();n++) {
				if (n==inode.getTransformId()) continue;
				pairs.add(new Pair<SubunitId>(
						new SubunitId(inode.getPdbChainCode(),n,null),
						new SubunitId(jnode.getPdbChainCode(),n,null)
						));
			}	
			// TODO does this apply also if the original relationship was intercell?
		}
		
		// 3 if they are different chains and different operators then the opposite relationship follows
		//   e.g. A0-B1 => B0-A1
		else {
			//TODO what if the original relationship was intercell?
			pairs.add(new Pair<SubunitId>(
				new SubunitId(jnode.getPdbChainCode(),inode.getTransformId(),null),
				new SubunitId(inode.getPdbChainCode(),jnode.getTransformId(),null)
				));
		}

		return pairs;
	}
	
	public int getOligomericState() {
		return graph.getVertexCount();
	}
	
	public boolean equals(Object o) {
		if (!(o instanceof InterfaceGraph)) return false;
		InterfaceGraph other = (InterfaceGraph) o;
		
		if (this.graph.getVertexCount()!=other.graph.getVertexCount()) return false;
		if (this.graph.getEdgeCount()!=other.graph.getEdgeCount()) return false;
		
		for (SubunitId v:graph.getVertices()) {
			if (!other.graph.containsVertex(v)) return false;
		}
		
		for (InterfaceGraphEdge e:graph.getEdges()) {
			if (!other.graph.containsEdge(e)) return false;
			if (graph.getEndpoints(e).equals(other.graph.getEndpoints(e))) return false;
		}
		return true;
	}

	public int hashCode() {
		//TODO implement a hash code function compatible with the equals above
		return 0;
	}
}

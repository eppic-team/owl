package embed;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.TreeMap;

import javax.vecmath.Vector3d;

import org.apache.commons.collections15.Transformer;

import Jama.Matrix;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;
import proteinstructure.AAinfo;
import proteinstructure.Pdb;
import proteinstructure.PdbasePdb;
import proteinstructure.RIGEdge;
import proteinstructure.RIGNode;
import proteinstructure.RIGraph;

/**
 * Implementation of the bounds smoothing part of the EMBED algorithm of Crippen and Havel
 * Given a sparse set of distance ranges between a set of atoms it finds distance bounds for
 * all pairs of atoms, using the triangle inequality.
 *  
 * Taken from "Distance Geometry: Theory, Algorithms, and Chemical Applications" (section 3.1) by T.F. Havel, 
 * in Encyclopedia of Computational Chemistry (Wiley, New York, 1998). 
 * See also "Distance Geometry and Molecular Conformation" (Chapter 5) by G.M. Crippen and T.F. Havel (Wiley)
 *   
 * @author duarte
 *
 */
public class BoundsSmoother {
	
	private static final double BB_CA_DIST = 3.8;
	private static final double HARD_SPHERES_BOUND = AAinfo.DIST_MIN_CA ;

	/*----------------- helper classes and transformers -----------------*/
	
	Transformer<SimpleEdge, Number> WeightTransformer = new Transformer<SimpleEdge, Number>() {
		public Number transform(SimpleEdge input) {
			return input.weight;
		}
	};
	
	private class SimpleEdge {
		public double weight;
		public SimpleEdge(double weight) {
			this.weight = weight;
		}
	}
	
	private class BoundsDigraphNode {
		public static final boolean LEFT = true;
		public static final boolean RIGHT = false;
		private boolean side; // true: left, false: right
		private int resSerial;
		public BoundsDigraphNode(int resSerial, boolean side) {
			this.resSerial = resSerial;
			this.side = side;
		}
		public boolean isRight() {
			return !side;
		}
		public boolean isLeft() {
			return side;
		}
		public int getSerial() {
			return resSerial;
		}
		public boolean equals(Object other) {
			if (! (other instanceof BoundsDigraphNode)) return false;
			BoundsDigraphNode otherNode = (BoundsDigraphNode)other;
			if (otherNode.resSerial==this.resSerial && otherNode.side == this.side) {
				return true;
			} else {
				return false;
			}
		}
		public String toString() {
			return resSerial+(side?"L":"R");
		}
	}

	private class Bound {
		public double lower;
		public double upper;
		public Bound(double lower, double upper) {
			this.lower = lower;
			this.upper = upper;
		}
		public String toString() {
			return String.format("[%4.1f %4.1f]", lower, upper);
		}
	}

	/*---------------------- member variables ----------------------------*/
	
	private RIGraph rig;
	private TreeMap<Integer,Integer> matIdx2Resser;
	private int conformationSize;
	private HashMap<Boolean, HashMap<Integer,BoundsDigraphNode>> nodesBoundsDigraph; // map of serial/side to nodes in the bounds digraph
	private double lmax; // maximum of the lower bounds: offset value for the boundsDigraph (not to have negative weights so that we can use Dijkstra's algo)
	
	/*------------------------ constructors ------------------------------*/
	
	/**
	 * Constructs a new BoundsSmoother object given a RIGraph.
	 * The RIGraph is converted into a set of distance ranges using the cutoff
	 * as upper limit and hard-spheres as lower limit. 
	 * @param graph
	 */
	public BoundsSmoother(RIGraph graph) {
		this.rig = graph;
		this.conformationSize = this.rig.getObsLength();
	}
	
	/*----------------------- public methods  ----------------------------*/
	
	/**
	 * Computes bounds for all pairs, based on the set of sparse distance ranges
	 * @return a 2-dimensional array with the bounds for all pairs of residues, the
	 * indices of the array can be mapped to residue serials through {@link #getResserFromIdx(int)}
	 * and are guaranteed to be in the same order as the residue serials.
	 */
	public Bound[][] getBoundsAllPairs() {
		SparseGraph<Integer,Bound> boundsGraph = convertRIGraphToDistRangeGraph(this.rig);
		return getBoundsAllPairs(boundsGraph);
	}
		
	/**
	 * Gets a random sample from a matrix of all pairs distance ranges
	 * @param bounds the matrix of all pairs distance ranges
	 * @return
	 */
	public Matrix sampleBounds(Bound[][] bounds) {
		double[][] matrix = new double[conformationSize][conformationSize];
		for (int i=0;i<conformationSize;i++) {
			for (int j=0;j<conformationSize;j++) {
				if (j>i) {
					Random rand = new Random();
					matrix[i][j]= bounds[i][j].lower+rand.nextDouble()*(bounds[i][j].upper-bounds[i][j].lower);
				} else if (j<i) {
					matrix[i][j]=matrix[j][i];
				}
			}
		}
		return new Matrix(matrix);
	}

	/**
	 * Maps from residue serials to indices of the matrices returned by {@link #getBoundsAllPairs()} and
	 * {@link #sampleBounds(Bound[][])}
	 * @param idx
	 * @return
	 */
	public int getResserFromIdx(int idx) {
		return matIdx2Resser.get(idx);
	}
	
	/*----------------------- private methods  ---------------------------*/
	
	/**
	 * Computes bounds for all pairs, given a graph containing a sparse set of lower/upper bounds
	 * 
	 * @param boundsGraph
	 * @return a 2-dimensional array with the bounds for all pairs of residues, the
	 * indices of the array can be mapped to residue serials through {@link #getResserFromIdx(int)}
	 * and are guaranteed to be in the same order as the residue serials.
	 */
	private Bound[][] getBoundsAllPairs(SparseGraph<Integer,Bound> boundsGraph) {
		Bound[][] bounds = new Bound[conformationSize][conformationSize];
		double[][] lowerBounds = getLowerBoundsAllPairs(boundsGraph);
		double[][] upperBounds = getUpperBoundsAllPairs(boundsGraph);
		for (int i=0;i<lowerBounds.length;i++) {
			for (int j=0;j<lowerBounds[i].length;j++) {
				// we fill in the lower half of the matrix which was missing from upperBounds/lowerBounds
				double upperBound = upperBounds[i][j];
				double lowerBound = lowerBounds[i][j];
				if (i>j) {
					upperBound = upperBounds[j][i];
					lowerBound = lowerBounds[j][i];
				}
				bounds[i][j]=new Bound(lowerBound,upperBound);
				// sanity check: lower bounds can be bigger than upper bounds!, i<j condition is only to do half of the matrix
				if (i<j && lowerBound>lowerBound) {
					System.err.println("Warning: lower bound ("+lowerBound+") for pair "+i+" "+j+" is bigger than upper bound ("+upperBound+")");
				}
			}
		}
		return bounds;
	}

	/**
	 * Computes upper bounds for all pairs from a sparse set of upper bounds based
	 * on the triangle inequality.
	 * The computation is actually just an all pairs shortest path using Dijkstra's 
	 * shortest path algorithm (as the set of distance coming from contact maps is
	 * very sparse this algorithm should be more efficient than Floyd's: Dijkstra's is
	 * o(nm logn) and Floyd's is o(n3))
	 * @return a 2-dimensional array with the upper bounds for all pairs of residues, the
	 * indices of the array can be mapped to residue serials through {@link #getResserFromIdx(int)}
	 * and are guaranteed to be in the same order as the residue serials.
	 */
	private double[][] getUpperBoundsAllPairs(SparseGraph<Integer,Bound> boundsGraph) {
		this.matIdx2Resser = new TreeMap<Integer,Integer>();
		double[][] upperBoundsMatrix = new double[conformationSize][conformationSize];
		SparseGraph<Integer,SimpleEdge> upperBoundGraph = convertBoundsGraphToUpperBoundGraph(boundsGraph);
		DijkstraDistance<Integer, SimpleEdge> dd = new DijkstraDistance<Integer, SimpleEdge>(upperBoundGraph,WeightTransformer);
		int iMatIdx = 0;
		for (int i:rig.getSerials()) {
			int jMatIdx = 0;
			for (int j:rig.getSerials()) {
				if (jMatIdx>iMatIdx) {
					upperBoundsMatrix[iMatIdx][jMatIdx] = dd.getDistance(i, j).doubleValue();
				}
				jMatIdx++;
			}
			this.matIdx2Resser.put(iMatIdx,i);
			iMatIdx++;
		}
		return upperBoundsMatrix;
	}
	
	/**
	 * Computes lower bounds for all pairs given a sparse set of upper/lower bounds based on the triangle inequality.
	 * 
	 * NOTE that because we use Dijkstra's algorithm for the computation of the shortest paths, we can't use negative 
	 * weights (see http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm). Because of this the boundsDigraph result 
	 * of calling {@link #convertBoundsGraphToBoundsDigraph(SparseGraph)}} have values offset with the maximum lower bound (lmax).
	 * Thus after computing the shortest paths we have to revert back that offset by counting the number of hops the shortest path has.
	 *  
	 * @return a 2-dimensional array with the lower bounds for all pairs of residues, the
	 * indices of the array can be mapped to residue serials through {@link #getResserFromIdx(int)}
	 * and are guaranteed to be in the same order as the residue serials.
	 * 
	 * @see {@link #convertBoundsGraphToBoundsDigraph(SparseGraph)} and {@link #getUpperBoundsAllPairs()}
	 */
	private double[][] getLowerBoundsAllPairs(SparseGraph<Integer,Bound> boundsGraph) {
		double[][] lowerBoundsMatrix = new double[conformationSize][conformationSize];		
		// this is the bounds digraph as described by Crippen and Havel
		SparseGraph<BoundsDigraphNode,SimpleEdge> boundsDigraph = convertBoundsGraphToBoundsDigraph(boundsGraph);
		DijkstraShortestPath<BoundsDigraphNode, SimpleEdge> dd = new DijkstraShortestPath<BoundsDigraphNode, SimpleEdge>(boundsDigraph,WeightTransformer);
		int iMatIdx = 0;
		for (int i:rig.getSerials()) {
			int jMatIdx = 0;
			for (int j:rig.getSerials()) {
				if (jMatIdx>iMatIdx) {
					int hops = dd.getPath(nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).get(i), nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).get(j)).size();
					double lower = Math.abs(
						(dd.getDistance(nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).get(i), 
									   nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).get(j)
							           ).doubleValue()) 
						- (hops*lmax)); // the lower limit for the triangle inequality is: Math.abs(shortestpath-(hops*lmax))
					lowerBoundsMatrix[iMatIdx][jMatIdx] = Math.max(lower, HARD_SPHERES_BOUND); // we only set the new lower bound to the one found if is above the HARD_SPHERES_BOUND
				}
				jMatIdx++;
			}
			iMatIdx++;
		}
		return lowerBoundsMatrix;
		
	}

	/**
	 * Converts the bounds graph to a graph with only the upper bounds: nodes residue 
	 * serials, edges upper bounds (in SimpleEdge objects containing the upper bound value 
	 * in their weight field).
	 * @param distanceGraph
	 * @return
	 */
	private SparseGraph<Integer,SimpleEdge> convertBoundsGraphToUpperBoundGraph(SparseGraph<Integer,Bound> distanceGraph) {
		SparseGraph<Integer,SimpleEdge> upperBoundGraph = new SparseGraph<Integer, SimpleEdge>();
		for (Bound bounds:distanceGraph.getEdges()) {
			Pair<Integer> pair = distanceGraph.getEndpoints(bounds);
			upperBoundGraph.addEdge(new SimpleEdge(bounds.upper), pair.getFirst(), pair.getSecond(), EdgeType.UNDIRECTED);
		}
		return upperBoundGraph;
	}

	/**
	 * Constructs a bounds digraph (as described by Crippen and Havel) to compute the triangle inequality limits
	 * for the lower bounds.
	 * The graph is composed by 2 subgraphs (we call them left and right) each of them containing the set of all atoms. 
	 * Within the subgraphs there is an undirected edge between atoms i,j with weight the upper bounds for i,j
	 * Between the subgraphs there is a directed edge from left to right between atoms i(left) to j(right) 
	 * with weight the negative of the lower bound i,j
	 * NOTE that because we use Dijkstra's algorithm for the computation of the shortest paths, we can't use negative 
	 * weights (see http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm). Thus we offset the values here to the maximum lower bound.
	 * After computing the shortest paths we have to revert back that offset by counting the number of hops the shortest path has.
	 * @param distanceGraph
	 * @return
	 */
	private SparseGraph<BoundsDigraphNode,SimpleEdge> convertBoundsGraphToBoundsDigraph(SparseGraph<Integer,Bound> distanceGraph) {
		// to do the offset thing (see docs above) we need to know first of all the max lower bound
		ArrayList<Double> lowerBounds = new ArrayList<Double>();
		for (Bound bounds:distanceGraph.getEdges()) {
			lowerBounds.add(bounds.lower);
		}
		lmax = Collections.max(lowerBounds); // this is the offset value
		
		SparseGraph<BoundsDigraphNode,SimpleEdge> boundsDigraph = new SparseGraph<BoundsDigraphNode, SimpleEdge>();
		// we have to store all nodes in a HashMap, so we can retrieve them by residue serial and side after 
		nodesBoundsDigraph = new HashMap<Boolean, HashMap<Integer,BoundsDigraphNode>>();
		nodesBoundsDigraph.put(BoundsDigraphNode.LEFT , new HashMap<Integer, BoundsDigraphNode>());
		nodesBoundsDigraph.put(BoundsDigraphNode.RIGHT, new HashMap<Integer, BoundsDigraphNode>());
		// first we create the nodes and store them into the HashMap
		for (int i:distanceGraph.getVertices()) {
			BoundsDigraphNode leftNode = new BoundsDigraphNode(i, BoundsDigraphNode.LEFT);
			BoundsDigraphNode rightNode = new BoundsDigraphNode(i, BoundsDigraphNode.RIGHT);
			boundsDigraph.addVertex(leftNode);
			boundsDigraph.addVertex(rightNode);
			nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).put(i,leftNode);
			nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).put(i,rightNode);
		}
		
		for (Bound bounds:distanceGraph.getEdges()) {
			Pair<Integer> pair = distanceGraph.getEndpoints(bounds);
			// first we add the upper bounds as undirected edges to the 2 subgraphs (left and right)
			boundsDigraph.addEdge(new SimpleEdge(lmax+bounds.upper), 
					nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).get(pair.getFirst()), 
					nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).get(pair.getSecond()), 
					EdgeType.UNDIRECTED);
			boundsDigraph.addEdge(new SimpleEdge(lmax+bounds.upper), 
					nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).get(pair.getFirst()), 
					nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).get(pair.getSecond()), 
					EdgeType.UNDIRECTED);
			// then we add the negative of the lower bounds as directed edges connecting nodes of subgraph left to subgraph right
			boundsDigraph.addEdge(new SimpleEdge(lmax-bounds.lower), 
					nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).get(pair.getFirst()), 
					nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).get(pair.getSecond()), 
					EdgeType.DIRECTED);			
		}
		return boundsDigraph;
	}

	/**
	 * Convert the given RIGraph to a bounds graph: residue serials as nodes and distance bounds as edges.
	 * Will only admit single atom contact type RIGraphs
	 * @param graph
	 * @return
	 * @throws IllegalArgumentException if contact type of given RIGraph is not a single atom contact type
	 */
	private SparseGraph<Integer,Bound> convertRIGraphToDistRangeGraph(RIGraph graph) {
		// code cloned from ConstraintsMaker.createDistanceConstraints with some modifications
		
		SparseGraph<Integer, Bound> distanceGraph = new SparseGraph<Integer, Bound>();
		
		double cutoff = graph.getCutoff();
		String ct = graph.getContactType();
		String i_ct = ct;
		String j_ct = ct;
		if (ct.contains("/")){
			i_ct = ct.split("/")[0];
			j_ct = ct.split("/")[1];
		}
		
		if (!AAinfo.isValidSingleAtomContactType(i_ct) || !AAinfo.isValidSingleAtomContactType(j_ct)){
			throw new IllegalArgumentException("Contact type "+i_ct+" or "+j_ct+" is not valid for reconstruction");
		}
		
		for (RIGEdge cont:graph.getEdges()){
			Pair<RIGNode> pair = graph.getEndpoints(cont);
			String i_res = pair.getFirst().getResidueType();
			String j_res = pair.getSecond().getResidueType();

			// as dist_min we take the average of the two dist mins, if i_ct and j_ct are the same then this will be the same as dist_min for ct
			double dist_min = (AAinfo.getLowerBoundDistance(i_ct,i_res,j_res)+AAinfo.getLowerBoundDistance(j_ct,i_res,j_res))/2;
			// for single atom contact types getUpperBoundDistance and getLowerBoundDistance will return 0 thus for those cases dist_max = cutoff
			double dist_max = AAinfo.getUpperBoundDistance(i_ct, i_res, j_res)/2+AAinfo.getUpperBoundDistance(i_ct, i_res, j_res)/2+cutoff;
			
			Bound bounds = new Bound(dist_min, dist_max);
			if (pair.getSecond().getResidueSerial()>pair.getFirst().getResidueSerial()+1) { //we don't add the first diagonal, we add it later as contiguous CA constraints 
				distanceGraph.addEdge(bounds, pair.getFirst().getResidueSerial(), pair.getSecond().getResidueSerial(), EdgeType.UNDIRECTED);
			}
		}
		// adding contiguous CA distance backbone constraints
		for (int i:graph.getSerials()) {
			if (i!=graph.getLastResidueSerial()) {
				distanceGraph.addEdge(new Bound(BB_CA_DIST,BB_CA_DIST), i, i+1, EdgeType.UNDIRECTED);
			}
		}
		return distanceGraph;
	}
	
	/*-------------------------- main  -------------------------------*/
	
	/**
	 * To test the class
	 */
	public static void main (String[] args) throws Exception {
		boolean debug = true;
		int numModels = 1;
		String pdbCode = "1bxy";
		String pdbChainCode = "A";
		
		Pdb pdb = new PdbasePdb(pdbCode);
		pdb.load(pdbChainCode);
		Pdb pdbEmbedded = new PdbasePdb(pdbCode);
		pdbEmbedded.load(pdbChainCode);

		RIGraph graph = pdb.get_graph("Ca", 8);
		BoundsSmoother bs = new BoundsSmoother(graph);
		Bound[][] bounds = bs.getBoundsAllPairs();
		if (debug) {
			for (int i=0;i<bounds.length;i++) {
				for (int j=0;j<bounds[i].length;j++) {
					System.out.print(bounds[i][j]);
				}
				System.out.println();
			}
			System.out.println();
		}

		System.out.printf("%6s\t%6s", "rmsd","rmsdm");
		System.out.println();
		for (int model=0;model<numModels;model++) {
			Matrix matrix = bs.sampleBounds(bounds);
			if (debug) {
				for (int i=0;i<matrix.getRowDimension();i++) {
					for (int j=0;j<matrix.getColumnDimension();j++) {
						System.out.printf("%4.1f ",matrix.get(i, j));
					}
					System.out.println();
				}
			}

			int size = pdb.get_length();
			Embedder embedder = new Embedder(matrix,Embedder.createTrivialVector(1.0, size), Embedder.createTrivialVector(1.0, size));
			Vector3d[] embedding = embedder.embed();
			pdbEmbedded.setAtomsCoords(embedding, "CA");

			double rmsd = pdb.rmsd(pdbEmbedded, "Ca");
			pdb.mirror();
			double rmsdm = pdb.rmsd(pdbEmbedded, "Ca");
			
			if (rmsd>rmsdm ) {
				pdbEmbedded.mirror();
			}
			pdbEmbedded.dump2pdbfile("/project/StruPPi/jose/embed_"+pdbCode+pdbChainCode+".pdb");
			
			System.out.printf("%6.3f\t%6.3f",rmsd,rmsdm);
			System.out.println();
			
			Matrix matrixEmbedded = pdbEmbedded.calculateDistMatrix("Ca");
			for (int i=0;i<matrixEmbedded.getRowDimension();i++) {
				for (int j=i+1;j<matrixEmbedded.getColumnDimension();j++) {
					if ((matrixEmbedded.get(i,j)<bounds[i][j].lower) || (matrixEmbedded.get(i,j)>bounds[i][j].upper)) {
						System.out.printf("%3d %3d %4.1f %s\n",i,j,matrixEmbedded.get(i,j),bounds[i][j].toString());
					}
				}
			}
//			IntPairSet violEdges = TinkerRunner.getViolatedEdges(graph, pdbEmbedded);
//			RIGraph graphEmbedded = pdbEmbedded.get_graph("Ca", 8);
//			for (Pair<Integer> violEdge:violEdges) {
//				System.out.println(violEdge+" "+graphEmbedded.getEdgeFromSerials(violEdge.getFirst(), violEdge.getSecond()).getDistance());
//			}
		}
	}
	
}

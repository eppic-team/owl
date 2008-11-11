package embed;

import java.util.Random;
import java.util.TreeMap;

import javax.vecmath.Vector3d;

import org.apache.commons.collections15.Transformer;

import Jama.Matrix;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraDistance;
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
	

	/*----------------- helper classes and transformers -----------------*/
	
	Transformer<SimpleEdge, Double> WeightTransformer = new Transformer<SimpleEdge, Double>() {
		public Double transform(SimpleEdge input) {
			return input.weight;
		}
	};
	
	private class SimpleEdge {
		public double weight;
		public SimpleEdge(double weight) {
			this.weight = weight;
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
	
	private SparseGraph<Integer,Bound> boundsGraph;
	private RIGraph rig;
	private TreeMap<Integer,Integer> matIdx2Resser;
	private int conformationSize;
	
	/*------------------------ constructors ------------------------------*/
	
	/**
	 * Constructs a new BoundsSmoother object given a RIGraph.
	 * The RIGraph is converted into a set of distance ranges using the cutoff
	 * as upper limit and hard-spheres as lower limit. 
	 * @param graph
	 */
	public BoundsSmoother(RIGraph graph) {
		this.rig = graph;
		this.boundsGraph = convertRIGraphToDistRangeGraph(graph);
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
		Bound[][] bounds = new Bound[conformationSize][conformationSize];
		double[][] lowerBounds = getLowerBoundsAllPairs();
		double[][] upperBounds = getUpperBoundsAllPairs();
		for (int i=0;i<lowerBounds.length;i++) {
			for (int j=0;j<lowerBounds[i].length;j++) {
				// we fill in the lower half of the upper bounds matrix which was missing from upperBounds
				double upperBound = upperBounds[i][j];
				if (i>j) upperBound = upperBounds[j][i];
				bounds[i][j]=new Bound(lowerBounds[i][j],upperBound);
			}
		}
		return bounds;
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
	private double[][] getUpperBoundsAllPairs() {
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
	 * Computes lower bounds for all pairs given a sparse set of upper/lower bounds
	 * At the moment it doesn't compute anything but just copies the first lower bound 
	 * and uses that value for all pairs.
	 * 
	 * @return a 2-dimensional array with the lower bounds for all pairs of residues, the
	 * indices of the array can be mapped to residue serials through {@link #getResserFromIdx(int)}
	 * and are guaranteed to be in the same order as the residue serials.
	 */
	private double[][] getLowerBoundsAllPairs() {
		// for now we implement the simplest approach: hard-spheres for all lower bounds
		// in the simplest way... we take the first lower bound from the boundsGraph as the lower bound for all pairs
		double lowerBound = 0;
		for (Bound bounds: boundsGraph.getEdges()){
			lowerBound = bounds.lower;
			break;
		}
		double[][] lowerBoundMatrix = new double[conformationSize][conformationSize];
		for (int iMatIdx=0;iMatIdx<conformationSize;iMatIdx++) {
			for (int jMatIdx=0;jMatIdx<conformationSize;jMatIdx++) {
				if (iMatIdx!=jMatIdx) {
					lowerBoundMatrix[iMatIdx][jMatIdx] = lowerBound;
				}
			}
		}
		return lowerBoundMatrix;
		
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
			distanceGraph.addEdge(bounds, pair.getFirst().getResidueSerial(), pair.getSecond().getResidueSerial(), EdgeType.UNDIRECTED);
		}
		return distanceGraph;
	}
	
	/*-------------------------- main  -------------------------------*/
	
	/**
	 * To test the class
	 */
	public static void main (String[] args) throws Exception {
		Pdb pdb = new PdbasePdb("1bxy");
		pdb.load("A");
		RIGraph graph = pdb.get_graph("Ca", 8);
		BoundsSmoother bs = new BoundsSmoother(graph);
		Bound[][] bounds = bs.getBoundsAllPairs();
		for (int i=0;i<bounds.length;i++) {
			for (int j=0;j<bounds[i].length;j++) {
				System.out.print(bounds[i][j]);
			}
			System.out.println();
		}
		System.out.println();
		
		Matrix matrix = bs.sampleBounds(bounds);
		for (int i=0;i<matrix.getRowDimension();i++) {
			for (int j=0;j<matrix.getColumnDimension();j++) {
				System.out.printf("%4.1f ",matrix.get(i, j));
			}
			System.out.println();
		}
		int size = pdb.get_length();
		Embedder embedder = new Embedder(matrix,Embedder.createTrivialVector(1.0, size), Embedder.createTrivialVector(1.0, size));
		Vector3d[] embedding = embedder.embed();
		Vector3d[] originalConformation = new Vector3d[size];
		for (int i=0;i<size;i++) {
			originalConformation[i]=new Vector3d(pdb.getAtomCoord(bs.getResserFromIdx(i), "CA"));
		}
		
		double rmsd = Pdb.calculate_rmsd(originalConformation, embedding);
		for (int i=0;i<originalConformation.length;i++){
			originalConformation[i].scale(-1);
		}
		double rmsdm = Pdb.calculate_rmsd(originalConformation, embedding);
		System.out.println("rmsd of embedded to original conformation: "+rmsd);
		System.out.println("rmsd of embedded to mirrored original conformation: "+rmsdm);

	}
	
}

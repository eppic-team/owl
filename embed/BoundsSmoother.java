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
	
	protected static final double BB_CA_DIST = 3.8;
	private static final double HARD_SPHERES_BOUND = AAinfo.DIST_MIN_CA ;
	
	private static final boolean DEBUG = false;
	private static final long DEBUG_SEED = 123456;
	
	private static final int NUM_ROOTS_PARTIAL_METRIZATION = 4;

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

	/*---------------------- member variables ----------------------------*/
	
	private RIGraph rig;
	private TreeMap<Integer,Integer> idx2resser;
	private int conformationSize;
	private HashMap<Boolean, HashMap<Integer,BoundsDigraphNode>> nodesBoundsDigraph; // map of serial/side to nodes in the bounds digraph
	private double lmax; // maximum of the lower bounds: offset value for the boundsDigraph (not to have negative weights so that we can use Dijkstra's algo)
	
	private Random rand; // the random generator for sampleBounds and metrize
	
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
		if (DEBUG) {
			rand = new Random(DEBUG_SEED);
		} else {
			rand = new Random();
		}
	}
	
	/*----------------------- public methods  ----------------------------*/
	
	/**
	 * Computes bounds for all pairs, based on the set of sparse distance ranges
	 * @return a 2-dimensional array with the bounds for all pairs of residues, the
	 * indices of the array can be mapped to residue serials through {@link #getResserFromIdx(int)}
	 * and are guaranteed to be in the same order as the residue serials.
	 */
	public Bound[][] getBoundsAllPairs() {
		Bound[][] bounds = convertRIGraphToBoundsMatrix(this.rig);
		return getBoundsAllPairs(bounds);
	}
		
	/**
	 * Gets a random sample from a matrix of all pairs distance ranges
	 * @param bounds the matrix of all pairs distance ranges
	 * @return a simmetric metric matrix (both sides filled)
	 */
	public Matrix sampleBounds(Bound[][] bounds) {
		double[][] matrix = new double[conformationSize][conformationSize];
		for (int i=0;i<conformationSize;i++) {
			for (int j=0;j<conformationSize;j++) {
				if (j>i) {
					matrix[i][j]= bounds[i][j].lower+rand.nextDouble()*(bounds[i][j].upper-bounds[i][j].lower);
				} else if (j<i) {
					matrix[i][j]=matrix[j][i];
				}
			}
		}
		return new Matrix(matrix);
	}

	/**
	 * Performs partial metrization for the given bounds 
	 * @param bounds
	 * @return
	 */
	public Matrix metrize(Bound[][] bounds) {

		ArrayList<Integer> roots = new ArrayList<Integer>();
		for (int count=1;count<=NUM_ROOTS_PARTIAL_METRIZATION;count++) {
			int root = rand.nextInt(conformationSize);
			if (roots.contains(root)) System.err.println("Warning: repeated root atom while doing metrization: "+root);
			if (DEBUG) System.out.println("Picked root: "+root);			
			roots.add(root);
			sampleBoundForRoot(bounds, root); // this alters directly the input bounds array
			bounds = getBoundsAllPairs(bounds);
		}
		// bounds after metrization
		if (DEBUG) printBounds(bounds);
		
		// finally pick a value at random for all the other bounds
		return sampleBounds(bounds);		
	}
	
	/**
	 * Maps from residue serials to indices of the matrices returned by {@link #getBoundsAllPairs()} and
	 * {@link #sampleBounds(Bound[][])}
	 * @param idx
	 * @return
	 */
	public int getResserFromIdx(int idx) {
		return idx2resser.get(idx);
	}
	
	/*----------------------- private methods  ---------------------------*/
	
	private void sampleBoundForRoot(Bound[][] bounds, int root) {
		for (int neighb=0;neighb<conformationSize;neighb++) {
			if (neighb==root) continue; // avoid the diagonal (which contains a null Bound)
			int i = root;
			int j = neighb;
			if (neighb<root) {
				i = neighb;
				j = root;
			}
			double sampledValue = bounds[i][j].lower+rand.nextDouble()*(bounds[i][j].upper-bounds[i][j].lower);
			bounds[i][j].lower = sampledValue;
			bounds[i][j].upper = sampledValue;
			if (DEBUG) System.out.print(bounds[i][j]);
		}
		if (DEBUG) System.out.println("\n");
	}
	
	/**
	 * Computes bounds for all pairs through triangle inequalities, given a Matrix containing a sparse set of lower/upper bounds
	 * 
	 * @param sparseBounds
	 * @return a 2-dimensional array with the bounds for all pairs of residues (only upper half filled), 
	 * the indices of the array can be mapped to residue serials through {@link #getResserFromIdx(int)}
	 * and are guaranteed to be in the same order as the residue serials.
	 */
	private Bound[][] getBoundsAllPairs(Bound[][] sparseBounds) { 
		double MARGIN = 0.0001; // for comparing doubles we need some tolerance value
		Bound[][] bounds = new Bound[conformationSize][conformationSize];
		double[][] lowerBounds = getLowerBoundsAllPairs(sparseBounds);
		double[][] upperBounds = getUpperBoundsAllPairs(sparseBounds);
		for (int i=0;i<conformationSize;i++) {
			for (int j=i+1;j<conformationSize;j++) {
				double upperBound = upperBounds[i][j];
				double lowerBound = lowerBounds[i][j];
				if (sparseBounds[i][j]!=null && lowerBound>upperBound+MARGIN) {
					//System.err.println("old: "+sparseBounds[i][j]+" new: "+new Bound(lowerBound,upperBound));
					
					// During metrization sometimes a new upper bound is found that is below the new lower bound 
					// (actually in these cases the new lower bound coincides with the old one i.e. nothing new was 
					// found through triangle inequality for the lower bound).
					// For some reason it doesn't happen the other way around: a new lower bound found that is 
					// above the new (coinciding with old) upper bound. I suppose this is because the triangle inequality 
					// "is a lot more effective at reducing the upper bounds than increasing the lower bounds" (quoting Havel) 
					// To correct this we set both lower and upper to the newly found upper, i.e. we assume that 
					// the new upper bound is better because is in accordance to the triangle inequality 
					lowerBound=upperBound;
					
					//if (upperBound<sparseBounds[i][j].upper-MARGIN) 
					//	System.err.printf("new upper bound (%4.1f) for pair %3d %3d is smaller than old upper bound (%4.1f)\n",upperBound,i,j,sparseBounds[i][j].upper);
					//if (lowerBound>sparseBounds[i][j].lower+MARGIN) 
					//	System.err.printf("new lower bound (%4.1f) for pair %3d %3d is bigger than old lower bound (%4.1f)\n",lowerBound,i,j,sparseBounds[i][j].lower);
				}
					
				bounds[i][j]=new Bound(lowerBound,upperBound);
				// sanity check: lower bounds can't be bigger than upper bounds!
				if (lowerBound>upperBound+MARGIN) {
					System.err.printf("Warning: lower bound (%4.1f) for pair "+i+" "+j+" is bigger than upper bound (%4.1f)\n",lowerBound,upperBound);
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
	private double[][] getUpperBoundsAllPairs(Bound[][] bounds) {
		double[][] upperBoundsMatrix = new double[conformationSize][conformationSize];
		SparseGraph<Integer,SimpleEdge> upperBoundGraph = convertBoundsMatrixToUpperBoundGraph(bounds);
		DijkstraDistance<Integer, SimpleEdge> dd = new DijkstraDistance<Integer, SimpleEdge>(upperBoundGraph,WeightTransformer);

		for (int i=0;i<conformationSize;i++) {
			for (int j=i+1;j<conformationSize;j++) {
				upperBoundsMatrix[i][j] = dd.getDistance(i, j).doubleValue();
			}
		}
		return upperBoundsMatrix;
	}
	
	/**
	 * Computes lower bounds for all pairs given a sparse set of upper/lower bounds based on the triangle inequality.
	 * 
	 * NOTE that because we use Dijkstra's algorithm for the computation of the shortest paths, we can't use negative 
	 * weights (see http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm). Because of this the boundsDigraph result 
	 * of calling {@link #convertBoundsMatrixToBoundsDigraph(Bound[][])}} have values offset with the maximum lower bound (lmax).
	 * Thus after computing the shortest paths we have to revert back that offset by counting the number of hops the shortest path has.
	 *  
	 * @return a 2-dimensional array with the lower bounds for all pairs of residues, the
	 * indices of the array can be mapped to residue serials through {@link #getResserFromIdx(int)}
	 * and are guaranteed to be in the same order as the residue serials.
	 * 
	 * @see {@link #convertBoundsMatrixToBoundsDigraph(Bound[][])} and {@link #getUpperBoundsAllPairs()}
	 */
	private double[][] getLowerBoundsAllPairs(Bound[][] bounds) {
		double[][] lowerBoundsMatrix = new double[conformationSize][conformationSize];		
		// this is the bounds digraph as described by Crippen and Havel
		SparseGraph<BoundsDigraphNode,SimpleEdge> boundsDigraph = convertBoundsMatrixToBoundsDigraph(bounds);
		DijkstraShortestPath<BoundsDigraphNode, SimpleEdge> dd = new DijkstraShortestPath<BoundsDigraphNode, SimpleEdge>(boundsDigraph,WeightTransformer);

		for (int i=0;i<conformationSize;i++) {
			for (int j=i+1;j<conformationSize;j++) {
				int hops = dd.getPath(nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).get(i), nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).get(j)).size();
				double lower = Math.abs(
						(dd.getDistance(nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).get(i), 
								nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).get(j)
						).doubleValue()) 
						- (hops*lmax)); // the lower limit for the triangle inequality is: Math.abs(shortestpath-(hops*lmax))
				lowerBoundsMatrix[i][j] = Math.max(lower, HARD_SPHERES_BOUND); // we only set the new lower bound to the one found if is above the HARD_SPHERES_BOUND
			}
		}
		return lowerBoundsMatrix;
		
	}

	/**
	 * Converts the bounds matrix to a graph with only the upper bounds: nodes indices of bounds 
	 * matrix, edges upper bounds (in SimpleEdge objects containing the upper bound value 
	 * in their weight field).
	 * @param bounds
	 * @return
	 */
	private SparseGraph<Integer,SimpleEdge> convertBoundsMatrixToUpperBoundGraph(Bound[][] bounds) {
		SparseGraph<Integer,SimpleEdge> upperBoundGraph = new SparseGraph<Integer, SimpleEdge>();
		for (int i=0;i<conformationSize;i++) {
			for (int j=0;j<conformationSize;j++) {
				if (bounds[i][j]!=null) {
					upperBoundGraph.addEdge(new SimpleEdge(bounds[i][j].upper), i, j, EdgeType.UNDIRECTED);
				}
			}
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
	 * @param bounds
	 * @return
	 */
	private SparseGraph<BoundsDigraphNode,SimpleEdge> convertBoundsMatrixToBoundsDigraph(Bound[][] bounds) {
		// to do the offset thing (see docs above) we need to know first of all the max lower bound
		ArrayList<Double> lowerBounds = new ArrayList<Double>();
		for (int i=0;i<conformationSize;i++) {
			for (int j=0;j<conformationSize;j++) {
				if (bounds[i][j]!=null) {
					lowerBounds.add(bounds[i][j].lower);
				}
			}
		}
		lmax = Collections.max(lowerBounds); // this is the offset value
		
		SparseGraph<BoundsDigraphNode,SimpleEdge> boundsDigraph = new SparseGraph<BoundsDigraphNode, SimpleEdge>();
		// we have to store all nodes in a HashMap, so we can retrieve them by residue serial and side after 
		nodesBoundsDigraph = new HashMap<Boolean, HashMap<Integer,BoundsDigraphNode>>();
		nodesBoundsDigraph.put(BoundsDigraphNode.LEFT , new HashMap<Integer, BoundsDigraphNode>());
		nodesBoundsDigraph.put(BoundsDigraphNode.RIGHT, new HashMap<Integer, BoundsDigraphNode>());
		// first we create the nodes and store them into the HashMap
		for (int i=0;i<conformationSize;i++) {
			BoundsDigraphNode leftNode = new BoundsDigraphNode(i, BoundsDigraphNode.LEFT);
			BoundsDigraphNode rightNode = new BoundsDigraphNode(i, BoundsDigraphNode.RIGHT);
			boundsDigraph.addVertex(leftNode);
			boundsDigraph.addVertex(rightNode);
			nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).put(i,leftNode);
			nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).put(i,rightNode);
		}
		
		for (int i=0;i<conformationSize;i++) {
			for (int j=0;j<conformationSize;j++) {
				if (bounds[i][j]!=null) {
					// first we add the upper bounds as undirected edges to the 2 subgraphs (left and right)
					boundsDigraph.addEdge(new SimpleEdge(lmax+bounds[i][j].upper), 
							nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).get(i), 
							nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).get(j), 
							EdgeType.UNDIRECTED);
					boundsDigraph.addEdge(new SimpleEdge(lmax+bounds[i][j].upper), 
							nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).get(i), 
							nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).get(j), 
							EdgeType.UNDIRECTED);
					// then we add the negative of the lower bounds as directed edges connecting nodes of subgraph left to subgraph right
					boundsDigraph.addEdge(new SimpleEdge(lmax-bounds[i][j].lower), 
							nodesBoundsDigraph.get(BoundsDigraphNode.LEFT).get(i), 
							nodesBoundsDigraph.get(BoundsDigraphNode.RIGHT).get(j), 
							EdgeType.DIRECTED);		
				}
			}
		}
		return boundsDigraph;
	}

	/**
	 * Convert the given RIGraph to a bounds matrix. The indices of the matrix can be mapped back to residue 
	 * serials through {@link #getResserFromIdx(int)}
	 * Will only admit single atom contact type RIGraphs
	 * @param graph
	 * @return
	 * @throws IllegalArgumentException if contact type of given RIGraph is not a single atom contact type
	 */
	private Bound[][] convertRIGraphToBoundsMatrix(RIGraph graph) {
		// code cloned from ConstraintsMaker.createDistanceConstraints with some modifications
		Bound[][] bounds = new Bound[conformationSize][conformationSize];
		TreeMap<Integer,Integer> resser2idx = new TreeMap<Integer, Integer>();
		this.idx2resser = new TreeMap<Integer, Integer>();
		int idx = 0;
		for (int resser:graph.getSerials()) {
			resser2idx.put(resser,idx);
			idx2resser.put(idx,resser);
			idx++;
		}
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
			
			if (pair.getSecond().getResidueSerial()>pair.getFirst().getResidueSerial()+1) { //we don't add the first diagonal, we add it later as contiguous CA constraints 
				bounds[resser2idx.get(pair.getFirst().getResidueSerial())][resser2idx.get(pair.getSecond().getResidueSerial())] = new Bound(dist_min, dist_max);
			}
		}
		// adding contiguous CA distance backbone constraints
		for (int i=0;i<conformationSize-1;i++) {
			bounds[i][i+1]=new Bound(BB_CA_DIST,BB_CA_DIST);
		}
		return bounds;
	}
	
	/**
	 * Deep copies given array of bounds
	 * @param bounds
	 * @return
	 */
	private static Bound[][] copyBounds(Bound[][] bounds) {
		Bound[][] newBounds = new Bound[bounds.length][bounds.length];
		for (int i=0;i<bounds.length;i++) {
			for (int j=0;j<bounds[i].length;j++) {
				if (bounds[i][j]!=null) {
					newBounds[i][j] = new Bound(bounds[i][j].lower,bounds[i][j].upper);
				}
			}
		}
		return newBounds;
	}
	
	/*------------------------ statics  ------------------------------*/
	
	private static void printBounds(Bound[][] bounds) {
		for (int i=0;i<bounds.length;i++) {
			for (int j=0;j<bounds[i].length;j++) {
				if (bounds[i][j]==null) {
					System.out.printf("%11s","null");
				} else {
					System.out.print(bounds[i][j]);
				}
			}
			System.out.println();
		}
		System.out.println();
	}
	
	private static void printViolations (Matrix matrixEmbedded, Bound[][] bounds) {
		int count = 0;
		for (int i=0;i<matrixEmbedded.getRowDimension();i++) {
			for (int j=i+1;j<matrixEmbedded.getColumnDimension();j++) {
				if ((matrixEmbedded.get(i,j)<bounds[i][j].lower) || (matrixEmbedded.get(i,j)>bounds[i][j].upper)) {
					System.out.printf("%3d %3d %4.1f %s\n",i,j,matrixEmbedded.get(i,j),bounds[i][j].toString());
					count++;
				}
			}
		}
		int cells = (matrixEmbedded.getRowDimension()*(matrixEmbedded.getRowDimension()-1))/2;
		System.out.println("Number of violations: "+count+" out of "+cells+" cells in half matrix");
	}
	
	private static int getNumberViolations(Matrix matrixEmbedded, Bound[][] bounds) {
		int count = 0;
		for (int i=0;i<matrixEmbedded.getRowDimension();i++) {
			for (int j=i+1;j<matrixEmbedded.getColumnDimension();j++) {
				if ((matrixEmbedded.get(i,j)<bounds[i][j].lower) || (matrixEmbedded.get(i,j)>bounds[i][j].upper)) {
					count++;
				}
			}
		}
		return count;
	}

	private static void printMatrix(Matrix matrix) {
		for (int i=0;i<matrix.getRowDimension();i++) {
			for (int j=0;j<matrix.getColumnDimension();j++) {
				System.out.printf("%4.1f ",matrix.get(i, j));
			}
			System.out.println();
		}
	}
	
	/*-------------------------- main  -------------------------------*/
	
	/**
	 * To test the class
	 */
	public static void main (String[] args) throws Exception {
		boolean debug = false;
		boolean writeFiles = false;
		int numModels = 5;
		Embedder.ScalingMethod scalingMethod = Embedder.ScalingMethod.AVRG_INTER_CA_DIST;
		
		String pdbCode = "1bxy";
		String pdbChainCode = "A";
		String ct = "Ca";
		double cutoff = 8.0;
		
		Pdb pdb = new PdbasePdb(pdbCode);
		pdb.load(pdbChainCode);
		Pdb pdbEmbedded = new PdbasePdb(pdbCode);
		pdbEmbedded.load(pdbChainCode);

		RIGraph graph = pdb.get_graph(ct, cutoff);
		BoundsSmoother bs = new BoundsSmoother(graph);
		Bound[][] bounds = bs.getBoundsAllPairs();

		Bound[][] initialAllPairs = copyBounds(bounds);
		if (debug) {
			// all pais bounds after triangle inequality
			printBounds(initialAllPairs);
		}
		
		
		System.out.printf("%6s\t%6s\t%6s", "rmsd","rmsdm","viols");
		System.out.println();
		for (int model=0;model<numModels;model++) {
			//Matrix matrix = bs.sampleBounds(bounds);
			
			Matrix matrix = bs.metrize(bounds);
			//if (debug) {
				// bounds after metrization
				//printBounds(bounds);
			//}
			
			if (debug) {
				printMatrix(matrix);
			}

			Embedder embedder = new Embedder(matrix);
			Vector3d[] embedding = embedder.embed(scalingMethod);
			pdbEmbedded.setAtomsCoords(embedding, "CA");

			double rmsd = pdb.rmsd(pdbEmbedded, "Ca");
			pdb.mirror();
			double rmsdm = pdb.rmsd(pdbEmbedded, "Ca");

			if (rmsd>rmsdm ) {
				pdbEmbedded.mirror();
			}

			if (writeFiles)
				pdbEmbedded.dump2pdbfile("/project/StruPPi/jose/embed_"+pdbCode+pdbChainCode+"_"+model+".pdb");

			Matrix matrixEmbedded = pdbEmbedded.calculateDistMatrix("Ca");
			
			System.out.printf("%6.3f\t%6.3f\t%6d",rmsd,rmsdm,getNumberViolations(matrixEmbedded, initialAllPairs));
			System.out.println();
			
			if (debug) {
				printViolations(matrixEmbedded, initialAllPairs);			
			}
		}
	}
	
}

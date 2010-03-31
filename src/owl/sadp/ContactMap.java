package owl.sadp;
import java.util.StringTokenizer;

import owl.core.structure.RIGEdge;
import owl.core.structure.RIGNode;
import owl.core.structure.RIGraph;

import edu.uci.ics.jung.graph.util.Pair;


/**
 * This class contains methods and fields to represent contact maps.  
 */
public class ContactMap {

	/**
	 * Adjacency matrix
	 */
	public boolean[][] A;

	/**
	 * Adjacency list
	 */
	public int[][] L;

	/**
	 * Number of nodes
	 */
	public int nNodes;

	/**
	 * Number of edges
	 */
	public int nEdges;

	/**
	 * Degree sequence
	 */
	public int[] deg;


	/**
	 * File name
	 */
	public String fn;


	/**
	 * Constructor.
	 **/
	public ContactMap(boolean[][] A) {

		this.A = A;
		this.nNodes = A.length;
		this.nEdges = 0;
		this.deg = new int[nNodes];
		this.fn = "NoName";

		// set number of edges and degree sequence
		for(int i = 0; i < nNodes; i++) {
			for(int j = i+1; j < nNodes; j++) {
				if(A[i][j]) {
					nEdges++;
					deg[i]++;
					deg[j]++;
				}
			}
		}

		// create adjacency list
		this.L = new int[nNodes][];
		for(int i = 0; i < nNodes; i++) {
			L[i] = new int[deg[i]];
			int k = 0; 
			for(int j = 0; j < nNodes; j++) {
				if(i == j) {
					continue;
				}
				if(A[i][j]) {
					L[i][k] = j;
					k++;
				}
			}
			if(k != deg[i]) {
				System.err.println("Error!");
				System.err.println("\t Mismatch of degrees in ContactMap()");
				System.exit(0);
			}
		}	
	}

	/**
	 * Constructs a contact map from a RIGraph. 
	 * This constructor is for CMView only. It does not provide the 
	 * System.exit(0) call if the given graph was inconsistent.
	 * Requirements to the format of the input graph:
	 * <ul>
	 * <li>smallest node index has to be larger than 0</li>
	 * <li>...</li>
	 * </ul>  
	 */
	public ContactMap(RIGraph g) throws ContactMapConstructorError {
		this.nNodes = g.getFullLength();//g.getNodes.size();
		this.nEdges = g.getEdgeCount();
		this.deg    = new int[nNodes];
		this.fn     = "NoName";
		this.A      = new boolean[nNodes][nNodes];
		this.L      = new int[nNodes][];

		// fill the adjacency matrix and the node degree array
		for( RIGEdge e : g.getEdges() ) {
			Pair<RIGNode> pair = g.getEndpoints(e);
			if( pair.getFirst().getResidueSerial() < 0 || pair.getSecond().getResidueSerial() < 0 ) {
				throw new ContactMapConstructorError("Graph contains negative node indices. I don't like it!");
			}

			A[pair.getFirst().getResidueSerial()-1][pair.getSecond().getResidueSerial()-1] = true;
			A[pair.getSecond().getResidueSerial()-1][pair.getFirst().getResidueSerial()-1] = true;

			++deg[pair.getFirst().getResidueSerial()-1];
			++deg[pair.getSecond().getResidueSerial()-1];
		}

		// create adjacency lists
		for(int i = 0; i < nNodes; ++i) {
			L[i] = new int[deg[i]];
			int k = 0; 
			for(int j = 0; j < nNodes; j++) {
				if(i == j) {
					continue;
				}
				if(A[i][j]) {
					L[i][k] = j;
					k++;
				}
			}
			if(k != deg[i]) {
				System.err.println("Error!");
				System.err.println("\t Mismatch of degrees in ContactMap(): def["+i+"]!="+k);
			}
		}
	}

	/**
	 * Returns adjacency matrix of this contact map.
	 * 
	 * @return adjacency matrix
	 */
	public boolean[][] getAdjacencyMatrix() {

		return A;
	}
	/**
	 * Returns adjacency list of this contact map.
	 * 
	 * @return adjacency list
	 */
	public int[][] getAdjacencyList() {

		return L;
	}

	/**
	 * Returns number of vertices of this graph.
	 * 
	 * @return number of vertices
	 */
	public int countNodes() {

		return nNodes;
	}

	/**
	 * Returns number of edges.
	 * 
	 * @return number of edges
	 */
	public int countEdges() {

		return nEdges;
	}

	public String getFileName() {

		return fn;
	}

	public void setFileName(String fn) {

		StringTokenizer st = new StringTokenizer(fn, ".");
		this.fn = st.nextToken();
	}

	public void show() {
		System.out.println(this.fn);
		System.out.println("\t Number of Nodes:" + this.nNodes);
		System.out.println("\t Number of Edges:" + this.nEdges);
	}
}
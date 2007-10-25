package sadp;
import java.util.StringTokenizer;

import proteinstructure.Edge;
import proteinstructure.Graph;

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
     * Constructs a contact map from a Graph. 
     * This constructor is for CMView only. It does not provide the 
     * System.exit(0) call if the given graph was inconsistent. 
     */
    public ContactMap(Graph g) throws ContactMapConstructorError {
	this.nNodes = g.getFullLength();//g.getNodes.size();
	this.nEdges = g.getNumContacts();
	this.deg    = new int[nNodes];
	this.fn     = "NoName";
	this.A      = new boolean[nNodes][nNodes];
	this.L      = new int[nNodes][];

	// fill the adjacency matrix and the node degree array
	for( Edge e : g.getContacts() ) {
	    if( e.i < 0 || e.j < 0 ) {
		throw new ContactMapConstructorError("Graph contains negative node indices. I don't like it!");
	    }
	    
	    A[e.i-1][e.j-1] = true;
	    A[e.j-1][e.i-1] = true;
	    
	    ++deg[e.i-1];
	    ++deg[e.j-1];
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
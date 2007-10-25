package sadp;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;
import java.util.logging.Logger;

import actionTools.Retriever;

import proteinstructure.Alignment;
import proteinstructure.Edge;
import proteinstructure.EdgeSet;
import proteinstructure.PairwiseAlignmentConverter;


/**
 * This class implements the softassign + dynamic programming algorithm.
 * 
 * Parameter setting:
 * 
 * The choice of parameters for softassign is problem dependent and crucial for
 * the time and quality performance. Parameters are set according to the Gold
 * and Rangarajan paper. Note that the annealing parameter is b = 1/T instead of
 * T, where T is the temperature (see paper).
 * 
 * @param b0
 *            initial value of annealing parameter b = 1/T.
 *            <p>
 *            Range : b0 > 0, vary within (0,2]
 *            <p>
 *            Default: b0 = 0.5
 * 
 * @param bf
 *            final value of annealing parameter b = 1/T.
 *            <p>
 *            Range : bf > b0, vary within (5, 20]
 *            <p>
 *            Default: b0 = 10
 * 
 * @param br
 *            factor by which the anealing parameter b = 1/T is increased.
 *            <p>
 *            Range : br > 1, vary within [1.075, 3.0]
 *            <p>
 *            Default: br = 1.075
 * @param I0
 *            number of iterations of assignment loop.
 *            <p>
 *            Range : I0 > 1, vary within {1, 2, ..., 10}
 *            <p>
 *            Default: I0 = 4
 * @param I1
 *            number of iterations of Sinkhorn loop.
 *            <p>
 *            Range : I1 > 1, vary within {1, 2, ..., 30}
 *            <p>
 *            Default: I1 = 30
 */
public class SADP {

    // control parameters of continuation method
    protected double b0 = 0.5;

    protected double bf = 10.0;

    protected double br = 1.075;

    // max # of iterations
    private int I0 = 4;

    private int I1 = 30;

    // precision
    private static double eps0 = 0.5;

    private static double eps1 = 0.05;

    // continuation parameter
    private double b = b0;

    /**
     * The reference protein represented as a graph to be matched.
     */
    protected ContactMap X;

    /**
     * The second protein represented as a graph to be matched.
     */
    protected ContactMap Y;

    /**
     * True if the input order of contact maps is preserved in the internal use
     * of the contact maps. This order might change as X always refers to the
     * smaller map and Y to the larger one. This variable has to be set
     * properly whenever a new contact map is set. 
     */
    protected boolean preservedInputOrder = true;
    
    /**
     * Holds sequences to the contact maps 
     */
    protected String[] sequences = new String[2];
        
    /**
     * The number of nodes of graph X.
     */
    protected int nNodes1;

    /**
     * The number of nodes of graph Y.
     */
    protected int nNodes2;

    /**
     * The maximum number of nodes
     */
    protected int maxNodes;

    /**
     * The mininmum number of nodes
     */
    protected int minNodes;

    /**
     * The maximum number of edges
     */
    protected int maxEdges;

    /**
     * The match matrix. M[i][j] is the match variable of node i from X and node
     * j from Y.
     */
    protected double[][] M;

    /**
     * The similarity score
     */
    protected double score;

    /**
     * The number of shared contacts
     */
    protected int ncc;

    /**
     * The state of the solution
     */
    protected boolean isFeasible;

    /**
     * The current number of iterations.
     */
    protected int nIterations;

    /**
     * Clock time needed to perform match
     */
    protected double time;
    
    /**
     * the logger 
     */
    protected Logger logger = null;
    
    /**
     * retrieves progress information
     * */
    protected Retriever retriever = null;
    
    /**
     * Default constructor. 
     * Do not use this constructor unless your are only interested in
     * retrieving the default values of this class.  
     */
    public SADP() {
	
    }
    
    /**
     * Constructor
     */
    public SADP(String fx, String fy) {

	ContactMap x = IOUtil.read(fx);
	ContactMap y = IOUtil.read(fy);

	if (x.countNodes() < y.countNodes()) {
	    this.X = x;
	    this.Y = y;
	    preservedInputOrder = true;
	} else {
	    this.X = y;
	    this.Y = x;
	    preservedInputOrder = false;
	}

	this.nNodes1 = X.countNodes();
	this.nNodes2 = Y.countNodes();
	this.maxNodes = Math.max(nNodes1, nNodes2);
	this.minNodes = Math.min(nNodes1, nNodes2);
	this.maxEdges = Math.max(X.countEdges(), Y.countEdges());
    }

    /**
     * Constructs an object instance based on ContactMaps.
     * Please note that the order of contact maps as passed to this constructor is not necessarily preserved for the internal use. 
     */
    public SADP(ContactMap x, ContactMap y) {
	if (x.countNodes() < y.countNodes()) {
	    this.X = x;
	    this.Y = y;
	    preservedInputOrder = true;
	} else {
	    this.X = y;
	    this.Y = x;
	    preservedInputOrder = false;
	}

	this.nNodes1 = X.countNodes();
	this.nNodes2 = Y.countNodes();
	this.maxNodes = Math.max(nNodes1, nNodes2);
	this.minNodes = Math.min(nNodes1, nNodes2);
	this.maxEdges = Math.max(X.countEdges(), Y.countEdges());
    }

    public Double getEps0() {
	return eps0;
    }

    public void setEps0( Double eps0 ) {
	SADP.eps0 = eps0;
    }

    public Double getEps1() {
	return eps1;
    }

    public void setEps1( Double eps1 ) {
	SADP.eps1 = eps1;
    }
    
    public Double getB0() {
	return b0;
    }

    public void setB0( Double b0 ) {
	this.b0 = b0;
    }
    
    public Double getBf() {
	return bf;
    }

    public void setBf( Double  bf ) {
	this.bf = bf;
    }
    
    public Double getBr() {
	return br;
    }

    public void setBr( Double  br ) {
	this.br = br;
    }

    public Integer getI0() {
	return I0;
    }

    public void setI0( Integer I0 ) {
	this.I0 = I0;
    }

    public Integer getI1() {
	return I1;
    }

    public void getI1( Integer  I1 ) {
	this.I1 = I1;
    }
    
    /**
     * Compares contact maps X and Y.
     * 
     * Use method runL() (= runLarge()) if runS() results in an out of memory
     * exception. In contrast to runS(), this method recomputes the
     * compatibility coefficients at each iteration. Therefore, runL() is slower
     * than runS().
     */
    public void run() {

	int[][] AX = X.getAdjacencyList();
	int[][] AY = Y.getAdjacencyList();
	int[] degX = X.deg;
	int[] degY = Y.deg;

	// set time
	this.time = System.currentTimeMillis();

	// initialize variables
	double[][] Q = new double[nNodes1][nNodes2];
	double[][] M0 = new double[nNodes1 + 1][nNodes2 + 1];
	double[][] M1 = new double[nNodes1 + 1][nNodes2 + 1];

	// initialize match matrix
	this.M = new double[nNodes1 + 1][nNodes2 + 1];
	for (int i = 0; i < nNodes1 + 1; i++) {
	    for (int j = 0; j < nNodes2 + 1; j++) {
		M[i][j] = 0.1;
	    }
	}

	// scale parameters
	double r = maxNodes / (double) nNodes1;

	// A loop iteration counter
	int countAloopIterations = 0;
	int maxAloopIterations   = 0;
	
	// initialize counter if logging is enabled
	if( retriever != null ) {
	    maxAloopIterations = getAloopIterations();
	}
		
	// A loop
	this.nIterations = 0;
	while (b < bf) {

	    // B loop
	    for (int t0 = 0; t0 < I0; t0++) {

		this.nIterations++;

		// copy
		for (int i = 0; i < nNodes1 + 1; i++) {
		    System.arraycopy(M[i], 0, M0[i], 0, nNodes2 + 1);
		}

		// softmax
		for (int i = 0; i < nNodes1; i++) {
		    for (int j = 0; j < nNodes2; j++) {
			Q[i][j] = 0;
			for (int k = 0; k < degX[i]; k++) {
			    int d1 = Math.abs(i - AX[i][k]);
			    for (int l = 0; l < degY[j]; l++) {
				if (i > AX[i][k] && j > AY[j][l]
				                              || i < AX[i][k] && j < AY[j][l]) {
				    int d2 = Math.abs(j - AY[j][l]);
				    double w = 1.0 / (1.0 + 0.1 * Math.abs(r
					    * d1 - d2)); // 0.1
				    Q[i][j] += w * M0[AX[i][k]][AY[j][l]];
				}
			    }
			}
			M[i][j] = Math.exp(b * Q[i][j]);
		    }
		}

		// C loop
		for (int t1 = 0; t1 < I1; t1++) {

		    // copy
		    for (int i = 0; i < nNodes1 + 1; i++) {
			System.arraycopy(M[i], 0, M1[i], 0, nNodes2 + 1);
		    }

		    /** * normalize across all rows ** */
		    for (int i = 0; i < nNodes1 + 1; i++) {
			double row_sum = 0.0;
			for (int j = 0; j < nNodes2 + 1; j++) {
			    row_sum += M[i][j];
			}
			for (int j = 0; j < nNodes2 + 1; j++) {
			    M[i][j] /= row_sum;
			}
		    }

		    /** * normalize across all columns ** */
		    double err1 = 0;
		    for (int j = 0; j < nNodes2 + 1; j++) {
			double col_sum = 0.0;
			for (int i = 0; i < nNodes1 + 1; i++) {
			    col_sum += M[i][j];
			}
			for (int i = 0; i < nNodes1 + 1; i++) {
			    M[i][j] /= col_sum;
			    err1 += Math.abs(M[i][j] - M1[i][j]);
			}
		    }

		    /** * check convergence ** */
		    if (err1 < eps1) {
			break;
		    }

		} // end C loop

		/** * check convergence ** */
		double err0 = 0;
		for (int i = 0; i < nNodes1; i++) {
		    for (int j = 0; j < nNodes2; j++) {
			err0 += Math.abs(M[i][j] - M0[i][j]);
		    }
		}
		if (err0 < eps0) {
		    break;
		}

	    } // end B loop

	    b *= br;
	    
	    // log progress if a logger is available
	    if( this.retriever != null ) {
		++countAloopIterations;
		try {
		    retriever.retrieve(100.0*((float) countAloopIterations/(float) maxAloopIterations));
		} catch (ClassCastException e) {
		    System.out.println(e.getMessage());
		}
//		logger.info("SADP-progress:"+100.0*((float) countAloopIterations/(float) maxAloopIterations));
	    }
	    
	} // end A loop

	cleanup();
	noncrossing();

	time = System.currentTimeMillis() - time;

	setScore();
    }
    
    private boolean isInputOrderPreserved() {
	return preservedInputOrder;
    }

    private void cleanup() {

	double[][] M2 = new double[nNodes1][nNodes2];

	boolean[] isSet = new boolean[nNodes2];
	for (int i = 0; i < nNodes1; i++) {
	    double max_val = -1.0;
	    int index = -1;
	    for (int j = 0; j < nNodes2; j++) {
		if (max_val < M[i][j] && !isSet[j]) {
		    index = j;
		    max_val = M[i][j];
		}
	    }
	    M2[i][index] = 1;
	    isSet[index] = true;
	}
	M = M2;
    }

    private void noncrossing() {

	double[][] S = new double[nNodes1][nNodes2];
	double[] sOpt = new double[nNodes2];
	for (int i = 0; i < nNodes1; i++) {
	    S[i][0] = M[i][0];
	    double max = 0;
	    for (int j = 1; j < nNodes2; j++) {
		max = Math.max(max, sOpt[j - 1]);
		S[i][j] = M[i][j] + max;
	    }
	    for (int j = 0; j < nNodes2; j++) {
		sOpt[j] = Math.max(sOpt[j], S[i][j]);
	    }
	}

	double[][] Q = new double[nNodes1][nNodes2];
	int k = nNodes2;
	for (int i = nNodes1 - 1; i > -1 && k > 0; i--) {
	    double max = -1;
	    int jOpt = k - 1;
	    for (int j = 0; j < k; j++) {
		if (max <= S[i][j]) {
		    max = S[i][j];
		    jOpt = j;
		}
	    }
	    Q[i][jOpt] = 1;
	    k = jOpt;
	}
	M = Q;
    }

    public void setScore() {

	isFeasible = true;
	score = 0;
	ncc = 0;

	if (M == null) {
	    return;
	}

	int[][] AX = X.getAdjacencyList();
	int[][] AY = Y.getAdjacencyList();

	for (int i = 0; i < nNodes1; i++) {
	    for (int j = 0; j < nNodes2; j++) {
		if (M[i][j] <= 0.0) {
		    continue;
		}
		for (int k = 0; k < AX[i].length; k++) {
		    for (int l = 0; l < AY[j].length; l++) {
			if (M[AX[i][k]][AY[j][l]] <= 0.0) {
			    continue;
			}
			if (((AX[i][k] < i) && (AY[j][l] < j))
				|| ((AX[i][k] > i) && (AY[j][l] > j))) {
			    score += 1.0;
			    ncc++;
			    continue;
			}
			score = -1.0;
			ncc = -1;
			isFeasible = false;
			return;
		    }
		}
	    }
	}
	score /= (2.0 * Math.min(X.countEdges(), Y.countEdges()));
	ncc /= 2;
    }

    /**
     * Returns the match matrix M.
     * 
     * @return match matrix
     */
    public double[][] getMatchMatrix() {

	return M;
    }

    /**
     * Compute the maximal iterations in the A loop.
     * The A loop denote the outer loop in function run(). This function is employed 
     * to estimate the running time of function run.
     */
    private int getAloopIterations() {
	
	// at the end we do a change of the logarithm with base 'br' to the logarithm with base 10
	// stop-criteria: b*br^x > bf
	// <=> x = \log_{br}(bf/b) = \frac{\log_{10}(bf/b)}{\log_{10}(br)}
	// x -> maxAloopIterations
	return (int) Math.floor(Math.log10(bf/b)/Math.log10(br)) + 1; 
    }
    
    /**
     * Returns similarity score of X and Y given the specified match matrix. The
     * score is given by lb/minEdges, where lb denotes the number of common
     * contacts found by softassign and minEdges is the minimum number of edges
     * of both contact maps X and Y.
     * 
     * @return the score of similarity between X and Y given the specified m.
     */
    public double getScore() {

	return Math.round(100.0 * score) / 100.0;
    }

    /**
     * Returns the number of common contacts of X and Y found by softassign.
     * 
     * @return number of common contacts.
     */
    public int getNumberOfCommonContacts() {

	return ncc;
    }

    /**
     * TODO: What is behind the term "solution is feasible"? 
     * suggestions: matching contains only noncrossing edges (which is guaranteed by the function noncrossing()), whatever ...
     * @return true if the matching is feasible
     */
    public boolean isFeasible() {

	return isFeasible;
    }

    /**
     * Returns the number of iterations.
     * 
     * @return Returns the number of Iterations.
     */
    public int getIterations() {

	return nIterations;
    }
    
    /**
     * @return computation time in milliseconds
     */
    public double getTime() {

	return time;
    }
    
    /**
     * This sequence refers to the first contact map assigned in any construtor.
     */
    public void setFirstSequence( String s ) {
	
	sequences[0] = s;
    }
    
    /**
     * This sequence refers to the second contact map assigned in any construtor.
     */
    public void setSecondSequence( String s ) {
	
	sequences[1] = s;
    }
        
    /**
     * Creates the alignment and the resulting protein structures graphs which
     * are based on the alignment and include nodes representing gaps. If there
     * is no sequence provided for either of the contact maps a sequence contait is used in the sequ
     * 
     * @param tag1  name of the first contact map to be taken as the sequence name tag resulting Alignment object
     * @param tag2  name of the second contact map 
     */
    public Alignment getAlignment( String tag1, String tag2 ) {
	
	ContactMap[] cm = new ContactMap[2];
	
	if( isInputOrderPreserved() ) {
	    cm[0] = X;
	    cm[1] = Y;	    
	} else {
	    cm[0] = X;
	    cm[1] = Y;
	}
	
	EdgeSet matching = getMatching();
	
	PairwiseAlignmentConverter pac = null;
	
	if( sequences[0] != null ) {
	    if( sequences[1] != null ) {
		pac = new PairwiseAlignmentConverter(matching.iterator(), sequences[0],       sequences[1],       tag1,tag2,0);
	    } else {
		pac = new PairwiseAlignmentConverter(matching.iterator(), sequences[0],       cm[1].countNodes(), tag1,tag2,0);
	    }
	} else {
	    if( sequences[1] != null ) {
		pac = new PairwiseAlignmentConverter(matching.iterator(), cm[0].countNodes(), sequences[1],       tag1,tag2,0);
	    } else {
		pac = new PairwiseAlignmentConverter(matching.iterator(), cm[0].countNodes(), cm[1].countNodes(), tag1,tag2,0);		
	    }
	}
	
	return pac.getAlignment();
    }

    /**
     * Retrieves the matching of nodes from the first contact map to the second
     * one in a map. The node indexing starts with 0 (not with 1 as in Graph).
     * @return edge set containing all noncrossing matching edges
     * @see #getMatchingAsEdgeList()
     */
    public EdgeSet getMatching() {

	EdgeSet matching = new EdgeSet();
	
	for (int i = 0; i < nNodes1; i++) {
	    for (int j = 0; j < nNodes2; j++) { // TODO: as the edge are non-crossing we do not need to start from zero ...
		if (M[i][j] > 0) {
		    if ( isInputOrderPreserved() ) {
			matching.add(new Edge(i,j));			
		    } else {
			matching.add(new Edge(j, i));
		    }
		}
	    }
	}
	
	return matching;
    }
    
    /**
     * Retrieves the matching of nodes from the first contact map to the second
     * one in a linked list.
     * @return linked list containg all noncrossing matching edges
     * @see #getMatching()  
     */
    public LinkedList<Edge> getMatchingAsEdgeList() {
	
	LinkedList<Edge> matching = new LinkedList<Edge>();
		
	for (int i = 0; i < nNodes1; i++) {
	    for (int j = 0; j < nNodes2; j++) {
		if (M[i][j] > 0) {
		    if ( isInputOrderPreserved() ) {
			matching.add(new Edge(i,j));			
		    } else {
			matching.add(new Edge(j, i));
		    }
		}
	    }
	}
	
	return matching;
    }
    
    /**
     * Sets logger.
     * Make use of this logging mechanism to send progress status informations to another thread.
     * @param logger  a logger
     * @see getLogger()
     * @deprecated
     */
    public void setLogger(Logger logger) {
	
	this.logger = logger;
    }
    
    /**
     * Gets the progress logger.
     * @return a logger if one has been set, else null.
     * @deprecated  
     */
    public Logger getLogger() {
	
	return this.logger;
    }
    
    /**
     * Sets the progress info retriever.
     * @param retr  a progress info retriever
     * */
    public void setProgressInfoRetriever(Retriever retr) {
	retriever = retr;
    }
    
    /**
     * Gets the progress info retriever.
     * @return the progress info retriever
     * */
    public Retriever getProgressInfoRetriever() {
	return retriever;
    }
    
    /**
     * Writes node mapping to standard output.
     */
    public void showMapping() {

	System.out.println("Match:");
	for (int i = 0; i < nNodes1; i++) {
	    for (int j = 0; j < nNodes2; j++) {
		if (M[i][j] <= 0) {
		    continue;
		}
		System.out.println("    " + i + " -> " + j);
	    }
	}
	System.out.println();
    }
    

    /**
     * Shows the result.
     */
    public void show() {
	X.show();
	Y.show();
	System.out.println();
	System.out.println("Result:");
	System.out.println("    #(shared contacts) : "
		+ this.getNumberOfCommonContacts());
	System.out.println("    Score              : " + this.getScore());
	System.out.println("    Time               : " + this.getTime()
		+ " msec");
	System.out.println();
	showMapping();
    }

    public static void main(String[] args) {

	Logger logger = Logger.getLogger("sasp.SADP");
	
	SADP sadp = new SADP(args[0], args[1]);
	sadp.setLogger(logger);
	sadp.run();
	sadp.show();
	
	// prints pseudo sequence alignment
	Alignment ali = sadp.getAlignment(args[0], args[1]);
	Set<String> tags = ali.getTags();
	Iterator<String> it = tags.iterator();
	while( it.hasNext() ) {
	    String tag = it.next();
	    System.out.println(tag+':');
	    System.out.println(ali.getAlignedSequence(tag));
	}
	
	System.out.println(ali.seq2al(args[0], 1));
	System.out.println(ali.seq2al(args[1], 1));
    }
}

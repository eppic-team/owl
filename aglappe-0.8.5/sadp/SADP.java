package sadp;
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
     * Constructor
     */
    public SADP(String fx, String fy) {

	ContactMap x = IOUtil.read(fx);
	ContactMap y = IOUtil.read(fy);

	if (x.countNodes() < y.countNodes()) {
	    this.X = x;
	    this.Y = y;
	} else {
	    this.X = y;
	    this.Y = x;
	}

	this.nNodes1 = X.countNodes();
	this.nNodes2 = Y.countNodes();
	this.maxNodes = Math.max(nNodes1, nNodes2);
	this.minNodes = Math.min(nNodes1, nNodes2);
	this.maxEdges = Math.max(X.countEdges(), Y.countEdges());
    }

    public SADP(ContactMap x, ContactMap y) {
	if (x.countNodes() < y.countNodes()) {
	    this.X = x;
	    this.Y = y;
	} else {
	    this.X = y;
	    this.Y = x;
	}

	this.nNodes1 = X.countNodes();
	this.nNodes2 = Y.countNodes();
	this.maxNodes = Math.max(nNodes1, nNodes2);
	this.minNodes = Math.min(nNodes1, nNodes2);
	this.maxEdges = Math.max(X.countEdges(), Y.countEdges());
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
	} // end A loop

	cleanup();
	noncrossing();

	time = System.currentTimeMillis() - time;

	setScore();
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
     * Returns the computation time.
     */
    public double getTime() {

	return time;
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

	SADP sadp = new SADP(args[0], args[1]);
	sadp.run();
	sadp.show();
    }
}

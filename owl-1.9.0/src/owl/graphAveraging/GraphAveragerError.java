package owl.graphAveraging;

/**
 * Exception to be thrown when GraphAverager construction fails becaus of inconsistent input.
 * @author duarte
 *
 */
public class GraphAveragerError extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public GraphAveragerError() {
	}

	public GraphAveragerError(String arg0) {
		super(arg0);
	}

	public GraphAveragerError(Throwable arg0) {
		super(arg0);
	}

	public GraphAveragerError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

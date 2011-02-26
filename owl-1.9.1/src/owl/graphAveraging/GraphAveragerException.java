package owl.graphAveraging;

/**
 * Exception to be thrown when GraphAverager construction fails becaus of inconsistent input.
 * @author duarte
 *
 */
public class GraphAveragerException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public GraphAveragerException() {
	}

	public GraphAveragerException(String arg0) {
		super(arg0);
	}

	public GraphAveragerException(Throwable arg0) {
		super(arg0);
	}

	public GraphAveragerException(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

package owl.core.runners.gromacs;

public class GromacsError extends Exception {

	/**
	 * To be thrown when a gromacs program exits abnormally 
	 */
	private static final long serialVersionUID = 1L;

	public GromacsError() {
	}

	public GromacsError(String arg0) {
		super(arg0);
	}

	public GromacsError(Throwable arg0) {
		super(arg0);
	}

	public GromacsError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

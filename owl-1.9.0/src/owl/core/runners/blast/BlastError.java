package owl.core.runners.blast;

public class BlastError extends Exception {

	/**
	 * To be thrown when a blast program exits abnormally 
	 */
	private static final long serialVersionUID = 1L;

	public BlastError() {
	}

	public BlastError(String arg0) {
		super(arg0);
	}

	public BlastError(Throwable arg0) {
		super(arg0);
	}

	public BlastError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

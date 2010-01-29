package sequence;

public class PsipredError extends Exception {

	/**
	 * To be thrown when a psipred program exits abnormally 
	 */
	private static final long serialVersionUID = 1L;

	public PsipredError() {
	}

	public PsipredError(String arg0) {
		super(arg0);
	}

	public PsipredError(Throwable arg0) {
		super(arg0);
	}

	public PsipredError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

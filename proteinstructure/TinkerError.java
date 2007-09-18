package proteinstructure;

public class TinkerError extends Exception {

	/**
	 * To be thrown when a tinker program exits abnormally 
	 */
	private static final long serialVersionUID = 1L;

	public TinkerError() {
	}

	public TinkerError(String arg0) {
		super(arg0);
	}

	public TinkerError(Throwable arg0) {
		super(arg0);
	}

	public TinkerError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

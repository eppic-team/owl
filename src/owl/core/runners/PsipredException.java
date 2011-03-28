package owl.core.runners;

public class PsipredException extends Exception {

	/**
	 * To be thrown when a psipred program exits abnormally 
	 */
	private static final long serialVersionUID = 1L;

	public PsipredException() {
	}

	public PsipredException(String arg0) {
		super(arg0);
	}

	public PsipredException(Throwable arg0) {
		super(arg0);
	}

	public PsipredException(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

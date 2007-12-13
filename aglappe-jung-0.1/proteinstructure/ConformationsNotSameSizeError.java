package proteinstructure;

public class ConformationsNotSameSizeError extends Exception {

	/**
	 * Exception to be thrown in rmsd calculation when given conformations are not of the same size
	 */
	private static final long serialVersionUID = 1L;

	public ConformationsNotSameSizeError() {
	}

	public ConformationsNotSameSizeError(String arg0) {
		super(arg0);
	}

	public ConformationsNotSameSizeError(Throwable arg0) {
		super(arg0);
	}

	public ConformationsNotSameSizeError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

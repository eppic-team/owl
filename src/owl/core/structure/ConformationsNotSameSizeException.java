package owl.core.structure;

public class ConformationsNotSameSizeException extends Exception {

	/**
	 * Exception to be thrown in rmsd calculation when given conformations are not of the same size
	 */
	private static final long serialVersionUID = 1L;

	public ConformationsNotSameSizeException() {
	}

	public ConformationsNotSameSizeException(String arg0) {
		super(arg0);
	}

	public ConformationsNotSameSizeException(Throwable arg0) {
		super(arg0);
	}

	public ConformationsNotSameSizeException(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

package proteinstructure;

public class CiffileFormatError extends Exception {

	/**
	 * Exception to be used when parsing a cif file and some formatting error is found
	 */
	private static final long serialVersionUID = 1L;

	public CiffileFormatError() {
	}

	public CiffileFormatError(String arg0) {
		super(arg0);
	}

	public CiffileFormatError(Throwable arg0) {
		super(arg0);
	}

	public CiffileFormatError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}
}

package proteinstructure;

public class PdbfileFormatError extends Exception {

	/**
	 * Exception to be used when a pdb file is found not to be in the right format
	 */
	private static final long serialVersionUID = 1L;

	public PdbfileFormatError() {
	}

	public PdbfileFormatError(String arg0) {
		super(arg0);
	}

	public PdbfileFormatError(Throwable arg0) {
		super(arg0);
	}

	public PdbfileFormatError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}
}

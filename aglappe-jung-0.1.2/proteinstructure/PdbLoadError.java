package proteinstructure;

/**
 * Exception to be thrown when loading pdb structures and catching IO/SQL/format exceptions 
 */
public class PdbLoadError extends Exception {

	
	private static final long serialVersionUID = 1L;

	public PdbLoadError() {
	}

	public PdbLoadError(String arg0) {
		super(arg0);
	}

	public PdbLoadError(Throwable arg0) {
		super(arg0);
	}

	public PdbLoadError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

package owl.core.structure;

/**
 * Exception to be thrown when loading pdb structures and catching IO/SQL/format exceptions 
 */
public class PdbLoadException extends Exception {

	
	private static final long serialVersionUID = 1L;

	public PdbLoadException() {
	}

	public PdbLoadException(String arg0) {
		super(arg0);
	}

	public PdbLoadException(Throwable arg0) {
		super(arg0);
	}

	public PdbLoadException(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

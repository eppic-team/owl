package proteinstructure;

public class PdbChainCodeNotFoundError extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public PdbChainCodeNotFoundError() {
	}

	public PdbChainCodeNotFoundError(String arg0) {
		super(arg0);
	}

	public PdbChainCodeNotFoundError(Throwable arg0) {
		super(arg0);
	}

	public PdbChainCodeNotFoundError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

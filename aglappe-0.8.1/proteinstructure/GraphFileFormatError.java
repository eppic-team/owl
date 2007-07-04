package proteinstructure;

public class GraphFileFormatError extends Exception {

	/**
	 * Exception to be used when reading a graph file and it is not in the right format
	 */
	private static final long serialVersionUID = 1L;

	public GraphFileFormatError() {
	}

	public GraphFileFormatError(String arg0) {
		super(arg0);
	}

	public GraphFileFormatError(Throwable arg0) {
		super(arg0);
	}

	public GraphFileFormatError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}
}

package proteinstructure;

public class GraphFileFormatError extends FileFormatError {

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
	
	public GraphFileFormatError(String arg0, String file, long line ) {
	    super(arg0,file,line);
	}

	public GraphFileFormatError(String arg0, Throwable arg1, String file, long line ) {
	    super(arg0,arg1,file,line);
	}

	public GraphFileFormatError(String arg0, String file, long startByte, long stopByte ) {
	    super(arg0,file,startByte,stopByte);
	}

	public GraphFileFormatError(String arg0, Throwable arg1, String file, long startByte, long stopByte ) {
	    super(arg0,arg1,file,startByte,stopByte);
	}
}

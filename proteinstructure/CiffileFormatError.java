package proteinstructure;

public class CiffileFormatError extends FileFormatError {

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

	public CiffileFormatError(String arg0, String file, long line ) {
	    super(arg0,file,line);
	}

	public CiffileFormatError(String arg0, Throwable arg1, String file, long line ) {
	    super(arg0,arg1,file,line);
	}

	public CiffileFormatError(String arg0, String file, long startByte, long stopByte ) {
	    super(arg0,file,startByte,stopByte);
	}

	public CiffileFormatError(String arg0, Throwable arg1, String file, long startByte, long stopByte ) {
	    super(arg0,arg1,file,startByte,stopByte);
	}
}

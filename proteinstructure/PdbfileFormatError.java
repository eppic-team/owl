package proteinstructure;

public class PdbfileFormatError extends FileFormatError {

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
	
	public PdbfileFormatError(String arg0, String file, long line ) {
	    super(arg0,file,line);
	}

	public PdbfileFormatError(String arg0, Throwable arg1, String file, long line ) {
	    super(arg0,arg1,file,line);
	}

	public PdbfileFormatError(String arg0, String file, long startByte, long stopByte ) {
	    super(arg0,file,startByte,stopByte);
	}

	public PdbfileFormatError(String arg0, Throwable arg1, String file, long startByte, long stopByte ) {
	    super(arg0,arg1,file,startByte,stopByte);
	}
}

package proteinstructure;

public class FileFormatError extends Exception {
    static final long serialVersionUID = 1L;

    public FileFormatError() {
    }

    public FileFormatError(String arg0) {
	super(arg0);
    }

    public FileFormatError(Throwable arg0) {
	super(arg0);
    }

    public FileFormatError(String arg0, Throwable arg1) {
	super(arg0, arg1);
    }
    
    public FileFormatError(String arg0, String file, long line ) {
	super(file+"["+line+"]: "+arg0);
    }
    
    public FileFormatError(String arg0, Throwable arg1, String file, long line ) {
	super(file+"["+line+"]: "+arg0, arg1);
    }
    
    public FileFormatError(String arg0, String file, long startByte, long stopByte ) {
	super(file+"["+startByte+"-"+stopByte+"]: "+arg0);
    }
    
    public FileFormatError(String arg0, Throwable arg1, String file, long startByte, long stopByte ) {
	super(file+"["+startByte+"-"+stopByte+"]: "+arg0, arg1);
    }
}
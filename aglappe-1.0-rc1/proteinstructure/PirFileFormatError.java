package proteinstructure;

/**
 * Exception to be thrown if Pir file seems to be corrupted or the 
 * considered file does not at all conform to the Pir file format, 
 * respectively. 
 *
 * @author Lars Petzold
 */
public class PirFileFormatError extends FileFormatError {
    private static final long serialVersionUID = 1L;

    /**
     * The empty constructor.
     */
    public PirFileFormatError() {
    }

    /**
     * Constructs an exception with a message.
     * @param arg0  the message
     */
    public PirFileFormatError(String arg0) {
	super(arg0);
    }

    /**
     * Constructs an exception with another throwable
     * @param arg0  the other throwable
     */
    public PirFileFormatError(Throwable arg0) {
	super(arg0);
    }

    /**
     * Constructs an exception with a message and another throwable
     * @param arg0  the message
     * @param arg1  the other throwable
     */
    public PirFileFormatError(String arg0, Throwable arg1) {
	super(arg0, arg1);
    }

    /**
     * Constructs an exception with a message, the name of the corrupted file 
     * and the line in the file where the error has been obtained.
     * @param arg0  a message describing the detected error
     * @param file  the name of the corrupted file
     * @param line  the corrupted line
     */
    public PirFileFormatError(String arg0, String file, long line ) {
	super(arg0,file,line);
    }

    /**
     * Constructs an exception with a message, another throwable, the name of 
     * the corrupted file and the line in the file where the error has been 
     * obtained.
     * @param arg0  a message describing the detected error
     * @param arg1  another throwable
     * @param file  the name of the corrupted file
     * @param line  the corrupted line
     */
    public PirFileFormatError(String arg0, Throwable arg1, String file, long line ) {
	super(arg0,arg1,file,line);
    }

    /**
     * Constructs an exception with a message, the name of the corrupted file 
     * and an interval indicating the erroneous region in the file. 
     * @param arg0  a message describing the detected error
     * @param file  the name of the corrupted file
     * @param startByte  first byte (inclusive) of the erroneous region
     * @param stopByte   last byte (exclusive) of the erroneous region
     */
    public PirFileFormatError(String arg0, String file, long startByte, long stopByte ) {
	super(arg0,file,startByte,stopByte);
    }

    /**
     * Constructs an exception with a message, the name of the corrupted file 
     * and an interval indicating the erroneous region in the file. 
     * @param arg0  a message describing the detected error
     * @param arg1  another throwable
     * @param file  the name of the corrupted file
     * @param startByte  first byte (inclusive) of the erroneous region
     * @param stopByte   last byte (exclusive) of the erroneous region
     */
    public PirFileFormatError(String arg0, Throwable arg1, String file, long startByte, long stopByte ) {
	super(arg0,arg1,file,startByte,stopByte);
    }
}
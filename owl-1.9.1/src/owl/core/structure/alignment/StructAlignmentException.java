package owl.core.structure.alignment;

/**
 * To be thrown when a Structural Alignment error occurs 
 */
public class StructAlignmentException extends Exception {

	private static final long serialVersionUID = 1L;

	public StructAlignmentException() {
	}

	public StructAlignmentException(String arg0) {
		super(arg0);
	}

	public StructAlignmentException(Throwable arg0) {
		super(arg0);
	}

	public StructAlignmentException(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

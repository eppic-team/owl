package owl.core.structure;

/**
 * To be thrown when a Structural Alignment error occurs 
 */
public class StructAlignmentError extends Exception {

	private static final long serialVersionUID = 1L;

	public StructAlignmentError() {
	}

	public StructAlignmentError(String arg0) {
		super(arg0);
	}

	public StructAlignmentError(Throwable arg0) {
		super(arg0);
	}

	public StructAlignmentError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

package owl.embed;
/**
 * A class indicating that some exceptional situation occurred during
 * computation of instances of the <code>{@link SparseMatrix}</code> class.
 * Typically, addition or multiplication of instances with non-matching dimensions
 * will cause such an exception, as well as <tt>double</tt> matrix instances, with
 * various row dimension or index pairs exceeding the dimensions.
 * @author gmueller
 *
 */
public class SparseMatrixException extends RuntimeException {
	private static final long serialVersionUID = 1L;
	/**
	 * Constructs a new sparse matrix exception indicating a serious error
	 * or exceptional situation occurred during computation, with a specified detail
	 * message.
	 * @param message the specified detail message
	 */
	public SparseMatrixException (String message){
		super(message);
	}
}

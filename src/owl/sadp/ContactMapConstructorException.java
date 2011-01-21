package owl.sadp;

/**
 * Exception to be used when constructing a contact map from a graph object.
 */
public class ContactMapConstructorException extends Exception {

    private static final long serialVersionUID = 1L;

    public ContactMapConstructorException() {
    }

    public ContactMapConstructorException(String arg0) {
	super(arg0);
    }

    public ContactMapConstructorException(Throwable arg0) {
	super(arg0);
    }

    public ContactMapConstructorException(String arg0, Throwable arg1) {
	super(arg0, arg1);
    }
}

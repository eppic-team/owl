package sadp;

/**
 * Exception to be used when constructing a contact map from a graph object.
 */
public class ContactMapConstructorError extends Exception {

    private static final long serialVersionUID = 1L;

    public ContactMapConstructorError() {
    }

    public ContactMapConstructorError(String arg0) {
	super(arg0);
    }

    public ContactMapConstructorError(Throwable arg0) {
	super(arg0);
    }

    public ContactMapConstructorError(String arg0, Throwable arg1) {
	super(arg0, arg1);
    }
}

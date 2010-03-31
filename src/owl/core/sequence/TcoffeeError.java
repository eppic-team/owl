package owl.core.sequence;

public class TcoffeeError extends Exception {

	/**
	 * To be thrown when tcoffee program exits abnormally 
	 */
	private static final long serialVersionUID = 1L;

	public TcoffeeError() {
	}

	public TcoffeeError(String arg0) {
		super(arg0);
	}

	public TcoffeeError(Throwable arg0) {
		super(arg0);
	}

	public TcoffeeError(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

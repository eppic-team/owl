package owl.core.runners;

public class TcoffeeException extends Exception {

	/**
	 * To be thrown when tcoffee program exits abnormally 
	 */
	private static final long serialVersionUID = 1L;

	public TcoffeeException() {
	}

	public TcoffeeException(String arg0) {
		super(arg0);
	}

	public TcoffeeException(Throwable arg0) {
		super(arg0);
	}

	public TcoffeeException(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}

}

package owl.core.sequence;

/**
 * A exception to be thrown when a sequence or codon can not be translated
 * @author duarte_j
 *
 */
public class TranslationException extends Exception {
	static final long serialVersionUID = 1L;

	public TranslationException() {
	}

	public TranslationException(String arg0) {
		super(arg0);
	}

	public TranslationException(Throwable arg0) {
		super(arg0);
	}

	public TranslationException(String arg0, Throwable arg1) {
		super(arg0, arg1);
	}
}

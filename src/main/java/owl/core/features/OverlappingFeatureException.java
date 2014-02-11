package owl.core.features;

/**
 * An exception which is thrown if on overlapping feature is being added to an object which implements the HasFeatures interface.
 * @author stehr
 */
public class OverlappingFeatureException extends Exception {
	private static final long serialVersionUID = 1L;
	
	public OverlappingFeatureException(String msg) {
		super(msg);
	}
}

package proteinstructure.features;

/**
 * Exception which is thrown by the HasFeatures.addFeature() method when attempting to add a feature where
 * the feature coordinates are out of bounds of the reference sequence (e.g. of the implementing Pdb object).
 * @author stehr
 */
public class InvalidFeatureCoordinatesException extends Exception {
	private static final long serialVersionUID = 1L;
	
	public InvalidFeatureCoordinatesException(String msg) {
		super(msg);
	}
}

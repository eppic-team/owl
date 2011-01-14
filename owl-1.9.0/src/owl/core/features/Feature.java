package owl.core.features;

import owl.core.util.IntervalSet;

/**
 * A feature in a sequence or structure.
 * Each subclass implementing this interface should have a corresponding instance of enum FeatureType.
 * Examples: Catalytic site, Scop domain.
 * 
 * This class represents local features which have a specified position in the reference sequence (which
 * implements the HasFeatures interface). The simplest type of features have just a beginning and end index
 * and corresponds to contiguous subsequences. A feature's position can also be a set of intervals
 * in which case it would be a non-contiguous feature. This interface allows different kinds of features
 * tp be visualized in a unified way. In a sequence this could be in a genome-browser like fashion. If the
 * reference sequence is associated with 3D coordinates the features could also be visualized in the structure
 * (e.g. with different coloured residues).
 * 
 * See also: @{@link FeatureType}, {@link HasFeatures}
 * 
 * @author stehr
 */
public interface Feature {
	
	/**
	 * Returns the position of this feature in the reference sequence.
	 * @return the position of this feature in the reference sequence
	 */
	public IntervalSet getIntervalSet();
	
	/**
	 * Returns the type of this feature. The type should be set by the
	 * constructor and needs to be defined in enum FeatureType. This is
	 * used to query a sequence for only a specific feature type using
	 * HasFeatures.getFeaturesOfType and HasFeatures.getFeaturesOfTypeForPosition.
	 * It could also be used to visualize different types of features in
	 * different ways (e.g. colors).
	 * @return the FeatureType
	 */
	public FeatureType getType();
	
	/**
	 * Returns a description for this particular instance of a feature. In a genome
	 * browser setting, this should be the label of the particular feature, not the
	 * stream (which corresponds to FeatureType).
	 * @return
	 */
	public String getDescription();
	
}

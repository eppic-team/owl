package owl.core.features;

import java.util.Collection;

/**
 * Classes implementing this interface provide methods to attach and retrieve sequence/structure features
 * complying to the {@link Feature} interface.
 * 
 * 2010-09-15: Changed the documented behaviour to always returning empty collections if no feature was found.
 * This still has to be enforced by the implementing classes (which should be currently only mutanom.gene).
 * 
 * @author stehr
 */
public interface HasFeatures {

		// minimal interface
	
		/**
		 * Returns all features attached to the implementing object (disregarding the FeatureType).
		 * @return a possibly empty collection of all features attached to the implementing object
		 */
		public Collection<Feature> getFeatures();
	
		/**
		 * Attach the given feature to the implementing object.
		 * @param feature
		 * @return true, if feature not present so it was added, false if feature already present
		 * @throws InvalidFeatureCoordinatesException if feature coordinates can not be mapped to the reference sequence (not yet implemented)
		 * @throws OverlappingFeatureException if the feature overlaps with an existing feature of the same type (not yet implemented)
		 */
		public boolean addFeature(Feature feature) throws InvalidFeatureCoordinatesException, OverlappingFeatureException;		
		
		// convenience methods; these could be implemented by an abstract class based on above methods
		
		/**
		 * Returns the FeatureTypes for which there are features attached to the implementing object.
		 * @return a possibly empty collection of FeatureTypes for which the implementing object has features
		 */
		public Collection<FeatureType> getFeatureTypes();
		
		/**
		 * Returns all features of a specific type which are attached to the implementing object.
		 * @param featureType the type of features which are to be retrieved
		 * @return a possibly empty collection of features if the given type
		 */
		public Collection<Feature> getFeaturesOfType(FeatureType featureType);
		
		/**
		 * Returns the features which annotate the given position.
		 * @param position the position in reference coordinates
		 * @return a possibly empty collection of features which annotate the given position
		 */		
		public Collection<Feature> getFeaturesForPositon(int position);
		
		/**
		 * Returns the features of the given type which annotate the given position.
		 * In most cases, this method should return at most one object, unless overlapping object of the same type exist.
		 * @param position the position in reference coordinates
		 * @param featureType the type which should be retrieved
		 * @return a possibly empty collection of features of the given type which annotate the given position
		 */
		public Collection<Feature> getFeaturesOfTypeForPosition(FeatureType featureType, int position);
		
}

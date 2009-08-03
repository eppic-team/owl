package proteinstructure.features;

import java.util.Collection;

/**
 * Classes implementing this interface provide methods to attach and retrieve sequence/structure features
 * complying to the {@link Feature} interface.
 * 
 * @author stehr
 */
public interface HasFeatures {

		/**
		 * Returns all features attached to the implementing object (disregarding the FeatureType).
		 * @return the features attached to the implementing object (disregarding the FeatureType)
		 */
		public Collection<Feature> getFeatures();
	
		/**
		 * Returns the FeatureTypes for which there are features attached to the implementing object.
		 * @return the FeatureTypes for which there are features attached to the implementing object
		 */
		public Collection<FeatureType> getFeatureTypes();
		
		/**
		 * Returns all features of a specific type which are attached to the implementing object.
		 * @param featureType the type of features which are to be retrieved
		 * @return the features of the given type or null if no such features exist
		 */
		public Collection<Feature> getFeaturesOfType(FeatureType featureType);
		
		/**
		 * Returns the features which annotate the given position or null if no such features exists.
		 * In most cases this method should return at most one object, unless overlapping object of the same type exist.
		 * @param position the position in reference coordinates
		 * @return the features of the given type which annotate the given position or null if no such features exists
		 */		
		public Collection<Feature> getFeaturesForPositon(int position);
		
		/**
		 * Returns the features of the given type which annotate the given position or null if no such features exists.
		 * In most cases this method should return at most one object, unless overlapping object of the same type exist.
		 * @param position the position in reference coordinates
		 * @param featureType the type which should be retrieved
		 * @return the features of the given type which annotate the given position or null if no such features exists
		 */
		public Collection<Feature> getFeaturesOfTypeForPosition(FeatureType featureType, int position);
		
		/**
		 * Attach the given feature to the implementing object.
		 * @param feature
		 */
		public void addFeature(Feature feature) throws InvalidFeatureCoordinatesException, OverlappingFeatureException;
}

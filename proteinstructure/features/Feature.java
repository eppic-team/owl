package proteinstructure.features;

import proteinstructure.IntervalSet;

/**
 * A feature in a sequence or structure.
 * Each subclass implementing this interface should have a corresponding instance of enum FeatureType.
 * Examples: Catalytic site, Scop domain.
 * 
 * See also: @{@link FeatureType}, {@link HasFeatures}
 * 
 * @author stehr
 */
public interface Feature {
	
	public IntervalSet getIntervalSet();
	public FeatureType getType();
	public String getDescription();
	
}

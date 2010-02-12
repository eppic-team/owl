package features;

import tools.Interval;
import tools.IntervalSet;

/**
 * A general sequence/structure feature.
 * @author stehr
 */
public class GeneralFeature implements Feature {

	/*--------------------------- member variables --------------------------*/

	FeatureType type;
	String description;
	IntervalSet position;
	
	/*----------------------------- constructors ----------------------------*/
	
	public GeneralFeature(IntervalSet position, String description) {
		this.type = FeatureType.GENERAL;
		this.position = position;
		this.description = description;
	}
	
	/*-------------------------- implemented methods ------------------------*/
	
	public String getDescription() {
		return this.description;
	}

	public IntervalSet getIntervalSet() {
		return this.position;
	}

	public FeatureType getType() {
		return this.type;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	public String toString() {
		return Interval.createSelectionString(this.getIntervalSet().getIntegerSet()) + " : " + this.description;
	}

}

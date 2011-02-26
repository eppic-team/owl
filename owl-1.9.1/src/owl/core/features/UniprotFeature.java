package owl.core.features;

import owl.core.util.Interval;
import owl.core.util.IntervalSet;

/**
 * A feature based on a Uniprot annotation.
 * @author stehr
 */
public class UniprotFeature implements Feature {

	/*--------------------------- member variables --------------------------*/
	IntervalSet position;
	String description;
	FeatureType type;
	
	// type specific data
	String uniprotTypeName;
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Creates a UniprotFeature based on a Feature object from the Uniprot Remote API.
	 */
	public UniprotFeature(int begRes, int endRes, String uniprotTypeName, String description) {
		this.type = FeatureType.UNIPROT;
		this.position = new IntervalSet();
		this.position.add(new Interval(begRes, endRes));
		this.description = description;
		this.uniprotTypeName = uniprotTypeName;
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
		return this.uniprotTypeName + " " + this.description + " : " + Interval.createSelectionString(this.getIntervalSet().getIntegerSet());
	}
	
	public String getUniprotTypeName() {
		return this.uniprotTypeName;
	}
}

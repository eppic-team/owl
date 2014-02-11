package owl.core.features;

import owl.core.util.IntervalSet;

public class StructuralDomainFeature implements Feature {

	/*--------------------------- member variables --------------------------*/
	IntervalSet position;
	String description;				// the name/id of this domain, e.g. 2RD0A1
	FeatureType type;
	
	// type specific data
	StructuralDomainType method;	// the method used to define this domain
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Createas a new StructuralDomainFeature.
	 */
	public StructuralDomainFeature(StructuralDomainType method, String name, IntervalSet position) {
		this.type = FeatureType.SDOMAIN;
		this.position = position;
		this.method = method;
		this.description = name;
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
	
	/**
	 * Returns a textual description of this domain
	 */
	public String toString() {
		return String.format("%s domain %s: %s", this.method.toString(), this.getName(), this.position.toString());
	}
	
	/**
	 * Returns the name/id of this domain, e.g. 2RD0A1
	 * @return
	 */
	public String getName() {
		return this.description;
	}
	
	/**
	 * Returns the method used to define this domain
	 * @return
	 */
	public StructuralDomainType getMethod() {
		return this.method;
	}	
}

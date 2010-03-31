package features;

/**
 * The type of a structurally defined domain.
 * Used by {@link StructuralDomainFeature}.
 * Based on the method available in the PDomains service (pdomains.sdsc.edu).
 */
public enum StructuralDomainType {
	
	/*------------------------------ instances ------------------------------*/
	
	CATH("CATH"),
	SCOP("SCOP"),
	DHCL("DHcL"),
	DP("dp"), 			
	PDP("pdp"), 		
	DDOMAIN("DDomain"), 
	NCBI("NCBI"),
	DUU("DUU");
	
	/*--------------------------- member variables --------------------------*/
	String parseString;
	
	/*----------------------------- constructors ----------------------------*/
	
	StructuralDomainType(String s) {
		this.parseString = s;			
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	public String toString() {
		return this.parseString;
	}
}

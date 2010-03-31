package owl.core.structure.decoyScoring;

/**
 * Enum to identify the different scoring methods
 * @author duarte
 *
 */
public enum ScoringMethod {

	RESCOUNT("residue count", "rescount"),
	ATOMCOUNT("atom count", "atomcount"),
	RESTYPE("residue type", "restype"),
	ATOMTYPE("atom type", "atomtype"),
	RESTRIPLET("residue triplet", "restriplet"),
	ATOMTRIPLET("atom triplet", "atomtriplet"),
	ATOMCOMBINED("atom combined", "atomcomb"),
	RESCOMBINED("residue combined", "rescomb"),
	ATOMDISTANCEDEP("atom distance dependent", "atomdist");
	
	
	private String id;
	private String description;
	
	private ScoringMethod(String description, String id) {
		this.description = description;
		this.id = id;
	}
	
	public String getDescription() {
		return description;
	}
	
	public String getId() {
		return id;
	}
	
	public static ScoringMethod getByDescription(String description) {
		for (ScoringMethod scMethod:values()) {
			if (scMethod.getDescription().equals(description)) {
				return scMethod;
			}
		}
		return null;
	}
}
package proteinstructure.DecoyScoring;

/**
 * Enum to identify the different scoring methods
 * @author duarte
 *
 */
public enum ScoringMethod {
	
	RESTYPE("residue type", "restype"),
	ATOMTYPE("atom type", "atomtype"),
	RESCOUNT("residue count", "rescount"),
	ATOMCOUNT("atom count", "atomcount"),
	ATOMCOMBINED("atom combined", "atomcomb"),
	RESCOMBINED("residue combined", "rescomb");
	
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
/**
 * 
 */
package owl.core.structure;

import java.util.HashMap;

/**
 * @author biyani_n
 *
 */
public enum BioUnitAssignmentType {
	authors(1,"authors"),
	pisa(2,"pisa"),
	pqs(3,"pqs"),
	eppic(4,"eppic"),
	none(5,"none");
	
	private int id;
	private String type;
	
	private static HashMap<String, BioUnitAssignmentType> string2type = initString2Type();
	
	//Default constructor
	private BioUnitAssignmentType(int id, String type){
		this.id = id;
		this.type = type;
	}
	
	/**
	 * initialize static map to get assignment type by its string
	 */
	private static HashMap<String, BioUnitAssignmentType> initString2Type() {
		HashMap<String, BioUnitAssignmentType> str2typeMap = new HashMap<String, BioUnitAssignmentType>();
		for(BioUnitAssignmentType type:BioUnitAssignmentType.values()) {
			str2typeMap.put(type.getType(), type);
		}
		return str2typeMap;
	}
	
	public static BioUnitAssignmentType getByString(String type) {
		if(string2type.containsKey(type)) return string2type.get(type.toLowerCase());
		else return none;
	}
	
	public String getType(){
		return this.type;
	}
	
	public int getId(){
		return this.id;
	}
}

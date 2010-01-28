package features;

import java.util.HashMap;
import java.util.Map;

/**
 * The type of a post-translational protein modification.
 * Used by PhosphoSitePlusFeature.
 */
public enum ProteinModificationType {
	
	/*------------------------------ instances ------------------------------*/
	
	PHOSPHORYLATION("p","Phosphorylation"),
	ACETYLATION("a","Acetylation"),
	METHYLATION("m","Methylation"),
	UBIQUITINATION("u","Ubiquitination"),
	SUMOYLATION("s","Sumoylation"),
	NEDDYLATYION("n","Neddylation"),
	GLYCOSYLATION("g","O-GlcNAc"),
	PALMITOYLATION("h","Palmitoylation");
	
	/*------------------------------ constants ------------------------------*/
	private static final Map<Character, ProteinModificationType> symb2type = initSymb2Type();
	
	/*--------------------------- member variables --------------------------*/
	
	private String symbol;
	private String description;
	
	/*----------------------------- constructors ----------------------------*/
	
	ProteinModificationType(String symbol, String description) {
		this.symbol = symbol;
		this.description = description;
	}
	
	/*---------------------------- private methods --------------------------*/
	private static Map<Character, ProteinModificationType> initSymb2Type() {
		Map<Character, ProteinModificationType> newMap = new HashMap<Character, ProteinModificationType>();
		for(ProteinModificationType t:ProteinModificationType.values()) {
			newMap.put(t.getSymbol(),t);
		}
		return newMap;
	}
	
	/*---------------------------- public methods ---------------------------*/
	public char getSymbol() { return this.symbol.charAt(0); }		
	public String getDescription() { return this.description; }
	public String toString() { return this.description; }
	
	/*---------------------------- static methods ---------------------------*/
	public static ProteinModificationType getBySymbol(char symb) {
		return symb2type.get(symb);
	}
}
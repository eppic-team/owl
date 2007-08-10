package proteinstructure;

/**
 * A particular secondary structure element within a protein structure
 */
public class SecStrucElement {
	
	/*--------------------------- member variables --------------------------*/
	
	char dsspType;
	Interval interval;
	
	/*----------------------------- constructors ----------------------------*/
	
	public SecStrucElement(char dsspType, int startRes, int endRes) {
		this.dsspType = dsspType;
		this.interval = new Interval(startRes, endRes);
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/** Returns the dssp type of this element. Valid values are H, B, E, G, I, T, S. */ 
	public char getDsspType() {
		return dsspType;
	}
	
	/** Returns the type of this ss element as one of the three states H, E, C. */
	public char getThreeStateType() {
		return dsspType;
	}
	
	/** Returns the range of this ss element in the sequence. */ 
	public Interval getInterval() {
		return interval;
	}
	
	public static char getSimpleTypeFromDsspType(char dsspType) {
		char type = dsspType;
		switch(dsspType) {
		case 'H':
		case 'G':
		case 'I': 
			type = 'H';
			break;
		case 'E':
		case 'B':
			type = 'S';
			break;
		case 'T':
			type = 'T';
			break;
		case 'S':
			type = ' ';
			break;
		}
		return type;
	}
}

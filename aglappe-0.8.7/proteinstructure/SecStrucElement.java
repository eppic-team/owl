package proteinstructure;

/**
 * A particular secondary structure element within a protein structure
 */
public class SecStrucElement {
	
	/*------------------------------ constants ------------------------------*/
	public enum ReducedState { THREESTATE, FOURSTATE, EIGHTSTATE };
	// three/four state secondary structure types (for 3 state, skip turn)
	public static final char HELIX = 'H';	// a helix
	public static final char STRAND = 'S';  // a beta strand
	public static final char TURN = 'T';    // a hydrogen bonded turn
	public static final char OTHER = 'O';   // all other states
	public static final char EXTENTED = 'E';// all extented 
	public static final char LOOP = 'L';  	// all other states 
	
	/*--------------------------- member variables --------------------------*/
	
	String secStrucId;		// legacy field for old ss ids (e.g. H1, S1, ...)
	char secStrucType;		// one of the above constants
	Interval interval;		// the location of this element in the sequence
	
	/*----------------------------- constructors ----------------------------*/
	
	public SecStrucElement(char secStrucType, int startRes, int endRes, String legacyId) {
		this.secStrucId = legacyId;
		this.secStrucType = secStrucType;
		this.interval = new Interval(startRes, endRes);
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	public SecStrucElement copy() {
		return new SecStrucElement(secStrucType, interval.beg, interval.end, secStrucId);
	}
	
	/** Returns the legacy ID of this element (e.g. H1, S1, ...) */
	public String getId() {
		return secStrucId;
	}
	
	/** Returns the dssp type of this element. Valid values are H, S, T, O */ 
	public char getType() {
		return secStrucType;
	}
	
	/** Returns the range of this ss element in the sequence. */ 
	public Interval getInterval() {
		return interval;
	}
	
	public char getSheetSerial() {
		if (isStrand()) {
			return Character.isLetter(secStrucId.charAt(1))?secStrucId.charAt(1):0;
		}
		return 0;
	}
	
	/** Returns true if this ss element is a helix */
	public boolean isHelix() {
		return secStrucType == HELIX;
	}
	
	/** Returns true if this ss element is a beta strand */
	public boolean isStrand() {
		return (secStrucType == STRAND || secStrucType == EXTENTED);
	}
	
	/** Returns true if this ss element is a hydrogen bonded turn */
	public boolean isTurn() {
		return secStrucType == TURN;
	}	
	
	/** Returns true if this ss element is other ... */
	public boolean isOther() {
		return (secStrucType == OTHER || secStrucType == LOOP);
	}	

	public boolean inSameSheet(SecStrucElement s) {
		boolean inSameSheet = false;
		if (s != null && this.isStrand() && s.isStrand() && (s.getSheetSerial() == this.getSheetSerial() && s.getSheetSerial() != 0)) {
			inSameSheet = true;
		}
		return inSameSheet;
	}
	
	/*---------------------------- static methods ---------------------------*/
	
	private static char getFourStateTypeFromDsspType(char dsspType) {
		char type = dsspType;
		switch(dsspType) {
		case 'H':
		case 'G':
		case 'I': 
			type = HELIX;
			break;
		case 'E':
			type = STRAND;
			break;
		case 'T':
			type = TURN;
			break;
		case 'S':
		case 'B':
			type = OTHER;
			break;
		default:
			type = OTHER;
		}
		return type;
	}
	
	private static char getThreeStateTypeFromDsspType(char dsspType) {
		char type = dsspType;
		switch(dsspType) {
		case 'H':
		case 'G':
		case 'I': 
			type = HELIX;
			break;
		case 'E':
		case 'B':
			type = EXTENTED;
			break;
		case 'T':
		case 'S':
			type = LOOP;
			break;
		default:
			type = LOOP;
		}
		return type;
	}
	
	public static char getReducedStateTypeFromDsspType(char dsspType, ReducedState state) {
		switch(state) {
		case THREESTATE:
			return getThreeStateTypeFromDsspType(dsspType);
		case FOURSTATE:
			return getFourStateTypeFromDsspType(dsspType);
		case EIGHTSTATE:
		default:
			return dsspType;
		}
	}
	
}

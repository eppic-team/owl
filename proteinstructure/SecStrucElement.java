package proteinstructure;

/**
 * A particular secondary structure element within a protein structure
 */
public class SecStrucElement {
	
	/*------------------------------ constants ------------------------------*/
	public enum ReducedState { THREESTATE, FOURSTATE, EIGHTSTATE };
	// four state secondary structure types
	public static final char HELIX = 'H';	// a helix
	public static final char STRAND = 'S';  // a beta strand
	public static final char TURN = 'T';    // a hydrogen bonded turn
	public static final char OTHER = 'O';   // all other states
	
	// three state secondary structure types, add also HELIX from above
	public static final char EXTENDED = 'E';// all extended 
	public static final char LOOP = 'L';  	// all other states 
	
	// eight state (dssp)
	// H, I, G, E, B, T, S, ' '
	
	/*--------------------------- member variables --------------------------*/
	
	String secStrucId;		// ss ids, i.e. the type + an element numerical identifier (e.g. H1, S1, ...)
	char secStrucType;		// ss types: one of the above char constants (3 in 3-state, 4 in 4-state, 8 in 8-state)
	Interval interval;		// the location of this element in the sequence 
	
	/*----------------------------- constructors ----------------------------*/
	
	public SecStrucElement(char secStrucType, int startRes, int endRes, String secStrucId) {
		this.secStrucId = secStrucId;
		this.secStrucType = secStrucType;
		this.interval = new Interval(startRes, endRes);
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Deep copies this ss element
	 * @return
	 */
	public SecStrucElement copy() {
		return new SecStrucElement(secStrucType, interval.beg, interval.end, secStrucId);
	}
	
	/** 
	 * Returns the ID of this element. The ID is the concatenation of the type 
	 * letter and the numerical element identifier (e.g. H1, S1, ...) 
	 * @return
	 */
	public String getId() {
		return secStrucId;
	}
	
	/** 
	 * Returns the type of this element. 
	 * The type depends on whether this secondary structure annotation has been 
	 * assigned with 3, 4 or 8 states (see constants above)
	 * @return
	 */ 
	public char getType() {
		return secStrucType;
	}
	
	/** 
	 * Returns the range of this ss element in the sequence.
	 * @return 
	 */ 
	public Interval getInterval() {
		return interval;
	}
	
	/**
	 * Returns the sheet serial of this ss element if it is a strand
	 * @return the sheet serial if this element is a strand and has a sheet 
	 * serial assigned, 0 otherwise
	 */
	public char getSheetSerial() {
		if (isStrand()) {
			return Character.isLetter(secStrucId.charAt(1))?secStrucId.charAt(1):0;
		}
		return 0;
	}
	
	/** 
	 * Returns true if this ss element is a helix 
	 * @return
	 */
	public boolean isHelix() {
		return secStrucType == HELIX;
	}
	
	/** 
	 * Returns true if this ss element is a beta strand 
	 * @return
	 */
	public boolean isStrand() {
		return (secStrucType == STRAND || secStrucType == EXTENDED);
	}
	
	/** 
	 * Returns true if this ss element is a hydrogen bonded turn 
	 * @return
	 */
	public boolean isTurn() {
		return secStrucType == TURN;
	}	
	
	/** 
	 * Returns true if this ss element is other 
	 * @return
	 */
	public boolean isOther() {
		return (secStrucType == OTHER || secStrucType == LOOP);
	}	

	/**
	 * Returns true if this ss element is in same sheet as given one
	 * @param s
	 * @return
	 */
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
			type = EXTENDED;
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
	
	public static char getThreeStateTypeFromPsiPredType(char psipredType) {
		switch (psipredType) {
		case 'C':
			return SecStrucElement.LOOP;
		case 'E':
			return SecStrucElement.EXTENDED;
		case 'H':
			return SecStrucElement.HELIX;
		}
		return 0; // if all fails we rather return nonsense
	}

	public boolean equals(Object o) {

		if (o==null) return false;
		if (!(o instanceof SecStrucElement)) return false;
		SecStrucElement other = (SecStrucElement) o;
		if (!this.interval.equals(other.interval)) 
			return false;
		if (!this.secStrucId.equals(other.secStrucId))
			return false;
		if (this.secStrucType!=other.secStrucType)
			return false;
		return true;
	}
}

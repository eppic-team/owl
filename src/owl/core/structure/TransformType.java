package owl.core.structure;

public enum TransformType {

	AU				(0,   1, false, "AU"), 
	XTALTRANSL		(1,   1, false, "XT"),  // translation
	CELLTRANSL		(2,   1, false, "FT"),  // fractional translation 
	
	TWOFOLD			(3,   2, false, "2" ), 
	TWOFOLDSCREW	(4,   2, true , "2S"), 
	
	THREEFOLD		(5,   3, false, "3" ), 
	THREEFOLDSCREW	(6,   3, true,  "3S"),
	
	FOURFOLD		(7,   4, false, "4" ), 
	FOURFOLDSCREW	(8,   4, true,  "4S"), 
	
	SIXFOLD       	(9,   6, false, "6" ), 
	SIXFOLDSCREW  	(10,  6, true,  "6S"),
	
	ONEBAR          (11, -1, false, "-1"),
	
	TWOBAR          (12, -2, false, "-2"),
	GLIDE           (13, -2, true,  "GL"),
	
	THREEBAR        (14, -3, false, "-3"),
	
	FOURBAR         (15, -4, false, "-4"),
	
	SIXBAR          (16, -6, false, "-6");
	
	
	
	private int id;
	private int foldType;
	private boolean isScrew;
	private String shortName;
	
	private TransformType(int id, int foldType, boolean isScrew, String shortName) {
		this.id = id;
		this.foldType = foldType;
		this.isScrew = isScrew;
		this.shortName = shortName;		
	}
	
	public int getId() {
		return id;
	}
	
	public int getFoldType() {
		return foldType;
	}
	
	/**
	 * Tells whether the transform is a screw or glide plane
	 * @return
	 */
	public boolean isScrew() {
		return isScrew;
	}
	 
	public String getShortName() {
		return shortName;
	}
}

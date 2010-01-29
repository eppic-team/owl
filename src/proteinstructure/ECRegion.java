package proteinstructure;

/**
 * A particular EC region within a protein structure
 */
public class ECRegion {
	
	/*--------------------------- member variables --------------------------*/
	
	String id;						// the EC number
	Interval interval;				// the location of this region in the sequence
	String startPdbRes, endPdbRes; 	// the starting and ending pdb residue serials of this region in the sequence
	
	/*----------------------------- constructors ----------------------------*/
	
	public ECRegion(String id, String startPdbRes, String endPdbRes, int startRes, int endRes) {
		this.id = id;
		this.interval = new Interval(startRes, endRes);
		this.startPdbRes = startPdbRes;
		this.endPdbRes = endPdbRes;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	public ECRegion copy() {
		return new ECRegion(id, startPdbRes, endPdbRes, interval.beg, interval.end);
	}
	
	/** Returns the ec number of this region */
	public String getECNum() {
		return id;
	}
	
	/** Returns the starting pdb residue serial of this region */ 
	public String getStartPdbRes() {
		return startPdbRes;
	}

	/** Returns the ending pdb residue serial of this region */ 
	public String getEndPdbRes() {
		return endPdbRes;
	}

	/** Returns the range of this region in the sequence. */ 
	public Interval getInterval() {
		return interval;
	}
	
	/*---------------------------- static methods ---------------------------*/
}	

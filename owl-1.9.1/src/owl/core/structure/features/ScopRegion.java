package owl.core.structure.features;

import owl.core.util.Interval;

/**
 * A particular region of a scop domain within a protein structure
 */
public class ScopRegion {
	
	/*------------------------------ constants ------------------------------*/
	public enum DomainType { UNKNOWN, WHOLECHAIN, SINGLEFRAGMENT, MULTIFRAGMENT, MULTICHAIN };

	/*--------------------------- member variables --------------------------*/
	
	String sid; 					// old SCOP identifier e.g. d1dlwa_
	String sccs;					// SCOP concise classification strings (sccs) e.g. a.1.1.1
	int sunid;						// SCOP unique identifiers e.g. 14982
	int orderIn;					// order of the region in the domain as in the original gene sequence (not the same as pdb)
	int numRegions;					// the number of regions of the domain
	Interval interval;				// the location of this region in the sequence
	String startPdbRes, endPdbRes; 	// the starting and ending pdb residue serials of this region in the sequence
	/*----------------------------- constructors ----------------------------*/
	
	public ScopRegion(String sid, String sccs, int sunid, int orderIn, int numRegions, String startPdbRes, String endPdbRes, int startRes, int endRes) {
		this.sid = sid;
		this.sccs = sccs;
		this.sunid = sunid;
		this.orderIn = orderIn;
		this.numRegions = numRegions;
		this.interval = new Interval(startRes, endRes);
		this.startPdbRes = startPdbRes;
		this.endPdbRes = endPdbRes;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	public ScopRegion copy() {
		return new ScopRegion(sid, sccs, sunid, orderIn, numRegions, startPdbRes, endPdbRes, interval.beg, interval.end);
	}
	
	/** Returns the old scop id of this region */
	public String getSId() {
		return sid;
	}
	
	/** Returns the sccs of this region */ 
	public String getSccs() {
		return sccs;
	}

	/** Returns the scop id of this region */ 
	public int getSunid() {
		return sunid;
	}
	
	/** Returns the order of this region */ 
	public int getOrder() {
		return orderIn;
	}
	
	/** Returns the number of regions of a domain*/ 
	public int getNumRegions() {
		return numRegions;
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
	
	/** Returns the domain type of this region. */	
	public DomainType getDomainType() {
		if (sid.endsWith("_")) {
			return DomainType.WHOLECHAIN;
		} else if (sid.charAt(sid.length()-2) == '.') {
			return DomainType.MULTICHAIN;
		} else if (numRegions == 1) {
			return DomainType.SINGLEFRAGMENT;
		} else if (numRegions > 1) {
			return DomainType.MULTIFRAGMENT;
		} else
			return DomainType.UNKNOWN;
	}
	
}	

package proteinstructure;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

/** 
 * This class encapsulates the scop annotation of a single protein chain. 
 */
public class Scop {

	/*------------------------------ constants ------------------------------*/
	
	public static final String LATEST_VERSION = "1.73"; // default version
	
	/*--------------------------- member variables --------------------------*/
	
	private HashMap<Integer,ScopRegion> resser2scopregion;  // residue serials to scop region
	private Vector<ScopRegion> scopRegions;			 		// the actual collection of scop regions
	private String version;									// the version 
	
	/*----------------------------- constructors ----------------------------*/
	
	/** 
	 * Create an empty scop object 
	 */
	public Scop() {
		this.scopRegions = new Vector<ScopRegion>();
		this.resser2scopregion = new HashMap<Integer,ScopRegion>();
		this.version = LATEST_VERSION;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Add the given scop region to this scop object.
	 */
	public void add(ScopRegion e) {
		this.scopRegions.add(e);
		Interval intv = e.getInterval();
		for(int i = intv.beg; i <= intv.end; i++) {
			this.resser2scopregion.put(i, e);
		}
	}
	
	public void remove(ScopRegion e) {
		Interval intv = e.getInterval();
		for(int i = intv.beg; i <= intv.end; i++) {
			this.resser2scopregion.remove(i);
		}
		this.scopRegions.remove(e);
	}
	
	/**
	 * Returns true iff this scop object contains no scop regions.
	 */
	public boolean isEmpty() {
		return this.scopRegions.isEmpty();
	}
	
	/** 
	 * Returns an iterator over all assigned scop regions in this object.
	 */
	public Iterator<ScopRegion> getIterator() {
		return scopRegions.iterator();
	}
	
	/** 
	 * For a given residue returns the scop region this residue participates in,
	 * or null if the residue is not in an assigned scop region.
	 */
	public ScopRegion getScopRegion(int resser){
		if(resser2scopregion.containsKey(resser)) {
			return this.resser2scopregion.get(resser);
		} else {
			return null;
		}
	}
	
	/** 
	 * Return a deep copy of this scop object 
	 */
	public Scop copy() {
		Scop newS = new Scop();
		for(ScopRegion e:scopRegions) {
			newS.add(e.copy());
		}
		return newS;
	}
		
	/** 
	 * Sets the optional version to v 
	 */
	public void setVersion(String v) {
		this.version = v;
	}
	
	/** 
	 * Returns the version 
	 */
	public String getVersion() {
		return this.version;
	}
	
	/** Returns the number of scop regions in this scop object */
	public int getNumRegions() {
		return this.scopRegions.size();
	}
	
	/**
	 * Returns a Vector of scop regions for the given sid
	 * @param sid
	 * @return
	 */
	public Vector<ScopRegion> getScopRegions (String sid) {
		Vector<ScopRegion> restrScopRegions = new Vector<ScopRegion>();
		Iterator<ScopRegion> it = scopRegions.iterator();
		while(it.hasNext()) {
			ScopRegion scopRegion = it.next();
			if (scopRegion.getSId().equals(sid)) {
				restrScopRegions.add(scopRegion);
			}
		}
		return restrScopRegions;
	}
	
	/**
	 * Returns a Vector of scop regions for the given sunid
	 * @param sunid
	 * @return
	 */
	public Vector<ScopRegion> getScopRegions (int sunid) {
		Vector<ScopRegion> restrScopRegions = new Vector<ScopRegion>();
		Iterator<ScopRegion> it = scopRegions.iterator();
		while(it.hasNext()) {
			ScopRegion scopRegion = it.next();
			if (scopRegion.getSunid()==(sunid)) {
				restrScopRegions.add(scopRegion);
			}
		}
		return restrScopRegions;
	}
}

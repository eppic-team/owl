package owl.core.structure;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import owl.core.util.Interval;


/** This class encapsulates the EC annotation of a single protein chain. */
public class EC {

	/*--------------------------- member variables --------------------------*/
	
	private HashMap<Integer,HashSet<ECRegion>> resser2ecregion;  	// residue serials to set of ec regions
	private Vector<ECRegion> ecRegions;			 					// the actual collection of ec regions
	
	/*----------------------------- constructors ----------------------------*/
	
	/** Create an empty ec object */
	public EC() {
		this.ecRegions = new Vector<ECRegion>();
		this.resser2ecregion = new HashMap<Integer,HashSet<ECRegion>>();
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Add the given ec region to this ec object.
	 */
	public void add(ECRegion e) {
		this.ecRegions.add(e);
		Interval intv = e.getInterval();
		for(int i = intv.beg; i <= intv.end; i++) {
			HashSet<ECRegion> values = resser2ecregion.get(i);
			if (values == null) {
				values = new HashSet<ECRegion>();
				resser2ecregion.put(i,values);
			}
			values.add(e);
		}
	}
	
	/**
	 * Returns true iff this ec object contains no scop regions.
	 */
	public boolean isEmpty() {
		return this.ecRegions.isEmpty();
	}
	
	/** 
	 * Returns an iterator over all assigned ec regions in this object.
	 */
	public Iterator<ECRegion> getIterator() {
		return ecRegions.iterator();
	}
	
	/** 
	 * For a given residue returns the ec regions this residue participates in,
	 * or null if the residue is not in an assigned ec region.
	 */
	public HashSet<ECRegion> getECRegion(int resser){
		if(resser2ecregion.containsKey(resser)) {
			return this.resser2ecregion.get(resser);
		} else {
			return null;
		}
	}
	
	public String getECNum(int resser) {
		String ecNums = null;
		if (resser2ecregion.containsKey(resser)) {
			ecNums = "";
			HashSet<ECRegion> ecs = getECRegion(resser);
			Iterator<ECRegion> it = ecs.iterator();
			while (it.hasNext()) {
				ECRegion er = (ECRegion)it.next();
				ecNums += er.getECNum();
				if (it.hasNext()) {
					ecNums += ",";
				}
			}
		}
		return ecNums;
	}
	
	/** Return a deep copy of this ec object */
	public EC copy() {
		EC newS = new EC();
		for(ECRegion e:ecRegions) {
			newS.add(e.copy());
		}
		return newS;
	}
		
	/** Returns the number of ec regions in this ec object */
	public int getNumRegions() {
		return this.ecRegions.size();
	}
}

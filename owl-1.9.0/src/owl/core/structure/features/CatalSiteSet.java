package owl.core.structure.features;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

/** This class encapsulates the catalytic site annotation of a single protein chain. */
public class CatalSiteSet {

	/*------------------------------ constants ------------------------------*/
	
	//private static final String LATEST_VERSION = "2.2.7"; // default version
	//public static final String LATEST_VERSION = "2.2.10"; // as of 2008-10-02
	public static final String LATEST_VERSION = "2.2.11"; // as of 2009-08-07
	
	/*--------------------------- member variables --------------------------*/
	
	private HashMap<Integer,HashSet<CatalyticSite>> resser2catalsite;   // residue serials to set of catalytic sites
	private Vector<CatalyticSite> catalyticSites;						// the actual collection of catalytic sites
	private String version;										 		// an optional comment describing the version of the catalytic site annotation
	
	/*----------------------------- constructors ----------------------------*/
	
	/** Create an empty catalytic site object */
	public CatalSiteSet() {
		this.catalyticSites = new Vector<CatalyticSite>();
		this.resser2catalsite = new HashMap<Integer,HashSet<CatalyticSite>>();
		this.version = LATEST_VERSION;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Add the given catalytic site to this set of catalytic sites object.
	 */
	public void add(CatalyticSite e) {
		this.catalyticSites.add(e);
		Set<Integer> ressers = e.getRes();
		Iterator<Integer> iter = ressers.iterator();
	    while (iter.hasNext()) {
	    	int resser = (Integer)iter.next();
			HashSet<CatalyticSite> values = resser2catalsite.get(resser);
			if (values == null) {
				values = new HashSet<CatalyticSite>();
				resser2catalsite.put(resser,values);
			}
			values.add(e);
	    }
	}
	
	/**
	 * Returns true iff this catalytic sites set object contains no catalytic sites.
	 */
	public boolean isEmpty() {
		return this.catalyticSites.isEmpty();
	}
	
	/** 
	 * Returns an iterator over all assigned catalytic sites in this object.
	 */
	public Iterator<CatalyticSite> getIterator() {
		return catalyticSites.iterator();
	}
		
	/** 
	 * For a given residue returns the catalytic sites this residue participates in,
	 * or null if the residue is not in an assigned catalytic site.
	 */
	public HashSet<CatalyticSite> getCatalSite(int resser){
		if(resser2catalsite.containsKey(resser)) {
			return this.resser2catalsite.get(resser);
		} else {
			return null;
		}
	}
	
	public String getCatalSiteNum(int resser) {
		String csNums = null;
		if (resser2catalsite.containsKey(resser)) {
			csNums = "";
			HashSet<CatalyticSite> css = getCatalSite(resser);
			Iterator<CatalyticSite> it = css.iterator();
			while (it.hasNext()) {
				CatalyticSite cs = (CatalyticSite)it.next();
				csNums += cs.getSerial();
				if (it.hasNext()) {
					csNums += ",";
				}
			}
		}
		return csNums;
	}
	
	public String getCatalSiteChemFunc(int resser) {
		String csChemFuncs = null;
		if (resser2catalsite.containsKey(resser)) {
			csChemFuncs = "";
			HashSet<CatalyticSite> css = getCatalSite(resser);
			Iterator<CatalyticSite> it = css.iterator();
			while (it.hasNext()) {
				CatalyticSite cs = (CatalyticSite)it.next();
				csChemFuncs += cs.getChemFuncFromResSerial(resser);
				if (it.hasNext()) {
					csChemFuncs += ",";
				}
			}
		}
		return csChemFuncs;
	}
	
	public String getCatalSiteEvid(int resser) {
		String csEvids = null;
		if (resser2catalsite.containsKey(resser)) {
			csEvids = "";
			HashSet<CatalyticSite> css = getCatalSite(resser);
			Iterator<CatalyticSite> it = css.iterator();
			while (it.hasNext()) {
				CatalyticSite cs = (CatalyticSite)it.next();
				csEvids += cs.getEvidence();
				if (it.hasNext()) {
					csEvids += ",";
				}
			}
		}
		return csEvids;
	}
	
	public void removeCatalSiteRes(int resser) {
		if (resser2catalsite.containsKey(resser)) {
			HashSet<CatalyticSite> css = getCatalSite(resser);
			Iterator<CatalyticSite> it = css.iterator();
			while (it.hasNext()) {
				CatalyticSite cs = (CatalyticSite)it.next();
				cs.remRes(resser);
				if (cs.getRes().size() == 0) {
					catalyticSites.remove(cs);
				}
			}
			resser2catalsite.remove(resser);
		}		
	}
	
	/** Return a deep copy of this catalytic site set object */
	public CatalSiteSet copy() {
		CatalSiteSet newS = new CatalSiteSet();
		for(CatalyticSite e:catalyticSites) {
			newS.add(e.copy());
		}
		return newS;
	}
		
	/** Sets the version to c */
	public void setVersion(String c) {
		this.version = c;
	}
	
	/** Returns the version */
	public String getVersion() {
		return this.version;
	}
	
	/** Returns the number of catalytic sites in this catalytic sites set object */
	public int getNumCatalSites() {
		return this.catalyticSites.size();
	}
	
	public void print() {
		Iterator<CatalyticSite> iterator = catalyticSites.iterator ();
        while (iterator.hasNext ()) {
        	CatalyticSite cs = (CatalyticSite)iterator.next();
        	cs.print();
        } 
	}
	
	public String toString() {
		String msg = "";
		Iterator<CatalyticSite> iterator = catalyticSites.iterator ();
        while (iterator.hasNext ()) {
        	CatalyticSite cs = (CatalyticSite)iterator.next();
        	msg += cs.toString();
        } 		
        return msg;
	}
}

package proteinstructure;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

/** This class encapsulates the secondary structure annotation of a single protein chain. */
public class SecondaryStructure {

	/*------------------------------ constants ------------------------------*/
	
	private static final String INITIAL_COMMENT = "None"; // default comment
	
	/*--------------------------- member variables --------------------------*/
	
	private HashMap<Integer,SecStrucElement> resser2secstruct;   // residue serials to secondary structure
	private Vector<SecStrucElement> secStructElements;			 // the actual collection of secondary structure elements
	private String comment;										 // an optional comment describing this secondary structure annotation
	
	/*----------------------------- constructors ----------------------------*/
	
	/** Create an empty secondary structure object */
	public SecondaryStructure() {
		this.secStructElements = new Vector<SecStrucElement>();
		this.resser2secstruct = new HashMap<Integer,SecStrucElement>();
		this.comment = INITIAL_COMMENT;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Add the given secondary structure element to this secondary structure object.
	 */
	public void add(SecStrucElement e) {
		this.secStructElements.add(e);
		Interval intv = e.getInterval();
		for(int i = intv.beg; i <= intv.end; i++) {
			this.resser2secstruct.put(i, e);
		}
	}
	
	/**
	 * Returns true iff this secondary structure object contains no secondary structure elements.
	 */
	public boolean isEmpty() {
		return this.secStructElements.isEmpty();
	}
	
	/** 
	 * Returns an iterator over all assigned secondary structure elements in this object.
	 */
	public Iterator<SecStrucElement> getIterator() {
		return secStructElements.iterator();
	}
	
	/** 
	 * For a given residue returns the secondary structure element this residue participates in,
	 * or null if the residue is not in an assigned secondary structure element.
	 */
	public SecStrucElement getSecStrucElement(int resser){
		if(resser2secstruct.containsKey(resser)) {
			return this.resser2secstruct.get(resser);
		} else {
			return null;
		}
	}
	
	/** Return a deep copy of this secondary structure object */
	public SecondaryStructure copy() {
		SecondaryStructure newSS = new SecondaryStructure();
		for(SecStrucElement e:secStructElements) {
			newSS.add(e.copy());
		}
		return newSS;
	}
		
	/** Sets the optional comment to c */
	public void setComment(String c) {
		this.comment = c;
	}
	
	/** Returns the comment */
	public String getComment() {
		return this.comment;
	}
	
	/** Returns the number of secondary structure elements in this SecondaryStructure object */
	public int getNumElements() {
		return this.secStructElements.size();
	}
	
	/**
	 * Returns true if the given SecStrucElement is already in this SecondaryStructure object 
	 * @param sselem
	 * @return
	 */
	public boolean contains(SecStrucElement sselem) {
		return secStructElements.contains(sselem);
	}
}

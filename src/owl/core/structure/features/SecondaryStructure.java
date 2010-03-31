package owl.core.structure.features;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Vector;

import owl.core.util.Interval;
import owl.graphAveraging.ConsensusInterval;
import owl.graphAveraging.ConsensusSquare;


/** This class encapsulates the secondary structure annotation of a single protein chain.
 * TODO eventually we should make this implement Iterable<SecStrucElement> (so then one can 
 * use foreach loops). The problem is we'd need to refactor getIterator() to iterator() and
 * it's difficult because there are many references all over the place some of them conflicting
 * with other local variable names 
 */
public class SecondaryStructure {

	/*------------------------------ constants ------------------------------*/
	
	private static final String INITIAL_COMMENT = "None"; // default comment
	
	/*--------------------------- member variables --------------------------*/
	
	private HashMap<Integer,SecStrucElement> resser2secstruct;   // residue serials to secondary structure
	private Vector<SecStrucElement> secStructElements;			 // the actual collection of secondary structure elements
	private String comment;										 // an optional comment describing this secondary structure annotation
	
	private HashMap<Integer,Double> resser2predConfidence; 		 // residue serials to secondary structure prediction confidence 
																 //(this is a per residue value, that's why we can't put it in SecStruElement)
																 // it will only be used for predictions, otherwise it will be null
	private String sequence;									 // the sequence for this secondary structure annotation
	
	
	/*----------------------------- constructors ----------------------------*/
	
	/** 
	 * Create an empty secondary structure object 
	 */
	public SecondaryStructure(String sequence) {
		//TODO we should check whether the given sequence is consistent with the other data (residue numbers used for sec. struct. intervals and so on)
		this.sequence = sequence; 
		this.secStructElements = new Vector<SecStrucElement>();
		this.resser2secstruct = new HashMap<Integer,SecStrucElement>();
		this.comment = INITIAL_COMMENT;
	}
	
	
	/**
	 * To create a secondary structure object by reading a psipred prediction horizontal file
	 * @param psipredHorizFile
	 * @throws IOException
	 */
	public SecondaryStructure(File psipredHorizFile) throws IOException{
		this.secStructElements = new Vector<SecStrucElement>();
		this.resser2secstruct = new HashMap<Integer,SecStrucElement>();
		this.comment = "Psipred prediction";
		this.resser2predConfidence = new HashMap<Integer, Double>();
		this.readPsiPredHorizFile(psipredHorizFile);
	}
	
	/*---------------------------- private methods ---------------------------*/
	
	/**
	 * Parses a psipred prediction horizontal file
	 * @param psipredHorizFile
	 * @throws IOException
	 */
	private void readPsiPredHorizFile(File psipredHorizFile) throws IOException {

		BufferedReader br = new BufferedReader(new FileReader(psipredHorizFile));
		String line;
		this.sequence = "";
		int resSerial = 0;
		char ssType = 0;
		char lastSsType = 0;
		int elemSerial = 0;
		int startRes = 1;
		// we keep a second residue serial counter for the confidence assignment
		int resSerial4Confidence = 0;
		
		while ((line=br.readLine())!=null) {
			if (line.length()==0) continue;
			if (line.startsWith("Conf: ")) {
				String confStr = line.substring(6);
				for (int i=0;i<confStr.length();i++) {
					resSerial4Confidence++;
					double confidence = Double.parseDouble("0."+String.valueOf(confStr.charAt(i)));
					this.resser2predConfidence.put(resSerial4Confidence, confidence);
				}
			} 
			else if (line.startsWith("Pred: ")) {
				String predStr = line.substring(6);
				for (int i=0;i<predStr.length();i++) {
					resSerial ++;
					
					ssType = predStr.charAt(i);
					if (lastSsType!=0 && ssType!=lastSsType){
						elemSerial++;
						char secStrucType = SecStrucElement.getThreeStateTypeFromPsiPredType(lastSsType);
						int endRes = resSerial-1; 
						String secStrucId = String.valueOf(secStrucType)+elemSerial;
						this.add(new SecStrucElement(secStrucType, startRes, endRes, secStrucId));
						
						startRes = resSerial; 
					}
					lastSsType = ssType;
	
				}
			}
			else if (line.startsWith("  AA: ")) {
				String aaStr = line.substring(6);
				for (int i=0;i<aaStr.length();i++) {
					sequence += aaStr.charAt(i);
				}				
			}
		}
		// adding the last element
		char secStrucType = SecStrucElement.getThreeStateTypeFromPsiPredType(lastSsType);
		this.add(new SecStrucElement(secStrucType, startRes, resSerial, String.valueOf(secStrucType)+(elemSerial+1)));
	}
	
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Add the given secondary structure element to this secondary structure object.
	 * @param e
	 */
	public void add(SecStrucElement e) {
		this.secStructElements.add(e);
		Interval intv = e.getInterval();
		for(int i = intv.beg; i <= intv.end; i++) {
			this.resser2secstruct.put(i, e);
		}
	}
	
	/**
	 * Returns true if this secondary structure object contains no secondary structure elements.
	 * @return
	 */
	public boolean isEmpty() {
		return this.secStructElements.isEmpty();
	}
	
	/** 
	 * Returns an iterator over all assigned secondary structure elements in this object.
	 * @return
	 */
	public Iterator<SecStrucElement> getIterator() {
		return secStructElements.iterator();
	}
	
	/** 
	 * For a given residue returns the secondary structure element this residue participates in,
	 * or null if the residue is not in an assigned secondary structure element.
	 * @param resser
	 * @return
	 */
	public SecStrucElement getSecStrucElement(int resser){
		if(resser2secstruct.containsKey(resser)) {
			return this.resser2secstruct.get(resser);
		} else {
			return null;
		}
	}
	
	/** 
	 * Return a deep copy of this secondary structure object
	 * @return 
	 */
	public SecondaryStructure copy() {
		SecondaryStructure newSS = new SecondaryStructure(this.getSequence());
		for(SecStrucElement e:secStructElements) {
			newSS.add(e.copy());
		}
		return newSS;
	}
		
	/** 
	 * Sets the optional comment to c
	 * @param c 
	 */
	public void setComment(String c) {
		this.comment = c;
	}
	
	/** 
	 * Returns the comment
	 * @return 
	 */
	public String getComment() {
		return this.comment;
	}
	
	/**
	 * Sets the sequence for this SecondaryStructure
	 * @param sequence
	 */
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	
	/**
	 * Returns the sequence for this SecondaryStructure
	 * @return
	 */
	public String getSequence() {
		return this.sequence;
	}
	
	/** 
	 * Returns the number of secondary structure elements in this SecondaryStructure object
	 * @return 
	 */
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
	
	/**
	 * Returns the confidence value (if this is a PsiPred sec. structure prediction)
	 * @param resser
	 * @return
	 * @throws NullPointerException if this is not a Psipred prediction and thus the confidence is not set
	 */
	public double getConfidence(int resser) {
		return this.resser2predConfidence.get(resser);
	}
	
	public TreeMap<Integer,ConsensusSquare> getPhiPsiConstraints() {
		TreeMap<Integer, ConsensusSquare> bounds = new TreeMap<Integer, ConsensusSquare>();
		for (int i=1; i <= sequence.length();i++) {
			if (resser2secstruct.get(i).getType() == SecStrucElement.HELIX) {
				bounds.put(i, new ConsensusSquare(new ConsensusInterval(-65,-55),new ConsensusInterval(-50,-40)));
			} else if (resser2secstruct.get(i).getType() == SecStrucElement.STRAND) {
				bounds.put(i, new ConsensusSquare(new ConsensusInterval(-135,-100),new ConsensusInterval(160,90)));
			}
		}
		return bounds;
	}

}

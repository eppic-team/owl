package owl.core.structure.features;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Vector;

import owl.core.util.Interval;
import owl.graphAveraging.ConsensusInterval;
import owl.graphAveraging.ConsensusSquare;


/** 
 * This class encapsulates the secondary structure annotation of a single protein chain.
 * Adapters for the following sources of secondary structure annotation have been implemented:
 * - PDB (see owl.core.structure.Pdb)
 * - DSSP (see owl.core.runners.DsspRunner)
 * - PsiPred (see owl.core.runners.PsipredRunner)
 * - JPred (see owl.core.connections.JPredConnection)
 * - Consensus (see getConsensusSecondaryStructure())
 * 
 */
public class SecondaryStructure implements Iterable<SecStrucElement> {

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
		this.resser2predConfidence = new HashMap<Integer, Double>();
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
	 * Parses a psipred prediction horizontal file.
	 * TODO: This should be in PsipredRunner.
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
	 * Sets the confidence value for the given residue.
	 * @param resNum the residue number
	 * @param conf the confidence value to assign
	 */
	public void setConfidence(int resNum, double conf) {
		this.resser2predConfidence.put(resNum, conf);
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
	public Iterator<SecStrucElement> iterator() {
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
	
	/**
	 * Gets phi/psi constraints for each position, e.g. for distance geometry.
	 * @return a map from position to ConsensusSquare with phi/psi constraints
	 */
	public TreeMap<Integer,ConsensusSquare> getPhiPsiConstraints() {
		TreeMap<Integer, ConsensusSquare> bounds = new TreeMap<Integer, ConsensusSquare>();
		for (int i=1; i <= sequence.length();i++) {
			SecStrucElement a = resser2secstruct.get(i);
			if (a == null) {
				continue;
			}
			if (resser2secstruct.get(i).getType() == SecStrucElement.HELIX) {
				bounds.put(i, new ConsensusSquare(new ConsensusInterval(-65,-55),new ConsensusInterval(-50,-40)));
			} else if (resser2secstruct.get(i).getType() == SecStrucElement.STRAND) {
				bounds.put(i, new ConsensusSquare(new ConsensusInterval(-135,-100),new ConsensusInterval(160,90)));
			}
		}
		return bounds;
	}

	/**
	 * Returns a new (three state) secondary structure object which represents the consensus of
	 * the given collection of secondary structure objects. It is assumed that all the
	 * objects in the list are three state and refer to the same sequence. For a given position
	 * i, the consensus will be assigned e.g. a helix if the fraction of individuals where i
	 * is in a helix is at least the threshold. Confidence values are currently ignored.
	 * @param the common underlying sequence
	 * @param ssList a collection of secondary structure objects
	 * @param consensusSSthresh the consensus threshold
	 * @return the consensus secondary structure object
	 */
	public static SecondaryStructure getConsensusSecondaryStructure(
			String sequence, Collection<SecondaryStructure> ssList, double thresh) {
		SecondaryStructure newSS = new SecondaryStructure(sequence);
		newSS.setComment("Consensus");
		newSS.setSequence(sequence);
		int numRes = sequence.length();
		int[] helix = new int[numRes+1];	// use sequence indices
		int[] extended = new int[numRes+1];	
		int[] total = new int[numRes+1];
		for(SecondaryStructure ss: ssList) {
			for (int i = 1; i <= numRes; i++) {
				total[i]++;	// conservative approach: take all templates
				SecStrucElement e = ss.getSecStrucElement(i);
				if(e != null) {
					//total[i]++ // alternative approach: only count templates where ss is at least 'loop' 
					if(e.isHelix()) helix[i]++;
					if(e.isStrand()) extended[i]++;
				}
			}
		}
		StringBuilder s = new StringBuilder(numRes);
		for (int i = 1; i <= numRes; i++) {
			char type = SecStrucElement.LOOP;
			if(extended[i] > helix[i]) {
				if(1.0 * extended[i] / total[i] >= thresh) {
					type = SecStrucElement.EXTENDED;
					newSS.resser2predConfidence.put(i, 1.0 * extended[i] / total[i]);
				}
			} else {
				if(1.0 * helix[i] / total[i] >= thresh) {
					type = SecStrucElement.HELIX;
					newSS.resser2predConfidence.put(i, 1.0 * helix[i] / total[i]);
				}
			}
			s.append(type);
		}
		parseSecondaryStructureString(newSS, s.toString());
		return newSS;
	}
	
	/**
	 * Helper function to parse secondary structure elements from a simple
	 * string representation.
	 * @param mapping
	 */
	private static void parseSecondaryStructureString(SecondaryStructure ss, String ssString) {
		int lastResSer = 1;
		char lastType = ssString.charAt(0);
		int start = 1;
		char thisType;
		int elementCount = 0;
		SecStrucElement ssElem;
		String ssId;
		int resSer;
		for (resSer = 1; resSer <= ssString.length(); resSer++) {
			thisType = ssString.charAt(resSer-1);
			if(thisType != lastType || resSer > lastResSer+1) {
				// finish previous element, start new one
				elementCount++;
				ssId = new Character(lastType).toString() + new Integer(elementCount).toString();
				ssElem = new SecStrucElement(lastType,start,lastResSer,ssId);
				ss.add(ssElem);
				start = resSer;
				lastType = thisType;
			}
			lastResSer = resSer;
		}
		// finish last element
		elementCount++;
		ssId = new Character(lastType).toString() + new Integer(elementCount).toString();
		ssElem = new SecStrucElement(lastType, start,ssString.length(),ssId);
		ss.add(ssElem);
	}


	/**
	 * Prints the secondary structure annotation to stdout.
	 */
	public void print() {
		String s = getSequence();
		System.out.println(s);
		for (int i = 1; i <= s.length(); i++) {
			System.out.print(resser2secstruct.containsKey(i)?resser2secstruct.get(i).secStrucType:"-");
		}
		System.out.println();
		if(resser2predConfidence != null) {
			for (int i = 1; i <= s.length(); i++) {
				if(resser2predConfidence.containsKey(i)) {
					double rawConf = resser2predConfidence.get(i);
					int conf = (int) (9.9999 * rawConf);
					System.out.print(conf);
				} else {
					System.out.print("-");
				}
			}
			System.out.println();
		}
	}

}

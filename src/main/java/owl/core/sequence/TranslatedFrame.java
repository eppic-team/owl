package owl.core.sequence;

import java.io.Serializable;

import owl.core.sequence.alignment.PairwiseSequenceAlignment;

/**
 * A class to encapsulate a sequence translation together with the reading frame used
 * for translation and the alignment to its reference protein sequence.
 * 
 * Comparable based on the sequence identities of the alignments.
 * 
 * @author duarte_j
 *
 */
public class TranslatedFrame implements Comparable<TranslatedFrame>, Serializable {

	private static final long serialVersionUID = 1L;
	
	private Sequence sequence; // the translated sequence (protein)
	private PairwiseSequenceAlignment psa; // the alignment between the the reference protein (1) sequence and the translated sequence (2)
	private ReadingFrame rf;

	public TranslatedFrame(Sequence sequence, ReadingFrame rf) {
		this.sequence = sequence;
		this.rf = rf;
	}
	
	public void setAln(PairwiseSequenceAlignment psa){
		this.psa = psa;
	}
	
	public float getPercentIdentity() {
		return psa.getPercentIdentity();
	}
	
	public int getMismatches() {
		return psa.getLength()-psa.getIdentity();
	}
	
	public PairwiseSequenceAlignment getAln(){
		return psa;
	}
	
	public int compareTo(TranslatedFrame o) {
		return Float.compare(this.getPercentIdentity(), o.getPercentIdentity());
	}
	
	public Sequence getSequence() {
		return sequence;
	}
	
	public ReadingFrame getReadingFrame() {
		return rf;
	}
	
	public int getNumMismatches() {
		return this.psa.getLength()-this.psa.getIdentity();
	}
	
	public int getNumGaps() {
		return this.psa.getGaps();
	}
	
	/**
	 * Given a reference protein sequence position, tells whether the translated sequence 
	 * matches (is an identity)
	 * @param i a reference protein sequence position (starting at 0)
	 * @return true if the position matches the translated sequence, false otherwise
	 */
	public boolean isMatch(int i) {
		return this.psa.isMatchingTo2(i);
	}
		
}

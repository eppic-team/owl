package owl.core.sequence;

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
public class TranslatedFrame implements Comparable<TranslatedFrame> {
	private Sequence sequence; // the translated sequence (protein)
	private PairwiseSequenceAlignment psa; // the alignment between the translated sequence and the reference protein sequence
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
	
	@Override
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
		
}

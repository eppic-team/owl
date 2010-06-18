package owl.core.sequence;

import owl.core.sequence.alignment.PairwiseSequenceAlignment;

/**
 * A class to encapsulate a sequence translation together with the reading frame used
 * for translation and the alignment to its reference protein sequence.
 * @author duarte_j
 *
 */
public class TranslatedSequence implements Comparable<TranslatedSequence> {
	private Sequence sequence;
	private PairwiseSequenceAlignment psa;
	private ReadingFrame rf;

	public TranslatedSequence(Sequence sequence, ReadingFrame rf) {
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
	public int compareTo(TranslatedSequence o) {
		return Float.compare(this.getPercentIdentity(), o.getPercentIdentity());
	}
	
	public Sequence getSequence() {
		return sequence;
	}
	
	public ReadingFrame getReadingFrame() {
		return rf;
	}
}

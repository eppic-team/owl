package owl.core.sequence;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import owl.core.sequence.alignment.PairwiseSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;

public class ProteinToCDSMatch {

	
	
	private class TranslatedSequence implements Comparable<TranslatedSequence> {
		Sequence sequence;
		PairwiseSequenceAlignment psa;
		ReadingFrame rf;
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
	}
	
	private Sequence protein;
	private Sequence cds;
	private GeneticCodeType gct;
	private Map<ReadingFrame, TranslatedSequence> translations;
	
	public ProteinToCDSMatch(Sequence protein, Sequence cds,GeneticCodeType gct){
		this.protein = protein;
		this.cds = cds;
		this.gct = gct;
		// we don't warn anymore about the nucleotide sequence not being multiple of 3
		// we do translations in all frames and alignments and if things don't match then we warn
		//if (cds.getLength()%3!=0) {
		//	System.err.println("Warning! Nucleotide sequence length is not multiple of 3 for sequence "+cds.getName());
		//}
		this.translations = get6FramesTranslations();
		this.align();
	}
	
	private Map<ReadingFrame,TranslatedSequence> get6FramesTranslations() {
		Map<ReadingFrame,TranslatedSequence> map = new HashMap<ReadingFrame, TranslatedSequence>();
		for (ReadingFrame rf:ReadingFrame.values()) {
			Sequence translated = Translator.translate(this.gct, this.cds, rf);
			translated.chopStopCodon();
			map.put(rf, new TranslatedSequence(translated,rf));
		}
		return map;
	}
	
	private void align() {
		for (ReadingFrame rf:translations.keySet()) {
			try {
				TranslatedSequence translated = translations.get(rf);
				PairwiseSequenceAlignment psa = new PairwiseSequenceAlignment(this.protein,translated.sequence);
				translated.setAln(psa);
			} catch (PairwiseSequenceAlignmentException e) {
				System.err.println("Problem aligning the mismatching sequences.");
				System.err.println(e.getMessage());
			}
		}
	}
	
	public TranslatedSequence getBestTranslation() {
		return Collections.max(translations.values());
	}
	
	public void printSummary(float cutoff) {
		TranslatedSequence best = getBestTranslation();
		//if (!best.rf.equals(ReadingFrame.ONE) || best.percentIdentity<cutoff) {
		if (best.getPercentIdentity()<cutoff) {
			System.out.printf("%s\t%s\tframe %2d\t%5.1f\t%d\n",protein.getName(),best.sequence.getSecondaryAccession(),best.rf.getNumber(),best.getPercentIdentity(),best.getMismatches());
			//best.getAln().printAlignment();
			printMatching(best);
		}
		
	}
	
	private void printMatching(TranslatedSequence translated) {
		int lineWidth = 80;
		int codonPrintWidth = 5;
		StringBuffer protLine = new StringBuffer();
		StringBuffer markLine = new StringBuffer();
		StringBuffer tranLine = new StringBuffer();
		StringBuffer codoLine = new StringBuffer();
		
		// protein sequence line
		for (int i=0;i<protein.getLength();i++) { 
			protLine.append("  "+protein.getSeq().charAt(i)+"  ");
		}
		
		// markup line
		char[] markupLine = translated.getAln().getMarkupLine();
		for (int i=0;i<protein.getLength();i++) { 
			markLine.append("  "+markupLine[i]+"  ");
		}
		
		// translated line
		for (int i=0;i<translated.sequence.getLength();i++) { 
			tranLine.append("  "+translated.sequence.getSeq().charAt(i)+"  ");
		}
		
		// codons 
		ReadingFrame rf = translated.rf;
		String scanningSeq = cds.getSeq();
		if (rf.isReverse()) {
			scanningSeq = new StringBuffer(cds.getSeq()).reverse().toString();
		}
		int startIndex = Math.abs(rf.getNumber())-1;
		for (int i=startIndex;i<cds.getLength() && i+3<=cds.getLength();i+=3) {
			codoLine.append(" "+scanningSeq.substring(i, i+3)+" ");
		}

		for (int i=0;i<protein.getLength()*codonPrintWidth;i+=lineWidth) {
			int endIndex = i+lineWidth;
			if (i+lineWidth>protLine.length()) {
				endIndex = protLine.length();
			}
			System.out.println(protLine.substring(i, endIndex));
			System.out.println(markLine.substring(i, endIndex));
			System.out.println(tranLine.substring(i, endIndex));
			System.out.println(codoLine.substring(i, endIndex));
			System.out.println();
		}
		System.out.println();


	}
	
	public ReadingFrame getBestTranslationFrame() {
		return getBestTranslation().rf;
	}
}

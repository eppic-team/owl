package owl.core.sequence;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import owl.core.sequence.alignment.PairwiseSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import owl.core.structure.AminoAcid;

/**
 * Class encapsulating an aminoacid sequence and its corresponding Uniprot-mapped CDS 
 * (nucleotide coding sequence). It contains methods to check the translation of the CDS 
 * sequence to the reference sequence and to obtain the best matching translation from
 * all the possible reading frames. 
 * 
 * @author duarte_j
 *
 */
public class ProteinToCDSMatch {

	private static final float PERFECT_MATCH_THRESHOLD = 99.99999f;
	
	
	private Sequence protein;
	private Sequence cds;
	private GeneticCodeType gct;
	private Map<ReadingFrame, TranslatedFrame> translations;
	private TranslatedFrame bestTranslation; // the cached best translation from the map above
	
	public ProteinToCDSMatch(Sequence protein, Sequence cds,GeneticCodeType gct) throws TranslationException {
		this.protein = protein;
		this.cds = cds;
		this.gct = gct;
		this.translations = get6FramesTranslations();
		this.align();
	}
	
	private Map<ReadingFrame,TranslatedFrame> get6FramesTranslations() throws TranslationException {
		Map<ReadingFrame,TranslatedFrame> map = new HashMap<ReadingFrame, TranslatedFrame>();
		for (ReadingFrame rf:ReadingFrame.values()) {
			Sequence translated = Translator.translate(this.gct, this.cds, rf);
			translated.chopStopCodon();
			map.put(rf, new TranslatedFrame(translated,rf));
		}
		return map;
	}
	
	private void align() {
		for (ReadingFrame rf:translations.keySet()) {
			try {
				TranslatedFrame translated = translations.get(rf);
				PairwiseSequenceAlignment psa = new PairwiseSequenceAlignment(this.protein,translated.getSequence());
				translated.setAln(psa);
			} catch (PairwiseSequenceAlignmentException e) {
				System.err.println("Problem aligning the mismatching sequences.");
				System.err.println(e.getMessage());
			}
		}
	}
	
	/**
	 * Gets the translation (sequence plus reading frame) of the CDS that matches with 
	 * the highest sequence identity the reference protein. 
	 * @return
	 */
	public TranslatedFrame getBestTranslation() {
		if (bestTranslation==null) {
			bestTranslation = Collections.max(translations.values()); 
		}
		return bestTranslation; 
	}
	
	/**
	 * Gets the nucleotide sequence of the best translation obtained from {@link #getBestTranslation()}
	 * by chopping off the non-translated tails. The output nucleotide sequence length will 
	 * be multiple of 3. 
	 * @return
	 */
	public String getNucleotideSeqForBestTranslation() {
		ReadingFrame rf = getBestTranslationFrame();
		StringBuffer sequenceSB = new StringBuffer(cds.getSeq());
		if (rf.isReverse()) {
			sequenceSB.reverse();
		}
		
		String seq = sequenceSB.toString();
		seq = seq.substring(Math.abs(rf.getNumber()-1));
		seq = seq.substring(0,seq.length()-seq.length()%3);
		if (seq.length()%3!=0) {
			System.err.println("Error! The nucleotide sequence is not multiple of 3! Please report a bug.");
			System.exit(1);
		}
		// chopping off the stop codon (if there is one)
		try {
			if (Translator.translate(this.gct, new Codon(seq.substring(seq.length()-3, seq.length()))).equals(AminoAcid.STP)) {
				seq = seq.substring(0,seq.length()-3);
			}
		} catch (TranslationException e) {
			System.err.println("Unexpected error while trying to chop off the STOP codon");
			System.err.println(e.getMessage());
			System.exit(1);
		}
		return seq;
	}
	
	/**
	 * Returns true if this ProteinToCDSMatch has a fully matching CDS (100% identity 
	 * to the reference protein sequence) in one of the reading frame. 
	 * @return
	 */
	public boolean hasFullMatch() {
		if (getBestTranslation().getPercentIdentity()>PERFECT_MATCH_THRESHOLD) {
			return true;
		} else {
			return false;
		}
	}
	
	public void printSummary(float cutoff) {
		TranslatedFrame best = getBestTranslation();
		//if (!best.rf.equals(ReadingFrame.ONE) || best.percentIdentity<cutoff) {
		if (best.getPercentIdentity()<cutoff) {
			System.out.printf("%s\t%s\tframe %2d\t%5.1f\t%d\n",protein.getName(),best.getSequence().getSecondaryAccession(),best.getReadingFrame().getNumber(),best.getPercentIdentity(),best.getMismatches());
			//best.getAln().printAlignment();
			printMatching(best);
		}
		
	}
	
	private void printMatching(TranslatedFrame translated) {
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
		for (int i=0;i<translated.getSequence().getLength();i++) { 
			tranLine.append("  "+translated.getSequence().getSeq().charAt(i)+"  ");
		}
		
		// codons 
		ReadingFrame rf = translated.getReadingFrame();
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
	
	/**
	 * Gets the translation frame for which the CDS translation matches with highest
	 * sequence identity to the reference protein sequence.
	 * @return
	 */
	public ReadingFrame getBestTranslationFrame() {
		return getBestTranslation().getReadingFrame();
	}
	
	public String getCDSName() {
		return this.cds.getName();
	}
	
	public Sequence getCDS() {
		return this.cds;
	}
	
	public Sequence getProteinSequence() {
		return this.protein;
	}
	
	public boolean hasStopCodonsInBestTranslation() {
		return (bestTranslation.getSequence().getSeq().contains("*"));
	}
}

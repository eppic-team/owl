package scripts;
import java.io.File;
import java.io.IOException;

import owl.core.connections.UniProtConnection;
import owl.core.sequence.Sequence;
import owl.core.sequence.alignment.PairwiseSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import owl.core.util.FileFormatException;




public class alignSeq {

	private static final float		DEFAULT_GAP_OPEN_SCORE =	10f;
	private static final float		DEFAULT_GAP_EXTEND_SCORE =	0.5f;
	private static final String		DEFAULT_MATRIX_NAME =		"BLOSUM62";

	
	/**
	 * Performs a pairwise alignment of two sequences using the NeedlemanWunsch algorithm implemented by JAligner.
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		if(args.length < 2) {
			System.out.println("Usage: alingSeq <seq1.fa or uniprot code 1> <seq2.fa or uniprot code 2> ");
			System.exit(1);
		}
		
		boolean fromFile = true;

		String uniId1 = null;
		String uniId2 = null;

		File file1 = new File(args[0]);
		File file2 = new File(args[1]);
		
		if (!file1.exists() || !file2.exists()) {
			uniId1 = args[0];
			uniId2 = args[1];
			fromFile = false;
		}

		
		Sequence seq1 = null;
		Sequence seq2 = null;
		
		if (fromFile) {
			seq1 = readSequence(file1);
			seq2 = readSequence(file2);
			
		} else {
			UniProtConnection uc = new UniProtConnection();
			seq1 = new Sequence(uniId1,uc.getEntry(uniId1).getSequence().getValue());
			seq2 = new Sequence(uniId2,uc.getEntry(uniId2).getSequence().getValue());
			
		}
		
				
		
		
		// setting the dummy file to trash all unwanted output from jaligner (done through a java Logger)
		File trashLogFile = null;
		try {
			trashLogFile = File.createTempFile("logger", ".trash");
		} catch (IOException e) {
			e.printStackTrace();
		}
		trashLogFile.deleteOnExit();
		System.setProperty("java.util.logging.config.file",trashLogFile.getAbsolutePath());
		
		// align
		PairwiseSequenceAlignment pa = null;
		try {
			pa = new PairwiseSequenceAlignment(seq1.getSeq(), seq2.getSeq(), seq1.getName(), seq2.getName(), DEFAULT_GAP_OPEN_SCORE, DEFAULT_GAP_EXTEND_SCORE, DEFAULT_MATRIX_NAME);
		} catch (PairwiseSequenceAlignmentException e) {
			System.err.println("Error creating alignment: " + e.getMessage());
			System.exit(4);
		}
		
		// print alignment
		pa.printReport();
		
	}
	
	private static Sequence readSequence (File file) {
		
		Sequence seq = new Sequence();
		try {
			seq.readFromFastaFile(file);
		} catch (FileFormatException e) {
			System.err.println("Error reading from file " + file + ": " + e.getMessage());
			System.exit(3);
		} catch (IOException e) {
			System.err.println("Error reading from file " + file + ": " + e.getMessage());
			System.exit(3);
		}
		return seq;
	}

}

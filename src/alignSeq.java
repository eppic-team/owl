import java.io.File;
import java.io.IOException;

import proteinstructure.FileFormatError;
import proteinstructure.PairwiseSequenceAlignment;
import proteinstructure.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;

import sequence.Sequence;


public class alignSeq {

	/**
	 * Performs a pairwise alignment of two sequences using the NeedlemanWunsch algorithm implemented by JAligner.
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length < 2) {
			System.out.println("Usage: alingPairwise seq1.fa seq2.fa");
			System.exit(1);
		}
		
		// read first sequence
		File fa1 = new File(args[0]);
		Sequence seq1 = new Sequence();
		try {
			seq1.readFromFastaFile(fa1);
		} catch (FileFormatError e) {
			System.err.println("Error reading from file " + fa1 + ": " + e.getMessage());
			System.exit(2);
		} catch (IOException e) {
			System.err.println("Error reading from file " + fa1 + ": " + e.getMessage());
			System.exit(2);
		}

		// read second sequence
		File fa2 = new File(args[1]);
		Sequence seq2 = new Sequence();
		try {
			seq2.readFromFastaFile(fa2);
		} catch (FileFormatError e) {
			System.err.println("Error reading from file " + fa2 + ": " + e.getMessage());
			System.exit(3);
		} catch (IOException e) {
			System.err.println("Error reading from file " + fa2 + ": " + e.getMessage());
			System.exit(3);
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
			pa = new PairwiseSequenceAlignment(seq1.getSeq(), seq2.getSeq(), seq1.getName(), seq2.getName());
		} catch (PairwiseSequenceAlignmentException e) {
			System.err.println("Error creating alignment: " + e.getMessage());
			System.exit(4);
		}
		
		// print alignment
		pa.printReport();
		
	}

}

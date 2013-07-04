package owl.tests.core.sequence;

import java.io.File;
import java.io.IOException;

import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;
import jaligner.util.SequenceParser;
import jaligner.util.SequenceParserException;
import jaligner.NeedlemanWunschGotoh;
import jaligner.Sequence;

import org.junit.Test;

import owl.core.util.FileFormatException;

public class JalignerTest {

	private static final String DATADIR = "src/owl/tests/core/sequence/data";
	
	private static final File TEST_LONG_SEQUENCE_FILE = new File(DATADIR+"/Q8WZ42.fasta");
	
	@Test
	public void testJalignerBug() throws IOException, FileFormatException, MatrixLoaderException, SequenceParserException {
			
		String name1 = "mytest";
		String name2 = "Q8WZ42";
		
		// subsequence of seq2 in the region under 32768: runs ok
		String seq1NoProblems = "PLSLGPNIEIIHEGLDYYALHIRDTLPEDT"; // works
		// subsequence of seq2 in the region above 32768: overflow bug, hangs forever (infinite loop)
		String seq1Problems = "GYYRVTATNTAGSTSCQAHLQVERLRYKKQ"; // doesn't work
		
		
		String seq2 = owl.core.sequence.Sequence.readSeqs(TEST_LONG_SEQUENCE_FILE, null).get(0).getSeq();


		Matrix matrix = MatrixLoader.load("BLOSUM50");

		System.out.println("Running sequence under 32768");
		Sequence s1 = SequenceParser.parse(seq1NoProblems);
		s1.setId(name1);
		Sequence s2 = SequenceParser.parse(seq2);
		s2.setId(name2);

		NeedlemanWunschGotoh.align(s1, s2, matrix, 10.0f, 0.5f);

		
		System.out.println("Running sequence above 32768");
		s1 = SequenceParser.parse(seq1Problems);
		s1.setId(name1);

		NeedlemanWunschGotoh.align(s1, s2, matrix, 10.0f, 0.5f);

		System.out.println("Finished successfully");
		
	}

}

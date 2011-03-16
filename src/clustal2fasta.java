import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;

import owl.core.sequence.alignment.AlignmentConstructionException;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.util.FileFormatException;

public class clustal2fasta {

	/*------------------------------ constants ------------------------------*/	
	
	/**
	 * To read a input file in clustalW format and write to a output file or stdout in fasta format 
	 * @param args
	 */
	public static void main(String[] args) {
		
		if(args.length < 1 || args.length > 2) {
			System.out.println("Usage: clustal2fasta input.aln [output.fa]");
			System.exit(1);
		}
		
		String inputfilename = args[0];
		String outputfilename = null;

		if(!new File(inputfilename).canRead()) {
			System.err.println("Error: Could not read from file " + inputfilename);
			System.exit(1);
		}
		
		PrintStream ps = System.out;
		
		if (args.length == 2){
			outputfilename = args[1];
			try {
				ps = new PrintStream(outputfilename);
			} catch (FileNotFoundException e1) {
				System.out.println("Failed to open file " + outputfilename + ": " + e1.getMessage());
				System.exit(1);
			}
		}
			
		
		//create a multiple sequence alignment from the input pir file
		MultipleSequenceAlignment ma = null;
		try {
			ma = new MultipleSequenceAlignment(inputfilename, MultipleSequenceAlignment.CLUSTALFORMAT);
			ma.writeFasta(ps, 80, true);
			if(outputfilename != null) {
				System.out.println(inputfilename + " -> " + outputfilename);
			}
		} catch (IOException e) {
			System.err.println("Error while processing " + inputfilename + ": " + e.getMessage());
		} catch (FileFormatException e) {
			System.err.println("Invalid format of " + inputfilename + ": " + e.getMessage());
		} catch (AlignmentConstructionException e) {
			System.err.println("Failed to create alignment from file " + inputfilename + ": " + e.getMessage());
		}
	}

}

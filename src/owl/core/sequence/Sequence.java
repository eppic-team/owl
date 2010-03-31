package owl.core.sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.FileFormatError;


/**
 * Presenting our sequence object(!). Provides method for reading and writing
 * sequences in Fasta format.
 */
public class Sequence {

	/*------------------------------ constants ------------------------------*/
	
	private static final String FASTAHEADER_REGEX = "^>\\s*([a-zA-Z0-9_|\\-.]+)"; // see also the REGEX in Alignment class
	
	/*--------------------------- member variables --------------------------*/
	
	private String name;
	private String seq;
	
	/*----------------------------- constructors ----------------------------*/
	
	/** 
	 * Creates a new empty sequence object
	 */
	public Sequence() {
		name = "";
		seq = "";
	}
	
	/**
	 * Creates a new sequence object with the given name and the given sequence string.
	 * @param name
	 * @param seq
	 */
	public Sequence(String name, String seq) {
		this.name = name;
		this.seq = seq;
	}

	/*-------------------------- getters and setters ------------------------*/
	
	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the sequence
	 */
	public String getSeq() {
		return seq;
	}

	/**
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * @param seq the sequence to set
	 */
	public void setSeq(String seq) {
		this.seq = seq;
	}
	
	/**
	 * Returns the length of this sequence.
	 * @return
	 */
	public int getLength() {
		return this.seq.length();
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * @return the string representation of this sequence (without header)
	 */
	public String toString() {
		return this.seq.toString();
	}
	
	/**
	 * Creates a new sequence object by parsing the given Fasta file.
	 * @param fastaFile
	 * @throws IOException
	 * @throws FileFormatError
	 */
	public void readFromFastaFile(File fastaFile) throws IOException, FileFormatError {
		BufferedReader fileIn = new BufferedReader(new FileReader(fastaFile));
		String nextLine;
		// read sequences
		while((nextLine = fileIn.readLine()) != null) {
			Pattern p = Pattern.compile(FASTAHEADER_REGEX);
			Matcher m = p.matcher(nextLine);
			if (m.find()) {
				name = m.group(1).trim();
			} else {
				seq += nextLine.trim();
			}
		}
		if (seq.contains(" ")) {
			throw new FileFormatError("The sequence in FASTA file "+fastaFile+" contains spaces.");
		}
	}
	
	/**
	 * Writes this sequence to fasta file
	 * Output line length is fixed at 80
	 * @param fastaFile
	 * @throws IOException
	 */
	public void writeToFastaFile(File fastaFile) throws IOException {
		String[] seqs = {this.getSeq()};
		String[] tags = {this.getName()};
		writeSeqs(fastaFile, seqs, tags);
	}
	
	/**
	 * Writes this sequence to the given PrintStream.
	 * Output line length is fixed at 80
	 * @param fastaFile
	 * @throws IOException
	 */
	public void writeToPrintStream(PrintStream out) throws IOException {
		String[] seqs = {this.getSeq()};
		String[] tags = {this.getName()};
		writeSeqs(out, seqs, tags);
	}

	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Writes given sequences and tags to given sequence file in FASTA format
	 * Output line length is fixed at 80
	 * @param seqFile
	 * @param seqs
	 * @param tags
	 * @throws FileNotFoundException 
	 */
	public static void writeSeqs(File seqFile, String[] seqs, String[] tags) throws FileNotFoundException {
		PrintStream Out = new PrintStream(new FileOutputStream(seqFile));
		writeSeqs(Out, seqs, tags);
		Out.close();
	} 
	
	/**
	 * Writes given sequences and tags to given PrintStream in FASTA format
	 * Output line length is fixed at 80
	 * @param seqFile
	 * @param seqs
	 * @param tags
	 * @throws FileNotFoundException 
	 */
	public static void writeSeqs(PrintStream Out, String[] seqs, String[] tags) {
		int len = 80;
		for (int seqIdx=0;seqIdx<seqs.length;seqIdx++) { 
			Out.println(">"+tags[seqIdx]);
			for(int i=0; i<seqs[seqIdx].length(); i+=len) {
				Out.println(seqs[seqIdx].substring(i, Math.min(i+len,seqs[seqIdx].length())));
			}		
		}
	} 
	
	/**
	 * Prints a ruler with column numbers in steps of 10 up to the given length.
	 * @param length
	 */
	public static void printSeqRuler(int length) {
		StringBuilder st = new StringBuilder(length);
		for (int i = 10; i <= length; i+=10) {
			st.append(String.format("%10d", i));
		}
		System.out.println(st);
		
		st = new StringBuilder(length);
		for (int i = 10; i <= length; i+=10) {
			st.append(String.format("%10s", "|"));
		}
		System.out.println(st);

	}
}
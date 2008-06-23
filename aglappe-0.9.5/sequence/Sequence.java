package sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.FastaFileFormatError;

/**
 * Presenting our sequence object(!).
 */
public class Sequence {
	
	private static final String FASTAHEADER_REGEX = "^>\\s*([a-zA-Z0-9_|\\-]+)"; // see also the REGEX in Alignment class
	
	private String name;
	private String seq;
	
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
	
	/**
	 * Creates a new sequence object by parsing the given Fasta file.
	 * @param fastaFile
	 * @throws IOException
	 * @throws FastaFileFormatError
	 */
	public void readFromFastaFile(File fastaFile) throws IOException, FastaFileFormatError {
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
			throw new FastaFileFormatError("The sequence in FASTA file "+fastaFile+" contains spaces.");
		}
	}

	/**
	 * @return the name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the seq
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
	 * @param seq the seq to set
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

	/**
	 * Writes given sequences and tags to given sequence file in FASTA format
	 * @param seqFile
	 * @param seqs
	 * @param tags
	 * @throws FileNotFoundException 
	 */
	public static void writeSeqs(File seqFile, String[] seqs, String[] tags) throws FileNotFoundException {
		PrintStream Out = new PrintStream(new FileOutputStream(seqFile));
		int len = 80;
		for (int seqIdx=0;seqIdx<seqs.length;seqIdx++) { 
			Out.println(">"+tags[seqIdx]);
			for(int i=0; i<seqs[seqIdx].length(); i+=len) {
				Out.println(seqs[seqIdx].substring(i, Math.min(i+len,seqs[seqIdx].length())));
			}		
		}
		Out.close();
	} 
}

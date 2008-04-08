package sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Presenting our sequence object(!).
 */
public class Sequence {
	
	String name;
	String seq;
	
	/** 
	 * Creates a new empty sequence object
	 */
	public Sequence() {
		name = new String("");
		seq = new String("");
	}
	
	/**
	 * Creates a new sequence object with the given name and the given sequence string.
	 * TODO: Inconsistency: Here we allow null values whereas in the constructor we do not.
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
	 */
	public void readFromFastaFile(File fastaFile) throws IOException {
			BufferedReader fileIn = new BufferedReader(new FileReader(fastaFile));
			String nextLine;
			// read sequences
			while((nextLine = fileIn.readLine()) != null) {
				Pattern p = Pattern.compile("^>(.*)$");
				Matcher m = p.matcher(nextLine);
				if (m.find()) {
					name = m.group(1).trim();
				} else {
					seq += nextLine.trim();
				}
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

}
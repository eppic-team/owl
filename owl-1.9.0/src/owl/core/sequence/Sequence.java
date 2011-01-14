package owl.core.sequence;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.AminoAcid;
import owl.core.util.FileFormatError;


/**
 * Presenting our sequence object(!). Provides method for reading and writing
 * sequences in Fasta format.
 */
public class Sequence {

	/*------------------------------ constants ------------------------------*/
	
	// see also the REGEX in Alignment class
	// TODO this is ignoring spaces, do we want to get them too and then  
	// parse other info from the whole FASTA header? (see getPrimaryAccession() and getSecondaryAccession)
	private static final Pattern FASTAHEADER_REGEX = Pattern.compile("^>\\s*([a-zA-Z0-9_|\\-.]+)");
	
	// in principle this works for emblcds and uniprot fasta headers (not tested anywhere else!) 
	// it also removes the version (the .1 suffix at the end of the identifier that sometimes is used for instance in emblcds)
	public static final Pattern DEFLINE_PRIM_ACCESSION_REGEX = Pattern.compile("^.*\\|([^.]+)(?:\\.\\d+)?\\|.*$");
	public static final Pattern DEFLINE_SEC_ACCESSION_REGEX = Pattern.compile("^.*\\|.*\\|([^. ]+)(?:\\.\\d+)?.*$");
	
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
	 * Returns the sequence's primary accession id (the one after the first pipe) extracted 
	 * from the name (FASTA header) of this Sequence if the FASTA header is from a Uniprot 
	 * or EMBL entry complying the FASTA Defline format. 
	 * See http://en.wikipedia.org/wiki/FASTA_format for Defline format 
	 * @return the accession code or null if none could be found
	 */
	public String getPrimaryAccession() {
		Matcher m = DEFLINE_PRIM_ACCESSION_REGEX.matcher(name);
		String acc = null;
		if (m.matches()) {
			acc = m.group(1);
		}
		return acc;
	}

	/**
	 * Returns the sequence's secondary accession id (the one after the second pipe) extracted 
	 * from the name (FASTA header) of this Sequence if the FASTA header is from a Uniprot 
	 * or EMBL entry complying the FASTA Defline format.
	 * See http://en.wikipedia.org/wiki/FASTA_format for Defline format 
	 * @return the accession code or null if none could be found
	 */
	public String getSecondaryAccession() {
		Matcher m = DEFLINE_SEC_ACCESSION_REGEX.matcher(name);
		String acc = null;
		if (m.matches()) {
			acc = m.group(1);
		}
		return acc;
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
			Matcher m = FASTAHEADER_REGEX.matcher(nextLine);
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
	 * Writes this sequence to given file in FASTA format
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
	 * Writes this sequence to the given PrintStream in FASTA format.
	 * Output line length is fixed at 80
	 * @param fastaFile
	 * @throws IOException
	 */
	public void writeToPrintStream(PrintStream out) throws IOException {
		String[] seqs = {this.getSeq()};
		String[] tags = {this.getName()};
		writeSeqs(out, seqs, tags);
	}

	/**
	 * Chops off the STOP symbol '*' from the end of this sequence if one is present
	 * @return true if the sequence ends with '*' and can be chopped off, false 
	 * otherwise
	 */
	public boolean chopStopCodon() {
		if (this.seq.charAt(this.seq.length()-1)==(AminoAcid.STP.getOneLetterCode())) {
			this.seq = this.seq.substring(0, this.seq.length()-1);
			return true;
		} else {
			return false;
		}
		
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
	 * Write given list of Sequence objects to given PrintStream in FASTA format
	 * Output line length is fixed at 80 
	 * @param ps
	 * @param sequences
	 */
	public static void writeSeqs(PrintStream ps, List<Sequence> sequences) {
		int len = 80;
		for (Sequence sequence:sequences) { 
			ps.println(">"+sequence.getName());
			for(int i=0; i<sequence.getSeq().length(); i+=len) {
				ps.println(sequence.getSeq().substring(i, Math.min(i+len,sequence.getSeq().length())));
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
	
	/**
	 * Reads a fasta file containing multiple sequences returning a 
	 * list of Sequence objects
	 * @param seqsFile
	 * @param fastaHeaderRegex a regex that specifies as its first capture group what
	 * will be read from the FASTA header as a sequence tag. If null then {@link #FASTAHEADER_REGEX}
	 * will be used
	 * @return
	 * @throws IOException
	 */
	public static List<Sequence> readSeqs(File seqsFile, Pattern fastaHeaderRegex) throws IOException, FileFormatError {
		if (fastaHeaderRegex==null) fastaHeaderRegex = FASTAHEADER_REGEX;
		List<Sequence> list = new ArrayList<Sequence>();
		BufferedReader br = new BufferedReader(new FileReader(seqsFile));
		String line;
		String lastTag = null;
		String currentTag = null;
		StringBuffer seq = null;
		while ((line=br.readLine())!=null){
			if (line.isEmpty()) continue;
			if (line.startsWith(">")) {
				Matcher m = fastaHeaderRegex.matcher(line);
				if (m.find()) {
					currentTag = m.group(1);
					if (lastTag!=null) {
						list.add(new Sequence(lastTag, seq.toString()));
					}
					seq = new StringBuffer();
					lastTag = currentTag;
				} else {
					throw new FileFormatError("FASTA file "+seqsFile+" does not seem to have proper FASTA headers");
				}
			} else {
				seq.append(line.trim());
			}
		}
		br.close();
		list.add(new Sequence(lastTag, seq.toString())); // adding the last sequence
		return list;
	}
}
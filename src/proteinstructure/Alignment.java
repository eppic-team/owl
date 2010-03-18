package proteinstructure;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import runners.DsspRunner;
import tools.Interval;
import tools.IntervalSet;
import tools.MySQLConnection;


/**
 * A multiple protein sequence alignment. This class represents a set of
 * protein sequences which are globally aligned and provides functions
 * to map between the original and the aligned sequences.
 * 
 * @author		Henning Stehr, Jose Duarte, Lars Petzold
 */
public class Alignment {
	
	/*------------------------------ constants ------------------------------*/	

	public static final String PIRFORMAT = "PIR";
	public static final String FASTAFORMAT = "FASTA";
	public static final String CLUSTALFORMAT = "CLUSTAL";
	public static final char GAPCHARACTER = '-';
	private static final String FASTAHEADER_REGEX = "^>\\s*([a-zA-Z0-9_|\\-.]+)";
	private static final String FASTAHEADER_CHAR = ">";
	/*--------------------------- member variables --------------------------*/		
	
	private String[] sequences;
	
	private TreeMap<String, Integer> tags2indices; 	// sequence tag to index in the sequences array (starting at 0)
	private TreeMap<Integer, String> indices2tags;	// sequence index to sequence tag

	private TreeMap<Integer,int[]> mapAlign2Seq; // map of seq index to arrays mapping alignment serials to sequence serials 
	private TreeMap<Integer,int[]> mapSeq2Align; // map of seq index to arrays mapping sequence serials to alignment serials
	
	private TreeMap<String, SecondaryStructure> secStructAnnotation;
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Creates an Alignment from a file in either FASTA, PIR or DALI format
	 * @param fileName
	 * @param format, one of {@link #PIRFORMAT}, {@link #FASTAFORMAT},  {@link #DALIFORMAT} or {@link #CLUSTALFORMAT}
	 * @throws IOException
	 * @throws FileFormatError
	 * @throws AlignmentConstructionError 
	 */
	public Alignment(String fileName, String format) throws IOException, FileFormatError, AlignmentConstructionError {
		if (format.equals(PIRFORMAT)){
			readFilePIRFormat(fileName);
		} else if (format.equals(FASTAFORMAT)){
			readFileFastaFormat(fileName);
		} else if (format.equals(CLUSTALFORMAT)) {
			readFileClustalFormat(fileName);
		} else {
			throw new IllegalArgumentException("Format "+format+" not supported by Alignment class");
		}
		
		// checking lengths, i.e. checking we read correctly from file
		checkLengths();
		// map sequence serials (starting at 1, no gaps) to alignment serials (starting at 1, possibly gaps)
		doMapping();
		
		// if indices2tags/tags2indices length don't match sequences length then there were duplicate tags in the file
		if (indices2tags.size()!=sequences.length || tags2indices.size()!=sequences.length) {
			throw new AlignmentConstructionError("There are duplicate tags in the file "+fileName);
		}

	}
		
	/**
	 * Creates a trivial alignment (i.e. without gaps) given a Map of tags to 
	 * sequences
	 * The sequences must have the same lengths. 
	 * @param sequences
	 * @throws AlignmentConstructionError if sequences lengths differ or if 
	 * size of given map is 0
	 */
	public Alignment(TreeMap<String, String> sequences) throws AlignmentConstructionError {

		if (sequences.size() == 0) {
			throw new AlignmentConstructionError("No sequences were passed for constructing the alignment.");
		}
		// check that sequences have the same length
		int length = sequences.get(sequences.firstKey()).length();
		for(String seqTag: sequences.keySet()) {
			if(sequences.get(seqTag).length() != length) {
				throw new AlignmentConstructionError("Cannot create trivial alignment. Sequence lenghts are not the same.");
			}
		}
		
		this.sequences = new String[sequences.size()];
		this.indices2tags = new TreeMap<Integer, String>();
		this.tags2indices = new TreeMap<String, Integer>();
		
		int i=0;
		for (String seqTag: sequences.keySet()) {
			this.sequences[i]=sequences.get(seqTag);
			this.indices2tags.put(i,seqTag);
			this.tags2indices.put(seqTag, i);
			i++;
		}
		doMapping();
		
	}
	
	/**
	 * Creates a trivial alignment (i.e. without gaps) given a set of tags and
	 * a set of sequences
	 * @param seqTags
	 * @param sequences
	 * @throws AlignmentConstructionError if different number of sequences and 
	 * tags given, or if size of given array is 0
	 */
	public Alignment(String[] seqTags, String[] sequences) throws AlignmentConstructionError {
		if (seqTags.length!=sequences.length) {
			throw new AlignmentConstructionError("Different number of sequences and tags given.");
		}
		if (sequences.length == 0) {
			throw new AlignmentConstructionError("No sequences were passed for constructing the alignment.");
		}
		// check that sequences have the same length
		int length = sequences[0].length();
		for(String sequence: sequences) {
			if(sequence.length() != length) {
				throw new AlignmentConstructionError("Cannot create trivial alignment. Sequence lenghts are not the same.");
			}
		}
		
		this.sequences = new String[sequences.length];
		this.indices2tags = new TreeMap<Integer, String>();
		this.tags2indices = new TreeMap<String, Integer>();
		
		for (int i=0;i<sequences.length;i++) {
			this.sequences[i]=sequences[i];
			this.indices2tags.put(i,seqTags[i]);
			this.tags2indices.put(seqTags[i], i);
		}
		doMapping();

	}
	
	/*---------------------------- private methods --------------------------*/
	
	/**
	 * Initializes the maps to map from sequence indices to alignment indices and vice versa.
	 * Both sequence and alignment indices start at 1
	 */
	private void doMapping() {
		this.mapAlign2Seq = new TreeMap<Integer, int[]>();
		this.mapSeq2Align = new TreeMap<Integer, int[]>();
				
		for (int i=0;i<sequences.length;i++){
			String seq = sequences[i];
			int[] mapAl2Seq = new int[seq.length()+1];
			int[] mapSeq2Al = new int[getSequenceNoGaps(indices2tags.get(i)).length()+1];
			int seqIndex = 1;
			for (int alignIndex=1;alignIndex<=seq.length();alignIndex++){
				if (seq.charAt(alignIndex-1)!=GAPCHARACTER) {
					mapAl2Seq[alignIndex] = seqIndex;
					mapSeq2Al[seqIndex] = alignIndex;
					seqIndex++;
				} else { // for gaps we assign a -1
					mapAl2Seq[alignIndex] = -1;
				}
			}
			mapAlign2Seq.put(i, mapAl2Seq);
			mapSeq2Align.put(i , mapSeq2Al);
		}
	}
	
	private void checkLengths() throws AlignmentConstructionError {
		if (sequences.length==0) return;
		
		int firstLength = 0;
		for (int i=0;i<sequences.length;i++) {
			if (i==0) {
				firstLength = sequences[i].length();
			} else {
				if (sequences[i].length()!=firstLength) {
					throw new AlignmentConstructionError("Error: Some sequences in alignment have different lengths.");
				}
			}
		}
	}
	
	private void readFilePIRFormat(String fileName) throws IOException, FileFormatError {
		String 	nextLine = "",
				currentSeq = "",
				currentSeqTag = "";
		boolean foundFastaHeader = false;
		int lineNum = 0;
		int seqIndex = 0;
		int nonEmptyLine = 0;

		// open file

		BufferedReader fileIn = new BufferedReader(new FileReader(fileName));

		// read file  	


		// otherwise initialize TreeMap of sequences and rewind file
		ArrayList<String> seqsAL = new ArrayList<String>();
		indices2tags = new TreeMap<Integer, String>();
		tags2indices = new TreeMap<String, Integer>();
		fileIn.reset();

		// read sequences
		while((nextLine = fileIn.readLine()) != null) {
		    ++lineNum;
			nextLine = nextLine.trim();					    // remove whitespace
			if(nextLine.length() > 0) {						// ignore empty lines
				nonEmptyLine++;
				if (nonEmptyLine==1 && !nextLine.startsWith(">")) { // quick check for PIR format
					throw new FileFormatError("First non-empty line of file "+fileName+" does not seem to be a FASTA header.");
				}
				if(nextLine.charAt(0) == '*') {				// finish last sequence
					seqsAL.add(currentSeq);
					indices2tags.put(seqIndex,currentSeqTag);
					tags2indices.put(currentSeqTag,seqIndex);
				} else {
					Pattern p = Pattern.compile(FASTAHEADER_REGEX);
					Matcher m = p.matcher(nextLine);
					if (m.find()){				// start new sequence
						currentSeq = "";						
						currentSeqTag=m.group(1);
						foundFastaHeader = true;
						seqIndex++;
					} else {
						currentSeq = currentSeq + nextLine;     // read sequence
					}
				}
			}
		} // end while

		sequences = new String[seqsAL.size()];
		seqsAL.toArray(sequences);
		
		fileIn.close();
		
		// if no fasta headers found, file format is wrong
		if(!foundFastaHeader) {
		    throw new FileFormatError("File does not conform with Pir file format (could not detect any fasta header in the file).",fileName,(long)lineNum);
		}
		
	}

	private void readFileFastaFormat(String fileName) throws IOException, FileFormatError {
		String 	nextLine = "",
				currentSeq = "",
				lastSeqTag = "";
		boolean foundFastaHeader = false;
		long lineNum = 0;
		int seqIndex = 0;
		int nonEmptyLine = 0;
		
		// open file

		BufferedReader fileIn = new BufferedReader(new FileReader(fileName));

		// read file  	

		// initialize TreeMap of sequences 
		ArrayList<String> seqsAL = new ArrayList<String>();
		tags2indices = new TreeMap<String, Integer>();
		indices2tags = new TreeMap<Integer, String>();

		// read sequences
		while((nextLine = fileIn.readLine()) != null) {
		    ++lineNum;
			nextLine = nextLine.trim();					    // remove whitespace
			if(nextLine.length() > 0) {						// ignore empty lines
				nonEmptyLine++;
				if (nonEmptyLine==1 && !nextLine.startsWith(">")) { // quick check for FASTA format
					throw new FileFormatError("First non-empty line of FASTA file "+fileName+" does not seem to be a FASTA header.");
				}
				Pattern p = Pattern.compile(FASTAHEADER_REGEX);
				Matcher m = p.matcher(nextLine);
				if (m.find()){
					if (!lastSeqTag.equals("")) {
						seqsAL.add(currentSeq);
						indices2tags.put(seqIndex, lastSeqTag);
						tags2indices.put(lastSeqTag, seqIndex);
						currentSeq = "";
						seqIndex++;
					}
					lastSeqTag=m.group(1);
					foundFastaHeader = true;
				} else {
					currentSeq += nextLine;
				}
			}
		} // end while
		// inserting last sequence
		seqsAL.add(currentSeq);
		indices2tags.put(seqIndex,lastSeqTag);
		tags2indices.put(lastSeqTag,seqIndex);

		sequences = new String[seqsAL.size()];
		seqsAL.toArray(sequences);
		
		fileIn.close();
		
		// if no fasta headers found, file format is wrong
		if(!foundFastaHeader) {
		    throw new FileFormatError("File does not conform with FASTA file format (could not find any FASTA header in the file).",fileName,lineNum);
		}
		
	}
    
	private void readFileClustalFormat(String fileName) throws IOException, FileFormatError {
		// open file
		BufferedReader fileIn = new BufferedReader(new FileReader(fileName));

		// initialize TreeMap of sequences 
		TreeMap<String,StringBuffer> tags2seqs = new TreeMap<String, StringBuffer>();
		tags2indices = new TreeMap<String, Integer>();
		indices2tags = new TreeMap<Integer, String>();

		// read sequences
		String line;
		int lineNum=0;
		int seqIdx = 0;
		Pattern p = Pattern.compile("^(\\S+)\\s+([a-zA-Z\\-]+).*"); // regex for the sequence lines
		while((line = fileIn.readLine()) != null) {
		    ++lineNum;
			if (lineNum == 1) {
				if (!line.startsWith("CLUSTAL"))
					throw new FileFormatError("File "+fileName+" does not conform with CLUSTAL format (first line does not start with CLUSTAL)");
				continue;
			}
			line = line.trim();
			if(line.length()==0) {
				continue;
			}
			
			Matcher m = p.matcher(line);
			if (m.matches()) {
				if (!tags2seqs.containsKey(m.group(1))) {
					tags2seqs.put(m.group(1),new StringBuffer(m.group(2)));
					tags2indices.put(m.group(1), seqIdx);
					indices2tags.put(seqIdx, m.group(1));
					seqIdx++;
				} else {
					tags2seqs.get(m.group(1)).append(m.group(2));
				}
			}
		} 

		
		sequences = new String[tags2seqs.size()];

		for (seqIdx=0;seqIdx<sequences.length;seqIdx++) {
			sequences[seqIdx]=tags2seqs.get(indices2tags.get(seqIdx)).toString().toUpperCase();
		}
		
		fileIn.close();
				
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * @return a deep copy of this alignment
	 */
	public Alignment copy() throws AlignmentConstructionError {
		String[] newSeqs = sequences.clone();
		String[] newTags = new String[newSeqs.length];
		for (int i = 0; i < newTags.length; i++) {
			newTags[i] = indices2tags.get(i);
		}
		return new Alignment(newTags, newSeqs);
	}
	
	/**
	 * Adds a sequence to the alignment. The sequence has to have the same length as the alignment
	 * (possibly containing gap characters) otherwise an AlignmentConstructionError is thrown.
	 * The internal mappings are regenerated for the new alignment.
	 * @param newTag the sequence name
	 * @param newSeq the new sequence
	 * @throws AlignmentConstructionError if given newSeq differs in length from this alignment or if 
	 * given newTag already exists in this alignment
	 */
	public void addSequence(String newTag, String newSeq) throws AlignmentConstructionError {
		int l = this.getAlignmentLength();
		// check length of new sequence
		if(newSeq.length() != l) {
			throw new AlignmentConstructionError("Cannot add sequence of length " + newSeq.length() + " to alignment of length " + l);
		}
		// make sure that tag does not exist yet
		if(tags2indices.containsKey(newTag)) {
			throw new AlignmentConstructionError("Cannot add sequence. Tag " + newTag + " exists in alignment.");
		}
		
		int oldNumSeqs = sequences.length;
		String[] newSequences = new String[oldNumSeqs+1];
		for (int i = 0; i < oldNumSeqs; i++) {
			newSequences[i] = this.sequences[i];
		}
		newSequences[oldNumSeqs] = newSeq;				// add to sequence array
		this.sequences = newSequences;
		
		this.indices2tags.put(oldNumSeqs, newTag);
		this.tags2indices.put(newTag, oldNumSeqs);
		
		// regenerate position mapppings
		doMapping();
		
	}
	
	/**
	 * Returns a deep copy of this alignment where the given sequence has been added.
	 * @param newTag the sequence name
	 * @param newSeq the new sequence (possibly containing gaps)
	 * @return a deep copy of this alignment where the given sequence has been added.
	 * @throws AlignmentConstructionError if given newSeq differs in length from this alignment or if 
	 * given newTag already exists in this alignment 
	 */
	public Alignment copyAndAdd(String newTag, String newSeq) throws AlignmentConstructionError {
		Alignment newAl = this.copy();
		newAl.addSequence(newTag, newSeq);
		return newAl;
	}
	
	/**
	 * Returns the gap character
	 * @return The gap character
	 */
    public static char getGapCharacter() { 
    	return GAPCHARACTER; 
    }
	
	/**
	 * Returns the sequence (with gaps) given a sequence tag
	 * @param seqTag
	 * @return
	 */
    public String getAlignedSequence(String seqTag) { 
    	return sequences[tags2indices.get(seqTag)]; 
    }
    
    /**
     * Returns the length of the alignment (including gaps) 
     * @return
     */
    public int getAlignmentLength() { 
    	return sequences[0].length(); 
    }
    
    /**
     * Returns the total number of sequences in the alignment
     * @return
     */
    public int getNumberOfSequences() { 
    	return sequences.length; 
    }
	
    /**
     * Gets the sequence tag from the sequence index
     * @param i
     * @return
     */
    public String getTagFromIndex(int i) {
    	return indices2tags.get(i);
    }
    
    /**
     * Gets the sequence index from the sequence tag
     * @param seqTag
     * @return
     */
    public int getIndexFromTag(int seqTag) {
    	return tags2indices.get(seqTag);
    }
    
    /**
     * Returns true if alignment contains the sequence identified by seqTag
     * @param seqTag
     * @return
     */
    public boolean hasTag(String seqTag){
    	return tags2indices.containsKey(seqTag);
    }
    
    /**
     * Reset the tags to the given set of tags
     * @param newTags
     * @throws IllegalArgumentException if length of array of tags is different
     *  from number of sequences in this Alignment
     */
   	public void resetTags(String[] newTags) {
   		if (newTags.length!=this.getNumberOfSequences()) {
   			throw new IllegalArgumentException("Number of tags given for resetting of tags differ from number of sequences in the alignment");
   		}
   		tags2indices = new TreeMap<String, Integer>();
   		indices2tags = new TreeMap<Integer, String>();
   		for (int i=0; i<newTags.length; i++) {
   			tags2indices.put(newTags[i],i);
   			indices2tags.put(i, newTags[i]);
   		}
    }
    
   	/**
   	 * Resets given existingTag to newTag
   	 * @param existingTag
   	 * @param newTag
   	 * @throws IllegalArgumentException if existingTag not present in this Alignment
   	 */
   	public void resetTag(String existingTag, String newTag) {
   		if (!tags2indices.containsKey(existingTag)) {
   			throw new IllegalArgumentException("Given tag "+existingTag+" doesn't exist in this Alignment");
   		} 
   		int i = tags2indices.get(existingTag);
   		indices2tags.put(i, newTag);
   		tags2indices.remove(existingTag);
   		tags2indices.put(newTag, i);
   		
   	}
   	
   	/**
   	 * Removes the sequence corresponding to the given tag.
   	 * @param tag
   	 */
   	public void removeSequence(String tag) {
   		if (!tags2indices.containsKey(tag)) {
   			throw new IllegalArgumentException("Given tag "+tag+" doesn't exist in this Alignment");
   		} 
 
   		int idxToRemove = tags2indices.get(tag);

   		String[] newSequences = new String[this.sequences.length-1];
		TreeMap<Integer,String> newIndices2tags = new TreeMap<Integer, String>();
		TreeMap<String,Integer> newTags2indices = new TreeMap<String, Integer>();
		
		int j = 0; // new indexing
		for (int i=0;i<this.sequences.length;i++) {
			if (i!=idxToRemove) {
				newSequences[j] = this.sequences[i];
				newIndices2tags.put(j,this.indices2tags.get(i));
				newTags2indices.put(this.indices2tags.get(i), j);
				j++;
			}
		}
		
		this.sequences = newSequences;
		this.indices2tags = newIndices2tags;
		this.tags2indices = newTags2indices;
		
		doMapping();

   	}
   	
    /**
     * Returns all sequence tags in a Collection<String>
     * Conserves the order of the sequences as they were added to the Alignment
     * @return
     */
    public Collection<String> getTags(){
    	return indices2tags.values();
    }
    
    /**
     * Returns sequence seqTag with no gaps
     * @param seqNumber
     * @return
     */
    public String getSequenceNoGaps(String seqTag){
    	String seq = "";
    	for (int i=0;i<getAlignmentLength();i++){
    		char letter = getAlignedSequence(seqTag).charAt(i);
    		if (letter!=GAPCHARACTER){
    			seq+=letter;
    		}
    	}
    	return seq;
    }
    
    /**
     * Given the alignment index (starting at 1, possibly gaps),
     * returns the sequence index (starting at 1, no gaps) of sequence seqTag
     * @param seqTag
     * @param alignIndex
     * @throws IndexOutOfBoundsException if 0 given as alignIndex or else if 
     * alignIndex is bigger than maximum stored index 
     * @return the sequence index, -1 if sequence is a gap at that position
     */
    public int al2seq(String seqTag, int alignIndex){
    	if (alignIndex==0) throw new IndexOutOfBoundsException("Disallowed alignment index (0) given");
    	return mapAlign2Seq.get(tags2indices.get(seqTag))[alignIndex];
    }
    
    /**
     * Given sequence index (starting at 1, no gaps) of sequence seqTag,
     * returns the alignment index (starting at 1, possibly gaps)
     * @param seqTag
     * @param seqIndex
     * @throws IndexOutOfBoundsException if 0 given as seqIndex or else if 
     * seqIndex is bigger than maximum stored index
     * @return the alignment index
     */
    public int seq2al(String seqTag, int seqIndex) {
    	if (seqIndex==0) throw new IndexOutOfBoundsException("Disallowed sequence index (0) given");
    	return mapSeq2Align.get(tags2indices.get(seqTag))[seqIndex];
    }
    
    /**
     * Gets column alignIndex of the alignment as a String
     * @param alignIndex
     * @return
     */
    public String getColumn(int alignIndex){
    	String col="";
    	for (String seq:sequences){
    		col+=seq.charAt(alignIndex-1);
    	}
    	return col;
    }
    
    /**
     * Prints alignment by columns in tab delimited format,
     * useful to import to MySQL
     */
    public void printTabDelimited(){
    	for (int alignIndex=1;alignIndex<getAlignmentLength();alignIndex++){
    		for (String seq:sequences){
    			System.out.print(seq.charAt(alignIndex-1)+"\t");
    		}
    		System.out.print(alignIndex+"\t");
    		for (int i=0; i<sequences.length;i++){
    			int seqIndex = al2seq(indices2tags.get(i), alignIndex); 
    			if (seqIndex!=-1){ // everything not gaps
    				System.out.print(seqIndex+"\t");
    			} else {  // gaps
    				System.out.print("\\N\t");
    			}
    		}
    		System.out.println();
    	}
    }
    
    /**
     * Prints the alignment in simple text format (without sequence tags) to stdout
     */
    public void printSimple() {
    	for(String sequence: sequences) {
    		System.out.println(sequence);
    	}
    }
    
    /**
     * Prints the alignment in fasta format to stdout
     */
    public void printFasta() {
   		writeFasta(System.out,80,true);
    }

	/**
	 * Writes alignment to the given output stream. The output format 
	 * conforms to the FASTA format.
	 * @param out  the PrintStream to print to
	 * @param lineLength  the maximal line length, setting this to null 
	 *  always results in 80 characters per line
	 * @param alignedSeqs  toggles the output of the aligned or ungapped 
	 *  sequences 
	 */
	public void writeFasta(PrintStream out, Integer lineLength, boolean alignedSeqs) {
		int len = 80;
		String seq = "";

		if( lineLength != null ) {
			len = lineLength;
		}

		for( String name : getTags() ) {
			seq = alignedSeqs ? getAlignedSequence(name) : getSequenceNoGaps(name);
			out.println(FASTAHEADER_CHAR + name);
			for(int i=0; i<seq.length(); i+=len) {
				out.println(seq.substring(i, Math.min(i+len,seq.length())));
			}
		}
	}
	
    /**
     * Gets list of consecutive non-gapped sub-sequences (by means of an interval set).
     * Example (X denotes any valid amino acid):
     * <p>
     * 
     * The data:<br>
     * <code>s1: XXX---XXX-X--X</code><br>
     * <code>s2: XXX---XXXX-XXX</code><br>
     * <code>s3: --XXXX--XX-XXX</code><br>
     * <p>
     * 
     * The function calls:<br>
     * <code>TreeMap m = new TreeMap();</code><br>
     * <code>m.put("s1","XXX---XXX-X--X");</code><br>
     * <code>m.put("s2","XXX---XXXX-XXX");</code><br>
     * <code>m.put("s3","--XXXX--XX-XXX");</code><br>
     * <code>Alignment ali = new Alignment(m);</code><br>
     * <code>String[] tagSet1 = new String[1];</code><br>
     * <code>String[] tagSet2 = new String[2];</code><br>
     * <code>tagSet1[0] = "s1";</code><br>
     * <code>tagSet2[0] = "s2";</code><br>
     * <code>tagSet2[1] = "s3";</code><br>
     * <code>System.out.println(ali.getMatchingBlocks("s2",tagSet1));</code><br>
     * <code>System.out.println(ali.getMatchingBlocks("s1",tagSet2));</code><br>
     * <p>
     * The output:<br>
     * <code>[0 6, 7 9]</code><br>
     * <code>[0 2, 3 5, 6 6, 7 7]</code><br>
     *
     * @param tag  tag of the sequence of which the chunk list is to be determined
     * @param projectionTags  list of tags of sequences in the alignment whose 
     *  projection along with sequence named tag is to be used as projection 
     *  from the whole alignment. Note, that invoking this function with 
     *  {@link #getTags()} set to this parameter, considers the whole alignment 
     *  matrix. 
     * @param tag
     * @param projectionTags
     * @param positions alignment columns for which we want to get the 2 matching interval sets
     * @param degOfConservation
     * @return interval set representing the sequence of non-gapped sequence 
     * chunks. 
     * @throws IndexOutOfBoundsException 
     */
    public IntervalSet getMatchingBlocks(String tag, Collection<String> projectionTags, TreeSet<Integer> positions, int degOfConservation) 
    throws IndexOutOfBoundsException {

    	/*
    	 * col        - current alignment column
    	 * prevCol    - previous alignment column
    	 * start      - start column for the next chunk to be added
    	 * foundStart - flag set whenever a start position for the next chunk 
    	 *               to be added has been encountered
    	 * c          - observed character in sequence 'tag' in column 'col'
    	 * limit      - maximal number of tolerated gap characters at a certain 
    	 *               alignment column with respect to the sequences 
    	 *               referencened in 'projectionTags'
    	 * chunks     - the list of consecutive chunks to be returned
    	 */
    	IntervalSet chunks = new IntervalSet();
    	int col;
    	int prevCol = 1;
    	int start = 1;
    	boolean foundStart = false;
    	char c = '-';
    	int limit =  Math.max(projectionTags.size() - degOfConservation,0);

    	if(positions.isEmpty()) return chunks;
    	col = positions.iterator().next();
    	
    	for(Iterator<Integer> it = positions.iterator(); it.hasNext(); ) {
    		prevCol = col;
    		col = it.next();
    		c = getAlignedSequence(tag).charAt(col-1);

    		if( c == getGapCharacter() ) {
    			if( foundStart ) {
    				// complete chunk
    				chunks.add(new Interval(al2seq(tag,start),al2seq(tag,prevCol)));
    				foundStart = false;
    			}
    		} else if ( limit >= count(projectionTags,col,getGapCharacter()) ) {
    			if( foundStart ) {
    				if( col - prevCol > 1 ) {
    					// we allow the in between positions only to consist 
    					// of gap characters. otherwise we have to complete the 
    					// current chunk as the in-between non-gap positions
    					// are not contained in 'positions'
    					if( isBlockOf(tag,prevCol,col,getGapCharacter()) ) {
    						for( String t : projectionTags) {
    							if( !isBlockOf(t,prevCol,col,getGapCharacter()) ) {
    								foundStart = false;
    								break;
    							}
    						}
    					} else {
    						foundStart = false;
    					}

    					// please note that the 'foundStart' variable is 
    					// abused in the preceding if-clause to avoid the 
    					// allocation of an additional boolean
    					if( !foundStart ) {
    						// complete chunk
    						chunks.add(new Interval(al2seq(tag,start),al2seq(tag,prevCol)));
    						foundStart = true;
    						start = col;
    					}
    				} // else: current chunk can easily be extended
    			} else {
    				foundStart = true;
    				start = col;
    			}
    		} else {
    			if( foundStart ) {
    				foundStart = false;
    				chunks.add(new Interval(al2seq(tag,start),al2seq(tag,prevCol)));
    			}
    		}
    	}

    	if( foundStart ) {
    		// complete last chunk
    		chunks.add(new Interval(al2seq(tag,start),al2seq(tag,col)));
    	}

    	return chunks;
    }

    /**
     * Extracts from the set of given alignment position those without gaps.
     * @param projectionTags  tags of the sequences to be considered 
     * @param positions  alignment positions, i.e. indices of some columns
     * @param extractInPlace  enable this flag to directly delete all nodes 
     *  pointing to "non-gapless" columns positions, set this parameter to 
     *  false to return a new node set, i.e., 'positions' remains unchanged!
     * @return a set of indices of alignment columns out of the set of 
     * considered columns ('positions'). Please note, that parameter 
     * 'extractInPlace' has an immense impact on the output generated.     
     */
    public TreeSet<Integer> getGaplessColumns(Collection<String> projectionTags, TreeSet<Integer> positions, boolean extractInPlace) {

    	// this node set will be filled and returned if the in place editing of 
    	// parameter 'positions' is disabled
    	TreeSet<Integer> output = null;
    	if( !extractInPlace ) {
    		output = new TreeSet<Integer>();
    	}

    	int col;

    	for( Iterator<Integer> it = positions.iterator(); it.hasNext(); ) {
    		col = it.next();
    		if(count(projectionTags, col, getGapCharacter()) > 0 ) {
    			// this column contains at least one gap
    			if( extractInPlace ) {
    				// remove corresponding item in 'positions'
    				it.remove();
    			}
    		} else if( !extractInPlace ) {
    			// gapless column found -> record this event in 'output' (as 
    			// 'positions' is not editable)
    			output.add(col);
    		}
    	}

    	// return the correct node set
    	if( extractInPlace ) {
    		return positions;
    	} else {
    		return output;
    	}
    }
    
    /**
     * Returns the set of gapless columns.
     * @return the set of gapless columns
     */
    public TreeSet<Integer> getGaplessColumns() {
    	TreeSet<Integer> cols = new TreeSet<Integer>();
    	Collection<String> tags = this.getTags();
    	for (int i = 1; i < this.getAlignmentLength(); i++) {
    		if(count(tags, i, getGapCharacter()) == 0) {
    			cols.add(i);
    		}
		}
    	return cols;
    }
    
    /**
     * Returns the set of gapless columns between start and end.
     * @return the set of gapless columns between start and end
     */
    public TreeSet<Integer> getGaplessColumns(int start, int end) {
    	TreeSet<Integer> cols = new TreeSet<Integer>();
    	Collection<String> tags = this.getTags();
    	for (int i = 1; i < this.getAlignmentLength(); i++) {
    		if(count(tags, i, getGapCharacter()) == 0) {
    			if(i >= start && i <= end ) {
    				cols.add(i);
    			}
    		}
		}
    	return cols;
    }
    
    /**
     * Returns the longest non-gapped region of columns in the alignments.
     * @return the interval with the longest non-gapped region
     */
    public TreeSet<Integer> getLongestNonGappedRegion() {
    	TreeSet<Integer> longestRegion = new TreeSet<Integer>();
    	int beg = 0; 
    	int bestBeg = 0;
    	int bestEnd = 0;
    	int l = getAlignmentLength();
    	Collection<String> tags = this.getTags(); 
    	
    	// find first matching column
    	int col = 1;
    	while(col <= l) {
    		if(count(tags, col, getGapCharacter()) == 0) {
    			beg = col;
    			bestBeg = col;
    			bestEnd = col;
    			break;
    		}
    		col++;
		}
    	if(col > l) return longestRegion; // no column without gaps found
    	
    	while(col < l) {
    		col++;
    		if(count(tags, col, getGapCharacter()) > 0) {
    			// end previous interval
    			if((col-1)-beg > bestEnd-bestBeg) {
    				bestBeg = beg;
    				bestEnd = col-1;
    			}
    			beg = col + 1;
    		}
    	}

    	for (int i = bestBeg; i <= bestEnd; i++) {
			longestRegion.add(i);
		}
    	return longestRegion;
    }

    /**
     * Counts the number of occurrences of the given character at the given 
     * alignment column. The sequences to be considered is limited to the 
     * given collection of alignment tags.
     * @param tags  tags of the sequences to be considered
     * @param col  
     * @param c  
     * @return
     */
    public int count(Collection<String> tags, int col, char c) throws IndexOutOfBoundsException {
    	int i=0;
    	for( String t : tags ) {
    		if( getAlignedSequence(t).charAt(col-1) == c ) {
    			++i;
    		}
    	}
    	return i;
    }

    /**
     * 
     * @param tag
     * @param begin
     * @param end
     * @param c
     * @return
     * @throws IndexOutOfBoundsException
     */
    public boolean isBlockOf( String tag, int begin, int end, char c ) throws IndexOutOfBoundsException {
    	for(int i=begin; i<end; ++i) {
    		if( getAlignedSequence(tag).charAt(i-1) != c ) {
    			return false;
    		}
    	}
    	return true;
    }
    
    /**
     * Associates a secondary structure annotation to each sequence of this Alignment 
     * whose tag is in the form pdbCode+pdbChainCode e.g. 1abcA
     * @param conn a db connection for getting the PDB data
     * @param pdbaseDb a pdbase database name
     * @param dsspExecutable
     */
    public void addSecStructAnnotation(MySQLConnection conn, String pdbaseDb, String dsspExecutable) {
    	this.secStructAnnotation = new TreeMap<String, SecondaryStructure>();
		for (String tag:this.getTags()) {
			SecondaryStructure secStruct = new SecondaryStructure("");
			Pattern p = Pattern.compile("(\\d\\w\\w\\w)(\\w)");
			Matcher m = p.matcher(tag);
			if (m.matches()) {
				PdbasePdb pdb = null;
				try {
					pdb = new PdbasePdb(m.group(1), pdbaseDb, conn);
					pdb.load(m.group(2));
					DsspRunner dsspRunner = new DsspRunner();
					pdb.setSecondaryStructure(dsspRunner.runDssp(pdb, dsspExecutable, "--", SecStrucElement.ReducedState.THREESTATE, SecStrucElement.ReducedState.THREESTATE));
					secStruct = pdb.getSecondaryStructure();
				} catch (PdbLoadError e) {
					System.err.println("Couldn't get secondary structure annotation for sequence "+tag+". Error: "+e.getMessage());
				} catch (SQLException e) {
					System.err.println("Couldn't get secondary structure annotation for sequence "+tag+". Error: "+e.getMessage());
				} catch (PdbCodeNotFoundError e) {
					System.err.println("Couldn't get secondary structure annotation for sequence "+tag+". Error: "+e.getMessage());					
				} catch (IOException e) {
					secStruct = pdb.getSecondaryStructure(); // we take author's assignment
					System.err.println("Couldn't run dssp for sequence "+tag+". Secondary structure will be the author's assignment for this sequence. Error: "+e.getMessage());
				}
			}
			secStructAnnotation.put(tag, secStruct);
		}
    }

    /**
     * Associates a secondary structure annotation to each sequence of this Alignment 
     * taking the PDB data from the given templates
     * @param templates a TemplateList containing templates with matching ids to tags in this Alignment
     * @param dsspExecutable
     */
    public void addSecStructAnnotation(TemplateList templates, String dsspExecutable) {
    	this.secStructAnnotation = new TreeMap<String, SecondaryStructure>();
    	int countTagsNotFound=0;
    	String tagsNotFound = "";
    	for (String tag:this.getTags()) {
    		SecondaryStructure secStructure = new SecondaryStructure("");
    		if (templates.contains(tag)) {
    			if (templates.getTemplate(tag).hasPdbData()) {
    				Pdb pdb =templates.getTemplate(tag).getPdb();
    				try {
    					DsspRunner dsspRunner = new DsspRunner();
    					pdb.setSecondaryStructure(dsspRunner.runDssp(pdb,dsspExecutable, "--", SecStrucElement.ReducedState.THREESTATE, SecStrucElement.ReducedState.THREESTATE));
    				} catch (IOException e) {
    					System.err.println("Couldn't run dssp for sequence "+tag+". Secondary structure will be the author's assignment for this sequence. Error: "+e.getMessage());
    				}
    				secStructure = pdb.getSecondaryStructure();
    			} 
    		} else {
    			countTagsNotFound++;
    			tagsNotFound+=tag+" ";
    		}
    		secStructAnnotation.put(tag,secStructure);
    	}
    	// we print a warning if more than 1 tag is not found in templates, the case of 1 tag (the target) 
    	// is not found is the usual 1, that's why we only warn for >1. If templates didn't contain a single one of the tags
    	// a call to writeWithSecStruct would print simply blank sec. struct. strings
		if (countTagsNotFound>1) {
			System.err.println("Warning: more than 1 tag was not found in given templates for secondary structure annotation: "+tagsNotFound);
		}

    }
    
    /**
     * Tells whether this Alignment has been assigned secondary structure annotations 
     * for its sequences
     * @return
     */
    public boolean hasSecStructAnnotation() {
    	return secStructAnnotation!=null;
    }
    
    /**
     * Writes to given PrintStream a "graphical" overview of the secondary structures 
     * from the given psipredFile and from the sequences in this alignment for which 
     * a secondary annotation exists. The secondary structure annotation must be assigned before
     * calling this method by using {@link #addSecStructAnnotation(TemplateList, String)} or
     * {@link #addSecStructAnnotation(MySQLConnection, String, String)}
     * This Alignment must contain the targetTag
     * @param Out
     * @param targetTag a tag of a sequence in this Alignment that corresponds to the 
     * sequence in the given psipredFile
     * @param psipredFile a file with a PsiPred sec. structure prediction for sequence 
     * of targetTag (horizontal format)
     * @param showSequence whether the sequences will be printed in addition to the secondary structure assignments
     * @throws IOException if sequence of target in alignment file and psipred file don't match 
     * or if targetTag not present in this alignment or if we can't read psipredFile 
     */
    public void writeWithSecStruct(PrintStream Out, String targetTag, File psipredFile, boolean showSequence) 
    throws IOException {
    	if (!this.hasSecStructAnnotation()) {
    		//TODO abusing here of IOException, we should use a more appropriate one
    		throw new IllegalArgumentException("This Alignment has no secondary structure annotation");
    	}
    	if (!this.hasTag(targetTag)) { 
    		//TODO abusing here of IOException, we should use a more appropriate one
    		throw new IOException("The alignment doesn't contain the tag "+targetTag);
    	}
    	
    	SecondaryStructure targetSecStruct = new SecondaryStructure(psipredFile);		
		if (!this.getSequenceNoGaps(targetTag).equals(targetSecStruct.getSequence())) {
    		//TODO abusing here of IOException, we should use a more appropriate one
			throw new IOException("Sequence in alignment file with tag "+targetTag+" doesn't match sequence in psipred file "+psipredFile);
		}

		this.secStructAnnotation.put(targetTag, targetSecStruct);
		
		// we get a new tags list reordered with targetTag as the first
		String[] tags = new String[this.getTags().size()];
		tags[0] = targetTag;
		int tagIdx = 1;
		for (String tag:getTags()) {
			if (!tag.equals(targetTag)) {
				tags[tagIdx] = tag;
				tagIdx++;
			}
		}

		// printing alignment positions
		Out.printf("%7s","");
		for(int i=10; i<=this.getAlignmentLength(); i+=10) {				
			Out.printf("%10d", i);
		}
		Out.println();
		
		//int len = 80;
		for( String tag : tags ) {
			SecondaryStructure secStruct = this.secStructAnnotation.get(tag);
			String seq = getAlignedSequence(tag);
			
			// printing sequence
			if(showSequence) {
				Out.printf("%7s",tag+": ");
				for(int i=0; i<this.getAlignmentLength(); i++) {
					Out.print(seq.charAt(i));
				}
				Out.println();
			}
			
			// printing psipred confidence values (only for target)
			if (tag.equals(targetTag)) {
				Out.printf("%7s",tag+": ");
				for(int i=0; i<this.getAlignmentLength(); i++) {
					if(this.al2seq(tag, i+1)==-1) {
						Out.print(" ");
					} else {
						SecStrucElement sselem = secStruct.getSecStrucElement(this.al2seq(tag, i+1));
						if (sselem!=null) {
							if (sselem.getType()==SecStrucElement.LOOP) Out.print(" ");
							else Out.printf("%1.0f",secStruct.getConfidence(this.al2seq(tag, i+1))*10);
						} else {
							Out.print(" ");
						}
					}
				}
				Out.println();
			}
			
			// printing sec. structure elements (as 'H' or 'E')
			Out.printf("%7s",tag+": ");
			for(int i=0; i<this.getAlignmentLength(); i++) {
				if(this.al2seq(tag, i+1)==-1) {
					Out.print(GAPCHARACTER);
				} else {
					SecStrucElement sselem = secStruct.getSecStrucElement(this.al2seq(tag, i+1));
					if (sselem!=null) {
						Out.print(sselem.getType()==SecStrucElement.LOOP?" ":sselem.getType());
					} else {
						Out.print(" ");
					}
				}
			}
			Out.println();
		}
		
    }
    

    /** 
     * to test the class 
     * @throws IOException 
     * @throws FileFormatError 
     * @throws AlignmentConstructionError */
    public static void main(String[] args) throws IOException, FileFormatError, AlignmentConstructionError {
    	if (args.length<1){
    		System.err.println("Must provide FASTA file name as argument");
    		System.exit(1);
    	}
    	String fileName=args[0];


    	Alignment al = new Alignment(fileName,Alignment.FASTAFORMAT);


    	// print columns
    	for (int i=1;i<=al.getAlignmentLength();i++){
    		System.out.println(al.getColumn(i));
    	}
    	// print all sequences tags and sequences
    	for (String seqTag:al.getTags()){
    		System.out.println(seqTag);
    		System.out.println(al.getAlignedSequence(seqTag));
    	}
    	// test of seq indices
    	for (int index:al.indices2tags.keySet()) {
    		System.out.println("index "+index+", tag: "+al.indices2tags.get(index));
    	}
    	// test of al2seq
    	for (int i=1;i<=al.getAlignmentLength();i++) {
    		System.out.println("alignment serial: "+i+", seq serial: "+al.al2seq(al.getTagFromIndex(0),i));
    	}
    	// test of seq2al 
    	for (int serial=1;serial<=al.getSequenceNoGaps(al.getTagFromIndex(0)).length();serial++){
    		System.out.println("seq serial: "+serial+", alignment serial: "+al.seq2al(al.getTagFromIndex(0), serial));
    	}
    	// print alignment by columns tab delimited
    	//al.printTabDelimited();
    }

}

package proteinstructure;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Package:		proteinstructure
 * Class: 		Alignment
 * Author:		Henning Stehr, Jose Duarte, Lars Petzold
 * 
 * A multiple protein sequence alignment. This class represents a set of
 * protein sequences which are globally aligned.
 * 
 */
public class Alignment {
	
	/*------------------------------ constants ------------------------------*/	
	
	private static final char GAPCHARACTER = '-';
	private static final String PIRFORMAT = "PIR";
	private static final String FASTAFORMAT = "FASTA";
	private static final String FASTAHEADER_REGEX = "^>\\s*([a-zA-Z0-9_|\\-]+)";
	private static final String FASTAHEADER_CHAR = ">";
	
	/*--------------------------- member variables --------------------------*/		
	
	private String[] sequences;
	
	private TreeMap<String, Integer> tags2indices; 	// sequence tag to index in the sequences array (starting at 0)
	private TreeMap<Integer, String> indices2tags;	// sequence index to sequence tag

	// we don't use arrays TreeMap<Integer,Integer>[] because generic arrays are not allowed by java 
	private TreeMap<Integer,TreeMap<Integer,Integer>> mapAlign2Seq; // map of seq index to maps of alignment serials to sequence serials 
	private TreeMap<Integer,TreeMap<Integer,Integer>> mapSeq2Align; // map of seq index to maps of sequence serials to alignment serials
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Creates an Alignment from a file in either FASTA or PIR format
	 * @param fileName
	 * @param format either PIR or FASTA
	 * @throws IOException
	 * @throws PirFileFormatError
	 * @throws FastaFileFormatError
	 * @throws AlignmentConstructionError 
	 */
	public Alignment(String fileName, String format) throws IOException, PirFileFormatError, FastaFileFormatError, AlignmentConstructionError {
		if (format.equals(PIRFORMAT)){
			readFilePIRFormat(fileName);
		} else if (format.equals(FASTAFORMAT)){
			readFileFastaFormat(fileName);
		} 
		
		// checking lengths, i.e. checking we read correctly from file
		checkLengths();
		// map sequence serials (starting at 1, no gaps) to alignment serials (starting at 0, possibly gaps)
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
	
	private void doMapping() {
		this.mapAlign2Seq = new TreeMap<Integer, TreeMap<Integer,Integer>>();
		this.mapSeq2Align = new TreeMap<Integer, TreeMap<Integer,Integer>>();
		for (int i=0;i<sequences.length;i++){
			this.mapAlign2Seq.put(i, new TreeMap<Integer,Integer>());
			this.mapSeq2Align.put(i, new TreeMap<Integer,Integer>());
		}
		for (int i=0;i<sequences.length;i++){
			String seq = sequences[i];
			int seqIndex = 1;
			for (int alignIndex=0;alignIndex<seq.length();alignIndex++){
				if (seq.charAt(alignIndex)!=GAPCHARACTER) {
					mapAlign2Seq.get(i).put(alignIndex,seqIndex);
					seqIndex++;
				} else { // for gaps we assing a -1
					mapAlign2Seq.get(i).put(alignIndex,-1);
				}
			}
		}
		for (int i=0;i<sequences.length;i++){
			for (int alignIndex:mapAlign2Seq.get(i).keySet()){
				int seqIndex = mapAlign2Seq.get(i).get(alignIndex);
				if (seqIndex!=-1){
					mapSeq2Align.get(i).put(seqIndex,alignIndex);
				}
			}
		}
		
	}
	
	private void checkLengths() {
		if (sequences.length==0) return;
		
		int firstLength = 0;
		for (int i=0;i<sequences.length;i++) {
			if (i==0) {
				firstLength = sequences[i].length();
			} else {
				if (sequences[i].length()!=firstLength) {
					System.err.println("Warning: Some sequences in alignment have different lengths.");
				}
			}
		}
	}
	
	private void readFilePIRFormat(String fileName) throws IOException, PirFileFormatError {
		String 	nextLine = "",
				currentSeq = "",
				currentSeqTag = "";
		boolean foundFastaHeader = false;
		int lineNum = 0;
		int seqIndex = 0;

		// open file

		BufferedReader fileIn = new BufferedReader(new FileReader(fileName));

		// read file  	


		// otherwise initalize TreeMap of sequences and rewind file
		ArrayList<String> seqsAL = new ArrayList<String>();
		indices2tags = new TreeMap<Integer, String>();
		tags2indices = new TreeMap<String, Integer>();
		fileIn.reset();

		// read sequences
		while((nextLine = fileIn.readLine()) != null) {
		    	++lineNum;
			nextLine = nextLine.trim();					    // remove whitespace
			if(nextLine.length() > 0) {						// ignore empty lines
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
		    throw new PirFileFormatError("File does not conform with Pir file format (could not detect any fasta header in the file).",fileName,(long)lineNum);
		}
		
	}

	private void readFileFastaFormat(String fileName) throws IOException, FastaFileFormatError {
		String 	nextLine = "",
				currentSeq = "",
				lastSeqTag = "";
		boolean foundFastaHeader = false;
		long lineNum = 0;
		int seqIndex = 0;
		
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
		    throw new FastaFileFormatError("File does not conform with Fasta file format (could not find any fasta header in the file).",fileName,lineNum);
		}
		
	}

	
	/*---------------------------- public methods ---------------------------*/
	
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
     * @throws IllegalArgumentException if length of array of tags is different
     *  from number of sequences in this Alignment
     * @param newTags
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
     * Given the alignment index (starting at 0, with gaps),
     * returns the sequence index (starting at 1, no gaps) of sequence seqTag
     * Returns -1 if sequence at that position is a gap
     * @param seqTag
     * @param alignIndex
     * @return
     */
    public int al2seq(String seqTag, int alignIndex){
    	return mapAlign2Seq.get(tags2indices.get(seqTag)).get(alignIndex);
    }
    
    /**
     * Given sequence index (starting at 1, no gaps) of sequence seqTag,
     * returns the alignment index (starting at 0, with gaps)
     * @param seqTag
     * @param seqIndex
     * @return the alignment index or -1 on error
     */
    public int seq2al(String seqTag, int seqIndex) {
    	int ret = -1;
    	TreeMap<Integer,Integer> map = mapSeq2Align.get(tags2indices.get(seqTag));
    	if(map.containsKey(seqIndex)) {
    		ret = map.get(seqIndex);
    	} else {
    		System.err.println("Sever error in seq2al: Index " + seqIndex + " out of bounds for sequence " + seqTag + " (max index=" + map.lastKey() + ")");
    	}
    	return ret;
    }
    
    /**
     * Gets column alignIndex of the alignment as a String
     * @param alignIndex
     * @return
     */
    public String getColumn(int alignIndex){
    	String col="";
    	for (String seq:sequences){
    		col+=seq.charAt(alignIndex);
    	}
    	return col;
    }
    
    /**
     * Prints alignment by columns in tab delimited format,
     * useful to import to MySQL
     */
    public void printTabDelimited(){
    	for (int alignIndex=0;alignIndex<getAlignmentLength();alignIndex++){
    		for (String seq:sequences){
    			System.out.print(seq.charAt(alignIndex)+"\t");
    		}
    		System.out.print((alignIndex+1)+"\t");
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
    	for(String seqTag:getTags()) {
    		System.out.println(FASTAHEADER_CHAR + seqTag);
    		System.out.println(getAlignedSequence(seqTag));
    	}
    }

	/**
	 * Writes alignment to the given output stream. The output format 
	 * conforms to the FASTA format.
	 * @param out  the output stream to be printed to
	 * @param lineLength  the maximal line length, setting this to null 
	 *  always results in 80 characters per line
	 * @param alignedSeqs  toggles the output of the aligned or ungapped 
	 *  sequences 
	 * @throws IOException
	 */
	public void writeFasta(OutputStream out, Integer lineLength, boolean alignedSeqs) throws IOException {
		int len = 80;
		String seq = "";

		if( lineLength != null ) {
			len = lineLength;
		}

		for( String name : getTags() ) {
			seq = alignedSeqs ? getAlignedSequence(name) : getSequenceNoGaps(name);
			out.write('>');
			out.write(name.getBytes());
			out.write(System.getProperty("line.separator").getBytes());
			for(int i=0; i<seq.length(); i+=len) {
				out.write(seq.substring(i, Math.min(i+len,seq.length())).getBytes());
				out.write(System.getProperty("line.separator").getBytes());
			}
		}
	}
    
    /**
     * Gets list of consecutive non-gapped sub-sequences (by means of an interval set).
     * Example (X-denote any valid atomic sequence component):
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
     * <code>System.out.println(ali.getConsecutiveChunks("s2",tagSet1));</code><br>
     * <code>System.out.println(ali.getConsecutiveChunks("s1",tagSet2));</code><br>
     * <p>
     * The output:<br>
     * <code>[0 6, 7 9]</code><br>
     * <code>[0 2, 3 5, 6 6, 7 7]</code><br>
     *
     * @param tag  tag of the sequence of which the chunk list is be determined
     * @param projectionTags  list of tags of sequences in the alignment whose 
     *  projection along with sequence named tag is to be used as projection 
     *  from the whole alignment. Note, that invoking this function with 
     *  {@link #getTags()} set to this parameter, considers the whole alignment 
     *  matrix. 
     * 
     * @return edge set representing the sequence of non-gapped sequence chunks, 
     *  where for each edge in the set Edge.i denotes the starting position and 
     *  Edge.j the ending position of a chunk, i.e., for a separated atomic 
     *  sequence component (e.g., between two gaps) this notion yields i == j. 
     *  The indexing respects the sequence indexing for this class, i.e., index 
     *  1 corresponds to the first position in the sequence.
     *  
     * @throws IndexOutOfBoundsException 
     */
    public IntervalSet getMatchingBlocks(String tag, Collection<String> projectionTags, int begin, int end, int degOfConservation) 
    throws IndexOutOfBoundsException {

    	if( end > getAlignmentLength() ) {
    		throw new IndexOutOfBoundsException("end position exceeds alignment length");
    	}

    	/*
    	 * col        - current alignment column
    	 * start      - start column for the next chunk to be added
    	 * foundStart - flag set whenever a start position for the next chunk 
    	 *               to be added has been encountered
    	 * c          - observed character in sequence 'tag' in column 'col'
    	 * limit      - maximal number of tolerated gap characters at a certain 
    	 *               alignment column with respect to the sequences 
    	 *               referencened in 'projectionTags'
    	 * chunks     - the list of consecutive chunks to be returned
    	 */
    	int col = begin;
    	int start = 0;
    	char c = '-';
    	boolean foundStart = false;
    	int limit =  Math.max(projectionTags.size() - degOfConservation,0);
    	IntervalSet chunks = new IntervalSet();

    	while( col<end ) {
    		c = getAlignedSequence(tag).charAt(col);

    		if( c == getGapCharacter() ) {
    			if( foundStart ) {
    				foundStart = false;
    				chunks.add(new Interval(al2seq(tag,start),al2seq(tag,Math.max(start, col-1))));
    			}
    		} else {
    			if( limit >= count(projectionTags,col,getGapCharacter()) ) {
    				if( !foundStart ) {
    					foundStart = true;
    					start = col;
    				}
    			} else {
    				if( foundStart ) {
    					foundStart = false;
    					chunks.add(new Interval(al2seq(tag,start),al2seq(tag,col-1)));
    				}
    			}
    		}
    		++col;
    	}

    	if( foundStart ) {
    		chunks.add(new Interval(al2seq(tag,start),al2seq(tag,Math.max(start,col-1))));
    	}

    	return chunks;
    }

    /**
     * 
     * @param tag
     * @param projectionTags
     * @param positions
     * @param degOfConservation
     * @return The indexing respects the sequence indexing for this class, i.e., index 1 corresponds to the first position in the sequence.
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
    	int col = positions.iterator().next();
    	int prevCol = 0;
    	int start = 0;
    	boolean foundStart = false;
    	char c = '-';
    	int limit =  Math.max(projectionTags.size() - degOfConservation,0);
    	IntervalSet chunks = new IntervalSet();

    	for(Iterator<Integer> it = positions.iterator(); it.hasNext(); ) {
    		prevCol = col;
    		col = it.next();
    		c = getAlignedSequence(tag).charAt(col);

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
     * @param positions  alignment positions, i.e., indices of some columns
     * @param extractInPlace  enable this flag to directly delete all nodes 
     *  pointing to "non-gapless" columns positions, set this parameter to 
     *  false to return a new node set, i.e., 'positions' remains unchanged!
     * @return a node set corresponding to indices of alignment columns out of 
     *  the set of considered columns ('positions'). Please note, that 
     *  parameter 'extractInPlace' has an immense impact on the output 
     *  generating.     
     */
    public TreeSet<Integer> getGaplessColumns(Collection<String> projectionTags, TreeSet<Integer> positions, boolean extractInPlace) {

    	// this node set will be filled and returned if the in place editing of 
    	// parameter 'positions' is disabled
    	TreeSet<Integer> output = null;
    	if( !extractInPlace ) {
    		output = new TreeSet<Integer>();
    	}

    	int numProjTags = projectionTags.size();
    	int col = 0;

    	for( Iterator<Integer> it = positions.iterator(); it.hasNext(); ) {
    		col = it.next();
    		if( numProjTags != count(projectionTags, col, getGapCharacter()) ) {
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
    		if( getAlignedSequence(t).charAt(col) == c ) {
    			++i;
    		}
    	}
    	return i;
    }

    public boolean isBlockOf( String tag, int begin, int end, char c ) throws IndexOutOfBoundsException {
    	for(int i=begin; i<end; ++i) {
    		if( getAlignedSequence(tag).charAt(i) != c ) {
    			return false;
    		}
    	}
    	return true;
    }

    /** to test the class 
     * @throws IOException
     * @throws FastaFileFormatError 
     * @throws PirFileFormatError 
     * @throws AlignmentConstructionError */
    public static void main(String[] args) throws IOException, PirFileFormatError, FastaFileFormatError, AlignmentConstructionError {
    	if (args.length<1){
    		System.err.println("Must provide FASTA file name as argument");
    		System.exit(1);
    	}
    	String fileName=args[0];


    	Alignment al = new Alignment(fileName,"FASTA");


    	// print columns
    	//for (int i=0;i<al.getSequenceLength();i++){
    	//	System.out.println(al.getColumn(i));
    	//}
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
    	for (int i=0;i<al.getAlignmentLength();i++) {
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

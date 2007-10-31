package proteinstructure;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Package:		proteinstructure
 * Class: 		Alignment
 * Author:		Henning Stehr, Jose Duarte
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
	
	private TreeMap<String,String> sequences;

	private TreeMap<String,TreeMap<Integer,Integer>> mapAlign2Seq; // key is sequence tag, the values are maps of alignment serials to sequence serials 
	private TreeMap<String,TreeMap<Integer,Integer>> mapSeq2Align; // key is sequence tag, the values are maps of sequence serials to alignment serials
	
	/*----------------------------- constructors ----------------------------*/
	
	public Alignment() {}
	
	/**
	 * Reads an alignment from a file in either FASTA or PIR format
	 */
	public Alignment(String fileName, String format) throws IOException, FileNotFoundException {
		if (format.equals(PIRFORMAT)){
			readFilePIRFormat(fileName);
		} else if (format.equals(FASTAFORMAT)){
			readFileFastaFormat(fileName);
		} 
		
		// checking lenghts, i.e. checking we read correctly from file
		checkLengths();
		// map sequence serials (starting at 1, no gaps) to alignment serials (starting at 0, possibly gaps)
		doMapping();		
	}
		
	/**
	 * Creates a trivial alignment (i.e. without gaps) for the given sequences. The sequences have to have the same lengths. 
	 * @param sequence
	 * @param numberOfCopies
	 */
	public Alignment(TreeMap<String, String> sequences) {
		
		// check that sequences have the same length
		int length = sequences.get(sequences.firstKey()).length();
		for(String seqTag: sequences.keySet()) {
			if(sequences.get(seqTag).length() != length) {
				System.err.println("Can not create trivial alignment. Sequence lenghts are not the same.");
				// TODO: throw exception
			}
		}
		
		this.sequences = sequences;
		doMapping();
		
	}

	/*---------------------------- private methods --------------------------*/
	
	private void doMapping() {
		this.mapAlign2Seq = new TreeMap<String, TreeMap<Integer,Integer>>();
		this.mapSeq2Align = new TreeMap<String, TreeMap<Integer,Integer>>();
		for (String seqTag:sequences.keySet()){
			this.mapAlign2Seq.put(seqTag, new TreeMap<Integer,Integer>());
			this.mapSeq2Align.put(seqTag, new TreeMap<Integer,Integer>());
		}
		for (String seqTag:sequences.keySet()){
			String seq = sequences.get(seqTag);
			int seqIndex = 1;
			for (int alignIndex=0;alignIndex<seq.length();alignIndex++){
				if (seq.charAt(alignIndex)!=GAPCHARACTER) {
					mapAlign2Seq.get(seqTag).put(alignIndex,seqIndex);
					seqIndex++;
				} else { // for gaps we assing a -1
					mapAlign2Seq.get(seqTag).put(alignIndex,-1);
				}
			}
		}
		for (String seqTag:sequences.keySet()){
			for (int alignIndex:mapAlign2Seq.get(seqTag).keySet()){
				int seqIndex = mapAlign2Seq.get(seqTag).get(alignIndex);
				if (seqIndex!=-1){
					mapSeq2Align.get(seqTag).put(seqIndex,alignIndex);
				}
			}
		}
		
	}
	
	private void checkLengths() {
		// checking all read sequences have same length
		TreeMap<String,Integer> seqLengths = new TreeMap<String, Integer>(); 
		for (String seqTag:sequences.keySet()){
			seqLengths.put(seqTag,sequences.get(seqTag).length());
		}
		int lastLength = 0;
		for (String seqTag:seqLengths.keySet()){
			if (lastLength!=0 && seqLengths.get(seqTag)!=lastLength) {
				System.err.println("Warning: Some sequences in alignment have different lengths.");
			}
			lastLength = seqLengths.get(seqTag);
		}		
	}
	
	private void readFilePIRFormat(String fileName) throws IOException, FileNotFoundException{
		String 	nextLine = "",
				currentSeq = "",
				currentSeqTag = "";
		boolean foundFastaHeader = false;

		// open file

		BufferedReader fileIn = new BufferedReader(new FileReader(fileName));

		// read file  	


		// otherwise initalize TreeMap of sequences and rewind file
		sequences = new TreeMap<String,String>();
		fileIn.reset();

		// read sequences
		while((nextLine = fileIn.readLine()) != null) {
			nextLine = nextLine.trim();					    // remove whitespace
			if(nextLine.length() > 0) {						// ignore empty lines
				if(nextLine.charAt(0) == '*') {				// finish last sequence
					sequences.put(currentSeqTag, currentSeq);
				} else {
					Pattern p = Pattern.compile(FASTAHEADER_REGEX);
					Matcher m = p.matcher(nextLine);
					if (m.find()){				// start new sequence
						currentSeq = "";						
						currentSeqTag=m.group(1);
						foundFastaHeader = true;
					} else {
						currentSeq = currentSeq + nextLine;     // read sequence
					}
				}
			}
		} // end while

		
		fileIn.close();
		
		// if no fasta headers found, file format is wrong
		if(!foundFastaHeader) {
			System.err.println("Error: " + fileName + " is not a "+PIRFORMAT+" file.");
			System.exit(1);	
			//TODO throw exception
		}
		
	}

	private void readFileFastaFormat(String fileName) throws IOException, FileNotFoundException{
		String 	nextLine = "",
				currentSeq = "",
				lastSeqTag = "";
		boolean foundFastaHeader = false;

		// open file

		BufferedReader fileIn = new BufferedReader(new FileReader(fileName));

		// read file  	

		// initalize TreeMap of sequences 
		sequences = new TreeMap<String,String>();

		// read sequences
		while((nextLine = fileIn.readLine()) != null) {
			nextLine = nextLine.trim();					    // remove whitespace
			if(nextLine.length() > 0) {						// ignore empty lines
				Pattern p = Pattern.compile(FASTAHEADER_REGEX);
				Matcher m = p.matcher(nextLine);
				if (m.find()){
					if (!lastSeqTag.equals("")) {
						sequences.put(lastSeqTag,currentSeq);
						currentSeq = "";
					}
					lastSeqTag=m.group(1);
					foundFastaHeader = true;
				} else {
					currentSeq += nextLine;
				}
			}
		} // end while
		// inserting last sequence
		sequences.put(lastSeqTag,currentSeq);

		fileIn.close();
		
		// if no fasta headers found, file format is wrong
		if(!foundFastaHeader) {
			System.err.println("Error: " + fileName + " is not a "+FASTAFORMAT+" file.");
			System.exit(1);
			//TODO throw exception
		}
		
	}

	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Returns the gap character
	 * @return The gap character
	 */
    public char getGapCharacter() { return GAPCHARACTER; }
    
    /**
     * Returns a TreeMap (keys: tags, values: sequences) with all sequences
     * @return
     */
	public TreeMap<String,String> getSequences() { return sequences; }
	
	/**
	 * Returns the sequence (with gaps) given a sequence tag
	 * @param seqTag
	 * @return
	 */
    public String getAlignedSequence(String seqTag) { return sequences.get(seqTag); }
    
    /**
     * Returns the length of the alignment (including gaps) 
     * @return
     */
    public int getAlignmentLength() { return sequences.get(sequences.firstKey()).length(); }
    
    /**
     * Returns the total number of sequences in the alignment
     * @return
     */
    public int getNumberOfSequences() { return sequences.size(); }
	
    /**
     * Returns true if alignment contains the sequence identified by seqTag
     * @param seqTag
     * @return
     */
    public boolean hasTag(String seqTag){
    	return sequences.containsKey(seqTag);
    }
    
    /**
     * Returns all sequence tags in a Set<String>
     * @return
     */
    public Set<String> getTags(){
    	return sequences.keySet();
    }
    
    /**
     * Returns sequence seqTag with no gaps
     * @param seqNumber
     * @return
     */
    public String getSequenceNoGaps(String seqTag){
    	String seq = "";
    	for (int i=0;i<getAlignmentLength();i++){
    		char letter = sequences.get(seqTag).charAt(i);
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
    	return mapAlign2Seq.get(seqTag).get(alignIndex);
    }
    
    /**
     * Given sequence index (starting at 1, no gaps) of sequence seqTag,
     * returns the alignment index (starting at 0, with gaps)
     * @param seqTag
     * @param seqIndex
     * @return
     */
    public int seq2al(String seqTag, int seqIndex) {
    	return mapSeq2Align.get(seqTag).get(seqIndex);
    }
    
    /**
     * Gets column alignIndex of the alignment as a String
     * @param alignIndex
     * @return
     */
    public String getColumn(int alignIndex){
    	String col="";
    	for (String seq:sequences.values()){
    		col+=seq.charAt(alignIndex);
    	}
    	return col;
    }
    
    /**
     * Writes alignment by columns in tab delimited format,
     * useful to import to MySQL
     *
     */
    public void writeTabDelimited(){
    	for (int alignIndex=0;alignIndex<getAlignmentLength();alignIndex++){
    		for (String seq:sequences.values()){
    			System.out.print(seq.charAt(alignIndex)+"\t");
    		}
    		System.out.print((alignIndex+1)+"\t");
    		for (String seqTag:sequences.keySet()){
    			int seqIndex = al2seq(seqTag, alignIndex); 
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
    	for(String seqTag:sequences.keySet()) {
    		System.out.println(getAlignedSequence(seqTag));
    	}
    }
    
    /**
     * Prints the alignment in fasta format to stdout
     */
    public void printFasta() {
    	for(String seqTag:sequences.keySet()) {
    		System.out.println(FASTAHEADER_CHAR + seqTag);
    		System.out.println(getAlignedSequence(seqTag));
    	}
    }

    /**
     * Gets list of consecutive non-gapped sub-sequences (by means of an edge set).
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
    public EdgeSet getMatchingBlocks(String tag, Collection<String> projectionTags, int begin, int end, int degOfConservation) 
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
	EdgeSet chunks = new EdgeSet();
		
	while( col<end ) {
	    c = getAlignedSequence(tag).charAt(col);
	    
	    if( c == getGapCharacter() ) {
		if( foundStart ) {
		    foundStart = false;
		    chunks.add(new Edge(al2seq(tag,start),al2seq(tag,Math.max(start, col-1))));
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
			chunks.add(new Edge(al2seq(tag,start),al2seq(tag,col-1)));
		    }
		}
	    }
	    ++col;
	}
	
	if( foundStart ) {
	    chunks.add(new Edge(al2seq(tag,start),al2seq(tag,Math.max(start,col-1))));
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
    public EdgeSet getMatchingBlocks(String tag, Collection<String> projectionTags, NodeSet positions, int degOfConservation) 
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
	int col = positions.iterator().next().num;
	int prevCol = 0;
	int start = 0;
	boolean foundStart = false;
	char c = '-';
	int limit =  Math.max(projectionTags.size() - degOfConservation,0);
	EdgeSet chunks = new EdgeSet();
	
	for(Iterator<Node> it = positions.iterator(); it.hasNext(); ) {
	    prevCol = col;
	    col = it.next().num;
	    c = getAlignedSequence(tag).charAt(col);
	    
	    if( c == getGapCharacter() ) {
		if( foundStart ) {
		    // complete chunk
		    chunks.add(new Edge(al2seq(tag,start),al2seq(tag,prevCol)));
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
			    chunks.add(new Edge(al2seq(tag,start),al2seq(tag,prevCol)));
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
		    chunks.add(new Edge(al2seq(tag,start),al2seq(tag,prevCol)));
		}
	    }
	}
	
	if( foundStart ) {
	    // complete last chunk
	    chunks.add(new Edge(al2seq(tag,start),al2seq(tag,col)));
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
    public NodeSet getGaplessColumns(Collection<String> projectionTags, NodeSet positions, boolean extractInPlace) {

	// this node set will be filled and returned if the in place editing of 
	// parameter 'positions' is disabled
	NodeSet output = null;
	if( !extractInPlace ) {
	    output = new NodeSet();
	}
	
	int numProjTags = projectionTags.size();
	int col = 0;
	
	for( Iterator<Node> it = positions.iterator(); it.hasNext(); ) {
	    col = it.next().num;
	    if( numProjTags != count(projectionTags, col, getGapCharacter()) ) {
		// this column contains at least one gap
		if( extractInPlace ) {
		    // remove corresponding item in 'positions'
		    it.remove();
		}
	    } else if( !extractInPlace ) {
		// gapless column found -> record this event in 'output' (as 
		// 'positions' is not editable)
		output.add(new Node(col));
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
    
    /** to test the class */
	public static void main(String[] args) throws FileNotFoundException, IOException {
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
		for (String seqTag:al.getSequences().keySet()){
			System.out.println(seqTag);
			System.out.println(al.getAlignedSequence(seqTag));
		}
		// test of al2seq
		//for (int i=0;i<al.getSequenceLength();i++) {
		//	System.out.println("alignment serial: "+i+", seq serial: "+al.al2seq(al.sequences.firstKey(),i));
		//}
		// test of seq2al 
		//for (int serial=1;serial<=al.getSequenceNoGaps(al.sequences.firstKey()).length();serial++){
		//	System.out.println("seq serial: "+serial+", alignment serial: "+al.seq2al(al.sequences.firstKey(), serial));
		//}
		// print alignment by columns tab delimited
		//al.writeTabDelimited();
	}

}

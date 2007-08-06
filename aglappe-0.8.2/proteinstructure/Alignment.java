package proteinstructure;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
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
	
	/*--------------------------- member variables --------------------------*/		
	
	private TreeMap<String,String> sequences;

	private TreeMap<String,TreeMap<Integer,Integer>> mapAlign2Seq; // key is sequence tag, the values are maps of alignment serials to sequence serials 
	private TreeMap<String,TreeMap<Integer,Integer>> mapSeq2Align; // key is sequence tag, the values are maps of sequence serials to alignment serials
	
	/*----------------------------- constructors ----------------------------*/
	
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
	
    public char getGapCharacter() { return GAPCHARACTER; }
	public TreeMap<String,String> getSequences() { return sequences; }
    public String getSequence(String seqTag) { return sequences.get(seqTag); }
    public int getSequenceLength() { return sequences.get(sequences.firstKey()).length(); }    
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
    	for (int i=0;i<getSequenceLength();i++){
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
     * Given sequence index of sequence seqTag,
     * returns the alignment index
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
    	for (int alignIndex=0;alignIndex<getSequenceLength();alignIndex++){
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
			System.out.println(al.getSequence(seqTag));
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

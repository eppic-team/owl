package proteinstructure;

import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
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
	
	/*--------------------------- member variables --------------------------*/		
	
	String[] sequences;
	
	String[] sequenceTags;

	int[][] map; // first index is the sequence number, second is the alignment serials (indices)
	
	/*----------------------------- constructors ----------------------------*/
	
	public Alignment(String fileName, String format) throws IOException, FileNotFoundException {
		if (format.equals(PIRFORMAT)){
			readFilePIRFormat(fileName);
		} else if (format.equals(FASTAFORMAT)){
			readFileFastaFormat(fileName);
		} 
		
		// map sequence serials (starting at 1, no gaps) to alignment serials (starting at 0, possibly gaps)
		doMapping();		
	}

	/*---------------------------- private methods --------------------------*/
	
	private void doMapping() {
		this.map = new int[getNumberOfSequences()][getSequenceLength()]; 
		for (int i=0;i<getNumberOfSequences();i++){
			String seq = getSequence(i);
			int serial = 1;
			for (int j=0;j<seq.length();j++){
				if (seq.charAt(j)!=GAPCHARACTER) {
					map[i][j]=serial;
					serial++;
				}
			}
		}
		
	}
	
	private void readFilePIRFormat(String fileName) throws IOException, FileNotFoundException{
		String 	nextLine = "",
				currentSeq = "";    	
		int 	numSeqs = 0,
				currentSeqNum = 0;

		// open file

		BufferedReader fileIn = new BufferedReader(new FileReader(fileName));

		// read file  	


		// count number of sequences in file, assume PIR format (sequences end with '\n*')
		fileIn.mark(100000);
		while((nextLine = fileIn.readLine()) != null) {
			if(nextLine.length() > 0 && nextLine.charAt(0) == '>') {
				numSeqs++;
			}
		}			

		// if no fasta headers found, file format is wrong
		if(numSeqs == 0) {
			System.err.println("Error: " + fileName + " is not a "+PIRFORMAT+" file.");
			System.exit(1);	
			//TODO throw exception
		}

		// otherwise initalize array of sequences and rewind file
		sequences = new String[numSeqs];
		sequenceTags = new String[numSeqs];
		fileIn.reset();

		// read sequences
		while((nextLine = fileIn.readLine()) != null) {
			nextLine = nextLine.trim();					    // remove whitespace
			if(nextLine.length() > 0) {						// ignore empty lines
				if(nextLine.charAt(0) == '*') {				// finish last sequence
					sequences[currentSeqNum] = currentSeq;
					currentSeqNum++;
				} else
					if(nextLine.charAt(0) == '>') {				// start new sequence
						currentSeq = "";
						Pattern p = Pattern.compile(">\\s*([a-zA-Z0-9_-]+)");
						Matcher m = p.matcher(nextLine);
						if (m.find()){
							sequenceTags[currentSeqNum]=m.group(1);
						}
					} else {
						currentSeq = currentSeq + nextLine;     // read sequence
					}
			}
		} // end while

		// verify number and lengths of sequences
		if(currentSeqNum != numSeqs) {
			System.err.println("Error: Could only read " + currentSeqNum + " out of " 
					+ numSeqs + " sequences in file " + fileName);
			System.exit(1);
			//TODO throw exception
		}

		for(currentSeqNum = 1; currentSeqNum < numSeqs; currentSeqNum++) {
			if(sequences[currentSeqNum].length() != sequences[currentSeqNum-1].length()) {
				System.err.println("Warning: Sequences " + (currentSeqNum-1) + " and " + currentSeqNum + " in alignment have different lengths.");
			}
		}

		// clean up
		fileIn.close();	
	}

	private void readFileFastaFormat(String fileName) throws IOException, FileNotFoundException{
		String 	nextLine = "",
				currentSeq = "";    	
		int 	numSeqs = 0,
				currentSeqNum = 0;

		// open file

		BufferedReader fileIn = new BufferedReader(new FileReader(fileName));

		// read file  	


		// count number of sequences in file
		fileIn.mark(100000);
		while((nextLine = fileIn.readLine()) != null) {
			if(nextLine.length() > 0 && nextLine.charAt(0) == '>') {
				numSeqs++;
			}
		}			

		// if no fasta headers found, file format is wrong
		if(numSeqs == 0) {
			System.err.println("Error: " + fileName + " is not a "+FASTAFORMAT+" file.");
			System.exit(1);
			//TODO throw exception
		}

		// otherwise initalize array of sequences and rewind file
		sequences = new String[numSeqs];
		sequenceTags = new String[numSeqs];
		fileIn.reset();

		// read sequences
		while((nextLine = fileIn.readLine()) != null) {
			nextLine = nextLine.trim();					    // remove whitespace
			if(nextLine.length() > 0) {						// ignore empty lines
				if(nextLine.charAt(0) == '>') {
					Pattern p = Pattern.compile(">\\s*([a-zA-Z0-9_-]+)");
					Matcher m = p.matcher(nextLine);
					if (m.find()){
						sequenceTags[currentSeqNum]=m.group(1);
					}
					if (currentSeqNum!=0) {
						sequences[currentSeqNum-1]=currentSeq;
						currentSeq = "";
					}
					currentSeqNum++;
				} else {
					currentSeq += nextLine;
				}
			}
		} // end while
		// inserting last sequence
		sequences[currentSeqNum-1]=currentSeq;

		// verify number and lengths of sequences
		if(currentSeqNum != numSeqs) {
			System.err.println("Error: Could only read " + currentSeqNum + " out of " 
					+ numSeqs + " sequences in file " + fileName);
			System.exit(1);
			//TODO throw exception
		}

		for(currentSeqNum = 1; currentSeqNum < numSeqs; currentSeqNum++) {
			if(sequences[currentSeqNum].length() != sequences[currentSeqNum-1].length()) {
				System.err.println("Warning: Sequences " + (currentSeqNum-1) + " and " + currentSeqNum + " in alignment have different lengths.");
			}
		}

		// clean up
		fileIn.close();	
	}

	
	/*---------------------------- public methods ---------------------------*/
	
    public char getGapCharacter() { return GAPCHARACTER; }
	public String[] getSequences() { return sequences; }
    public String getSequence(int i) { return sequences[i]; }
    public int getSequenceLength() { return sequences[0].length(); }
    //public int getSequenceLength(int i) { return sequences[i].length(); }    
    public int getNumberOfSequences() { return sequences.length; }
	
    /**
     * Get sequence's tag as present in fasta file
     * @param i
     * @return
     */
    public String getSequenceTag(int i){
    	return sequenceTags[i];
    }
    
    /**
     * Given the alignment serial (starting at 0, with gaps),
     * returns the sequence serial (starting at 1, no gaps) of sequence seqNumber
     * @param seqNumber
     * @param alignmentSerial
     * @return
     */
    public int getSerialInSequence(int seqNumber, int alignmentSerial){
    	return map[seqNumber][alignmentSerial];
    }
    
    /**
     * Gets column i of the alignment as a String
     * @param i
     * @return
     */
    public String getColumn(int i){
    	String col="";
    	for (String seq:sequences){
    		col+=seq.charAt(i);
    	}
    	return col;
    }
    
    /**
     * Writes alignment by columns in tab delimited format,
     * useful to import to MySQL
     *
     */
    public void writeTabDelimited(){
    	for (int i=0;i<getSequenceLength();i++){
    		for (String seq:sequences){
    			System.out.print(seq.charAt(i)+"\t");
    		}
    		System.out.print((i+1)+"\t");
    		for (int seqNumber=0;seqNumber<getNumberOfSequences();seqNumber++){
    			int serial = getSerialInSequence(seqNumber, i); 
    			if (serial!=0){
    				System.out.print(serial+"\t");
    			} else {
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
//		for (int i=0;i<al.getSequenceLength();i++){
//			System.out.println(al.getColumn(i));
//		}
		
		for (int i=0;i<al.getNumberOfSequences();i++){
			System.out.println(al.getSequenceTag(i));
			System.out.println(al.getSequence(i));
		}
//		for (int i=0;i<al.getSequenceLength();i++) {
//			System.out.println("alignment serial: "+i+", seq serial: "+al.getSerialInSequence(0,i));
//		}
		//al.writeTabDelimited();
	}

}
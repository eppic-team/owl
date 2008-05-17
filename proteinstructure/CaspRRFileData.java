package proteinstructure;

import java.io.*;
import java.util.*;
import java.util.regex.*;

/**
 * This class represents the data contained in a Casp RR (contact prediction) file.
 * It also provides methods for reading from and writing to RR files and
 * for verifying that a given file is of the desired type.
 * @author stehr
 * @date 2007-12-19
 */
public class CaspRRFileData {
	
	/*------------------------------ constants ------------------------------*/
	
	// default values (as used in Casp)
	public static final double DEFAULT_MIN_DIST = 0.0;
	public static final double DEFAULT_MAX_DIST = 8.0;
	public static final double DEFAULT_WEIGHT = 1.0;
	public static final String DEFAULT_CONT_TYPE = "Cb";
	
	// output line templates (printf style, for writing to file)
	public static final String OUT_HEADER = "PFRMAT RR\n";
	public static final String OUT_TARGET = "TARGET T%04d\n";
	public static final String OUT_MODEL = "MODEL %d\n";
	public static final String OUT_AUTHOR = "AUTHOR %s\n";
	public static final String OUT_METHOD = "METHOD %s\n";
	public static final String OUT_SEQ = "%s\n";
	public static final String OUT_CONTACT = "%3d %3d %2.0f %2.0f %3.2f";
	public static final String OUT_END = "END\n";
	private static final String OUT_FILENAME = "T%04dRR%03d_%d";
	
	public static final int MAX_CHARS_PER_SEQ_LINE = 50;
	
	// regular expressions (for parsing from file)
	private static final String REGEX_END = "END";
	private static final String REGEX_CONTACT = "(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d(.\\d+)?(e[+-]\\d\\d)?)";
	private static final String REGEX_SEQUENCE = "([A-Y ]+)";
	private static final String REGEX_MODEL = "MODEL (.*)";
	private static final String REGEX_METHOD = "METHOD (.*)";
	private static final String REGEX_REMARK = "REMARK (.*)";
	private static final String REGEX_AUTHOR = "AUTHOR (.*)";
	private static final String REGEX_TARGET = "TARGET\\s+T(\\d\\d\\d\\d)";
	private static final String REGEX_HEADER = "PFRMAT\\s+RR";
	private static final String REGEX_FILENAME = "T\\d{4}RR(\\d{3})_\\d";

	/*--------------------------- type definitions --------------------------*/
	
	/**
	 * A contact in a casp RR file data object
	 */
	public class RRContact {
		int i;
		int j;
		double minDist;
		double maxDist;
		double weight;
		
		public RRContact(int i, int j) {
			this(i, j, DEFAULT_MIN_DIST, DEFAULT_MAX_DIST, DEFAULT_WEIGHT);
		}
		
		public RRContact(int i, int j, double minDist, double maxDist, double weight) {
			this.i = i;
			this.j = j;
			this.minDist = minDist;
			this.maxDist = maxDist;
			this.weight = weight;
		}
		
		public String toString() {
			return String.format(OUT_CONTACT, i, j, minDist, maxDist, weight);
		}
		
		public boolean equals(RRContact c2) {
			if(this.i != c2.i) return false;
			if(this.j != c2.j) return false;
			if(this.minDist != c2.minDist) return false;
			if(this.maxDist != c2.maxDist) return false;
			if(Math.abs(this.weight - c2.weight) > 0.005) {
				return false;
			}
			return true;
		}
	}
	
	/*--------------------------- member variables --------------------------*/
	int targetNum;
	int groupNum;
	int modelNum;
	String sequence;
	String author;
	String method;
	ArrayList<RRContact> contacts;
	String contactType;
	
	/*----------------------------- constructors ----------------------------*/
	public CaspRRFileData() {
		targetNum = 0;
		groupNum = 0;
		modelNum = 1;
		sequence = "";
		author = null;
		method = null;
		contacts = new ArrayList<RRContact>();
		this.contactType = DEFAULT_CONT_TYPE;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	public boolean equals(CaspRRFileData d2) {
		if(this.targetNum != d2.targetNum) return false;
		if(this.groupNum != d2.groupNum) return false;
		if(this.modelNum != d2.modelNum) return false;
		if(!this.sequence.equals(d2.sequence)) return false;
		if(this.contacts.size() != d2.contacts.size()) return false;
		if(!this.contactType.equals(d2.contactType)) return false;
		for(int i = 0; i < contacts.size(); i++) {
			if(!contacts.get(i).equals(d2.contacts.get(i))) return false;
		}
		return true;
	}
	
	// getters
	
	public int getTargetNum() {
		return this.targetNum;
	}
	
	public int getGroupNum() {
		return this.groupNum;
	}
	
	public int getModelNum() {
		return this.modelNum;
	}
	
	public String getSequence() {
		return this.sequence;
	}
	
	public ArrayList<RRContact> getContacts() {
		return this.contacts;
	}	
	
	public String getContactType() {
		return this.contactType;
	}
	
	/**
	 * Returns the name of the file this object should be saved to according to the default Casp naming convention.
	 */
	public String getDefaultFileName() {
		String name = String.format(OUT_FILENAME, targetNum, groupNum, modelNum);
		return name;
	}
	
	// setters
	
	public void setTargetNum(int targetNum) {
		this.targetNum = targetNum;
	}
	
	public void setGroupNum(int groupNum) {
		this.groupNum = groupNum;
	}
	
	public void setModelNum(int modelNum) {
		this.modelNum = modelNum;
	}
	
	public void setSequence(String seq) {
		this.sequence = seq;
	}
	
	/**
	 * Sets the author field needed for Casp submissions, set to null to supress output.
	 * @param author the author string
	 */
	public void setAuthor(String author) {
		this.author = author;
	}
	
	/**
	 * Sets the method field needed for Casp submissions, set to null to supress output.
	 * @param method the method string
	 */
	public void setMethod(String method) {
		this.method = method;
	}
	
	public void addContact(RRContact contact) {
		this.contacts.add(contact);
	}
	
	// input/output
	
	/**
	 * Prints a summary of this object to screen.
	 */
	public void printSummary() {
		System.out.println("Target:   " + targetNum);
		System.out.println("Group:    " + groupNum);
		System.out.println("Model:    " + modelNum);
		System.out.println("Contacts: " + contacts.size());
		System.out.println("Sequence: " + sequence);
	}
	
	/**
	 * Fills this object with data read from a Casp RR contact prediction file.
	 */
	public void readFromFile(File inFile) throws FileNotFoundException, IOException{
		String line;
		String match;
		
		// Try to take group from file name
		match = getFirstGroup(inFile.getName(), REGEX_FILENAME);
		if(match != null) {
			try {
				groupNum = Integer.parseInt(match);
			} catch (NumberFormatException e) {
				// could not convert, ignore
			}
		}
		
		BufferedReader in = new BufferedReader(new FileReader(inFile));
		line = nextLine(in);
		if(!lineMatches(line, REGEX_HEADER)) throw new IOException("PFRMAT RR expected, found " + line);	// obligatory header
		line = nextLine(in);
		if(!lineMatches(line, REGEX_TARGET)) throw new IOException("TARGET expected, found " + line);		// obligatory target number
		match = getFirstGroup(line, REGEX_TARGET);
		this.targetNum = Integer.parseInt(match);
		line = nextLine(in);
		if(lineMatches(line, REGEX_AUTHOR)) line = nextLine(in);											// optional author line (ignore)
		while(lineMatches(line, REGEX_REMARK)) line = nextLine(in);											// optional remark lines (ignore)
		while(lineMatches(line, REGEX_METHOD)) line = nextLine(in);											// optional method lines (ignore)
		if(!lineMatches(line, REGEX_MODEL)) throw new IOException("MODEL expected, found " + line);			// obligatory model line
		match = getFirstGroup(line, REGEX_MODEL);
		this.modelNum = Integer.parseInt(match.trim());		
		line = nextLine(in);
		if(!lineMatches(line, REGEX_SEQUENCE)) throw new IOException("Sequence expected, found " + line);	// obligatory sequence line
		match = getFirstGroup(line, REGEX_SEQUENCE);
		this.sequence = match;
		line = nextLine(in);
		while(lineMatches(line, REGEX_SEQUENCE)) {
			match = getFirstGroup(line, REGEX_SEQUENCE);
			this.sequence += match;
			line = nextLine(in);										// optional further sequence lines
		}
		if(!lineMatches(line, REGEX_CONTACT)) throw new IOException("Contact line expected, found " + line);
		RRContact c = parseContactLine(line);
		contacts.add(c);
		line = nextLine(in);
		while(lineMatches(line, REGEX_CONTACT)) {
			c = parseContactLine(line);
			contacts.add(c);
			line = nextLine(in);
		}
		if(!lineMatches(line, REGEX_END)) throw new IOException("END expected, found " + line);				// obligatory end tag
	}
	
	/**
	 * Writes this data object to an output stream (i.e. a FileOutputStream).
	 */
	public void writeToStream(OutputStream out) {
		PrintStream bout = new PrintStream(out);
		StringBuffer seq = new StringBuffer(this.sequence);
		bout.format(OUT_HEADER);
		bout.format(OUT_TARGET, targetNum);
		if(author != null) bout.format(OUT_AUTHOR, author);
		if(method != null) bout.format(OUT_METHOD, method);
		bout.format(OUT_MODEL, modelNum);
		while(seq.length() > MAX_CHARS_PER_SEQ_LINE) {
			bout.format(OUT_SEQ, seq.substring(0,MAX_CHARS_PER_SEQ_LINE));
			seq.delete(0, MAX_CHARS_PER_SEQ_LINE);
		}
		bout.format(OUT_SEQ, seq);
		for(RRContact c:contacts) {
			bout.println(c.toString());
		}
		bout.format(OUT_END);
		
	}
	
	/**
	 * Writes this object to a file.
	 */
	public void writeToFile(File outFile) throws FileNotFoundException {
		PrintStream out = new PrintStream(new FileOutputStream(outFile));
		this.writeToStream(out);
	}
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Parses the file and return true iff the file format is correct.
	 * @param inFile the file to read from
	 */
	public static boolean isValidCaspRRFile(File inFile) throws FileNotFoundException {		
		String line;
		BufferedReader in = new BufferedReader(new FileReader(inFile));
		line = nextLine(in);
		if(!lineMatches(line, REGEX_HEADER)) return false;	// obligatory header
		line = nextLine(in);
		if(!lineMatches(line, REGEX_TARGET)) return false;		// obligatory target number
		line = nextLine(in);
		if(lineMatches(line, REGEX_AUTHOR)) line = nextLine(in);											// optional author line
		while(lineMatches(line, REGEX_REMARK)) line = nextLine(in);											// optional remark lines
		while(lineMatches(line, REGEX_METHOD)) line = nextLine(in);											// optional method lines
		if(!lineMatches(line, REGEX_MODEL)) return false;			// obligatory model line
		line = nextLine(in);
		if(!lineMatches(line, REGEX_SEQUENCE)) return false;	// obligatory sequence line
		line = nextLine(in);
		while(lineMatches(line, REGEX_SEQUENCE)) line = nextLine(in);										// optional further sequence lines
		if(!lineMatches(line, REGEX_CONTACT)) return false;
		line = nextLine(in);
		while(lineMatches(line, REGEX_CONTACT)) line = nextLine(in);
		if(!lineMatches(line, REGEX_END)) return false;				// obligatory end tag
		return true;
	}
		
	/*---------------------------- private methods --------------------------*/
	
	/**
	 * Returns true if the current line in the stream matches the regular expression, false if it does not match
	 * or something went wrong.
	 * @param line
	 * @param regex
	 * @return true if regex matches the next line, false otherwise
	 */
	private static boolean lineMatches(String line, String regex) {
			Pattern pattern = 
	            Pattern.compile(regex);
			Matcher matcher = 
	            pattern.matcher(line.trim());
			return matcher.matches();
	}
	
	/** Reads the next line from the input stream
	 * @param in
	 * @return the line read or null if something went wrong
	 */
	private static String nextLine(BufferedReader in) {
		try {
			return in.readLine();
		} catch (IOException e) {
			return null;
		}	
	}
	
	/**
	 * Returns the first match group of regex in line or null if no match is found.
	 * @param line
	 * @param regex
	 * @return
	 */
	private String getFirstGroup(String line, String regex) {
		Pattern pattern = 
            Pattern.compile(regex);
		Matcher matcher = 
            pattern.matcher(line.trim());
		if(matcher.matches() && matcher.groupCount() > 0) {
			return matcher.group(1);
		} else {
			return null;
		}
	}
	
	/**
	 * Parses a contact line and returns the resulting RRContact object or null if the line does not match.
	 * @param line
	 * @return
	 */
	private RRContact parseContactLine(String line) {
		Pattern pattern = 
            Pattern.compile(REGEX_CONTACT);
		Matcher matcher = 
            pattern.matcher(line.trim());
		if(matcher.matches() && matcher.groupCount() >= 5) {
			int i = Integer.parseInt(matcher.group(1));
			int j = Integer.parseInt(matcher.group(2));
			double minDist = Double.parseDouble(matcher.group(3));
			double maxDist = Double.parseDouble(matcher.group(4));
			double weight = Double.parseDouble(matcher.group(5));
			RRContact c = new RRContact(i,j,minDist, maxDist, weight);
			return c;
		} else {
			return null;
		}
	}
	
	/*--------------------------------- main --------------------------------*/
		
	private static void testReadWrite(File inFile) {
		
		String inFileName = inFile.getName();
		String tmpFileName = inFile.getName();
		String tmpDir = System.getProperty("java.io.tmpdir");
		File tmpFile = new File(tmpDir, tmpFileName);
		tmpFile.deleteOnExit();
		
		// read from inFile
		CaspRRFileData d = new CaspRRFileData();
		//System.out.println(inFile.getName());
		try {
			d.readFromFile(inFile);
		} catch(IOException e) {
			System.err.println("Parsing error in " + inFile + ":" + e.getMessage());
		}
		//d.printSummary();
		
		// check file name
		if(!d.getDefaultFileName().equals(inFileName)) {
			System.err.println("File name " + inFileName + " does not match naming convention. Expected: " + d.getDefaultFileName());
		}
		
		// write to tmpFile
		try {
			d.writeToFile(tmpFile);
		} catch (FileNotFoundException e) {
			System.err.println("Error processing " + inFile.getName() + ": " + e.getMessage());
		}
		
		// reload from tmpFile
		CaspRRFileData d2 = new CaspRRFileData();
		try {
			d2.readFromFile(tmpFile);
		} catch(IOException e) {
			System.err.println("Parsing error in " + tmpFile + ":" + e.getMessage());
		}
		
		// compare with original data
		if(!d.equals(d2)) {
			System.err.println("Error while testing " + inFileName + ". Files do not match.");
			//d.printSummary();
			//d2.printSummary();
		}
	}
	
	public static void main(String args[]) throws FileNotFoundException {
		
		// counters
//		int valid = 0;		// number of valid prediction files
//		int invalid = 0;	// number of invalid prediction files
//		TreeMap<Integer, Integer> modelCount = new TreeMap<Integer, Integer>();   // predictions per target
//		TreeMap<Integer, Integer> targetCount = new TreeMap<Integer, Integer>();  // predicted targets per group
		
		// process command line
		if(args.length != 1) {
			System.out.println("Usage: CaspRRFileParser directory");
			System.exit(1);
		}
		String dirName = args[0];
		File inDir = new File(dirName);
			
		// check input parameters
		if(!inDir.canRead()) {
			System.err.println("Can not read from directory " + dirName);
			System.exit(1);
		}
		if(!inDir.isDirectory()) {
			System.err.println(dirName + " is not a directory");
			System.exit(1);			
		}
		File[] inFiles = inDir.listFiles();
		if(inFiles.length == 0) {
			System.err.println(dirName + " does not contain any files");
			System.exit(1);				
		}
				
		//testReadWrite(new File("/project/StruPPi/CASP7/predictions/RR/RR283-386/T0364RR763_1"));
		
		for(File inFile:inFiles) {
			testReadWrite(inFile);
		}
		System.out.println("done.");
		
		// process input files
//		for(File inFile:inFiles) {
//			int target = Integer.parseInt(inFile.getName().substring(2,5));
//			if(modelCount.containsKey(target)) {
//				int count = modelCount.get(target);
//				modelCount.put(target, count+1);
//			} else {
//				modelCount.put(target, 1);
//			}
//			try {
//				parseFile(inFile);
//				valid++;
//			} catch(IOException e) {
//				System.out.println(inFile.getName() + " is wrong: " + e.getMessage());
//				invalid++;				
//			}
//		}
//		
//		// print results
//		System.out.println();
//		System.out.println("Summary:");
//		System.out.println("Total files: " + inFiles.length);
//		System.out.println("Valid files: " + valid);
//		System.out.println("Wrong files: " + invalid);
//		for(int target:modelCount.keySet()) {
//			int count = modelCount.get(target);
//			if(target % 10 == 0) System.out.println();
//			System.out.print(target + ":" + count + ", ");
//		}
//		System.out.println();
	}
	
}


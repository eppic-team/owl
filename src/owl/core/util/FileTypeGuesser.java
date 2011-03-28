package owl.core.util;

import java.io.*;
import java.util.regex.*;

/**
 * Provides static methods for guessing the type of several protein related file formats.
 * Currently supported:
 * - PDB files
 * - PDB atom line files
 * - Casp TS (3D prediction) files
 * - Casp RR (contact prediction) files
 * - mmCIF files
 * - CMView contact map files
 * @author stehr
 * @date 2007-12-20
 */
public class FileTypeGuesser {

	/*------------------------------ constants ------------------------------*/
	//TODO: Create enum
	public static final int UNKNOWN_FILE = 0;		// any unknown file type
	public static final int PDB_FILE = 1;			// pdb file with header
	public static final int RAW_PDB_FILE = 2;		// file with pdb atom lines
	public static final int CASP_TS_FILE = 3;		// casp 3d prediction file
	public static final int CASP_RR_FILE = 4;		// casp contact prediction file
	public static final int OWL_CM_FILE = 5;			// contact map file
	public static final int CIF_FILE = 6;			// mmCIF file from PDB	
	
	// names of the file types as above
	private static final String[] FILE_TYPE_NAMES  = {
		"Unknown file",
		"PDB file",
		"File with PDB atom lines",
		"Casp 3D prediction file",
		"Casp contact prediction file",
		"Contact map file",
		"PDB mmCIF file"
	};
	
	// signatures for the files in order as above
	private static final String[] FILE_SIGNATURES  = {
		"(?:HEADER|SEQRES|CRYST1).*",
		"ATOM\\s+.*",
		"PFRMAT\\s+TS.*",
		"PFRMAT\\s+RR.*",
		"#(?:AGLAPPE|CMVIEW|OWL) GRAPH FILE.*",
		"data_\\d\\w\\w\\w"
	};
	
	/*---------------------------- public methods ---------------------------*/
	public static int guessFileType(File file) throws FileNotFoundException, IOException {
		BufferedReader in = new BufferedReader(new FileReader(file));
		String firstLine = in.readLine();
		if(firstLine.trim().length() == 0) firstLine = in.readLine();
		for(int i=1; i <= FILE_SIGNATURES.length; i++) {
			if(lineMatches(firstLine, FILE_SIGNATURES[i-1])) {
				in.close();
				return i;
			}
		}
		in.close();
		return UNKNOWN_FILE;
	}
	
	/**
	 * Returns a string with the name of the given file type constant or null if no name is defined.
	 * The file types are defined as public constants in this class.
	 * @param fileType the file type number
	 * @return The name of the file type or null if no name is defined for the given constant
	 */
	public static String getFileTypeName(int fileType) {
		if(fileType > FILE_TYPE_NAMES.length) return null;
		return FILE_TYPE_NAMES[fileType];
	}
	
	/*---------------------------- private methods --------------------------*/
	/**
	 * Returns true if the line matches the regular expression, false if it does not match
	 * or something went wrong.
	 * @param line the line of text
	 * @param regex the regular expression to match
	 * @return true if regex matches the line, false otherwise
	 */
	private static boolean lineMatches(String line, String regex) {
			Pattern pattern = 
	            Pattern.compile(regex);
			Matcher matcher = 
	            pattern.matcher(line.trim());
			return matcher.matches();
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		if(args.length < 1) {
			System.out.println("Usage: FileTypeGuesser.java file [file2 ...]");
			System.exit(1);
		}
		
		int fileType;
		String fileTypeName;
	
		for(int f=0; f < args.length; f++) {
			String fileName = args[f];
			try {
				fileType = guessFileType(new File(fileName));
				fileTypeName = getFileTypeName(fileType);
				if(fileType > 0) {
					System.out.println("I believe " + fileName + " is a " + fileTypeName);
				} else {
					System.out.println("I really can't figure out the file type of " + fileName);
				}
			} catch(IOException e) {
				System.err.println("Error: " + e.getMessage());
			}
		}

	}

}


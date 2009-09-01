package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import tools.MySQLConnection;


public abstract class Scorer {
	
	// to compile a scoring matrix from the given pdb ids we won't 
	// take anything that has observed length below this value
	private static final int MIN_VALID_CHAIN_LENGTH = 30;

	protected static final String DB = "pdbase";
	
	/**
	 * Enum to identify the different scoring methods
	 * @author duarte
	 *
	 */
	public enum ScoringMethod {
		RESTYPE("residue type", "restype"),
		ATOMTYPE("atom type", "atomtype"),
		RESCOUNT("residue count", "rescount"),
		ATOMCOUNT("atom count", "atomcount");
		
		private String id;
		private String description;
		private ScoringMethod(String description, String id) {
			this.description = description;
			this.id = id;
		}
		public String getDescription() {
			return description;
		}
		public String getId() {
			return id;
		}
		public static ScoringMethod getByDescription(String description) {
			for (ScoringMethod scMethod:values()) {
				if (scMethod.getDescription().equals(description)) {
					return scMethod;
				}
			}
			return null;
		}
	}

	protected File listFile;							// the file containing the training set list (list of pdbCode+pdbChainCode)
	protected ArrayList<String> structureIds;			// ids (pdbCode+pdbChainCode) of structures used for scoring matrix (only the valide ones from the listFile)
	protected int totalStructures;						// the number of valid structures used for scoring matrix (i.e. size of structureIds above). 
														// We use a separate variable to store the value for the case we read scoring matrix from file
	
	protected String ct;								// the contact type, only used in the case of residue-based scoring
	protected double cutoff;							// the distance cutoff
	protected int minSeqSep;							// the minimum sequence separation for which contacts are considered 
														// when compiling the scoring matrix
	
	protected MySQLConnection conn;

	
	/**
	 * Assigns a score to the given Pdb structure based on a scoring matrix read from file
	 * through the appropriate constructor. 
	 * Only contacts with a minimum of minSeqSep sequence separation will be considered.
	 * @param pdb the structure to be scored
	 * @param minSeqSep the minimum sequence separation for which contacts will be considered
	 * @return the score
	 */
	public abstract double scoreIt(Pdb pdb, int minSeqSep);

	/**
	 * Returns the list file, i.e. the files with the list of pdbCode+pdbChainCode used for 
	 * compiling the scoring matrix
	 * @return
	 */
	public File getListFile() {
		return listFile;
	}
	
	/**
	 * Returns the number of structures used for compiling the scoring matrix (only the
	 * ones that were actually used, i.e. that passed the quality checks. See {@link #isValidPdb(Pdb)})
	 * @return
	 */
	public int sizeOfTrainingSet() {
		return totalStructures;
	}
	
	/**
	 * Gets the contact type used for defining contacts (only makes sense for residue-based scoring).
	 * @return the contact type used for defining contacts or null if scoring method is atom-based
	 */
	public String getContactType() {
		return ct;
	}
	
	/**
	 * Gets the cutoff used for defining contacts
	 * @return
	 */
	public double getCutoff() {
		return cutoff;
	}
	
	/**
	 * Gets the minimum sequence separation used for filtering contacts when counting to compile
	 * the scoring matrix.
	 * @return
	 */
	public int getMinSeqSep() {
		return minSeqSep;
	}
	
	public ScoringMethod getScoringMethod() {
		if (this instanceof ResTypeScorer) {
			return ScoringMethod.RESTYPE;
		} else if (this instanceof AtomTypeScorer) {
			return ScoringMethod.ATOMTYPE;
		} else if (this instanceof ResCountScorer) {
			return ScoringMethod.RESCOUNT;
		} else if (this instanceof AtomCountScorer) {
			return ScoringMethod.ATOMCOUNT;
		}
		return null;
	}
	
	/**
	 * Returns a Scorer object of the appropriate type given a scoring matrix file.
	 * It first detects the type of scoring matrix by reading its header and based on 
	 * that returns the Scorer object of the appropriate subclass.
	 * @param scMatFile
	 * @return
	 * @throws IOException
	 * @throws FileFormatError if no header found in file or type of scoring matrix file
	 * not recognised
	 */
	public static Scorer readScoreMatFromFile(File scMatFile) throws IOException, FileFormatError {
		BufferedReader br = new BufferedReader(new FileReader(scMatFile));
		String line;
		while ((line=br.readLine())!=null) {
			Pattern p = Pattern.compile("^# SCORE METHOD: (.*)$");
			Matcher m = p.matcher(line);
			if (m.matches()) {
				String method = m.group(1);
				if (method.equals(ScoringMethod.ATOMCOUNT.getDescription())) {
					return new AtomCountScorer(scMatFile);
				} else if (method.equals(ScoringMethod.RESCOUNT.getDescription())) {
					return new ResCountScorer(scMatFile);
				} else if (method.equals(ScoringMethod.ATOMTYPE.getDescription())) {
					return new AtomTypeScorer(scMatFile);
				} else if (method.equals(ScoringMethod.RESTYPE.getDescription())) {
					return new ResTypeScorer(scMatFile);
				} else {
					throw new FileFormatError("Header #SCORE METHOD doesn't contain one of the recognized method descriptions");
				}
			} else {
				throw new FileFormatError("File "+scMatFile+" doesn't contain a #SCORE METHOD header");
			}
			
		}
		br.close();
		
		return new AtomCountScorer(scMatFile);
	}
	
	/**
	 * Tells whether given Pdb passes the minimal checks to be valid
	 * as a training set member (i.e. for compiling a scoring matrix). The checks are: 
	 * observed length above {@value #MIN_VALID_CHAIN_LENGTH} and given Pdb is an
	 * all-atom pdb (see Pdb.isAllAtom()) 
	 * @param pdb
	 * @return
	 */
	public static boolean isValidPdb(Pdb pdb) {
		if (pdb.get_length()<=MIN_VALID_CHAIN_LENGTH) {
			return false;
		}
		return pdb.isAllAtom();
	}
}

package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public abstract class Scorer {
	
	// to compile a scoring matrix from the given pdb ids we won't 
	// take anything that has observed length below this value
	private static final int MIN_VALID_CHAIN_LENGTH = 30;


	
	public enum ScoringMethod {
		RESTYPE("residue type"),
		ATOMTYPE("atom type"),
		RESCOUNT("residue count"),
		ATOMCOUNT("atom count");
		
		private String description;
		private ScoringMethod(String description) {
			this.description = description;
		}
		public String getDescription() {
			return description;
		}
	}
	
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

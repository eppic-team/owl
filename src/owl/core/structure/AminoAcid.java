package owl.core.structure;

import java.util.Comparator;
import java.util.HashMap;

/**
 * Package:		proteinstructure
 * Class: 		AminoAcid
 * Author:		Henning Stehr, stehr@molgen.mpg.de
 * Date:		6/Feb/2006, updated 5/Jan/2009
 * 
 * An amino acid in a protein sequence. Each of the twenty naturally occuring 
 * amino acids is an instance of this class. Additionally, a placeholder 'X' for
 * an unknown amino acid and a stop codon are defined. This class also provides
 * static methods for obtaining amino acid objects by their number, one letter
 * code or three letter code and to convert between these representations. 
 * 
 * Note: Currently, most classes in the package still use a simple string representation of
 * amino acids. Helper functions for these are provided in the class AAinfo.
 * 
 * See also: {@link AAinfo}
 * 
 * Changelog:
 * 2006/02/06 first created by HS
 * 2009/01/05 moved to package proteinstructure
 * 2009/03/18 adding stop codon
 */
public enum AminoAcid {
		
	/*---------------------- member variables --------------------------*/

	//                                         hydro  hydro  arom   aliph  polar  charg  pos    neg    small  tiny 
	 ALA ( 1, "Alanine",       'A', "ALA",  1, -0.20, true,  false, false, false, false, false, false, true,  true ), 
	 ARG ( 2, "Arginine",      'R', "ARG",  7,  1.43, false, false, false, true,  true,  true,  false, false, false), 
	 ASN ( 3, "Asparagine",    'N', "ASN",  4,  0.69, false, false, false, true,  false, false, false, true,  false),
	 ASP ( 4, "Aspartic acid", 'D', "ASP",  4,  0.72, false, false, false, true,  true,  false, true,  true,  false),
	 CYS ( 5, "Cysteine",      'C', "CYS",  2, -0.67, true,  false, false, true,  false, false, false, true,  false),
	 GLN ( 6, "Glutamine",     'Q', "GLN",  5,  0.74, false, false, false, true,  false, false, false, false, false),
	 GLU ( 7, "Glutamic acid", 'E', "GLU",  5,  1.09, false, false, false, true,  true,  false, true,  false, false),
	 GLY ( 8, "Glycine",       'G', "GLY",  0, -0.06, true,  false, false, false, false, false, false, true,  true ),
	 HIS ( 9, "Histidine",     'H', "HIS",  6, -0.04, true,  true,  false, true,  true,  true,  false, false, false),
	 ILE (10, "Isoleucine",    'I', "ILE",  4, -0.74, true,  false, true,  false, false, false, false, false, false),
	 LEU (11, "Leucine",       'L', "LEU",  4, -0.65, true,  false, true,  false, false, false, false, false, false),
	 LYS (12, "Lysine",        'K', "LYS",  5,  2.00, true,  false, false, true,  true,  true,  false, false, false),
	 MET (13, "Methionine",    'M', "MET",  4, -0.71, true,  false, false, false, false, false, false, false, false),
	 PHE (14, "Phenylalanine", 'F', "PHE",  7, -0.67, true,  true,  false, false, false, false, false, false, false),
	 PRO (15, "Proline",       'P', "PRO",  3, -0.44, false, false, false, false, false, false, false, true , false),
	 SER (16, "Serine",        'S', "SER",  2,  0.34, false, false, false, true,  false, false, false, true,  true ),
	 THR (17, "Threonine",     'T', "THR",  3,  0.26, false, false, false, true,  false, false, false, true,  false),
	 TRP (18, "Tryptophan",    'W', "TRP", 10, -0.45, true,  true,  false, true,  false, false, false, false, false),
	 TYR (19, "Tyrosine",      'Y', "TYR",  8,  0.22, true,  true,  false, true,  false, false, false, false, false),
	 VAL (20, "Valine",        'V', "VAL",  3, -0.61, true,  false, true,  false, false, false, false, true , false),
	 XXX ( 0, "Unknown",       'X', "XXX", -1,  Double.NaN, false, false, false, false, false, false, false, false, false),
	 STP (-1, "Stop codon",    '*', "STP", -1,  Double.NaN, false, false, false, false, false, false, false, false, false);	 
		
	private int number;
	private String name;			
	private char oneLetterCode;
	private String threeLetterCode;
	private int numberOfAtoms;		// number of heavy (non-Hydrogen) side chain atoms
	private double hydrophobicity;
	private boolean hydrophobic;
	private boolean aromatic;	
	private boolean aliphatic;		
	private boolean polar;
	private boolean charged;	
	private boolean positive; 		// = basic (TODO: is H in this or not?)
	private boolean negative; 		// = acidic	
	private boolean small;	
	private boolean tiny;
	
	/*------------------------- constants ------------------------------*/
	public static final char 	INVALID_ONE_LETTER_CODE 	= '?';
	public static final String 	INVALID_THREE_LETTER_CODE 	= null;
	public static final int 	INVALID_AA_NUMBER 			= -2;
	
	/*---------------------- static variables --------------------------*/
	
	// These variables are being initialized once, when they are first
	// accessed. 
	
	// Retrieval maps: This allows retrieval of amino acids from different
	// representations in O(1) time.
	private static HashMap<Character, AminoAcid> one2aa = initOne2aa();
	private static HashMap<String, AminoAcid> three2aa = initThree2aa();
	private static HashMap<Integer, AminoAcid> num2aa = initNum2aa();
	
	/* ---------------------- constructors -----------------------------*/
	
	AminoAcid(int number, String name, char one, String three, int atoms,
			  double hydrophobicity, // empirical hydrophibicity scale by Miller in kcal/mol
			  boolean hydrophobic, boolean aromatic, boolean aliphatic,
			  boolean polar,       boolean charged,  boolean positive,
			  boolean negative,    boolean small,    boolean tiny) {
		
		this.number = number;
		this.name = name;
		this.oneLetterCode = one;
		this.threeLetterCode = three;
		this.numberOfAtoms = atoms;
		this.hydrophobicity = hydrophobicity;
		this.aromatic = aromatic;
		this.hydrophobic = hydrophobic;
		this.aliphatic = aliphatic;
		this.small = small;
		this.tiny = tiny;
		this.positive = positive;
		this.polar = polar;
		this.charged = charged;
		this.negative = negative;
	}
	
	/*---------------------- standard methods --------------------------*/
	public int getNumber() { return this.number; }
	public String getName() { return this.name; }
	public char getOneLetterCode() { return this.oneLetterCode; }
	public String getThreeLetterCode() { return this.threeLetterCode; }
	
	/**
	 * Returns the number of side chain heavy (non-Hydrogen) atoms for this
	 * AminoAcid 
	 * @return number of side chain heavy atoms
	 */
	public int getNumberOfAtoms() {return this.numberOfAtoms; }
	
	public double getHydrophobicity() {return this.hydrophobicity; }
	public boolean isAromatic() { return this.aromatic; }
	public boolean isHydrophobic() { return this.hydrophobic; }
	public boolean isAliphatic() { return this.aliphatic; }
	public boolean isSmall() { return this.small; }
	public boolean isTiny() { return this.tiny; }
	public boolean isPositive() { return this.positive; }
	public boolean isPolar() { return this.polar; }
	public boolean isCharged() { return this.charged; }
	public boolean isNegative() { return this.negative; }
	
	/**
	 * Returns true if this AminoAcid is one of the 20 standard amino acids
	 * or false otherwise
	 * @return true if this AminoAcid is one of the 20 standard amino acids
	 * or false otherwise
	 */
	public boolean isStandardAA() {
		return (this.getNumber()<21 && this.getNumber()>0);
	}
	
	/*----------------------- static methods ---------------------------*/
	
	/**
	 * Get amino acid object by number
	 * @param num amino acid number (between 0 and 20) (e.g. 0 for Unknown, 1 for Alanine)
	 * @return An amino acid object of the given type 
	 * or null, if num is invalid
	 */
	public static AminoAcid getByNumber(int num) {
		return num2aa.get(num);
	}	
	
	/**
	 * Get amino acid object by one letter code
	 * @param one amino acid one letter code (e.g. "A" for Alanine)
	 * @return An amino acid object of the given type 
	 * or null, if one letter code is invalid
	 */
	public static AminoAcid getByOneLetterCode(char one) {
		return one2aa.get(Character.toUpperCase(one));
	}
	
	/**
	 * Get amino acid object by three letter code
	 * @param three amino acid three letter code (e.g. "ALA" for Alanine)
	 * @return An amino acid object of the given type 
	 * or null, if three letter code is invalid
	 */
	public static AminoAcid getByThreeLetterCode(String three) {
		return three2aa.get(three.toUpperCase());
	}
	
	/**
	 * Get amino acids in ascending order of hydrophobicity.
	 * @return an array containing the amino acids in order of hydrophobicity
	 */
	public static AminoAcid[] valuesInOrderOfHydrophobicity() {
		AminoAcid[] values = AminoAcid.values();
		AminoAcid[] newValues = new AminoAcid[values.length];
		System.arraycopy(values, 0, newValues, 0, values.length);
		java.util.Arrays.sort(newValues, new Comparator<AminoAcid>() {
			public int compare(AminoAcid a, AminoAcid b) {
				return new Double(a.hydrophobicity).compareTo(new Double(b.hydrophobicity));
			}
		});		
		return newValues;
	}
	
	// conversion methods
	
	/**
	 * convert amino acid one letter code to three letter code
	 * @param one amino acid one letter code (e.g. "A" for Alanine)
	 * @return amino acid three letter code or null if input is invalid
	 */
	public static String one2three(char one) {
		AminoAcid aa = getByOneLetterCode(one);
		return aa==null?INVALID_THREE_LETTER_CODE:aa.getThreeLetterCode();
	}
	
	/**
	 * convert amino acid three letter code to one letter code
	 * @param three amino acid three letter code (e.g. "ALA" for Alanine)
	 * @return amino acid one letter code or '?' if input is invalid
	 */
	public static char three2one(String three) {
		AminoAcid aa = getByThreeLetterCode(three);
		return aa==null?INVALID_ONE_LETTER_CODE:aa.getOneLetterCode();
	}
	
	/**
	 * convert amino acid three letter code to number
	 * @param three amino acid three letter code (e.g. "ALA" for Alanine)
	 * @return amino acid number (between 1 and 20) or -2 if input is invalid
	 */
	public static int three2num(String three) {
		AminoAcid aa = getByThreeLetterCode(three);
		return aa==null?INVALID_AA_NUMBER:aa.getNumber();		
	}
	
	/**
	 * convert amino acid number to three letter code 
	 * @param num amino acid number (between -1 and 20) (e.g. 1 for Alanine, -1 for stop codon)
	 * @return amino acid three letter code or null if input is invalid  
	 */
	public static String num2three(int num) {
		AminoAcid aa = getByNumber(num);
		return aa==null?INVALID_THREE_LETTER_CODE:aa.getThreeLetterCode();
	}

	/**
	 * convert amino acid one letter code to number
	 * @param one amino acid one letter code (e.g. "ALA" for Alanine)
	 * @return amino acid number (between 1 and 20 or 0 for unknown type, -1 for stop codon) 
	 * or -2 if input is invalid
	 */
	public static int one2num(char one) {
		AminoAcid aa = getByOneLetterCode(one);
		return aa==null?INVALID_AA_NUMBER:aa.getNumber();
	}
	
	/**
	 * convert amino acid number to one letter code 
	 * @param num amino acid number (between -1 and 20) (e.g. 1 for Alanine, -1 for stop codon, 0 for unknown)
	 * @return amino acid one letter code or '?' if input is invalid  
	 */
	public static char num2one(int num) {
		AminoAcid aa = getByNumber(num);
		return aa==null?INVALID_ONE_LETTER_CODE:aa.getOneLetterCode();		
	}	
	
	/**
	 * Returns true if given string is the three-letter code of a standard aminoacid,
	 * false otherwise
	 * @param three
	 * @return true if given string is the three-letter code of a standard aminoacid, 
	 * false otherwise
	 */
	public static boolean isStandardAA(String three) {
		AminoAcid aa = getByThreeLetterCode(three);
		return aa==null?false:aa.isStandardAA();
	}

	/**
	 * Returns true if given char is the one-letter code of a standard aminoacid,
	 * false otherwise
	 * @param one
	 * @return true if given char is the one-letter code of a standard aminoacid,
	 * false otherwise
	 */
	public static boolean isStandardAA(char one) {
		AminoAcid aa = getByOneLetterCode(one);
		return aa==null?false:aa.isStandardAA();		
	}

	/*----------------------- private methods --------------------------*/
	/**
	 * initialize static map to get amino acid by its ordinal number
	 */
	private static HashMap<Integer, AminoAcid> initNum2aa() {
		HashMap<Integer, AminoAcid> num2aa = new HashMap<Integer, AminoAcid>();
		for(AminoAcid aa:AminoAcid.values()) {
			num2aa.put(aa.getNumber(), aa);
		}
		return num2aa;
	}
	
	/**
	 * initialize static map to get amino acid by its one letter code
	 */
	private static HashMap<Character, AminoAcid> initOne2aa() {
		HashMap<Character, AminoAcid> one2aa = new HashMap<Character, AminoAcid>();
		for(AminoAcid aa:AminoAcid.values()) {
			one2aa.put(aa.getOneLetterCode(), aa);
		}
		return one2aa;
	}
	
	/**
	 * initialize static map to get amino acid by its three letter code
	 */
	private static HashMap<String, AminoAcid> initThree2aa() {
		HashMap<String, AminoAcid> three2aa = new HashMap<String, AminoAcid>();
		for(AminoAcid aa:AminoAcid.values()) {
			three2aa.put(aa.getThreeLetterCode(), aa);
		}
		return three2aa;
	}	
	
	/*--------------------------- main ---------------------------------*/
	
    /*
     * some test for class AminoAcid
     */
	
	public static void main(String[] args) {
		
    	// iterate over all 20 amino acids
		// perform circular conversion and verify result
		// make sure that wrong input will lead to desired result

		String three; char one; int num; AminoAcid result; boolean ok, allright = true;
		
		System.out.println("Testing conversion functions of class AminoAcid:");

		// Step 1
		System.out.print("Step 1");
		ok = true;
		for (AminoAcid aa : AminoAcid.values()) {
			System.out.print(".");
			num = aa.getNumber();
			result = AminoAcid.getByNumber(AminoAcid.three2num(AminoAcid.one2three(AminoAcid.num2one(num))));
			if(!aa.equals(result)) { ok = false; break;}
		}
		if(ok) System.out.println("passed.");
		else System.out.println("failed.");
		allright = allright && ok;
		
		// Step 2
		System.out.print("Step 2");
		ok = true;
		for (AminoAcid aa : AminoAcid.values()) {
			System.out.print(".");
			one = aa.getOneLetterCode();
			result = AminoAcid.getByOneLetterCode(AminoAcid.three2one(AminoAcid.num2three(AminoAcid.one2num(one))));
			if(!aa.equals(result)) { ok = false; break;}
		}
		if(ok) System.out.println("passed.");
		else System.out.println("failed.");
		allright = allright && ok;
		
		// Step 3
		System.out.print("Step 3");
		ok = true;
		for (AminoAcid aa : AminoAcid.values()) {	
			System.out.print(".");
			three = aa.getThreeLetterCode();
			result = AminoAcid.getByThreeLetterCode(AminoAcid.one2three(AminoAcid.num2one(AminoAcid.three2num(three))));
			if(!aa.equals(result)) { ok = false; break;}
		}
		if(ok) System.out.println("passed.");
		else System.out.println("failed.");
		allright = allright && ok;
		
		// Step 4
		System.out.print("Step 4");
		ok = true;	
	
		System.out.print("."); if(AminoAcid.getByNumber(-2) != null) { ok = false; } 
		System.out.print("."); if(AminoAcid.num2one(-2) != '?') { ok = false; }
		System.out.print("."); if(AminoAcid.num2three(-2) != null) { ok = false;  }

		System.out.print("."); if(AminoAcid.getByOneLetterCode('?') != null) { ok = false;  }
		System.out.print("."); if(AminoAcid.one2num('?') != -2) { ok = false;  }
		System.out.print("."); if(AminoAcid.one2three('?') != null) { ok = false;  }
		
		System.out.print("."); if(AminoAcid.getByThreeLetterCode("") != null) { ok = false;  }
		System.out.print("."); if(AminoAcid.three2one("") != '?') { ok = false;  }
		System.out.print("."); if(AminoAcid.three2num("") != -2) { ok = false;  }	
		
		if(ok) System.out.println("passed.");
		else System.out.println("failed.");
		allright = allright && ok;		
		
		// Step 5 : compare number of atoms with data in AAinfo
		System.out.print("Step 5");
		ok = true;
		
		for(AminoAcid aa : AminoAcid.values()) {
			System.out.print(".");
			// for unknown amino acids, the reported number of atoms should be -1
			if(aa == XXX || aa == STP) {
				if(aa.getNumberOfAtoms() != -1) {
					System.err.printf("Unexpected number of atoms reported for unknown amino acid type or stop codon");
					ok = false; 
					break; 					
				}
			} else {	
				if(aa.getNumberOfAtoms() + 4 != AAinfo.getNumberAtoms(aa.getThreeLetterCode())) {
					System.err.printf("Number of atoms reported for %s do not agree\n", aa.getName());
					ok = false; 
					break; 
				}
			}
		}
		
		if(ok) System.out.println("passed.");
		else System.out.println("failed.");
		allright = allright && ok;			
		
		// Summary
		
		if(allright) System.out.println("All tests passed.");
		else System.out.println("Some tests failed.");
		
    }
	
}


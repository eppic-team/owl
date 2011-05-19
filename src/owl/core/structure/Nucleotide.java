package owl.core.structure;

import java.util.HashMap;

public enum Nucleotide {

	A ( 1, "Adenine",  'A', "DA"),
	C ( 2, "Cytosine", 'C', "DC"), 
	G ( 3, "Guanine",  'G', "DG"), 
	U ( 4, "Uracil",   'U', "DU"),
	T ( 5, "Thymine",  'T', "DT");

	private int number;
	private String name;			
	private char oneLetterCode;
	private String twoLetterCode; 

	private Nucleotide(int number, String name, char oneLetterCode, String twoLetterCode) {
		this.number = number;
		this.name = name;
		this.oneLetterCode = oneLetterCode;
		this.twoLetterCode = twoLetterCode;
	}
	
	public int getNumber() {
		return number;
	}

	public String getName() {
		return name;
	}

	public char getOneLetterCode() {
		return oneLetterCode;
	}

	public String getTwoLetterCode() {
		return twoLetterCode;
	}
	
	/**
	 * Returns true if this Nucleotide is one of the 5 standard nucleotides
	 * or false otherwise
	 * @return true if this Nucleotide is one of the 5 standard nucleotides
	 * or false otherwise
	 */
	public boolean isStandardNuc() {
		return (this.getNumber()<=5 && this.getNumber()>0);
	}

	
	private static HashMap<Character, Nucleotide> one2aa = initOne2aa();
	private static HashMap<String, Nucleotide> two2aa = initTwo2aa();
	private static HashMap<String, Nucleotide> full2aa = initFull2aa();
	
	/**
	 * Returns true if given string is the 1-letter or 2-letter code of a standard nucleotide,
	 * false otherwise
	 * @param two
	 * @return true if given string is the 1-letter or 2-letter code of a standard nucleotide, 
	 * false otherwise
	 */
	public static boolean isStandardNuc(String code) {
		if (code.length()==2) {
			Nucleotide nuc = getByTwoLetterCode(code);
			return nuc==null?false:nuc.isStandardNuc();
		} else if (code.length()==1) {
			return isStandardNuc(code.charAt(0));
		} else {
			return false;
		}
	}
	
	/**
	 * Returns true if given char is the one-letter code of a standard nucleotide,
	 * false otherwise
	 * @param one
	 * @return true if given char is the one-letter code of a standard nucleotide,
	 * false otherwise
	 */
	public static boolean isStandardNuc(char one) {
		Nucleotide nuc = getByOneLetterCode(one);
		return nuc==null?false:nuc.isStandardNuc();		
	}
	
	/**
	 * Get Nucleotide object by one letter code
	 * @param one Nucleotide one letter code (e.g. "C" for Cytosine)
	 * @return a Nucleotide object of the given type 
	 * or null, if one letter code is invalid
	 */
	public static Nucleotide getByOneLetterCode(char one) {
		return one2aa.get(Character.toUpperCase(one));
	}
	
	/**
	 * Get Nucleotide object by 2-letter code 
	 * @param two Nucleotide 2-letter code (e.g. "DC" for Cytosine)
	 * @return a Nucleotide object of the given type 
	 * or null, if 2-letter code is invalid
	 */
	public static Nucleotide getByTwoLetterCode(String two) {
		return two2aa.get(two.toUpperCase());
	}
	
	/**
	 * Get Nucleotide object by 1-letter or 2-letter code 
	 * @param two Nucleotide 1-letter or 2-letter code (e.g. "C" or "DC" for Cytosine)
	 * @return a Nucleotide object of the given type 
	 * or null, if 1-letter or 2-letter code is invalid
	 */	
	public static Nucleotide getByCode(String code) {
		if (code.length()==1) {
			return one2aa.get(Character.toUpperCase(code.charAt(0)));
		} else if (code.length()==2){
			return two2aa.get(code.toUpperCase());
		}
		return null;
	}
	
	/**
	 * Get Nucleotide object by full nucleotide name
	 * @param fullName the full name of a nucleotide with first letter capitalised (e.g. "Adenine") 
	 * @return a Nucleotide object of the given type or null if full name is invalid
	 */
	public static Nucleotide getByFullName(String fullName) {
		return full2aa.get(fullName);
	}
	
	public static char two2one(String two) {
		return two2aa.get(two.toUpperCase()).oneLetterCode;
	}
	
	public static String one2two(char one) {
		return one2aa.get(one).twoLetterCode;
	}
	
	/**
	 * initialize static map to get Nucleotide by its one letter code
	 */
	private static HashMap<Character, Nucleotide> initOne2aa() {
		HashMap<Character, Nucleotide> one2aa = new HashMap<Character, Nucleotide>();
		for(Nucleotide aa:Nucleotide.values()) {
			one2aa.put(aa.getOneLetterCode(), aa);
		}
		return one2aa;
	}
	
	/**
	 * initialize static map to get Nucleotide by its two letter code
	 */
	private static HashMap<String, Nucleotide> initTwo2aa() {
		HashMap<String, Nucleotide> two2aa = new HashMap<String, Nucleotide>();
		for(Nucleotide aa:Nucleotide.values()) {
			two2aa.put(aa.getTwoLetterCode(), aa);
		}
		return two2aa;
	}	
	
	private static HashMap<String, Nucleotide> initFull2aa()	{
		HashMap<String, Nucleotide> full2aa = new HashMap<String, Nucleotide>();
		for (Nucleotide aa:Nucleotide.values()) {
			full2aa.put(aa.getName(),aa);
		}
		return full2aa;
	}
	

}

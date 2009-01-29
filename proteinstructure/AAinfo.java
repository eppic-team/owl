package proteinstructure;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * Class with static methods to get aminoacids and contact type information
 * example usage:
 *  AAinfo.isValidContactType("Cg"); // returns true
 *  AAinfo.three2oneletter("ALA"); // returns "A"
 *  
 * The contact types and aas2atoms data are defined in separate text file contactTypes.dat
 * New contact types can be added simply by editing the file 
 * 
 * Beware that everything is static in this file. The JVM will initialise the static
 * variables when they are first called and keep them as if if the static class itself 
 * was a global instantiated object
 * 
 * See also: {@link AminoAcid}
 */
public class AAinfo {

	/*--------------------------- constants ------------------------------*/
	// file with contact type definitions (refers to root of the aglappe package)
	private static final String RESOURCES_DIR = "/proteinstructure/"; // path to data files to be retrieved with getResourceAsStream
	private static final String CT_DEFS_FILE = "contactTypes.dat";
	// file with aminoacid pairs distance bounds definitions (refers to root of aglappe package)
	private static final String AAPAIRSBOUNDS_FILE = "aapairsBounds.dat";
		
	// lower bound distances used for our ConstraintsMaker class
	// from our "empirical" calculations
	public static final double BB_DIAMETER_GYRATION=4.6;
	public static final double DIST_MIN_CA=2.8;
	// "guessed" general min distance from hydrogen the hydrogen bond length (we used it in Cb and Cg)
	public static final double DIST_MIN=2.6;
	
	// if two adjacent CAs are more than this distance apart, we assume a chain break (value taken from Casp assessment)
	public static final double DIST_CHAIN_BREAK=4.5;

	// 1/3 letter code we assign to nonstandard aas to use in sequence
	public static final String NONSTANDARD_AA_ONE_LETTER="X";
	public static final String NONSTANDARD_AA_THREE_LETTER="XXX";
	
	// 1/3 letter code for unknown unobserved residues (used when reading from pdb file with no SEQRES and we have to introduce gaps)
	public static final String UNKNOWN_UNOBSERVED_RES_ONE_LETTER = "X"; 
	public static final String UNKNOWN_UNOBSERVED_RES_THREE_LETTER = "XXX";
	
	
	/*----------------------- "member" variables ---------------------------*/
	private final static Map<String,Double> lowerBoundDistances = initialiseLowerBoundDistances(); 
	
	private final static Map<String,String> one2threeletter = initialiseOne2threeletter();
	private final static Map<String,String> three2oneletter = initialiseThree2oneletter();
	private final static Set<String> aas = initialiseAAs();
	
	private final static Map<String,ContactType> cts = initialiseCTsFromFile();
	
	private final static Map<String,Set<String>> aas2atoms = initialiseAas2atoms(); // depends on cts 
	
	private final static Map<String,String> fullname2threeletter = initialiseFullNames2Threeletter();
	
	private final static Map<String,double[]> aapairs2bounds = initialiseAapairs2BoundsFromFile();
	
	/*----------------------- private methods ----------------------------*/
	private static Map<String,Double> initialiseLowerBoundDistances() {
		Map<String,Double> lowerBoundDistances = new HashMap<String, Double>();
		lowerBoundDistances.put("Ca", DIST_MIN_CA);
		lowerBoundDistances.put("Cb", DIST_MIN);
		lowerBoundDistances.put("Cg", DIST_MIN);
		lowerBoundDistances.put("C", DIST_MIN_CA);
		return lowerBoundDistances;
	}
	
	private static Map<String,String> initialiseOne2threeletter() {
		Map<String,String> one2threeletter = new HashMap<String,String>();
		one2threeletter.put("C", "CYS");
		one2threeletter.put("D", "ASP");
		one2threeletter.put("S", "SER");
		one2threeletter.put("Q", "GLN");
		one2threeletter.put("K", "LYS");
		one2threeletter.put("I", "ILE");
		one2threeletter.put("P", "PRO");
		one2threeletter.put("T", "THR");
		one2threeletter.put("F", "PHE");
		one2threeletter.put("A", "ALA");
		one2threeletter.put("G", "GLY");
		one2threeletter.put("H", "HIS");
		one2threeletter.put("E", "GLU");
		one2threeletter.put("L", "LEU");
		one2threeletter.put("R", "ARG");
		one2threeletter.put("W", "TRP");
		one2threeletter.put("V", "VAL");
		one2threeletter.put("N", "ASN");
		one2threeletter.put("Y", "TYR");
		one2threeletter.put("M", "MET");
		one2threeletter.put(NONSTANDARD_AA_ONE_LETTER, NONSTANDARD_AA_THREE_LETTER);
		one2threeletter.put(UNKNOWN_UNOBSERVED_RES_ONE_LETTER, UNKNOWN_UNOBSERVED_RES_THREE_LETTER);
		return one2threeletter;
	}
	
	private static Map<String,String> initialiseThree2oneletter() {
		Map<String,String> three2oneletter = new HashMap<String,String>();
		three2oneletter.put("CYS", "C");
		three2oneletter.put("ASP", "D");
		three2oneletter.put("SER", "S");
		three2oneletter.put("GLN", "Q");
		three2oneletter.put("LYS", "K");
		three2oneletter.put("ILE", "I");
		three2oneletter.put("PRO", "P");
		three2oneletter.put("THR", "T");
		three2oneletter.put("PHE", "F");
		three2oneletter.put("ALA", "A");
		three2oneletter.put("GLY", "G");
		three2oneletter.put("HIS", "H");
		three2oneletter.put("GLU", "E");
		three2oneletter.put("LEU", "L");
		three2oneletter.put("ARG", "R");
		three2oneletter.put("TRP", "W");
		three2oneletter.put("VAL", "V");
		three2oneletter.put("ASN", "N");
		three2oneletter.put("TYR", "Y");
		three2oneletter.put("MET", "M");
		three2oneletter.put(NONSTANDARD_AA_THREE_LETTER, NONSTANDARD_AA_ONE_LETTER);
		three2oneletter.put(UNKNOWN_UNOBSERVED_RES_THREE_LETTER, UNKNOWN_UNOBSERVED_RES_ONE_LETTER);
		return three2oneletter;
	}
	
	private static Map<String,String> initialiseFullNames2Threeletter() {
		Map<String,String> fullnames2threeletter = new HashMap<String,String>();
		fullnames2threeletter.put("Alanine","ALA");
		fullnames2threeletter.put("Arginine","ARG");
		fullnames2threeletter.put("Asparagine","ASN");
		fullnames2threeletter.put("Aspartic Acid","ASP");
		fullnames2threeletter.put("Cysteine","CYS");
		fullnames2threeletter.put("Glutamic Acid","GLU");
		fullnames2threeletter.put("Glutamine","GLN");
		fullnames2threeletter.put("Glycine","GLY");
		fullnames2threeletter.put("Histidine","HIS");
		fullnames2threeletter.put("Isoleucine","ILE");
		fullnames2threeletter.put("Leucine","LEU");
		fullnames2threeletter.put("Lysine","LYS");
		fullnames2threeletter.put("Methionine","MET");
		fullnames2threeletter.put("Phenylalanine","PHE");
		fullnames2threeletter.put("Proline","PRO");
		fullnames2threeletter.put("Serine","SER");
		fullnames2threeletter.put("Threonine","THR");
		fullnames2threeletter.put("Tryptophan","TRP");
		fullnames2threeletter.put("Tyrosine","TYR");
		fullnames2threeletter.put("Valine","VAL");
		return fullnames2threeletter;
	}
	
	private static Set<String> initialiseAAs() {
		Set<String> aas = new TreeSet<String>();
		aas.add("TRP");
		aas.add("CYS");
		aas.add("GLN");
		aas.add("ALA");
		aas.add("VAL");
		aas.add("LEU");
		aas.add("ASP");
		aas.add("SER");
		aas.add("PRO");
		aas.add("THR");
		aas.add("PHE");
		aas.add("ARG");
		aas.add("LYS");
		aas.add("MET");
		aas.add("HIS");
		aas.add("GLY");
		aas.add("ILE");
		aas.add("ASN");
		aas.add("GLU");
		aas.add("TYR");
		return aas;
	}
	
	private static Map<String,Set<String>> initialiseAas2atoms() {
		Map<String,Set<String>> aas2atoms = new HashMap<String, Set<String>>();
		aas2atoms = cts.get("ALL");
		return aas2atoms;
	}
	
	private static Map<String,ContactType> initialiseCTsFromFile() {
		Map<String,ContactType> cts = new TreeMap<String,ContactType>();
		
		InputStream inp = AAinfo.class.getResourceAsStream(RESOURCES_DIR+CT_DEFS_FILE);
		BufferedReader br = new BufferedReader(new InputStreamReader(inp));
		String line;
		try {
			ContactType contactType = null;
			String ct = "";
			boolean multiAtom = false;
			while ((line = br.readLine())!= null) {
				// skip comments and empty lines
				if (line.startsWith("#")) continue;
				if (line.trim().equals("")) continue;
				if (line.startsWith(">")){ 
					if (!ct.equals("")) { // except for first ct put last res2atoms HashMap for the last ct
						cts.put(ct, contactType);
					}
					Pattern p = Pattern.compile("^>\\s(\\w+)\\s(\\w+)$");
					Matcher m = p.matcher(line);
					if (m.matches()){
						ct = m.group(1);
						String type = m.group(2);
						if (type.equals("multi")) {
							multiAtom = true ;
						} else {
							multiAtom = false;
						}
					}
					contactType = new ContactType(ct,multiAtom);
				} else { // for all other lines
					String aa = line.substring(0,3);
					String atomsStr = line.substring(4).trim();
					String[] atomsArray = new String[0]; // initialisation to empty array
					if (!atomsStr.equals("")) { // if not atomsArray stays empty (for cases of no atoms for a given residue)
						atomsArray = atomsStr.split("\\s");
					}
					Set<String> atoms = new TreeSet<String>();
					for (String atom: atomsArray) {
						atoms.add(atom); // if atomsArray was empty then atoms will be an empty (not null) Set
					}
					contactType.put(aa, atoms);
				}
			}
			cts.put(ct, contactType);
			br.close();
		} catch (IOException e) {
			System.err.println("IO error while reading contact types definition file: "+CT_DEFS_FILE+". Exiting.");
			System.err.println("Error was: "+e.getMessage());
			System.exit(1);
		}
		return cts;
	}
	
	private static Map<String,double[]> initialiseAapairs2BoundsFromFile (){
		Map<String,double[]> aapairs2bounds = new TreeMap<String,double[]>();
		InputStream inp = AAinfo.class.getResourceAsStream(RESOURCES_DIR+AAPAIRSBOUNDS_FILE);
		BufferedReader br = new BufferedReader(new InputStreamReader(inp));
		String line;
		try {
			while ((line = br.readLine())!= null) {
				String[] tokens = line.split("\\s+");
				String pair = tokens[0];
				double min = Double.parseDouble(tokens[1]);
				double max = Double.parseDouble(tokens[2]);
				double[] bounds = {min,max};
				aapairs2bounds.put(pair, bounds);
			}
			br.close();
		} 
		catch (IOException e) {
			System.err.println("IO error while reading aminoacid pairs distance bounds definition file: "+AAPAIRSBOUNDS_FILE+". Exiting.");
			System.err.println("Error was: "+e.getMessage());
			System.exit(1);
		}
		return aapairs2bounds;
	}
	
	/*----------------------- public methods ---------------------------*/
	
	/**
	 * Given a three letter code returns true if is a standard aminoacid
	 */
	public static boolean isValidAA(String three) {
		return aas.contains(three);
	}
	
	/**
	 * Gets all three letter code standard aminoacids in a Set
	 * @return
	 */
	public static Set<String> getAAs() {
		return aas;
	}
	
	/**
	 * Gets all contact type names in a Set
	 * @return
	 */
	public static Set<String> getAllContactTypes() {
		return cts.keySet();
	}
	
	/**
	 * Gets all single atom contact types in a Set
	 * @return
	 */
	public static Set<String> getSingleAtomContactTypes() {
		Set<String> singleAtomCts = new TreeSet<String>();
		for (ContactType contactType:cts.values()) {
			if (!contactType.isMultiAtom()) singleAtomCts.add(contactType.getName());
		}
		return singleAtomCts;
	}
	
	/**
	 * Gets all multiple atom contact types in a Set
	 * @return
	 */
	public static Set<String> getMultiAtomContactTypes() {
		Set<String> multiAtomCts = new TreeSet<String>();
		for (ContactType contactType:cts.values()) {
			if (contactType.isMultiAtom()) multiAtomCts.add(contactType.getName());
		}
		return multiAtomCts;

	}
	
	/**
	 * Returns true if ct is a valid contact type name
	 * Crossed contacts (e.g. BB/SC or Ca/Cg) will also be valid
	 * @param ct
	 * @return
	 */
	public static boolean isValidContactType(String ct){
		Set<String> allCts = getAllContactTypes(); // depends on cts being initialised
		String[] cts1 = ct.split("\\+");
		for(int i = 0; i < cts1.length; i++) {
			String[] cts2 = cts1[i].split("/");
			for(int j = 0; i < cts2.length; i++) {
				if (!allCts.contains(cts2[j])) {
					return false;
				}
			}
		}
		return true;
	}
	
	/**
	 * Returns true if ct is a valid single atom contact type name
	 * Crossed contacts (e.g. Ca/Cg) will also be valid
	 * @param ct
	 * @return
	 */
	// TODO: contact type object
	public static boolean isValidSingleAtomContactType(String ct, boolean directed){
		Set<String> singleAtomCts = getSingleAtomContactTypes(); // depends on cts being initialised
		if (ct.contains("+")) {
			return false;
		}
		if (ct.contains("/")) {
			if (!directed) {
				return false;
			} else {
				String[] cts = ct.split("/");
				if (singleAtomCts.contains(cts[0]) && singleAtomCts.contains(cts[1])) {
					return true;
				} else {
					return false;
				}
			}
		}
		return singleAtomCts.contains(ct);
	}
	
	public static boolean isValidSingleAtomContactType(String ct){
		boolean crossed = false;
		if (ct.contains("/")) {
			crossed = true;
		}		
		return isValidSingleAtomContactType(ct, crossed);
	}
	
	/**
	 * Returns true if ct is a valid multiple atom contact type name
	 * Crossed contacts (e.g. BB/SC) will also be valid
	 * @param ct
	 * @return
	 */
	public static boolean isValidMultiAtomContactType(String ct, boolean directed){
		Set<String> multiAtomCts = getMultiAtomContactTypes(); // depends on cts being initialised
		if (ct.contains("+")) {
			return isValidContactType(ct);
		}
		if (ct.contains("/")) {
			if (!directed) {
				return isValidContactType(ct);
			} else {
				String[] cts = ct.split("/");
				if (multiAtomCts.contains(cts[0]) && multiAtomCts.contains(cts[1])) {
					return true;
				} else {
					return false;
				}			
			}
		}
		return multiAtomCts.contains(ct);
	}
	
	public static boolean isValidMultiAtomContactType(String ct){
		boolean crossed = false;
		if (ct.contains("/")) {
			crossed = true;
		}		
		return isValidMultiAtomContactType(ct, crossed);
	}

	public static boolean isOverlapping(String ct){
		String[] inputCts = ct.split("[+/]");
		for(int i=0;i<(inputCts.length-1);i++) {
			for(int j=(i+1);j<inputCts.length;j++) {
				for(String aa : aas) {
					for (String atom : cts.get(inputCts[j]).get(aa)) {
						if (cts.get(inputCts[i]).get(aa).contains(atom)) return true;
					}
				}
			}
		}
		return false;
	}
	
	/**
	 * Returns the lower bound distance for assigning distance restraints 
	 * to contacts given a contact type
	 * Returns -1.0 if contact type is not one of the contact types for which we can assign constraints
	 * @param ct
	 * @param aa1
	 * @param aa2
	 * @return
	 */
	public static double getLowerBoundDistance(String ct, String aa1, String aa2) {
		if (isValidSingleAtomContactType(ct)) {
			return lowerBoundDistances.get(ct);
		} else if (ct.equals("BB")){
			return DIST_MIN_CA;
		} else if (ct.equals("SC")) {
			// in aapairs2bounds we have the pairs only in one direction (order defined by alphabetical order)
			String minAA = aa1;
			String maxAA = aa2;
			if (aa1.compareTo(aa2)>0) {
				minAA = aa2;
				maxAA = aa1;
			}
			return aapairs2bounds.get(minAA+"_"+maxAA)[0];
		}
		return -1.0;
	}
	
	/**
	 * Returns the upper bound distance for assigning distance restraints
	 * to contacts given a contact type
	 * Returns -1.0 if contact type is not one of the contact types for which we can assign constraints
	 * @param ct
	 * @param aa1
	 * @param aa2
	 * @return
	 */
	public static double getUpperBoundDistance(String ct, String aa1, String aa2) {
		if (isValidSingleAtomContactType(ct)) {
			return 0.0;
		} else if (ct.equals("BB")) {
			return BB_DIAMETER_GYRATION;
		} else if (ct.equals("SC")) {
			// in aapairs2bounds we have the pairs only in one direction (order defined by alphabetical order)
			String minAA = aa1;
			String maxAA = aa2;
			if (aa1.compareTo(aa2)>0) {
				minAA = aa2;
				maxAA = aa1;
			}
			return aapairs2bounds.get(minAA+"_"+maxAA)[1];
		}
		return -1.0;
	}
	
	/**
	 * Converts from one letter aminoacid codes to three letter codes
	 * If invalid input returns null
	 * @param one
	 * @return
	 */
	public static String oneletter2threeletter(String one) {
		return one2threeletter.get(one);
	}
	
	/**
	 * Converts from three letter aminoacid codes to one letter codes
	 * If invalid input returns null
	 * @param three
	 * @return
	 */
	public static String threeletter2oneletter(String three) {
		return three2oneletter.get(three);
	}
	
	/**
	 * Converts from aminoacid full names (capitalised first letter, rest lower case) 
	 * to three letter codes
	 * @param full
	 * @return
	 */
	public static String fullname2threeletter(String full){
		return fullname2threeletter.get(full);
	}
	
	/**
	 * Returns true if given String is a valid aminoacid name
	 * (first letter capitalised, rest lower case)
	 * @param full
	 * @return
	 */
	public static boolean isValidFullName(String full) {
		return fullname2threeletter.keySet().contains(full);
	}
	
	/**
	 * Returns all aminoacid full names in a Set
	 * @return
	 */
	public static Set<String> getAAFullNames(){
		return fullname2threeletter.keySet();
	}
	
	/**
	 * Given a three letter code aminoacid and an atom name say whether 
	 * the atom is a valid atom for that aminoacid 
	 * Doesn't consider OXT to be a valid atom for any aminoacid
	 * @param aa
	 * @param atom
	 * @return
	 */
	public static boolean isValidAtom(String aa, String atom) {
		return aas2atoms.get(aa).contains(atom);
	}

	/**
	 * Given a three letter code aminoacid and an atom name say whether 
	 * the atom is a valid atom for that aminoacid 
	 * Considers OXT to be a valid atom for all aminoacids
	 * @param aa
	 * @param atom
	 * @return
	 */
	public static boolean isValidAtomWithOXT(String aa, String atom) {
		if (atom.equals("OXT")) return true;
		return aas2atoms.get(aa).contains(atom);
	}

	/**
	 * Gets all (non-Hydrogen) atoms for an aminoacid (three letter code)
	 * @param aa
	 * @return
	 */
	public static Set<String> getAtoms(String aa) {
		return aas2atoms.get(aa);
	}
	
	/**
	 * Gets the number of non-hydrogen atoms for an aminoacid (three letter code)
	 * @param aa
	 * @return
	 */
	public static int getNumberAtoms(String aa) {
		return aas2atoms.get(aa).size();
	}
	
	/**
	 * Returns a Set of all atom names given an aminoacid and a contact type
	 * e.g. for aa="SER" and ct="SC" returns ["CB", "CG"] 
	 * @param ct
	 * @param aa
	 * @return
	 */
	public static Set<String> getAtomsForCTAndRes(String ct, String aa) {
		return cts.get(ct).get(aa);
	} 
	
	/**
	 * Gets the identifier for gap-characters in three-letter-code-sequences. 
	 * @return a three letter code representation of the gap character.
	 */
	public static String getGapCharacterThreeLetter() {
	    return "GAP";
	}
	
	/**
	 * Returns the one letter representation of a gap in a protein sequence
	 * @return a one letter representation of a gap in a protein sequence
	 */
	public static char getGapCharacterOneLetter() {
	    return '-';
	}
	
}

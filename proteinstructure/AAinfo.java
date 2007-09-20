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
 */
public class AAinfo {

	/*--------------------------- constants ------------------------------*/
	// file with contact type definitions
	// refers to root of the aglappe package
	private static final String CT_DEFS_FILE = "/proteinstructure/contactTypes.dat";
	
	// lower bound distances used for our ConstraintsMaker class
	// from our "empirical" calculations
	public static final double BB_DIAMETER_GYRATION=4.6;
	public static final double DIST_MIN_CA=2.8;
	// "guessed" general min distance from hydrogen the hydrogen bond length (we used it in Cb and Cg)
	public static final double DIST_MIN=2.6;

	
	/*----------------------- member variables ---------------------------*/
	private final static Map<String,Double> lowerBoundDistances = initialiseLowerBoundDistances(); 
	
	private final static Map<String,String> one2threeletter = initialiseOne2threeletter();
	private final static Map<String,String> three2oneletter = initialiseThree2oneletter();
	private final static Set<String> aas = initialiseAAs();
	
	private final static Map<String,ContactType> cts = initialiseCTsFromFile();
	
	private final static Map<String,Set<String>> aas2atoms = initialiseAas2atoms(); // depends on cts 
	
	
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
		return three2oneletter;
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
		
		InputStream inp = Runtime.getRuntime().getClass().getResourceAsStream(CT_DEFS_FILE);
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
		} catch (IOException e) {
			System.err.println("IO error while reading contact types definition file: "+CT_DEFS_FILE+". Exiting.");
			System.err.println("Error was: "+e.getMessage());
			System.exit(1);
		}
		return cts;
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
		if (ct.contains("/")){
			String[] cts = ct.split("/");
			if (allCts.contains(cts[0]) && allCts.contains(cts[1])) {
				return true;
			} else {
				return false;
			}
		}
		return allCts.contains(ct);
	}
	
	/**
	 * Returns true if ct is a valid single atom contact type name
	 * Crossed contacts (e.g. Ca/Cg) will also be valid
	 * @param ct
	 * @return
	 */
	public static boolean isValidSingleAtomContactType(String ct){
		Set<String> singleAtomCts = getSingleAtomContactTypes(); // depends on cts being initialised
		if (ct.contains("/")){
			String[] cts = ct.split("/");
			if (singleAtomCts.contains(cts[0]) && singleAtomCts.contains(cts[1])) {
				return true;
			} else {
				return false;
			}
		}
		return singleAtomCts.contains(ct);
	}
	
	/**
	 * Returns true if ct is a valid multiple atom contact type name
	 * Crossed contacts (e.g. BB/SC) will also be valid
	 * @param ct
	 * @return
	 */
	public static boolean isValidMultiAtomContactType(String ct){
		Set<String> multiAtomCts = getMultiAtomContactTypes(); // depends on cts being initialised
		if (ct.contains("/")){
			String[] cts = ct.split("/");
			if (multiAtomCts.contains(cts[0]) && multiAtomCts.contains(cts[1])) {
				return true;
			} else {
				return false;
			}
		}
		return multiAtomCts.contains(ct);
	}
	
	/**
	 * Gets the lower bound distance for assigning distance restraints 
	 * to contacts given a contact type
	 * @param ct
	 * @return
	 */
	public static double getLowerBoundDistance(String ct) {
		return lowerBoundDistances.get(ct);
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
	 * Returns a Set of all atom names given an aminoacid and a contact type
	 * e.g. for aa="SER" and ct="SC" returns ["CB", "CG"] 
	 * @param ct
	 * @param aa
	 * @return
	 */
	public static Set<String> getAtomsForCTAndRes(String ct, String aa) {
		return cts.get(ct).get(aa);
	}

}

package owl.core.structure;

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
 * Class with static methods to get contact type information
 * example usage:
 *  AAinfo.isValidContactType("Cg"); // returns true
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
	// file with contact type definitions (refers to root of the OWL library)
	private static final String RESOURCES_DIR = "/owl/core/structure/"; // path to data files to be retrieved with getResourceAsStream
	private static final String CT_DEFS_FILE = "contactTypes.dat";
	// file with aminoacid pairs distance bounds definitions (refers to root of OWL library)
	private static final String AAPAIRSBOUNDS_FILE = "aapairsBounds.dat";
		
	// lower bound distances used for our ConstraintsMaker class
	// from our "empirical" calculations
	public static final double BB_DIAMETER_GYRATION=4.6;
	public static final double DIST_MIN_CA=2.8;
	// "guessed" general min distance from hydrogen the hydrogen bond length (we used it in Cb and Cg)
	public static final double DIST_MIN=2.6;
	
	// if two adjacent CAs are more than this distance apart, we assume a chain break (value taken from Casp assessment)
	public static final double DIST_CHAIN_BREAK=4.5;

	/*----------------------- "member" variables ---------------------------*/
	private final static Map<String,Double> lowerBoundDistances = initialiseLowerBoundDistances(); 
	
	private final static Map<String,ContactType> cts = initialiseCTsFromFile();
	
	private final static Map<String,Set<String>> aas2heavyAtoms = initialiseAas2heavyAtoms(); // grabs the data from the cts map above 
	private final static Map<String,Set<String>> aas2atoms = initialiseAas2atoms(); // grabs the data from the cts map above
	
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
	
	private static Map<String,Set<String>> initialiseAas2heavyAtoms() {
		Map<String,Set<String>> aas2atoms = new HashMap<String, Set<String>>();
		aas2atoms = cts.get("ALL");
		return aas2atoms;
	}
	
	private static Map<String,Set<String>> initialiseAas2atoms() {
		Map<String,Set<String>> aas2atoms = new HashMap<String, Set<String>>();
		aas2atoms = cts.get("ALL_H");
		return aas2atoms;
	}
	
	private static Map<String,ContactType> initialiseCTsFromFile() {
		Map<String,ContactType> cts = new TreeMap<String,ContactType>();
		
		InputStream inp = null;
		inp = AAinfo.class.getResourceAsStream(RESOURCES_DIR+CT_DEFS_FILE);
		if(inp == null) {
			System.err.println("Severe Error: Resource " + RESOURCES_DIR+CT_DEFS_FILE + " not found. Could not initialize contact types. Exiting.");
			System.exit(1);
		}
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
		InputStream inp = null;
		inp = AAinfo.class.getResourceAsStream(RESOURCES_DIR+AAPAIRSBOUNDS_FILE);
		if(inp == null) {
			System.err.println("Severe Error: Resource " + RESOURCES_DIR+AAPAIRSBOUNDS_FILE + " not found. Could not initialize AA pair bounds. Exiting.");
			System.exit(1);
		}
		BufferedReader br = new BufferedReader(new InputStreamReader(inp));
		String line;
		try {
			while ((line = br.readLine())!= null) {
				// skip comments and empty lines
				if (line.startsWith("#")) continue;
				if (line.trim().equals("")) continue;
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
				for(AminoAcid aa : AminoAcid.getAllStandardAAs()) {
					for (String atom : cts.get(inputCts[j]).get(aa.getThreeLetterCode())) {
						if (cts.get(inputCts[i]).get(aa.getThreeLetterCode()).contains(atom)) return true;
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
	 * Given a three letter code aminoacid and an atom name say whether 
	 * the atom is a valid (non-Hydrogen) atom for that aminoacid 
	 * Doesn't consider OXT to be a valid atom for any aminoacid
	 * @param aa
	 * @param atom
	 * @return
	 */
	public static boolean isValidHeavyAtom(String aa, String atom) {
		return aas2heavyAtoms.get(aa).contains(atom);
	}

	/**
	 * Given a three letter code aminoacid and an atom name say whether 
	 * the atom is a valid (non-Hydrogen) atom for that aminoacid 
	 * Considers OXT to be a valid atom for all aminoacids
	 * @param aa
	 * @param atom
	 * @return
	 */
	public static boolean isValidHeavyAtomWithOXT(String aa, String atom) {
		if (atom.equals("OXT")) return true;
		return aas2heavyAtoms.get(aa).contains(atom);
	}

	/**
	 * Given a three letter code aminoacid and an atom name say whether 
	 * the atom is a valid atom (including Hydrogens) for that aminoacid 
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
	 * the atom is a valid atom (including Hydrogens) for that aminoacid 
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
	 * Gets all (non-Hydrogen) atoms given a three letter code aminoacid
	 * @param aa
	 * @return
	 */
	public static Set<String> getAtoms(String aa) {
		return aas2heavyAtoms.get(aa);
	}
	
	/**
	 * Gets the number of non-Hydrogen atoms given a three letter code aminoacid
	 * @param aa
	 * @return
	 * @throws NullPointerException if aa not a valid 3 letter code standard aminoacid
	 */
	public static int getNumberAtoms(String aa) {
		return aas2heavyAtoms.get(aa).size();
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

package proteinstructure;
import java.util.HashMap;
import java.util.ArrayList;

public class AA {
	
	public static final String CONTACT_TYPE_C_ALPHA = "Ca";
	public static final String CONTACT_TYPE_C_BETA = "Cb";
	
	public static String[] contactTypes() {
		// NOTE: whenever a new contact type is added it needs to be added here as well as in ct2atoms
		String[] cts ={"ALL","BB","SC","Ca","Cb","Cg","C"};
		return cts;
	}
	
	public static String[] singleAtomContactTypes() {
		// NOTE: whenever a new contact type is added it needs to be added here as well as in ct2atoms
		String[] cts ={"Ca","Cb","Cg","C"};
		return cts;
	}

	public static String[] multiAtomContactTypes() {
		// NOTE: whenever a new contact type is added it needs to be added here as well as in ct2atoms
		String[] cts ={"ALL","BB","SC"};
		return cts;
	}
	
	/**
	 * Returns true if given contact type is a valid one
	 * @param ct
	 * @return
	 */
	public static boolean isValidCT(String ct){
		if (!ct.contains("/")){
			for (String validCt:contactTypes()){
				if (ct.equals(validCt)) return true;
			}
		} else {
			String[] cts = ct.split("/");
			String i_ct = cts[0];
			String j_ct = cts[1];
			if (isValidCT(i_ct) && isValidCT(j_ct)) return true;
		}
		return false;
	}

	/**
	 * Returns true if given contact type is a valid single atom contact type
	 * @param ct
	 * @return
	 */
	public static boolean isValidSingleAtomCT(String ct) {
		if (!ct.contains("/")){
			for (String validCt:singleAtomContactTypes()){
				if (ct.equals(validCt)) return true;
			}
		} else {
			String[] cts = ct.split("/");
			String i_ct = cts[0];
			String j_ct = cts[1];
			if (isValidSingleAtomCT(i_ct) && isValidSingleAtomCT(j_ct)) return true;			
		}
		return false;		
	}

	/**
	 * Returns true if given contact type is a valid multiple atom contact type
	 * @param ct
	 * @return
	 */
	public static boolean isValidMultiAtomCT(String ct) {
		if (!ct.contains("/")){
			for (String validCt:multiAtomContactTypes()){
				if (ct.equals(validCt)) return true;
			}
		} else {
			String[] cts = ct.split("/");
			String i_ct = cts[0];
			String j_ct = cts[1];
			if (isValidMultiAtomCT(i_ct) && isValidMultiAtomCT(j_ct)) return true;						
		}
		return false;		
	}

	private static HashMap<String,String> getThreeletter2oneletter() {
		HashMap<String,String> three2oneletter = new HashMap<String,String>();
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
		three2oneletter.put("TYR", "Y") ;
		three2oneletter.put("MET", "M");
		return three2oneletter;
	}
	
	public static String threeletter2oneletter(String three) {
		HashMap<String,String> three2oneletter = getThreeletter2oneletter();
		return three2oneletter.get(three);
	}
	
	private static HashMap<String,String> getOneletter2Threeletter(){
		HashMap<String,String> one2threeletter = new HashMap<String,String>();
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
	
	public static String oneletter2threeletter(String one) {
		HashMap<String,String> one2threeletter = getOneletter2Threeletter();
		return one2threeletter.get(one);
	}
	
	public static ArrayList<String> aas() {
		HashMap<String,String> three2oneletter = getThreeletter2oneletter();
		ArrayList<String> aas = new ArrayList<String>();
		for (String aa:three2oneletter.keySet()) {
			aas.add(aa);
		}
		return aas;
	}
	
	public static HashMap<String,ArrayList<String>> getaas2atoms() {
		HashMap<String,ArrayList<String>> aas2atoms = new HashMap<String,ArrayList<String>>();
		aas2atoms.put("CYS", new ArrayList<String>());
		aas2atoms.get("CYS").add("SG");
		aas2atoms.get("CYS").add("CB");
		aas2atoms.get("CYS").add("O");
		aas2atoms.get("CYS").add("CA");
		aas2atoms.get("CYS").add("C");
		aas2atoms.get("CYS").add("N");
		aas2atoms.put("ASP", new ArrayList<String>());
		aas2atoms.get("ASP").add("N");
		aas2atoms.get("ASP").add("CB");
		aas2atoms.get("ASP").add("OD2");
		aas2atoms.get("ASP").add("CG");
		aas2atoms.get("ASP").add("CA");
		aas2atoms.get("ASP").add("C");
		aas2atoms.get("ASP").add("O");
		aas2atoms.get("ASP").add("OD1");
		aas2atoms.put("SER", new ArrayList<String>());
		aas2atoms.get("SER").add("C");
		aas2atoms.get("SER").add("N");
		aas2atoms.get("SER").add("CA");
		aas2atoms.get("SER").add("O");
		aas2atoms.get("SER").add("CB");
		aas2atoms.get("SER").add("OG");
		aas2atoms.put("GLN", new ArrayList<String>());
		aas2atoms.get("GLN").add("CA");
		aas2atoms.get("GLN").add("CG");
		aas2atoms.get("GLN").add("N");
		aas2atoms.get("GLN").add("NE2");
		aas2atoms.get("GLN").add("OE1");
		aas2atoms.get("GLN").add("CD");
		aas2atoms.get("GLN").add("C");
		aas2atoms.get("GLN").add("O");
		aas2atoms.get("GLN").add("CB");
		aas2atoms.put("LYS", new ArrayList<String>());
		aas2atoms.get("LYS").add("C");
		aas2atoms.get("LYS").add("CA");
		aas2atoms.get("LYS").add("N");
		aas2atoms.get("LYS").add("O");
		aas2atoms.get("LYS").add("CB");
		aas2atoms.get("LYS").add("CG");
		aas2atoms.get("LYS").add("CD");
		aas2atoms.get("LYS").add("CE");
		aas2atoms.get("LYS").add("NZ");
		aas2atoms.put("ASN", new ArrayList<String>());
		aas2atoms.get("ASN").add("ND2");
		aas2atoms.get("ASN").add("OD1");
		aas2atoms.get("ASN").add("CB");
		aas2atoms.get("ASN").add("O");
		aas2atoms.get("ASN").add("CG");
		aas2atoms.get("ASN").add("C");
		aas2atoms.get("ASN").add("CA");
		aas2atoms.get("ASN").add("N");
		aas2atoms.put("PRO", new ArrayList<String>());
		aas2atoms.get("PRO").add("O");
		aas2atoms.get("PRO").add("C");
		aas2atoms.get("PRO").add("CB");
		aas2atoms.get("PRO").add("CG");
		aas2atoms.get("PRO").add("CA");
		aas2atoms.get("PRO").add("CD");
		aas2atoms.get("PRO").add("N");
		aas2atoms.put("THR", new ArrayList<String>());
		aas2atoms.get("THR").add("CA");
		aas2atoms.get("THR").add("N");
		aas2atoms.get("THR").add("C");
		aas2atoms.get("THR").add("CG2");
		aas2atoms.get("THR").add("OG1");
		aas2atoms.get("THR").add("CB");
		aas2atoms.get("THR").add("O");
		aas2atoms.put("PHE", new ArrayList<String>());
		aas2atoms.get("PHE").add("O");
		aas2atoms.get("PHE").add("CE2");
		aas2atoms.get("PHE").add("CE1");
		aas2atoms.get("PHE").add("CG");
		aas2atoms.get("PHE").add("C");
		aas2atoms.get("PHE").add("N");
		aas2atoms.get("PHE").add("CA");
		aas2atoms.get("PHE").add("CB");
		aas2atoms.get("PHE").add("CD2");
		aas2atoms.get("PHE").add("CD1");
		aas2atoms.get("PHE").add("CZ");
		aas2atoms.put("ALA", new ArrayList<String>());
		aas2atoms.get("ALA").add("CA");
		aas2atoms.get("ALA").add("C");
		aas2atoms.get("ALA").add("N");
		aas2atoms.get("ALA").add("CB");
		aas2atoms.get("ALA").add("O");
		aas2atoms.put("HIS", new ArrayList<String>());
		aas2atoms.get("HIS").add("CE1");
		aas2atoms.get("HIS").add("CD2");
		aas2atoms.get("HIS").add("ND1");
		aas2atoms.get("HIS").add("CG");
		aas2atoms.get("HIS").add("CB");
		aas2atoms.get("HIS").add("O");
		aas2atoms.get("HIS").add("C");
		aas2atoms.get("HIS").add("CA");
		aas2atoms.get("HIS").add("N");
		aas2atoms.get("HIS").add("NE2");
		aas2atoms.put("GLY", new ArrayList<String>());
		aas2atoms.get("GLY").add("N");
		aas2atoms.get("GLY").add("CA");
		aas2atoms.get("GLY").add("C");
		aas2atoms.get("GLY").add("O");
		aas2atoms.put("ILE", new ArrayList<String>());
		aas2atoms.get("ILE").add("CG2");
		aas2atoms.get("ILE").add("CD1");
		aas2atoms.get("ILE").add("O");
		aas2atoms.get("ILE").add("N");
		aas2atoms.get("ILE").add("CA");
		aas2atoms.get("ILE").add("C");
		aas2atoms.get("ILE").add("CB");
		aas2atoms.get("ILE").add("CG1");
		aas2atoms.put("LEU", new ArrayList<String>());
		aas2atoms.get("LEU").add("N");
		aas2atoms.get("LEU").add("CA");
		aas2atoms.get("LEU").add("C");
		aas2atoms.get("LEU").add("O");
		aas2atoms.get("LEU").add("CB");
		aas2atoms.get("LEU").add("CG");
		aas2atoms.get("LEU").add("CD2");
		aas2atoms.get("LEU").add("CD1");
		aas2atoms.put("ARG", new ArrayList<String>());
		aas2atoms.get("ARG").add("NH2");
		aas2atoms.get("ARG").add("CZ");
		aas2atoms.get("ARG").add("NE");
		aas2atoms.get("ARG").add("CD");
		aas2atoms.get("ARG").add("CG");
		aas2atoms.get("ARG").add("C");
		aas2atoms.get("ARG").add("CA");
		aas2atoms.get("ARG").add("N");
		aas2atoms.get("ARG").add("NH1");
		aas2atoms.get("ARG").add("CB");
		aas2atoms.get("ARG").add("O");
		aas2atoms.put("TRP", new ArrayList<String>());
		aas2atoms.get("TRP").add("CE2");
		aas2atoms.get("TRP").add("CA");
		aas2atoms.get("TRP").add("C");
		aas2atoms.get("TRP").add("O");
		aas2atoms.get("TRP").add("CB");
		aas2atoms.get("TRP").add("CG");
		aas2atoms.get("TRP").add("CD1");
		aas2atoms.get("TRP").add("CD2");
		aas2atoms.get("TRP").add("CH2");
		aas2atoms.get("TRP").add("CZ3");
		aas2atoms.get("TRP").add("CZ2");
		aas2atoms.get("TRP").add("CE3");
		aas2atoms.get("TRP").add("NE1");
		aas2atoms.get("TRP").add("N");
		aas2atoms.put("VAL", new ArrayList<String>());
		aas2atoms.get("VAL").add("CG1");
		aas2atoms.get("VAL").add("CB");
		aas2atoms.get("VAL").add("O");
		aas2atoms.get("VAL").add("C");
		aas2atoms.get("VAL").add("CA");
		aas2atoms.get("VAL").add("N");
		aas2atoms.get("VAL").add("CG2");
		aas2atoms.put("GLU", new ArrayList<String>());
		aas2atoms.get("GLU").add("C");
		aas2atoms.get("GLU").add("N");
		aas2atoms.get("GLU").add("CA");
		aas2atoms.get("GLU").add("OE2");
		aas2atoms.get("GLU").add("OE1");
		aas2atoms.get("GLU").add("CD");
		aas2atoms.get("GLU").add("CG");
		aas2atoms.get("GLU").add("CB");
		aas2atoms.get("GLU").add("O");
		aas2atoms.put("TYR", new ArrayList<String>());
		aas2atoms.get("TYR").add("CD2");
		aas2atoms.get("TYR").add("CD1");
		aas2atoms.get("TYR").add("OH");
		aas2atoms.get("TYR").add("CZ");
		aas2atoms.get("TYR").add("CE2");
		aas2atoms.get("TYR").add("N");
		aas2atoms.get("TYR").add("CA");
		aas2atoms.get("TYR").add("C");
		aas2atoms.get("TYR").add("O");
		aas2atoms.get("TYR").add("CB");
		aas2atoms.get("TYR").add("CG");
		aas2atoms.get("TYR").add("CE1");
		aas2atoms.put("MET", new ArrayList<String>());
		aas2atoms.get("MET").add("CA");
		aas2atoms.get("MET").add("O");
		aas2atoms.get("MET").add("C");
		aas2atoms.get("MET").add("CB");
		aas2atoms.get("MET").add("CE");
		aas2atoms.get("MET").add("N");
		aas2atoms.get("MET").add("CG");
		aas2atoms.get("MET").add("SD");
		return aas2atoms;
		
	}
	
	public static HashMap<String,String[]> ct2atoms(String ct) {
		ArrayList<String> aas = aas();
		HashMap<String,String[]> ct2atoms = new HashMap<String,String[]>();
		if (ct.equals("Ca")){
			String[] atoms = {"CA"};
			for (String aa:aas) {
				ct2atoms.put(aa, atoms);
			}
		} 
		else if (ct.equals("Cb")){
			String[] atoms = {"CB"};
			for (String aa:aas) {
				ct2atoms.put(aa, atoms);
			}
			atoms = new String[1];
			atoms[0]="CA";
			ct2atoms.put("GLY", atoms);
		}
		else if (ct.equals("C")){
			String[] atoms = {"C"};
			for (String aa:aas) {				
				ct2atoms.put(aa, atoms);
			}			
		}
		else if (ct.equals("Cg")){
			String[] atoms = {"SG"};
			ct2atoms.put("CYS", atoms);
			atoms = new String[1];
			atoms[0]= "CG";
			ct2atoms.put("ASP", atoms);
			atoms = new String[1];
			atoms[0]="OG";
			ct2atoms.put("SER", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("GLN", atoms);
			atoms = new String[1];
			atoms[0]="CD";
			ct2atoms.put("LYS", atoms);
			atoms = new String[1];
			atoms[0]="CG1";
			ct2atoms.put("ILE", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("PRO", atoms);
			atoms = new String[1];
			atoms[0]="OG1";
			ct2atoms.put("THR", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("PHE", atoms);
			ct2atoms.put("ALA", new String[0]);
			ct2atoms.put("GLY", new String[0]);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("HIS", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("GLU", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("LEU", atoms);
			atoms = new String[1];
			atoms[0]="NE";
			ct2atoms.put("ARG", atoms);
			atoms = new String[1];
			atoms[0]="CD2";
			ct2atoms.put("TRP", atoms);
			atoms = new String[1];
			atoms[0]="CG1";
			ct2atoms.put("VAL", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("ASN", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("TYR", atoms);
			atoms = new String[1];
			atoms[0]="CG";
			ct2atoms.put("MET", atoms);			
		}
		else if (ct.equals("BB")){
			String[] atoms = {"CA", "N", "C", "O"};
			for (String aa:aas) {
				ct2atoms.put(aa, atoms);
			}
		}
		else if (ct.equals("SC")){
			HashMap<String,ArrayList<String>> aas2atoms = getaas2atoms();			
			for (String aa:aas) {
				ArrayList<String> SCatoms =aas2atoms.get(aa);
				SCatoms.remove("CA");
				SCatoms.remove("N");
				SCatoms.remove("C");
				SCatoms.remove("O");
				String[] SCatomsar= new String[SCatoms.size()];
				SCatomsar=SCatoms.toArray(SCatomsar);
				ct2atoms.put(aa, SCatomsar);
			}
		}
		else if (ct.equals("ALL")){
			HashMap<String,ArrayList<String>> aas2atoms = getaas2atoms();
			for (String aa:aas) {
				ArrayList<String> ALLatoms = aas2atoms.get(aa);
				String[] ALLatomsar= new String[ALLatoms.size()];
				ALLatomsar=ALLatoms.toArray(ALLatomsar);				
				ct2atoms.put(aa,ALLatomsar);
			}
		}		
		return ct2atoms;
	}
}

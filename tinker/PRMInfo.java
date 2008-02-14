package tinker;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.AAinfo;

/**
 * Class containing information (atom identities) parsed from 
 * tinker's prm files (force field parameter sets)
 *
 */
public class PRMInfo {

	public enum PRMType {amber};	// currently only amber supported
	
	String prmFileName;
	
	TreeMap<Integer,Integer> prmid2tinkerid;
	TreeMap<Integer,String> prmid2prmname;
	TreeMap<Integer,String> prmid2fullname;
	
	TreeMap<Integer,String> prmid2pdbname;
	
	TreeMap<Integer,String> prmid2res_atom;
	
	/**
	 * Constructs a PRMInfo object, mapping the PRM atom names to PDB atom names
	 * @param prmFileName
	 * @param type
	 * @throws IOException
	 */
	public PRMInfo(String prmFileName, PRMType type) throws IOException{
		this.prmFileName = prmFileName;		
		
		readPrmFile();
		
		if (type==PRMType.amber){
			mapPdbAtomNamesAmber();
		} else {
			System.err.println("prm file type "+type+" not supported");
		}
	}
	
	private void readPrmFile() throws IOException{
		prmid2tinkerid = new TreeMap<Integer, Integer>();
		prmid2prmname = new TreeMap<Integer, String>();
		prmid2fullname = new TreeMap<Integer, String>();
		
		
		BufferedReader fprm = new BufferedReader(new FileReader(prmFileName));
		String line;
		while((line = fprm.readLine()) != null ) {
			//atom      1    14    N       "Glycine N"                 7     14.010     3
			Pattern p = Pattern.compile("^atom\\s+(\\d+)\\s+(\\d+)\\s+([a-zA-Z0-9*]+)\\s+\"([^\"]+)\""); // note for atom names we use [a-zA-Z0-9*] instead of \w (because * also appears in some atom names)
			Matcher m = p.matcher(line);
			if (m.find()){
				int prmId = Integer.parseInt(m.group(1));
				int tinkerId = Integer.parseInt(m.group(2));
				String prmName = m.group(3);
				String fullName = m.group(4);
				prmid2tinkerid.put(prmId, tinkerId);
				prmid2prmname.put(prmId,prmName);
				prmid2fullname.put(prmId,fullName);
			}
		}
		fprm.close();
	}
	
	private void mapPdbAtomNamesAmber() {
		//NOTE amber uses some special atom names as compared to pdb:
		//					pdb			amber
		// ARG:				NH1,NH2		NH  	atoms are indistinguishable
		// GLU:				OE1,OE2		OE		atoms are indistinguishable
		// PHE: 			CD1,CD2		CD		atoms are indistinguishable
		// 					CE1,CE2		CE		atoms are indistinguishable
		// TYR: 			CD1,CD2		CD		atoms are indistinguishable
		//					CE1,CE2		CE		atoms are indistinguishable
		// ILE:				CD1			CD		change of nomenclature
		// ASP: 			OD1,OD2		OD		they don't seem to be indistiguishable, why do they use the same name??
		// all c-term aas:	O, OXT		OXT		they don't seem to be indistiguishable, why do they use the same name??

		prmid2res_atom = new TreeMap<Integer, String>();
		
		for (int prmid:prmid2fullname.keySet()){
			String fullname = prmid2fullname.get(prmid);

			String molName = fullname.substring(0, fullname.lastIndexOf(" ")).trim();
			molName = molName.replace("C-Term ", "");
			molName = molName.replace("N-Term ", "");
			String atomName = fullname.substring(fullname.lastIndexOf(" ")).trim();
			// amber uses OXT as the atom name for BOTH Oxygens in the mainchain of the C-Term aminoacid
			// we don't want OXTs so here we replace it by O which is the symbol for the usual backbone Oxygen
			// After when we map xyz to pdb atom serials, the first amber OXT type atom encountered in the c-term aminoacid of the xyz file will be considered as the backbone Oxygen
			if (atomName.equals("OXT")) {
				atomName = "O";
			}
			
			for (String name:AAinfo.getAAFullNames()) {
				// here we can't just call AAinfo.isValidFullName as we test only if the string STARTS like a full amino acid name
				// this is because some aminoacid names are 1 word and some 2 words (e.g. Aspartic Acid) 
				if(molName.startsWith(name)) molName=AAinfo.fullname2threeletter(name);
			}
			
			if (AAinfo.isValidAA(molName.substring(0, 3))) {
				// Hydrogens will also be in the mapping, in case we need them later
				prmid2res_atom.put(prmid, molName.substring(0,3)+"_"+atomName);
			}

		}		
		
	}
	
	public String getRes_AtomFromPrmid(int prmid) {
		return prmid2res_atom.get(prmid);
	}
	
	/*
	"Alanine"
	"Arginine"
	"Asparagine"
	"Aspartic Acid"
	"Cysteine (-SH)"
	"Glutamic Acid"
	"Glutamine"
	"Glycine"
	"Histidine (+)"
	"Histidine (HD)"
	"Histidine (HE)"
	"Isoleucine"
	"Leucine"
	"Lysine"
	"Methionine"
	"Phenylalanine"
	"Proline"
	"Serine"
	"Threonine"
	"Tryptophan"
	"Tyrosine"
	"Valine"
*/
}

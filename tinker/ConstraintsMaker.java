package tinker;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Formatter;
import java.util.Locale;
import java.util.Set;
import java.util.TreeMap;

import proteinstructure.AAinfo;
import proteinstructure.Edge;
import proteinstructure.Graph;
import proteinstructure.Pdb;
import proteinstructure.PdbChainCodeNotFoundError;
import proteinstructure.PdbfileFormatError;
import proteinstructure.PdbfilePdb;

/**
 * Reads tinker's xyz file and pdb file (result of converting the xyz file using xyzpdb program) and
 * maps the xyz atom serials to the pdb atom serials
 * The mapping is done through the PRMInfo class that reads prm files and map pdb atom names to prm atom identifiers
 * 
 * Method createConstraints takes a Graph object and writes to file atom distance constraints in tinker key file format              
 *  
 */
public class ConstraintsMaker {

	private static final double DEFAULT_FORCECONSTANT = 100.0;
	
	private File xyzFile;
	private Pdb pdb;
	private PrintWriter fkey;
	
	private TreeMap<Integer,Integer> pdb2xyz;
	
	private PRMInfo prminfo;
	private String type; // amber, charmm, ...
	
	private String lastPdbResSerial_Atom;
	
	public ConstraintsMaker(File pdbFile, File xyzFile, File prmFile, String type, File keyFile) throws FileNotFoundException, IOException, PdbfileFormatError {
		this.xyzFile = xyzFile;
		this.fkey = new PrintWriter(new FileOutputStream(keyFile));
		try {
			this.pdb = new PdbfilePdb(pdbFile.getAbsolutePath(),"NULL");
		} catch (PdbChainCodeNotFoundError e){
			// this shouldn't happen as the pdb chain code is hard coded, we print stack trace and exit
			System.err.println(e.getMessage());
			e.printStackTrace();
		}
		this.type = type;
		
		this.lastPdbResSerial_Atom = "";
		
		prminfo = new PRMInfo(prmFile.getAbsolutePath(),type);
		
		this.mapAtomSerials();
	}
	
	private void mapAtomSerials() throws IOException {
		pdb2xyz = new TreeMap<Integer, Integer>();
		
		int pdbResSerial = 1; // our pointer to the current pdb residue serial as we read the xyz file
		
		String sequence = pdb.getSequence();
		int numAtoms = 0; // count of atoms per residue
		String resFromSeq = ""; // 3 letter residue code that we take from the sequence
				
		// reading xyz file
		BufferedReader fxyz = new BufferedReader(new FileReader(xyzFile));
		String line;
		fxyz.readLine(); // we skip first line which contains the title
		while((line = fxyz.readLine()) != null ) {
			int xyzAtomSer = Integer.parseInt(line.substring(0,6).trim());
			int prmId = Integer.parseInt(line.substring(48,53).trim());
			
			// from the prmId we can get residue type and atom from our prminfo object
			String res_atom = prminfo.getRes_AtomFromPrmid(prmId);
			String res = res_atom.split("_")[0];
			String atom = res_atom.split("_")[1];
		
			resFromSeq = AAinfo.oneletter2threeletter(String.valueOf((sequence.charAt(pdbResSerial-1))));
			int totalNumAtoms = AAinfo.getNumberAtoms(resFromSeq);
			
			// when current atom counts coincides with the totalNumAtoms this residue type should have (Hs excluded) then we increment residue serial and reset atom count
			if (numAtoms==totalNumAtoms) {
				pdbResSerial++;
				numAtoms = 0;
			}
			
			// if we are in last aminoacid and we've counted all non-Hydrogen atoms for it then we don't want to continue reading lines in the xyz file
			if (pdbResSerial>sequence.length()) { 
				break;
			}
			
			if (!atom.startsWith("H")) { // we don't want Hydrogens as we don't have them in the pdb object
				numAtoms++;
				atom = getCorrectedPdbAtomName(res,atom,pdbResSerial);
				int pdbAtomSer = pdb.getAtomSerFromResSerAndAtom(pdbResSerial, atom);
				if (!pdb.getResTypeFromResSerial(pdbResSerial).equals(res)){
					// sanity check, shouldn't happen at all
					System.err.println("error! res types don't match for res serial "+pdbResSerial+" res type "+res);
				}
				// we assign the atom serial mapping
				pdb2xyz.put(pdbAtomSer, xyzAtomSer);
			}
		}
		fxyz.close();
	}
	
	private String getCorrectedPdbAtomName(String res, String atom, int pdbResSerial) {
		boolean first; // if true it is the first time we have this atom and residue serial
		String thisPdbResSerial_Atom = pdbResSerial+"_"+atom;
		if (thisPdbResSerial_Atom.equals(lastPdbResSerial_Atom)) {
			first = true;
		} else {
			first = false;
		}
		lastPdbResSerial_Atom = thisPdbResSerial_Atom;
		
		if (type.equals("amber")) {
			// amber uses some special atom names as compared to pdb:
			//		pdb			amber
			// ARG:	NH1,NH2		NH  	atoms are indistinguishable
			// GLU:	OE1,OE2		OE		atoms are indistinguishable
			// PHE: CD1,CD2		CD		atoms are indistinguishable
			// 		CE1,CE2		CE		atoms are indistinguishable
			// TYR: CD1,CD2		CD		atoms are indistinguishable
			//		CE1,CE2		CE		atoms are indistinguishable
			// ILE:	CD1			CD		change of nomenclature
			// ASP: OD1,OD2		OD		they don't seem to be indistiguishable, why do they use the same name??
			if (res.equals("ILE") && atom.equals("CD")) {
				return "CD1";
			}			
			if (res.equals("ARG") && atom.equals("NH")) {
				if (first)	return "NH1";
				else 		return "NH2";
			}
			if (res.equals("GLU") && atom.equals("OE")) {
				if (first) 	return "OE1";
				else 		return "OE2";
			}
			if ((res.equals("PHE") || res.equals("TYR")) && atom.equals("CD")) {
				if (first) 	return "CD1";
				else 		return "CD2";
			}
			if ((res.equals("PHE") || res.equals("TYR")) && atom.equals("CE")) {
				if (first) 	return "CE1";
				else 		return "CE2";
			}			
			if (res.equals("ASP") && atom.equals("OD")) {
				if (first) 	return "OD1"; 
				else 		return "OD2";
			}
		}
		// if nothing else returned then we return the same atom name as we passed
		return atom;
	}
	
	/**
	 * Writes to keyFile distance constraints in tinker's format (RESTRAIN-DISTANCE lines) b
	 * based on the mapping of pdb atom serials to xyz atom serials done in mapAtomSerials
	 * Valid contact types to use are: 
	 * 	- all single atom types and any cross between them
	 *  - BB, SC, BB/SC
	 *  - all crosses between BB, SC and single atom contact types, e.g. BB/Cg
	 * That means: 'ALL' CANNOT be used
	 * 
	 * @param graph
	 * @throws Exception If the graph's contact type is ALL
	 */
	public void createConstraints(Graph graph) throws Exception {
		
		for (Edge cont:graph.contacts){
			String i_res = graph.getResType(cont.i);
			String j_res = graph.getResType(cont.j);
			String ct = graph.getContactType();
			String i_ct = ct;
			String j_ct = ct;
			if (ct.contains("/")){
				i_ct = ct.split("/")[0];
				j_ct = ct.split("/")[1];
			}
			if (i_ct.equals("ALL") || j_ct.equals("ALL")) {
				throw new Exception("ALL is not a valid contact type for creating constraints");
			}
			double cutoff = graph.getCutoff();
			double forceConstant = DEFAULT_FORCECONSTANT;
			//TODO get force constants from weights

			Set<String> i_atoms = AAinfo.getAtomsForCTAndRes(i_ct, i_res);
			Set<String> j_atoms = AAinfo.getAtomsForCTAndRes(j_ct, j_res);

			// as dist_min we take the average of the two dist mins, if i_ct and j_ct are the same then this will be the same as dist_min for ct
			double dist_min = (AAinfo.getLowerBoundDistance(i_ct,i_res,j_res)+AAinfo.getLowerBoundDistance(j_ct,i_res,j_res))/2;
			// for single atom contact types getUpperBoundDistance and getLowerBoundDistance will return 0 thus for those cases dist_max = cutoff
			double dist_max = AAinfo.getUpperBoundDistance(i_ct, i_res, j_res)/2+AAinfo.getUpperBoundDistance(i_ct, i_res, j_res)/2+cutoff;
			
			for (String i_atom:i_atoms) {
				for (String j_atom:j_atoms) {
					int i_pdb = pdb.getAtomSerFromResSerAndAtom(cont.i, i_atom);
					int i_xyz = pdb2xyz.get(i_pdb);
					int j_pdb = pdb.getAtomSerFromResSerAndAtom(cont.j, j_atom);
					int j_xyz = pdb2xyz.get(j_pdb);
					fkey.println(new Formatter().format(Locale.US,"RESTRAIN-DISTANCE %s %s %5.1f %2.1f %2.1f",i_xyz,j_xyz,forceConstant,dist_min,dist_max).toString());
				}
			}
			
		}
	}
	
	public void closeKeyFile() {
		fkey.close();
	}
}

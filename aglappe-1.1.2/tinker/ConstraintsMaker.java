package tinker;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Locale;
import java.util.Set;
import java.util.TreeMap;

import edu.uci.ics.jung.graph.util.Pair;
import graphAveraging.ConsensusSquare;

import proteinstructure.AAinfo;
import proteinstructure.AminoAcid;
import proteinstructure.Pdb;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbfilePdb;
import proteinstructure.RIGEdge;
import proteinstructure.RIGNode;
import proteinstructure.RIGraph;

/**
 * Reads tinker's xyz file and pdb file (result of converting the xyz file using xyzpdb program) and
 * maps the xyz atom serials to the pdb atom serials
 * The mapping is done through the PRMInfo class that reads prm files and map pdb atom names to prm atom identifiers
 * 
 * Method createConstraints takes a Graph object and writes to file atom distance constraints in tinker key file format              
 *  
 */
public class ConstraintsMaker {

	public static final int OMEGA_LOWER_BOUND = 178 - 5; // taken from Whatcheck output: omega angles for trans-peptides
	public static final int OMEGA_UPPER_BOUND = 178 + 5; // are gaussian distributed with mean 178 and stddev 5.5
	
	private File xyzFile;
	private Pdb pdb;
	private PrintWriter fkey;
	
	private TreeMap<Integer,Integer> pdb2xyz;
	
	private PRMInfo prminfo;
	private PRMInfo.PRMType type; // amber, charmm, ...
	
	private String lastPdbResSerial_Atom;
	
	/**
	 * Creates a new constraints maker.
	 * @param pdbFile the output file of Tinker's protein program converted to PDB format
	 * @param xyzFile the output file of Tinker's protein program
	 * @param prmFile a Tinker compliant force field parameter file
	 * @param type the type of the force field parameter file (currently only 'amber' is supported)
	 * @param keyFile output tinker restraints file
	 * @throws FileNotFoundException if one of the files was not found
	 * @throws IOException if something went wrong while reading from or writing to files
	 * @throws PdbLoadError if the PDB file could not be read
	 */
	public ConstraintsMaker(File pdbFile, File xyzFile, File prmFile, PRMInfo.PRMType type, File keyFile) throws FileNotFoundException, IOException, PdbLoadError {
		this.xyzFile = xyzFile;
		this.fkey = new PrintWriter(new FileOutputStream(keyFile));

		this.pdb = new PdbfilePdb(pdbFile.getAbsolutePath());
		this.pdb.load(Pdb.NULL_CHAIN_CODE);
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
		
			resFromSeq = AminoAcid.one2three(sequence.charAt(pdbResSerial-1));
			int totalNumAtoms = 0;
			boolean inUNKres = false;
			if (!AminoAcid.isStandardAA(resFromSeq)) { // i.e. a non-standard aa (X)
				// the 'protein' program interprets the Xs in the sequence as 'UNK', but actually writes
				// an xyz file treating the residue as if it were a 'GLY' (just the seq file keeps the UNK 
				// in the sequence). Then 'xyzpdb' writes a pdb file that contains coordinates for that 
				// residue (again as if it were a GLY), but names the residue 'UNK'.
				// Because of this we don't want to consider these atoms at all
				// we set a flag and use GLY for the totalNumAtoms (we still need the totalNumAtoms to be
				// able to keep track of on which residue we are on)
				inUNKres = true;
				totalNumAtoms = AAinfo.getNumberAtoms("GLY");
			} else {
				totalNumAtoms = AAinfo.getNumberAtoms(resFromSeq);
			}
			
			// when current atom counts coincides with the totalNumAtoms this residue type should have (Hydrogens excluded) 
			// then we increment residue serial and reset atom count
			if (numAtoms==totalNumAtoms) {
				pdbResSerial++;
				numAtoms = 0;
			}
			
			// if we are in last aminoacid and we've counted all non-Hydrogen atoms for it then we don't want to continue 
			// reading lines in the xyz file
			if (pdbResSerial>sequence.length()) { 
				break;
			}
			
			if (!atom.startsWith("H")) { // we are not interested in Hydrogens as we don't use them for constraints
				numAtoms++;
				if (!inUNKres) {
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
		
		if (type==PRMInfo.PRMType.amber) {
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
	 * 
	 * @param graph
	 * @param defaultForceConstant a global force constant used for all distance restraints
	 * @throws IllegalArgumentException if contact type of graph is 'ALL' or if it's a non valid contact type
	 */
	public void createDistanceConstraints(RIGraph graph, double defaultForceConstant) {

		double cutoff = graph.getCutoff();
		String ct = graph.getContactType();
		String i_ct = ct;
		String j_ct = ct;
		if (ct.contains("/")){
			i_ct = ct.split("/")[0];
			j_ct = ct.split("/")[1];
		}
		
		// ALL is not a valid contact type for creating constraints!
		if (i_ct.equals("ALL") || j_ct.equals("ALL")) {
			throw new IllegalArgumentException("ALL is not a valid contact type for creating constraints.");
		}
		if (!AAinfo.isValidContactType(i_ct) || !AAinfo.isValidContactType(j_ct)) {
			throw new IllegalArgumentException("Either "+i_ct+" or "+j_ct+" are not valid contact types");
		}
		
		for (RIGEdge cont:graph.getEdges()){
			Pair<RIGNode> pair = graph.getEndpoints(cont);
			String i_res = pair.getFirst().getResidueType();
			String j_res = pair.getSecond().getResidueType();

			//TODO get force constants from weights

			if (!AminoAcid.isStandardAA(i_res) || !AminoAcid.isStandardAA(j_res)) {
				// we have to skip contacts that involve non standard aminoacids
				continue;
			}
			Set<String> i_atoms = AAinfo.getAtomsForCTAndRes(i_ct, i_res);
			Set<String> j_atoms = AAinfo.getAtomsForCTAndRes(j_ct, j_res);

			// as dist_min we take the average of the two dist mins, if i_ct and j_ct are the same then this will be the same as dist_min for ct
			double dist_min = (AAinfo.getLowerBoundDistance(i_ct,i_res,j_res)+AAinfo.getLowerBoundDistance(j_ct,i_res,j_res))/2;
			// for single atom contact types getUpperBoundDistance and getLowerBoundDistance will return 0 thus for those cases dist_max = cutoff
			double dist_max = AAinfo.getUpperBoundDistance(i_ct, i_res, j_res)/2+AAinfo.getUpperBoundDistance(i_ct, i_res, j_res)/2+cutoff;
			
			for (String i_atom:i_atoms) {
				for (String j_atom:j_atoms) {
					int i_pdb = pdb.getAtomSerFromResSerAndAtom(pair.getFirst().getResidueSerial(), i_atom);
					int i_xyz = pdb2xyz.get(i_pdb);
					int j_pdb = pdb.getAtomSerFromResSerAndAtom(pair.getSecond().getResidueSerial(), j_atom);
					int j_xyz = pdb2xyz.get(j_pdb);
					fkey.printf(Locale.US,"RESTRAIN-DISTANCE %s %s %5.1f %2.1f %2.1f\n",
							i_xyz,j_xyz,defaultForceConstant,dist_min,dist_max);
				}
			}
			
		}
	}
	
	public void closeKeyFile() {
		fkey.close();
	}
	
	public void createPhiPsiConstraints(TreeMap<Integer,ConsensusSquare> phiPsiConsensus, double defaultForceConstantPhiPsi) {
		for (int resser:phiPsiConsensus.keySet()) {
			// we can't assign a phi angle for residue 1 or a psi for the last residue so we have to skip those
			if (resser==1) continue;
			if (resser==pdb.getFullLength()) continue;
			if (!pdb.hasCoordinates(resser-1) || !pdb.hasCoordinates(resser) || !pdb.hasCoordinates(resser+1)) {
				// Here we don't have to take care of unobserved as the pdb file we read is output of tinker's protein and 
				// thus with coordinates for all residues
				// BUT! if there were non-standard aas in the sequence, then our pdb object will not have coordinates for 
				// those (the file does have the coordinates but we don't read non-standard aas).
				// We thus have to check if resser or +1 or -1 are not there
				continue;
			}
			
			// get all atoms necessary for the phi/psi angles 
			int Ciminus1 = pdb2xyz.get(pdb.getAtomSerFromResSerAndAtom(resser-1, "C"));
			int Ni       = pdb2xyz.get(pdb.getAtomSerFromResSerAndAtom(resser, "N"));
			int CAi      = pdb2xyz.get(pdb.getAtomSerFromResSerAndAtom(resser, "CA"));
			int Ci       = pdb2xyz.get(pdb.getAtomSerFromResSerAndAtom(resser, "C"));
			int Niplus1  = pdb2xyz.get(pdb.getAtomSerFromResSerAndAtom(resser+1, "N"));
			
			// phi restraint
			fkey.printf("RESTRAIN-TORSION %d %d %d %d %5.1f %d %d\n",
					Ciminus1, Ni, CAi, Ci, defaultForceConstantPhiPsi, 
					phiPsiConsensus.get(resser).getConsInterval1stDim().beg, phiPsiConsensus.get(resser).getConsInterval1stDim().end);
			// psi restraint
			fkey.printf("RESTRAIN-TORSION %d %d %d %d %5.1f %d %d\n", 
					Ni, CAi, Ci, Niplus1, defaultForceConstantPhiPsi, 
					phiPsiConsensus.get(resser).getConsInterval2ndDim().beg, phiPsiConsensus.get(resser).getConsInterval2ndDim().end);
		}
	}
	
	public void createOmegaConstraints(double defaultForceConstantOmega) {
		for (int resser:pdb.getAllSortedResSerials()) {
			// we can't assign an omega angle for the last residue so we have to skip those
			if (resser==pdb.getFullLength()) continue;
			if (!pdb.hasCoordinates(resser) || !pdb.hasCoordinates(resser+1)) {
				// Here we don't have to take care of unobserved as the pdb file we read is output of tinker's protein and 
				// thus with coordinates for all residues
				// BUT! if there were non-standard aas in the sequence, then our pdb object will not have coordinates for 
				// those (the file does have the coordinates but we don't read non-standard aas).
				// We thus have to check if resser or resser+1 are not there
				continue;
			}
			
			// get all atoms necessary for the phi/psi angles 
			int CAi = pdb2xyz.get(pdb.getAtomSerFromResSerAndAtom(resser, "CA"));
			int Ci       = pdb2xyz.get(pdb.getAtomSerFromResSerAndAtom(resser, "C"));
			int Niplus1  = pdb2xyz.get(pdb.getAtomSerFromResSerAndAtom(resser+1, "N"));
			int CAiplus1 = pdb2xyz.get(pdb.getAtomSerFromResSerAndAtom(resser+1, "CA"));
			
			// phi restraint
			fkey.printf("RESTRAIN-TORSION %d %d %d %d %5.1f %d %d\n",
					CAi, Ci, Niplus1, CAiplus1, defaultForceConstantOmega, 
					OMEGA_LOWER_BOUND, OMEGA_UPPER_BOUND);
		}
	}
}

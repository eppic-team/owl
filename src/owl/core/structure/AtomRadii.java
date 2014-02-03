package owl.core.structure;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Class containing methods to parse the vdw.radii resource file from NACCESS.
 * The radii values are from Chothia (1976) J.Mol.Biol.105,1-14
 * @author duarte_j
 *
 */
public class AtomRadii {

	// the van der waals radii file from NACCESS
	private static final String RADFILE = "/owl/core/structure/vdw.radii";

	private static final InputStream vdwradIS = AtomRadii.class.getResourceAsStream(RADFILE);
	private static final HashMap<String, HashMap<String,Double>> radii = parseRadFile();
	
	public static final double TRIGONAL_CARBON_VDW = 1.76;
	public static final double TETRAHEDRAL_CARBON_VDW = 1.87;
	public static final double TRIGONAL_NITROGEN_VDW = 1.65;
	public static final double TETRAHEDRAL_NITROGEN_VDW = 1.50;
	public static final double SULFUR_VDW = 1.85;
	public static final double OXIGEN_VDW = 1.40;

	/**
	 * Gets the radius for given amino acid and atom (PDB convention)
	 * @param aa
	 * @param atom
	 * @return
	 */
	public static double getRadius(AminoAcid aa, Atom atom) {
		if (atom.getType()==null) {
			System.err.println("Unrecognised atom "+atom.getCode()+" in residue "+atom.getParentResSerial()+"-"+aa.getThreeLetterCode()+
					", setting its radius to radius of default unknown atom (Nitrogen).");
			return AtomType.X.getRadius();
		}
		if (atom.getType().equals(AtomType.H)) return AtomType.H.getRadius();
		// some unusual entries (e.g. 1tes) contain Deuterium atoms in standard aminoacids
		if (atom.getType().equals(AtomType.D)) return AtomType.D.getRadius();
		String atomCode = atom.getCode();
		if (atomCode.equals("OXT")) atomCode="O";
		if (!radii.get(aa.getThreeLetterCode()).containsKey(atomCode)) {
			System.err.println("Unexpected atom "+atomCode+" in a "+aa.getThreeLetterCode()+" standard residue. Will use its generic vdw radius.");
			return atom.getType().getRadius();
		}
		return radii.get(aa.getThreeLetterCode()).get(atomCode);
	}
	
	public static double getRadiusNew(AminoAcid aa, Atom atom) {
		if (atom.getType()==null) {
			System.err.println("Warning: unrecognised atom "+atom.getCode()+" in residue "+atom.getParentResSerial()+"-"+aa.getThreeLetterCode()+
					", setting its vdw radius to radius of default unknown atom (Nitrogen).");
			return AtomType.X.getRadius();
		}
		if (atom.getType().equals(AtomType.H)) return AtomType.H.getRadius();
		// some unusual entries (e.g. 1tes) contain Deuterium atoms in standard aminoacids
		if (atom.getType().equals(AtomType.D)) return AtomType.D.getRadius();

		String atomCode = atom.getCode();
		
		// here we use the values that Chothia gives in his paper (as NACCESS does)
		if (atom.getType()==AtomType.O) {
			return OXIGEN_VDW;
		} 
		else if (atom.getType()==AtomType.S) {
			return SULFUR_VDW;
		}
		else if (atom.getType()==AtomType.N) {
			if (atomCode.equals("NZ")) return TETRAHEDRAL_NITROGEN_VDW; // tetrahedral Nitrogen
			return TRIGONAL_NITROGEN_VDW;								// trigonal Nitrogen
		}
		else if (atom.getType()==AtomType.C) { // it must be a carbon
			if (atomCode.equals("C") || 
					atomCode.equals("CE1") || atomCode.equals("CE2") || atomCode.equals("CE3") ||
					atomCode.equals("CH2") || 
					atomCode.equals("CZ") || atomCode.equals("CZ2") || atomCode.equals("CZ3")) {
				return TRIGONAL_CARBON_VDW; 							// trigonal Carbon
			}
			else if (atomCode.equals("CA") || atomCode.equals("CB") || 
					atomCode.equals("CE") ||
					atomCode.equals("CG1") || atomCode.equals("CG2")) {
				return TETRAHEDRAL_CARBON_VDW;							// tetrahedral Carbon
			}
			// left cases depend on amino acid: CD, CD1, CD2, CG
			else {
				switch (aa) {
					case PHE:						
					case TRP:
					case TYR:
					case HIS:
					case ASP:
					case ASN:
						return TRIGONAL_CARBON_VDW;
						
					case PRO:
					case LYS:
					case ARG:
					case MET:
					case ILE:
					case LEU:
						return TETRAHEDRAL_CARBON_VDW;
						
					case GLN:
					case GLU:
						if (atomCode.equals("CD")) return TRIGONAL_CARBON_VDW;
						else if (atomCode.equals("CG")) return TETRAHEDRAL_CARBON_VDW;
						
					default:
						System.err.println("Warning: unexpected carbon atom "+atomCode+" for aminoacid "+aa.getThreeLetterCode()+", assigning its standard vdw radius");
						return AtomType.C.getRadius();
				}
			}
			
		// not any of the expected atoms
		} else {
			System.err.println("Warning: unexpected atom "+atomCode+" for aminoacid "+aa.getThreeLetterCode()+", assigning its standard vdw radius");
			return atom.getType().getRadius();
		}
	}

	
	/**
	 * Gets the radius for given nucleotide and atom (PDB convention)
	 * @param nuc
	 * @param atom
	 * @return
	 */
	public static double getRadius(Nucleotide nuc, Atom atom) {
		if (atom.getType()==null) {
			System.err.println("Unrecognised atom "+atom.getCode()+" in residue "+atom.getParentResSerial()+"-"+nuc.getTwoLetterCode()+
					", setting its radius to radius of default unknown atom (Nitrogen).");
			return AtomType.X.getRadius();
		}
		if (atom.getType().equals(AtomType.H)) return AtomType.H.getRadius();
		if (atom.getType().equals(AtomType.D)) return AtomType.D.getRadius();
		String atomCode = atom.getCode();
		// in RNA there's an additional O2' that is not in DNA (which is what vdw.radii has)
		if (atomCode.equals("O2'")) return AtomType.O.getRadius();
		// in T of DNA there's C7, for some reason not present in vdw.radii
		if (atomCode.equals("C7")) atomCode = "C6";
		return radii.get(nuc.getTwoLetterCode()).get(atomCode);
	}
	
	/**
	 * Gets the radius for given HET residue 3 letter code and atom (PDB convention)
	 * The NACCESS vdw.radii file contains radii values for only a few HET residues.
	 * For all others the value returned is the default vdw radius from the {@link AtomType} enum. 
	 * @param mol3lettercode
	 * @param atom
	 * @return
	 */
	public static double getRadius(String mol3lettercode, Atom atom) {
		if (atom.getType()==null) {
			System.err.println("Unrecognised atom "+atom.getCode()+" in residue "+atom.getParentResSerial()+"-"+mol3lettercode+
					", setting its radius to radius of default unknown atom (Nitrogen).");
			return AtomType.X.getRadius();
		}
		if (atom.getType().equals(AtomType.H)) return AtomType.H.getRadius();
		if (!radii.containsKey(mol3lettercode) || 
			!radii.get(mol3lettercode).containsKey(atom.getCode())) { // this can happen if for a het aa there is an unknown atom (X atom)
			return atom.getType().getRadius();
		}
		return radii.get(mol3lettercode).get(atom.getCode());
	}
	
	private static HashMap<String, HashMap<String,Double>> parseRadFile() {
		HashMap<String, HashMap<String,Double>> radii = new HashMap<String, HashMap<String,Double>>();
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(vdwradIS));
			String line;
			HashMap<String,Double> perResRad = null;
			boolean isNucleotide = false;
			while ((line=br.readLine())!=null) {
				if (line.startsWith("#")) continue;
				
				if (line.startsWith("RESIDUE ")) {
					isNucleotide = false;
					String res = line.substring(13,16);
					if (res.startsWith("__")) {
						// the nucleotides start with a __ that we strip here and convert to 2-letter code
						res = Nucleotide.getByOneLetterCode(res.charAt(2)).getTwoLetterCode();
						isNucleotide = true;
					}
					
					perResRad = new HashMap<String, Double>();
					radii.put(res, perResRad);
				}
				if (line.startsWith("ATOM")) {
					String atomCode = line.substring(5,9).trim();
					// in HET residue HEM the NA, NB, NC, ND have a space in between (no idea why). This is to fix that
					if (atomCode.contains(" ")) atomCode = atomCode.replaceAll(" ", "");
					if (isNucleotide) {
						// in nucleotide atoms * are used instead of ', this is to fix that:
						if (atomCode.contains("*")) atomCode = atomCode.replaceAll("\\*","'");
						// and the oxigens of the P are not in standard PDB nomenclature either
						if (atomCode.equals("O1P")) atomCode = "OP1";
						if (atomCode.equals("O2P")) atomCode = "OP2";
					}
					double radius = Double.parseDouble(line.substring(10,14));
					perResRad.put(atomCode,radius);
				}
			}
			br.close();
		} catch (IOException e) {
			System.err.println("Fatal error! Can't read resource file "+RADFILE+". Error: "+e.getMessage()+". Exiting.");
			System.exit(1);
		}
		return radii;
	}
	
	
	// a tester method to print all aas from a pdb and show the differences between generic 
	// element vdw radii and those taken from the naccess vdw radii file  
	public static void main (String[] args) throws Exception {
		
		PdbAsymUnit pdb = new PdbAsymUnit(new File("/home/duarte_j/3hbx.cif"));
		
		TreeMap<String, Residue> uniqueAas = new TreeMap<String,Residue>();
		
		for (PdbChain chain:pdb.getPolyChains()) {
			
			chain.setAtomRadii();
			
			for (Residue res:chain) {
				if (res instanceof AaResidue) {
				
					uniqueAas.put(res.getLongCode(),res);
				}
			}
		}

		
		System.out.println(uniqueAas.size()+" aminoacids");
		
		TreeMap<String,String> uniqueValues = new TreeMap<String,String>();
		TreeMap<String,TreeSet<String>> uniqueValues2List = new TreeMap<String,TreeSet<String>>();
		
		for (Residue res:uniqueAas.values()) {
			
			
			
			System.out.println("## "+res.getLongCode());
			for (Atom atom:res) {
				String valueStr = String.format("%4.2f",atom.getRadius());
				uniqueValues.put(valueStr, atom.getType().getSymbol());
				if (!uniqueValues2List.containsKey(valueStr)) {
					TreeSet<String> list = new TreeSet<String>();
					list.add(atom.getCode());
					uniqueValues2List.put(valueStr, list);
				} else {
					uniqueValues2List.get(valueStr).add(atom.getCode());
				}
					
				
				
				if (Math.abs(atom.getType().getRadius()-atom.getRadius())>0.001) {
					System.out.printf("%s\t%5.2f\t%5.2f\n",atom.getCode(),atom.getType().getRadius(),atom.getRadius());				
				
				}
			}
		}
		
		System.out.println("Unique radii values (different from default AtomType vdw radius): ");
		for (String value:uniqueValues.keySet()) {
			System.out.println(uniqueValues.get(value)+" "+value);
		}
		
		System.out.println("Unique values and atom classes belonging to them: ");
		for (String value:uniqueValues2List.keySet()) {
			System.out.print(value);
			for (String atomCode:uniqueValues2List.get(value)) {
				System.out.print(" "+atomCode);
			}
			System.out.println();
		}
		
		// comparing old and new implementation

		for (PdbChain chain:pdb.getPolyChains()) {

			chain.setAtomRadii();

			for (Residue res:chain) {
				if (res instanceof AaResidue) {
					for (Atom atom:res) {
						double radiusOld = getRadius(((AaResidue) res).getAaType(), atom);
						double radiusNew = getRadiusNew(((AaResidue) res).getAaType(), atom);
						if (Math.abs(radiusOld-radiusNew)>0.0000001) 
							System.err.println("Mismatch for "+res.getLongCode()+" atom "+atom.getCode());
					}
					
				}
			}
		}

	}
}

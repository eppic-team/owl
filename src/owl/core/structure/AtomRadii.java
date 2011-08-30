package owl.core.structure;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

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

	/**
	 * Gets the radius for given amino acid and atom (PDB convention)
	 * @param aa
	 * @param atom
	 * @return
	 */
	public static double getRadius(AminoAcid aa, Atom atom) {
		if (atom.getType().equals(AtomType.H)) return AtomType.H.getRadius();
		String atomCode = atom.getCode();
		if (atomCode.equals("OXT")) atomCode="O";
		return radii.get(aa.getThreeLetterCode()).get(atomCode);
	}
	
	/**
	 * Gets the radius for given nucleotide and atom (PDB convention)
	 * @param nuc
	 * @param atom
	 * @return
	 */
	public static double getRadius(Nucleotide nuc, Atom atom) {
		if (atom.getType().equals(AtomType.H)) return AtomType.H.getRadius();
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
}

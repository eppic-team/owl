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

	private static final InputStream vdwradIS = SymoplibParser.class.getResourceAsStream(RADFILE);
	private static final HashMap<String, HashMap<String,Double>> radii = parseRadFile();

	/**
	 * Gets the radius for given amino acid and atom code (PDB convention)
	 * @param aa
	 * @param atomCode
	 * @return
	 */
	public static double getRadius(AminoAcid aa, Atom atom) {
		if (atom.getType().equals(AtomType.H)) return AtomType.H.getRadius();
		String atomCode = atom.getCode();
		if (atomCode.equals("OXT")) atomCode="O";
		return radii.get(aa.getThreeLetterCode()).get(atomCode);
	}
		
	private static HashMap<String, HashMap<String,Double>> parseRadFile() {
		HashMap<String, HashMap<String,Double>> radii = new HashMap<String, HashMap<String,Double>>();
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(vdwradIS));
			String line;
			HashMap<String,Double> perResRad = null;
			while ((line=br.readLine())!=null) {
				if (line.startsWith("#")) continue;
				
				if (line.startsWith("RESIDUE ATOM")) {
					String res = line.substring(13,16);
					if (res.equals("ASX")) { 
						break; // so that we don't continue reading atoms (for the moment not interested in non-standard aas)
					}
					perResRad = new HashMap<String, Double>();
					radii.put(res, perResRad);
				}
				if (line.startsWith("ATOM")) {
					String atomCode = line.substring(6,9).trim();
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

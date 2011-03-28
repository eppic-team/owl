import java.io.File;

import owl.core.structure.Pdb;
import owl.core.structure.PdbfilePdb;
import owl.core.structure.Residue;


/**
 * Executable class to check the chirality of amino acids of a given PDB file
 * @author duarte
 *
 */
public class checkChirality {

	public static void main(String[] args) throws Exception {
		if (args.length==0) {
			System.err.println("Usage: checkChirality <pdb_file_name>");
			System.exit(1);
		}
		
		File file = new File(args[0]);
		Pdb pdb = new PdbfilePdb(file.getAbsolutePath());
		String[] chains = pdb.getChains();
		
		for (String chain:chains) {
			System.out.println("Chain "+chain+". D-form residues:");
			pdb.load(chain);
			int dcount=0;
			int lcount=0;
			int ucount=0;
			for (int resser:pdb.getAllSortedResSerials()) {
				Residue res = pdb.getResidue(resser);
				Residue.Chirality chir = res.getChirality();
				if (chir == Residue.Chirality.D) {
					//System.out.println(res+" "+chir.getAbbrev());
					dcount++;
				} else if (chir == Residue.Chirality.L) {
					lcount++;
				} else if (chir == Residue.Chirality.U) {
					ucount++;
				}
			}
			System.out.println("L-form: "+lcount+" D-form: "+dcount+" Undetermined: "+ucount);
		}

	}

}

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbasePdb;
import proteinstructure.PdbfilePdb;
import proteinstructure.SecStrucElement;

/**
 * Given a structure and a list of mutations, visualize the mutated sites on the structure
 * and calculate their surface accessibility and secondary structure state.
 * @author stehr
 */
public class mapMutations {

	private static final String NACCESS_EXECUTABLE = "/project/StruPPi/bin/naccess";
	private static final String NACCESS_PARAMETERS = "";
	private static final String DSSP_EXECUTABLE = "/project/StruPPi/Software/dssp/dsspcmbi";
	private static final String DSSP_PARAMETERS = "--";
	private static final double EXPOSURE_CUTOFF = 5.0; // everything above this cutoff is considered exposed
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length < 2) {
			System.out.println(mapMutations.class.getName() + " <pdbcode or filename> <mutated positions>");
			System.exit(1);
		}
		
		Pdb pdb = readStructureOrExit(args[0]);		
		
		int[] mutations = new int[args.length-1];
		for (int i = 1; i < args.length; i++) {
			int pos = Integer.parseInt(args[i]);
			if(pos > pdb.get_length()) {
				System.err.println("Error: Position " + pos + " is bigger than length of protein (" + pdb.get_length() + "). Skipping.");
			} else {
				mutations[i-1] = pos;
			}
		}
		
		try {
			pdb.runNaccess(NACCESS_EXECUTABLE, NACCESS_PARAMETERS);
		} catch (IOException e) {
			System.err.println("Error running NACCESS: " + e.getMessage());
			System.exit(1);
		}
			
		try {
			pdb.runDssp(DSSP_EXECUTABLE, DSSP_PARAMETERS);
		} catch (IOException e) {
			System.err.println("Error running DSSP: " + e.getMessage());
			System.exit(1);
		}
		
		for (int i = 0; i < mutations.length; i++) {
			int pos = mutations[i];
			boolean exposed = pdb.getAllRsaFromResSerial(pos) > EXPOSURE_CUTOFF;
			char ssState = pdb.getSecondaryStructure().getSecStrucElement(pos).getType();
			System.out.printf("Position %d :                     %s\n", pos, pdb.getResTypeFromResSerial(pos));
			System.out.printf("Relative surface accessibility :  %2.0f%% (%s)\n", pdb.getAllRsaFromResSerial(pos), exposed?"exposed":"buried");
			System.out.printf("Secondary structure state:        %3s (%s)\n", ssState, SecStrucElement.getTypeDescription(ssState));
			System.out.println();
		}
	}

	/**
	 * Loads a pdb structure where arg can be a pdbcode+chaincode or a pdb file name.
	 * If something goes wrong, prints an error message and exits.
	 * @param arg a pdbcode+chaincode (e.g. 1tdrB) or a pdb file name
	 * @return the structure object
	 */
	public static Pdb readStructureOrExit(String arg) {
		Pdb pdb = null;
		
		// check if argument is a filename
		File inFile = new File(arg);
		if(inFile.canRead()) {
			System.out.println("Reading file " + arg);
			pdb = new PdbfilePdb(arg);
			try {
				String[] chains = pdb.getChains();
				System.out.println("Loading chain " + chains[0]);
				pdb.load(chains[0]);
			} catch (PdbLoadError e) {
				System.err.println("Error loading file " + arg + ":" + e.getMessage());
			}
		} else {
			// check if argument is a pdb code
			if(arg.length() < 4 || arg.length() > 5) {
				System.err.println(arg + "is neither a valid file name nor a valid pdb code");
				System.exit(1);
			} else {
				String pdbCode = arg.substring(0,4);
				String chainCode = arg.substring(4,5);
				try {
					System.out.println("Loading pdb code " + pdbCode);
					pdb = new PdbasePdb(pdbCode);
					if(chainCode.length() == 0) {
						try {
							chainCode = pdb.getChains()[0];
						} catch (PdbLoadError e) {
							System.err.println("Error loading pdb structure:" + e.getMessage());
							System.exit(1);
						}
					}
					try {
						System.out.println("Loading chain " + chainCode);
						pdb.load(pdb.getChains()[0]);
					} catch (PdbLoadError e) {
						System.err.println("Error loading pdb structure:" + e.getMessage());
						System.exit(1);
					}
				} catch (SQLException e) {
					System.err.println("Database error: " + e.getMessage());
					System.exit(1);
				} catch (PdbCodeNotFoundError e) {
					System.err.println("Pdb code " + pdbCode + " not found in database.");
					System.exit(1);
				}

			}
		}
		return pdb;
	}
	
}

import java.io.*;
import java.sql.SQLException;
import java.util.HashMap;

import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbasePdb;
import proteinstructure.PdbfilePdb;
import runners.NaccessRunner;

/**
 * Loads a structure from a file or from PDBase and saves it as a pdb file where the b-factor column is the
 * relative per-residue solvent accessible surface area (SASA) as calculated by NACCESS. These values can be
 * displayed in e.g. Pymol with the command 'spectrum b'.
 * If a cutoff is given, the b-factor is set to 1 (exposed) or 0 (buried) depending on whether the SASA exceeds
 * the given cutoff or not. Otherwise the SASA value itself (in percent but can sometimes exceed 100) is written.
 * @author stehr
 */
public class writeSASAtoBFactor {

	private static final String NACCESS_EXECUTABLE = "/project/StruPPi/bin/naccess";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length < 2) {
			System.out.println(writeSASAtoBFactor.class.getName() + " <pdbcode or filename> <outfile> [<cutoff for buried/exposed>]");
			System.exit(1);
		}
		String arg = args[0];
		String outFile = args[1];
		double cutoff = -1;
		if(args.length > 2) {
			cutoff = Double.parseDouble(args[2]);
		}

		Pdb pdb = readStructureOrExit(arg);
				
		// now we have a pdb structure
		try {
			NaccessRunner naccRunner = new NaccessRunner(new File(NACCESS_EXECUTABLE), "");
			naccRunner.runNaccess(pdb);
			HashMap<Integer, Double> sasas = pdb.getSurfaceAccessibilities();
			if(sasas == null) {
				System.err.println("Error: no surface accessibility values found.");
				System.exit(1);
			}
			if(cutoff >= 0) {
				System.out.println("Discretizing values");
				sasas = discretizeValues(sasas, cutoff);
			}
			pdb.setBFactorsPerResidue(sasas);
			try {
				System.out.println("Writing " + outFile);
				pdb.writeToPDBFile(outFile);
			} catch (IOException e) {
				System.err.println("Error writing to file " + outFile + ": " + e.getMessage());
				System.exit(1);
			}
		} catch (IOException e) {
			System.err.println("Error running NACCESS: " + e.getMessage());
			System.exit(1);
		}
		
		System.out.println("done.");
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
	
	/**
	 * Transform all values above the threshold to 1 and all below or equal to 0. 
	 * @param map
	 * @param cutoff
	 * @return the map with the new values
	 */
	public static HashMap<Integer, Double> discretizeValues(HashMap<Integer, Double> map, double cutoff) {
		HashMap<Integer, Double> map2 = new HashMap<Integer, Double>();
		for(int key:map.keySet()) {
			if(map.get(key) > cutoff) {
				map2.put(key, 1.0);
			} else {
				map2.put(key, 0.0);
			}
		}
		return map2; 
	}
	
}

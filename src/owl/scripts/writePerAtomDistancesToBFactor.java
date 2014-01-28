package owl.scripts;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.TreeSet;

import owl.core.structure.PdbChain;



public class writePerAtomDistancesToBFactor {

	/**
	 * Loads two structures from file or PDBase and outputs the first structure as a pdb file where the b-factor for
	 * each atom is set to the distance between the atom in structure 1 and the corresponding atom in structure 2.
	 * If there is no such atom in structure 2, the default b-factor 0 is used. These values can be displayed in Pymol
	 * with the command 'spectrum b'.
	 * @author stehr
	 */
	public static void main(String[] args) {
		if(args.length < 2) {
			System.out.println("Usage: " + writePerAtomDistancesToBFactor.class.getName() + " <structure1> <structure2> <outfile>");
			System.out.println("Structures can be given as Pdb+chain code or pdb file name. For pdb files, always the");
			System.out.println("first chain in the file is used.");
			System.exit(1);
		}
		String arg1 = args[0];
		String arg2 = args[1];
		String outFileName = args[2];
		
		if(!new File(arg1).canRead()) {
			System.err.println("Can not read file " + arg1);
		}
		if(!new File(arg2).canRead()) {
			System.err.println("Can not read file " + arg2);
		}
		PdbChain pdb1 = PdbChain.readStructureOrExit(arg1);
		PdbChain pdb2 = PdbChain.readStructureOrExit(arg2);
		
		HashMap<Integer, Double> distances = pdb1.getPerAtomDistances(pdb2);
		pdb1.setBFactorsPerAtom(distances);
		TreeSet<Double> sortedDistances = new TreeSet<Double>();
		sortedDistances.addAll(distances.values());
		double min = sortedDistances.first();
		double max = sortedDistances.last();
		try {
			pdb1.writeToPDBFile(new File(outFileName));
			System.out.println("File " + outFileName + " written.");
			System.out.println("Minimum distance: " + min);
			System.out.println("Maximum distance: " + max);
		} catch (IOException e) {
			System.err.println("Error writing to " + outFileName + ": " + e.getMessage());
		}
	}

}

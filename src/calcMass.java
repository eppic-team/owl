import java.io.File;

import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbChain;


public class calcMass {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		File inFile = new File(args[0]);
		PdbAsymUnit pdb = new PdbAsymUnit(inFile);
		System.out.println("Per chain:");
		for (PdbChain chain:pdb.getAllChains()) {
			System.out.printf("%s %10.2fDa\n",chain.getPdbChainCode(),chain.getMass());
		}
		System.out.println("Total:");
		System.out.printf("%10.2fDa\n",pdb.getMass());

	}

}

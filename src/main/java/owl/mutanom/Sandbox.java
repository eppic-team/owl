package owl.mutanom;

import java.io.File;
import java.io.IOException;

import owl.core.structure.PdbChain;
import owl.core.structure.AaResidue;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.util.FileFormatException;

/**
 * Class for trying out things in the main method.
 * @author stehr
 *
 */
public class Sandbox {

	public static void main(String[] args) throws IOException, PdbLoadException,FileFormatException {
		
		// pdb to cif
		PdbAsymUnit pdb = new PdbAsymUnit(new File("3na3.cif"));
		PdbChain chain = pdb.getChain("A");
		String[] pdbResNums = {"3","208","209","321"};
		for(String pdbResNum:pdbResNums) {
			int cifResNum = chain.getResSerFromPdbResSer(pdbResNum);
			String resType = ((AaResidue)chain.getResidue(cifResNum)).getAaType().toString();
			System.out.printf("pdb=%s cif=%d type=%s\n", pdbResNum, cifResNum, resType);			
		}
				
	}
	
}

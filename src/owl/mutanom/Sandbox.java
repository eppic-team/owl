package owl.mutanom;

import java.sql.SQLException;

import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbasePdb;

/**
 * Class for trying out things in the main method.
 * @author stehr
 *
 */
public class Sandbox {

	public static void main(String[] args) throws SQLException, PdbCodeNotFoundException, PdbLoadException {
		
		// pdb to cif
		Pdb pdb = new PdbasePdb("3na3");
		pdb.load("A");
		String[] pdbResNums = {"3","208","209","321"};
		for(String pdbResNum:pdbResNums) {
			int cifResNum = pdb.getResSerFromPdbResSer(pdbResNum);
			String resType = pdb.getResidue(cifResNum).getAaType().toString();
			System.out.printf("pdb=%s cif=%d type=%s\n", pdbResNum, cifResNum, resType);			
		}
				
	}
	
}

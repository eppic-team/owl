package owl.mutanom;

import java.sql.SQLException;

import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadException;
import owl.core.util.MySQLConnection;

/**
 * Class for trying out things in the main method.
 * @author stehr
 *
 */
public class Sandbox {

	public static void main(String[] args) throws SQLException, PdbCodeNotFoundException, PdbLoadException {
		
		// pdb to cif
		PdbAsymUnit pdb = new PdbAsymUnit("3na3", new MySQLConnection(),"pdbase");
		PdbChain chain = pdb.getChain("A");
		String[] pdbResNums = {"3","208","209","321"};
		for(String pdbResNum:pdbResNums) {
			int cifResNum = chain.getResSerFromPdbResSer(pdbResNum);
			String resType = chain.getResidue(cifResNum).getAaType().toString();
			System.out.printf("pdb=%s cif=%d type=%s\n", pdbResNum, cifResNum, resType);			
		}
				
	}
	
}

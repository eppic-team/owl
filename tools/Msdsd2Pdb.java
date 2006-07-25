package tools;

import java.io.*;

/**
 * Package:		tools
 * Class: 		PyMol
 * Author:		Ioannis Filippis, filippis@molgen.mpg.de
 * Date:		21/03/2006
 *
 * Msdsd2Pdb's static export2File method creates a pdb file by exporting the atom
 * lines directly from the msdsd. This is needed for the visualization of the
 * biological units since pdb files contain the ASUs. Moreover, contact graphs
 * are defined using msdsd-custom fields (like chain code and residue serial)
 * and mapping to pdb fields would be necessary for the graph visualisation 
 * if the original pdb files were preferred.
 * 
 * Notes:
 * - Hetatoms are excluded (pdb_group = "A") and in case of multiple locations of
 * 	 amino acids, only the default location is considered (alt_code = "A") 
 * 	 (VendruscoloM_00_PSFG.pdf currently found in LitNet/incomingPDF/LAST_ROUND/)
 * - The filename is either accessionCode_assemblyId_modelId.pdb (biological unit)
 *   or accessionCode.pdb (asu)
 * - There is also the oportunity to send the atom lines directly to PyMol and 
 * 	 loading the structure without intermediate files. Look at PyMol class and
 *   sendAtomLines method.
 * 
 * Changelog:
 * 21/03/06 first created by IF
 */

public class Msdsd2Pdb {
	
    /**
     * exports to file the atom lines of a model (modelId) of a biological unit (assemblyId) 
     * of a protein (accessionCode) directly from msdsd. The filename is returned.
     * 
     * Notes:
     * - Hetatoms are excluded (pdb_group = "A") and in case of multiple locations of
     * 	 amino acids, only the default location is considered (alt_code = "A") 
     * 	 (VendruscoloM_00_PSFG.pdf currently found in LitNet/incomingPDF/LAST_ROUND/)
     * - The filename is accessionCode_assemblyId_modelId.pdb (biological unit).
     * - The chain_pdb_code is used in the chainID field in the atom line, while the chain_code is used 
     * 	 the segID field (due to its length). Therefore, "segi" and not "chain" must be used in pymol
     * 	 selections.
     * - There are two versions of export2File. One that takes the atomic coordinates from the 
     * 	 partial atom_data tables (needs the table number e.g. 1 for atom_data_1, but is faster), 
     *	 while the other uses the merged table (really slow - should be avoided)
     */ 
	public static String export2File(String accessionCode, int assemblyId, int modelId, String pdbDir) {

		String pdbFileName = accessionCode+"_"+assemblyId+"_"+modelId+".pdb";
		String query = "SELECT CONCAT("+
		    "RPAD(\"ATOM\", 6, \" \"), "+
		    "LPAD(serial, 5, \" \"), "+
		    "\" \", "+
		    "LPAD(chem_atom_name, 4, \" \"), "+
		    "IF(alt_code IS NULL, \" \", alt_code), "+
		    "code_3_letter, "+
		    "\" \", "+
		    "IF(chain_pdb_code IS NULL, \" \", chain_pdb_code), "+
		    "LPAD(residue_serial, 4, \" \"), "+
		    "IF(residue_pdb_insert_code IS NULL, \" \", residue_pdb_insert_code), "+
		    "REPEAT(\" \", 3), "+
		    "LPAD(x, 8, \" \"), "+
		    "LPAD(y, 8, \" \"), "+
		    "LPAD(z, 8, \" \"), "+
		    "LPAD(occupancy, 6, \" \"), "+
		    "REPEAT(\" \", 6), "+
		    "REPEAT(\" \", 6), "+
		    "RPAD(chain_code, 4, \" \") "+
		    ") AS atom_lines FROM msdsd.atom_data WHERE "+
		    "(assembly_id = "+assemblyId+") AND "+
		    "(model_id = "+modelId+") AND "+
		    "((alt_code = \"A\") OR (alt_code IS NULL)) AND "+
		    "(pdb_group = \"A\") "+
		    "ORDER BY chain_code, residue_serial, serial;";
		
		gen(pdbDir+"/"+pdbFileName, query);
	
		return pdbFileName;

    }

    /**
     * exports to file the atom lines of a model (modelId) of a biological unit (assemblyId) 
     * of a protein (accessionCode) directly from msdsd. The filename is returned.
     * 
     * Notes:
     * - Hetatoms are excluded (pdb_group = "A") and in case of multiple locations of
     * 	 amino acids, only the default location is considered (alt_code = "A") 
     * 	 (VendruscoloM_00_PSFG.pdf currently found in LitNet/incomingPDF/LAST_ROUND/)
     * - The filename is accessionCode_assemblyId_modelId.pdb (biological unit).
     * - The chain_pdb_code is used in the chainID field in the atom line, while the chain_code is used 
     * 	 the segID field (due to its length). Therefore, "segi" and not "chain" must be used in pymol
     * 	 selections.
     * - There are two versions of export2File. One that takes the atomic coordinates from the 
     * 	 partial atom_data tables (needs the table number e.g. 1 for atom_data_1, but is faster), 
     *	 while the other uses the merged table (really slow - should be avoided)
     */ 	
    public static String export2File(String accessionCode, int assemblyId, int modelId, int atomDataTblNum, String pdbDir) {

		String pdbFileName = accessionCode+"_"+assemblyId+"_"+modelId+".pdb";
		String query = "SELECT CONCAT("+
		    "RPAD(\"ATOM\", 6, \" \"), "+
		    "LPAD(serial, 5, \" \"), "+
		    "\" \", "+
		    "LPAD(chem_atom_name, 4, \" \"), "+
		    "IF(alt_code IS NULL, \" \", alt_code), "+
		    "code_3_letter, "+
		    "\" \", "+
		    "IF(chain_pdb_code IS NULL, \" \", chain_pdb_code), "+
		    "LPAD(residue_serial, 4, \" \"), "+
		    "IF(residue_pdb_insert_code IS NULL, \" \", residue_pdb_insert_code), "+
		    "REPEAT(\" \", 3), "+
		    "LPAD(x, 8, \" \"), "+
		    "LPAD(y, 8, \" \"), "+
		    "LPAD(z, 8, \" \"), "+
		    "LPAD(occupancy, 6, \" \"), "+
		    "REPEAT(\" \", 6), "+
		    "REPEAT(\" \", 6), "+
		    "RPAD(chain_code, 4, \" \") "+
		    ") AS atom_lines FROM msdsd.atom_data_"+atomDataTblNum+" WHERE "+
		    "(assembly_id = "+assemblyId+") AND "+
		    "(model_id = "+modelId+") AND "+
		    "((alt_code = \"A\") OR (alt_code IS NULL)) AND "+
		    "(pdb_group = \"A\") "+
		    "ORDER BY chain_code, residue_serial, serial;";
		
		gen(pdbDir+"/"+pdbFileName, query);

		return pdbFileName;

    }
    
    /**
     * exports to file the atom lines of the assymetric unit of a protein (accessionCode) 
     * directly from msdsd. The filename is returned.
     * 
     * Notes:
     * - Hetatoms are excluded (pdb_group = "A") and in case of multiple locations of
     * 	 amino acids, only the default location is considered (alt_code = "A") 
     * 	 (VendruscoloM_00_PSFG.pdf currently found in LitNet/incomingPDF/LAST_ROUND/)
     * - The filename is accessionCode.pdb (asu).
     * - The chain_pdb_code is used in the chainID field in the atom line, while the chain_code is used 
     * 	 the segID field (due to its length). Therefore, "segi" and not "chain" must be used in pymol
     * 	 selections.
     * - There are two versions of export2File. One that takes the atomic coordinates from the 
     * 	 partial atom_data tables (needs the table number e.g. 1 for atom_data_1, but is faster), 
     *	 while the other uses the merged table (really slow - should be avoided)
     */ 
	public static String export2File(String accessionCode, String pdbDir) {

		String pdbFileName = accessionCode+".pdb";
		String query = "SELECT CONCAT("+
		    "RPAD(\"ATOM\", 6, \" \"), "+
		    "LPAD(serial, 5, \" \"), "+
		    "\" \", "+
		    "LPAD(chem_atom_name, 4, \" \"), "+
		    "IF(alt_code IS NULL, \" \", alt_code), "+
		    "code_3_letter, "+
		    "\" \", "+
		    "IF(chain_pdb_code IS NULL, \" \", chain_pdb_code), "+
		    "LPAD(residue_serial, 4, \" \"), "+
		    "IF(residue_pdb_insert_code IS NULL, \" \", residue_pdb_insert_code), "+
		    "REPEAT(\" \", 3), "+
		    "LPAD(x, 8, \" \"), "+
		    "LPAD(y, 8, \" \"), "+
		    "LPAD(z, 8, \" \"), "+
		    "LPAD(occupancy, 6, \" \"), "+
		    "REPEAT(\" \", 6), "+
		    "REPEAT(\" \", 6), "+
		    "RPAD(chain_code, 4, \" \") "+
		    ") AS atom_lines FROM msdsd.atom_data WHERE "+
		    "(accession_code = \""+accessionCode+"\") AND "+
		    "(non_assembly_valid = \"Y\") AND "+
		    "((alt_code = \"A\") OR (alt_code IS NULL)) AND "+
		    "(pdb_group = \"A\") "+
		    "ORDER BY chain_code, residue_serial, serial;";
		
		gen(pdbDir+"/"+pdbFileName, query);
	
		return pdbFileName;

    }

    /**
     * exports to file the atom lines of the assymetric unit of a protein (accessionCode) 
     * directly from msdsd. The filename is returned.
     * 
     * Notes:
     * - Hetatoms are excluded (pdb_group = "A") and in case of multiple locations of
     * 	 amino acids, only the default location is considered (alt_code = "A") 
     * 	 (VendruscoloM_00_PSFG.pdf currently found in LitNet/incomingPDF/LAST_ROUND/)
     * - The filename is accessionCode.pdb (asu).
     * - The chain_pdb_code is used in the chainID field in the atom line, while the chain_code is used 
     * 	 the segID field (due to its length). Therefore, "segi" and not "chain" must be used in pymol
     * 	 selections.
     * - There are two versions of export2File. One that takes the atomic coordinates from the 
     * 	 partial atom_data tables (needs the table number e.g. 1 for atom_data_1, but is faster), 
     *	 while the other uses the merged table (really slow - should be avoided)
     */ 
	public static String export2File(String accessionCode, int atomDataTblNum, String pdbDir) {

		String pdbFileName = accessionCode+".pdb";
		String query = "SELECT CONCAT("+
		    "RPAD(\"ATOM\", 6, \" \"), "+
		    "LPAD(serial, 5, \" \"), "+
		    "\" \", "+
		    "LPAD(chem_atom_name, 4, \" \"), "+
		    "IF(alt_code IS NULL, \" \", alt_code), "+
		    "code_3_letter, "+
		    "\" \", "+
		    "IF(chain_pdb_code IS NULL, \" \", chain_pdb_code), "+
		    "LPAD(residue_serial, 4, \" \"), "+
		    "IF(residue_pdb_insert_code IS NULL, \" \", residue_pdb_insert_code), "+
		    "REPEAT(\" \", 3), "+
		    "LPAD(x, 8, \" \"), "+
		    "LPAD(y, 8, \" \"), "+
		    "LPAD(z, 8, \" \"), "+
		    "LPAD(occupancy, 6, \" \"), "+
		    "REPEAT(\" \", 6), "+
		    "REPEAT(\" \", 6), "+
		    "RPAD(chain_code, 4, \" \") "+
		    ") AS atom_lines FROM msdsd.atom_data_"+atomDataTblNum+" WHERE "+
		    "(accession_code = \""+accessionCode+"\") AND "+
		    "(non_assembly_valid = \"Y\") AND "+
		    "((alt_code = \"A\") OR (alt_code IS NULL)) AND "+
		    "(pdb_group = \"A\") "+
		    "ORDER BY chain_code, residue_serial, serial;";
		
		gen(pdbDir+"/"+pdbFileName, query);
	
		return pdbFileName;

    }

	/**
     * creates a temporary sql script called export.sql with the necessary sql query to dump the atom lines,
     * executes the script redirecting the output to a file and deletes the sql script
     */ 	
    private static void gen(String pdbFileName, String query) {
	
        try {

		    File sqlScript = new File("export.sql");
		    PrintWriter scriptOut = new PrintWriter(new FileWriter(sqlScript));
		    scriptOut.println(query);
		    if (scriptOut != null) { scriptOut.close(); }
		    
		    System.out.println(SystemCmd.exec(new String[] {"/bin/sh", "-c", "my_lila < export.sql > "+pdbFileName}));
	
		    sqlScript.delete();

		} catch (Exception e) { 
			System.out.println(e); 
		}
		
	}

} // end of class Msdsd2Pdb

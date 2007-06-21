package tools;

import java.io.*;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

/**
 * Package:		tools
 * Class: 		Msdsd2Pdb
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
 * 	 amino acids, only the default location is considered (graph_alt_code_used = 1) 
 * 	 (VendruscoloM_00_PSFG.pdf currently found in LitNet/incomingPDF/LAST_ROUND/)
 * - There is also the possibility to send the atom lines directly to PyMol and 
 * 	 loading the structure without intermediate files. Look at PyMol class and
 *   sendAtomLines method.
 * 
 * Changelog:
 * 21/03/06 first created by IF
 * 02/03/07 JD, Major changes: adapted to msdsd_00_07_a and using my_msdsd_00_07_a. 
 * 			Using MySQLConnection and file PrintStream instead of shell mysql client for output
 */

public class Msdsd2Pdb {

	
	public static String MSDSDDB="msdsd_00_07_a";
	public static String INFODB="my_msdsd_00_07_a";
	public static String HOST="white";
	public static String PWD="nieve";
	
	
    /**
     * Exports to file in pdb format the atom lines of a model (modelId) of a biological unit (assemblyId) 
     * of a protein (accessionCode) directly from msdsd.
     * 
     * Notes:
     * - Hetatoms are excluded (pdb_group = "A") and in case of multiple locations of
     * 	 amino acids, only the default location is considered (graph_alt_code_used = 1) 
     * 	 (VendruscoloM_00_PSFG.pdf currently found in LitNet/incomingPDF/LAST_ROUND/)
     * - The chain_pdb_code is used in the chainID field in the atom line, while the chain_code is used 
     * 	 the segID field (due to its length). Therefore, "segi" and not "chain" must be used in pymol
     * 	 selections.
     * 
     * @param accessionCode
     * @param assemblyId
     * @param modelId
     * @param pdbFile
     * @param user
     * @throws SQLException 
     */ 
	public static void export2File(String accessionCode, int assemblyId, int modelId, String pdbFile, String user) throws FileNotFoundException, SQLException{
		PrintStream Pdb = new PrintStream(new FileOutputStream(pdbFile));
		MySQLConnection conn;
		conn = new MySQLConnection(HOST,user,PWD,MSDSDDB);

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
		    ") AS atom_lines FROM "+MSDSDDB+".atom_data WHERE "+
		    "(assembly_id = "+assemblyId+") AND "+
		    "(model_id = "+modelId+") AND "+
		    "(graph_alt_code_used = 1) AND "+
		    "(pdb_group = \"A\") "+
		    "ORDER BY chain_code, residue_serial, serial;";
		
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(query);
			while (rsst.next()) {
				Pdb.println(rsst.getString(1));
			}
			stmt.close();
			rsst.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		Pdb.close();
		conn.close();
    }

 
    /**
     * Exports to file in pdb format the atom lines of the assymetric unit of a protein (accessionCode) 
     * (all chains) directly from msdsd. 
     * 
     * Notes:
     * - Hetatoms are excluded (pdb_group = "A") and in case of multiple locations of
     * 	 amino acids, only the default location is considered (graph_alt_code_used = 1) 
     * 	 (VendruscoloM_00_PSFG.pdf currently found in LitNet/incomingPDF/LAST_ROUND/)
     * - The chain_pdb_code is used in the chainID field in the atom line, while the chain_code is used 
     * 	 the segID field (due to its length). Therefore, "segi" and not "chain" must be used in pymol
     * 	 selections.
     * 
     * @param accessionCode
     * @param pdbFile
     * @param user
     * @throws SQLException 
     */ 
	public static void export2File(String accessionCode, String pdbFile, String user) throws FileNotFoundException, SQLException {
		PrintStream Pdb = new PrintStream(new FileOutputStream(pdbFile));
		MySQLConnection conn = new MySQLConnection(HOST,user,PWD,MSDSDDB);

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
		    ") AS atom_lines FROM "+MSDSDDB+".atom_data WHERE "+
		    "(accession_code = \""+accessionCode+"\") AND "+
		    "(non_assembly_valid = \"Y\") AND "+
		    "(graph_alt_code_used = 1) AND "+
		    "(pdb_group = \"A\") "+
		    "ORDER BY chain_code, residue_serial, serial;";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(query);
			while (rsst.next()) {
				Pdb.println(rsst.getString(1));
			}
			stmt.close();
			rsst.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		Pdb.close();
		conn.close();
    }
	
	/**
	 * Exports to file in pdb format the atom coordinates for the assymetric unit of a protein given model_id and chain_id
	 * @param chainId
	 * @param modelId
	 * @param pdbFile
	 * @param user
	 * @return
	 * @throws SQLException 
	 */
	public static void export2File(int chainId, int modelId, String pdbFile, String user) throws FileNotFoundException, SQLException {
		PrintStream Pdb = new PrintStream(new FileOutputStream(pdbFile));
		MySQLConnection conn = new MySQLConnection(HOST,user,PWD,MSDSDDB);

		String query = "SELECT CONCAT("+
		    "RPAD(\"ATOM\", 6, \" \"), "+
		    "LPAD(serial, 5, \" \"), "+
		    "\" \", "+
		    "LPAD(chem_atom_name, 4, \" \"), "+
		    "IF(alt_code IS NULL, \" \", alt_code), "+
		    "code_3_letter, "+
		    "\" \", "+
		    "IF(chain_pdb_code IS NULL, \" \", chain_pdb_code), "+
		    "LPAD(residue_serial, 4, \" \"), "+// check if this is msd or pdb residue serials, do we care?
		    "IF(residue_pdb_insert_code IS NULL, \" \", residue_pdb_insert_code), "+
		    "REPEAT(\" \", 3), "+
		    "LPAD(x, 8, \" \"), "+
		    "LPAD(y, 8, \" \"), "+
		    "LPAD(z, 8, \" \"), "+
		    "LPAD(occupancy, 6, \" \"), "+
		    "REPEAT(\" \", 6), "+
		    "REPEAT(\" \", 6), "+
		    "RPAD(chain_code, 4, \" \") "+
		    ") AS atom_lines FROM "+MSDSDDB+".atom_data WHERE "+
		    "(model_id = "+modelId+") AND "+
		    "(chain_id = "+chainId+") AND "+
		    "(graph_alt_code_used = 1) AND "+
		    "(pdb_group = \"A\") "+
		    "ORDER BY chain_code, residue_serial, serial;";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(query);
			while (rsst.next()) {
				Pdb.println(rsst.getString(1));
			}
			stmt.close();
			rsst.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		Pdb.close();
		conn.close();		
    }

	/**
	 * Exports to file in pdb format the atom coordinates for the assymetric unit of a protein given an accession_code and chain_pdb_code 
	 * (if NMR just the model with model_serial=1)
	 * @param accessionCode
	 * @param chainPdbCode
	 * @param pdbFile
	 * @param user
	 * @return
	 * @throws SQLException 
	 */
	public static void export2File(String accessionCode, String chainPdbCode, String pdbFile, String user) throws FileNotFoundException, SQLException{
		MySQLConnection conn = new MySQLConnection(HOST,user,PWD,MSDSDDB);
		int chainId=0;
		int modelId=0;
		String chainStr="='"+chainPdbCode+"'";
		if (chainPdbCode.equals("NULL")) {
			chainStr="IS NULL";
		}
		String query = "SELECT chain_id, model_id " +
				"FROM "+INFODB+".mmol_chain_info " +
				"WHERE accession_code='"+accessionCode+"' " +
				"AND chain_pdb_code " + chainStr +
				" AND chain_type='C' " +
				"AND asu_chain=1 " +
				"AND model_serial=1;";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(query);
			while (rsst.next()) {
				chainId=rsst.getInt(1);
				modelId=rsst.getInt(2);
			}
			stmt.close();
			rsst.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		conn.close();
		export2File(chainId,modelId,pdbFile,user);
    }




} // end of class Msdsd2Pdb

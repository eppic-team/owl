import java.sql.*;
import java.util.*;
import proteinstructure.*;
import tools.MySQLConnection;


public class calculateGridDensity {

	/*------------------------------ constants ------------------------------*/
	
	// database with a list of pdb codes and chain codes to process
	public static String 			DB_NAME = "pdb_reps";
	public static String 			DB_TABLE = "reps";
	public static String 			DB_COL_PDB = "accession_code";
	public static String 			DB_COL_CHAIN = "chain_pdb_code";	
	
	public static String			PDB_DB = "pdbase";
	public static String			DB_HOST = "white";								
	public static String			DB_USER = getUserName();
	public static String			DB_PWD = "nieve";

	public static String			PDB_CODE = "1tdr";
	public static String			CHAIN_CODE = "B";
	public static double			cutoff = 4.0;
	public static String			edgeType = "Ca";
	
	/*---------------------------- private methods --------------------------*/
	/** 
	 * Get user name from operating system (for use as database username). 
	 * */
	private static String getUserName() {
		String user = null;
		user = System.getProperty("user.name");
		if(user == null) {
			System.err.println("Could not get user name from operating system.");
		}
		return user;
	}
	
	private static void calcDensity(String pdbCode, String chainCode, double cutoff, String egdeType, MySQLConnection conn, Map<Integer, Integer> densityCount) {
		Pdb pdb = null;
		try {
			pdb = new PdbasePdb(pdbCode, chainCode, PDB_DB, conn);
			// add to density count vector
			pdb.calcGridDensity(edgeType, cutoff, densityCount);
			
		} catch (PdbaseInconsistencyError e) {
			System.out.println("Inconsistency in " + pdbCode + chainCode);
		} catch (PdbCodeNotFoundError e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (PdbChainCodeNotFoundError e) {
			e.printStackTrace();
		}
			
	}
	
	public static void printValues(Map<Integer, Integer> v) {
		int atoms = 0;
		for(int size:v.keySet()) {
			System.out.println(size + ": " + v.get(size));
			atoms += size*v.get(size);
		}
		System.out.println("Atoms: " + atoms);
	}
	
	public static void main(String[] args) {
		MySQLConnection conn = null;
		Map<Integer, Integer> densityCount = new TreeMap<Integer,Integer>();		
		String pdbCode, chainCode;
		int numPdbs = 0;
		
		// read command line parameters
		
		// read structures from database
		try{
			conn = new MySQLConnection(DB_HOST, DB_USER, DB_PWD);
		} catch (Exception e) {
			System.err.println("Error opening database connection");
		}
		String query = "SELECT DISTINCT " + DB_COL_PDB + "," + DB_COL_CHAIN + " FROM " + DB_NAME + "." + DB_TABLE + ";" ;
		Statement stmt;
		try {
			stmt = conn.createStatement();
			ResultSet rs = stmt.executeQuery(query);
			while(rs.next()) {
				pdbCode = rs.getString(1);
				chainCode = rs.getString(2);
				
				if(chainCode == null) {
					chainCode = "NULL";
				}
				
				numPdbs++;
				// calculate statistics
				calcDensity(pdbCode, chainCode, cutoff, edgeType, conn, densityCount); // will add to densityCount
				System.out.print(".");
				
				// for each protein write to db: pdb, chain, num_res, volume, max_density
			}
			System.out.println();
			
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		// output results
		System.out.println("Number of structures: " + numPdbs);
		printValues(densityCount);

	}

}

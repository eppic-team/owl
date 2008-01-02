import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import proteinstructure.PdbLoadError;
import proteinstructure.RIGraph;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbasePdb;
import tools.MySQLConnection;


public class benchmarkGraphAlgorithm {
	/*------------------------------ constants ------------------------------*/
	
	// database with a list of pdb codes and chain codes to process
	public static String 			DB_NAME = "pdb_reps";
	public static String 			DB_TABLE = "reps";
	public static String 			DB_COL_PDB = "accession_code";
	public static String 			DB_COL_CHAIN = "chain_pdb_code";	
	
	public static String			PDB_DB = "pdbase_test";
	public static String			DB_HOST = "white";								
	public static String			DB_USER = getUserName();
	public static String			DB_PWD = "nieve";

	public static String			PDB_CODE = "1tdr";
	public static String			CHAIN_CODE = "B";
	public static double			cutoff = 4.2;
	public static String			edgeType = "ALL";
	
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
	
	public static void main(String[] args) throws IOException {
		MySQLConnection conn = null;		
		String pdbCode, chainCode;
		int numPdbs = 0;
		
		// read command line parameters
		
		// read structures from database
		try{
			conn = new MySQLConnection(DB_HOST, DB_USER, DB_PWD);
		} catch (Exception e) {
			System.err.println("Error opening database connection");
		}
		String query = "SELECT DISTINCT " + DB_COL_PDB + "," + DB_COL_CHAIN + " FROM " + DB_NAME + "." + DB_TABLE + " LIMIT 1000;" ;
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
				// get graphs
			
				Pdb pdb = null;
				try {
					pdb = new PdbasePdb(pdbCode, PDB_DB, conn);
					pdb.load(chainCode);
					int length = pdb.get_length();
					int atoms = pdb.getNumAtoms();

					// get graph
					long start = System.currentTimeMillis();
					RIGraph graph = pdb.get_graph(edgeType, cutoff);
					long end = System.currentTimeMillis();

					graph.write_graph_to_file(pdbCode+chainCode+"_"+edgeType+"_"+cutoff);
					
					System.out.print(pdbCode+"_"+chainCode);
					System.out.print("\t"+length+"\t"+atoms);
					System.out.printf("\t%4.3f",(double) (end-start)/1000);
					System.out.println();

					
				} catch (PdbLoadError e) {
					System.out.println("pdb load error in " + pdbCode + chainCode+", specific error: "+e.getMessage());
				} catch (PdbCodeNotFoundError e) {
					e.printStackTrace();
				} catch (SQLException e) {
					e.printStackTrace();
				}					
				
				
				
			
			}
			System.out.println();
			
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		// output results
		System.out.println("Number of structures: " + numPdbs);


	}

}

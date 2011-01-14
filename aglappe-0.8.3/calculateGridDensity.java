import java.io.*;
import java.sql.*;
import java.util.*;
import proteinstructure.*;
import tools.MySQLConnection;

// TODO:
// - how will the distribution look like for randomly drawing from the unit square? (use Ioannis' random graphs?)
// - what is the maximum and expected number of points per square? how does that translate into runtime
// Problems:
// - surface effect (distinguish surface/core)
// - translation of grid (-> kdtrees)
// Space of grid:
// - when assuming spherical proteins and constant density, grid size grows linearly
// - in the worst case of completely linear proteins, grid size growth is cubic
// - for the protein case, do benchmark (to prove linear runtime, estimate runtime constants and space exponent)
// - test whether rotating (by using principle components?) improves runtime/space

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
	public static String			edgeType = "Ca";
	public static double			cutoff_from = 4.0;
	public static double			cutoff_to = 15.0;
	public static double			cutoff_step = 1.0;
	public static int 				limit = 100;
	
	public static String			outFileName = "grid_nbs_pdbreps100";
	
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
	
	public static void printValues(Map<Integer, Integer> v, PrintStream out) {
		int atoms = 0;
		for(int size:v.keySet()) {
			out.println(size + "\t" + v.get(size));
			atoms += size*v.get(size);
		}
		//out.println("Atoms: " + atoms);
	}
	
	
	public static void writeResultToFile(Map<Integer, Integer> v, String baseName, String edgeType, double cutoff) {
		try {
			File outFile = new File(baseName + "_" + edgeType + "_" + cutoff + ".out");
			if(outFile.exists()) {
				outFile.delete();
			}
			outFile.createNewFile();
			printValues(v, new PrintStream(outFile));
			System.out.println("Results written to file " + outFile.getName());
		
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public static void writeResultToDb(Map<Integer, Integer> v) {
		// insert into runs(run_id, edgeType, cutoff, timestamp, proteins, points, cells, avg_pts_per_cell, max_pts_per_cell, avg_pts_per_area, max_pts_per_area)
		// insert into density_distr(run_id, points_per_cell, num_cells)
		// do another run for values per protein (run_id, pdb_id, chain_id, num_res, num_atoms, num_cells, surface_cells, core_cells, surface_area, volume)

	}
	
	public static Map<Integer, Integer> calcDensity(MySQLConnection conn, String edgeType, double cutoff, boolean verbose) {
		Map<Integer, Integer> densityCount = new TreeMap<Integer,Integer>();
		String pdbCode, chainCode;
		int numPdbs = 0;
		
		if(verbose) {
			System.out.print(edgeType + " " + cutoff + ": ");
		}
		
		// read structures from database
		String query = "SELECT DISTINCT " + DB_COL_PDB + "," + DB_COL_CHAIN + " FROM " + DB_NAME + "." + DB_TABLE + " LIMIT " + limit + ";" ;
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
				if(verbose) {
					if(numPdbs %2 == 0) {
						System.out.print('\b');
					} else {
						System.out.print(".");
					}
					if(numPdbs % 500 == 0) System.out.print(numPdbs + " ");
				}
				
				// for each protein write to db: pdb, chain, num_res, volume, max_density
			}
			rs.close();
			stmt.close();
			if(verbose) System.out.println(".");
			
		} catch (SQLException e) {
			e.printStackTrace();
		}		
		return densityCount;
	}
	
	public static void main(String[] args) {
		MySQLConnection conn = null;
		Map<Integer, Integer> densityCount = null;
		double cutoff;

		// opening db connection
		try{
			conn = new MySQLConnection(DB_HOST, DB_USER, DB_PWD);
		} catch (Exception e) {
			System.err.println("Error opening database connection. Exiting");
			System.exit(1);
		}

		for(cutoff = cutoff_from; cutoff <= cutoff_to; cutoff += cutoff_step) {
			// run calculation
			densityCount = calcDensity(conn, edgeType, cutoff, true);

			// output results
			//printValues(densityCount, System.out);
			writeResultToFile(densityCount, outFileName, edgeType, cutoff);
		}		
		System.out.println("Done.");

	}

}
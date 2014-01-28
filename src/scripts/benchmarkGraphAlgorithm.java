package scripts;
import java.io.File;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatException;
import owl.core.util.MySQLConnection;



public class benchmarkGraphAlgorithm {
	/*------------------------------ constants ------------------------------*/
	
	// database with a list of pdb codes and chain codes to process
	public static String 			DB_NAME = "pdb_reps";
	public static String 			DB_TABLE = "reps";
	public static String 			DB_COL_PDB = "accession_code";
	public static String 			DB_COL_CHAIN = "chain_pdb_code";	
	
	public static final String		CIFREPODIR = "/path/to/mmCIF/gz/all/repo/dir";

	public static String			PDB_CODE = "1tdr";
	public static String			CHAIN_CODE = "B";
	public static double			cutoff = 4.2;
	public static String			edgeType = "ALL";
	
	/*---------------------------- private methods --------------------------*/
	
	public static void main(String[] args) throws IOException {
		MySQLConnection conn = null;		
		String pdbCode, chainCode;
		int numPdbs = 0;
		
		// read command line parameters
		
		// read structures from database
		try{
			conn = new MySQLConnection();
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
					chainCode = "A";
				}
				
				numPdbs++;
				// get graphs
			
				PdbChain pdb = null;
				try {
					
					File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".cif");
					cifFile.deleteOnExit();
					PdbAsymUnit.grabCifFile(CIFREPODIR, null, pdbCode, cifFile, false);				
					PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile);
					
					pdb = fullpdb.getChain(chainCode);
					int length = pdb.getStdAaObsLength();
					int atoms = pdb.getNumStdAaHeavyAtoms();

					// get graph
					long start = System.currentTimeMillis();
					RIGraph graph = pdb.getRIGraph(edgeType, cutoff);
					long end = System.currentTimeMillis();

					graph.writeToFile(pdbCode+chainCode+"_"+edgeType+"_"+cutoff);
					
					System.out.print(pdbCode+"_"+chainCode);
					System.out.print("\t"+length+"\t"+atoms);
					System.out.printf("\t%4.3f",(double) (end-start)/1000);
					System.out.println();

					
				} catch (PdbLoadException e) {
					System.out.println("pdb load error in " + pdbCode + chainCode+", specific error: "+e.getMessage());
				} catch (FileFormatException e) {
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

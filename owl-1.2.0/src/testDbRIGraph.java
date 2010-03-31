import gnu.getopt.Getopt;

import java.io.IOException;
import java.sql.SQLException;

//import proteinstructure.CiffilePdb;
import proteinstructure.DbRIGraph;
import proteinstructure.GraphIdNotFoundError;
import tools.MySQLConnection;


public class testDbRIGraph {
	/*------------------------------ constants ------------------------------*/
	
	public static final String			PDB_DB = "pdbase";
	public static final String			DB_HOST = "talyn";								
	public static final String			DB_USER = getUserName();
	public static final String			DB_PWD = "nieve";
	
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
		
		
		String help = "Usage, 3 options:\n" +
				"1)  testDbRIGraph -s <sid> -d <distance_cutoff> -t <contact_type> -r directed -w weighted -F <from_db> -T <to_db> \n" +
				"2)  testDbRIGraph -p <pdb_code> -c <chain_pdb_code> -d <distance_cutoff> -t <contact_type> -r directed -w weighted -F <from_db> -T <to_db> \n" +
				"3)  testDbRIGraph -g <graph_id> -F <from_db> -T <to_db> \n"; 

		String pdbCode = null;
		String pdbChainCode = null;
		String sid = null;
		int graphId = 0;
		String fromDb = null;
		String toDb = null;
		String edgeType = null;
		double cutoff = 0;
		boolean directed = false; 
		boolean weighted = false; 
		
		Getopt g = new Getopt("testDbRIGraph", args, "d:t:r:w:p:c:s:g:F:T:h?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'd':
				cutoff = Double.valueOf(g.getOptarg());
				break;
			case 't':
				edgeType = g.getOptarg();
				break;
			case 'r':
				directed = Boolean.valueOf(g.getOptarg());
				break;
			case 'w':
				weighted = Boolean.valueOf(g.getOptarg());
				break;
			case 'p':
				pdbCode = g.getOptarg();
				break;
			case 'c':
				pdbChainCode = g.getOptarg();
				break;
			case 's':
				sid = g.getOptarg();
				break;
			case 'g':
				graphId = Integer.valueOf(g.getOptarg());
				break;
			case 'F':
				fromDb = g.getOptarg();
				break;
			case 'T':
				toDb = g.getOptarg();
				break;
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		MySQLConnection conn = null;		

		try{
			conn = new MySQLConnection(DB_HOST, DB_USER, DB_PWD);
			conn.setSqlMode("NO_UNSIGNED_SUBTRACTION,TRADITIONAL");
		} catch (Exception e) {
			System.err.println("Error opening database connection. Exiting");
			System.exit(1);
		}
		
		if (fromDb == null || toDb == null) {
			System.exit(0);
		} else if (((pdbCode == null) && (pdbChainCode != null)) || ((pdbCode != null) && (pdbChainCode == null))) {
			System.exit(0);
		} else if ((pdbCode == null) && (sid == null) && (graphId == 0)) {
			System.exit(0);
		}
		
		DbRIGraph graph = null;
		try {
			if (pdbCode != null && pdbChainCode != null) {
				graph = new DbRIGraph(fromDb, conn, pdbCode, pdbChainCode, cutoff, edgeType, directed, weighted);
				graph.setSingleModelsDb(fromDb);
				graph.writeToDbFast(conn,toDb, weighted);
			} else if (sid != null) {
				graph = new DbRIGraph(fromDb, conn, sid, cutoff, edgeType, directed, weighted);
				graph.setSingleModelsDb(fromDb);
				graph.writeToDbFast(conn,toDb, weighted);
			} else if (graphId != 0) {
				graph = new DbRIGraph(fromDb, conn, graphId);
				graph.setSingleModelsDb(fromDb);
				graph.writeToDbFast(conn,toDb);
			}		
		} catch (GraphIdNotFoundError e) {
			System.err.println("Couldn't find such graph!");
		} catch (SQLException e) {
			System.err.println("SQL error for structure "+pdbCode+"_"+pdbChainCode+", error: "+e.getMessage());
		} 
	
		// closing db connection
		conn.close();
	}

}

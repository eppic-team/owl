package owl.casp.benchmarking;

import java.sql.SQLException;

import owl.core.util.MySQLConnection;

/**
 * Scores for a single Casp server model compared to the native structure.
 * Corresponds to one line in the result table.
 * See also: scoreAllServerModels
 * @author stehr
 *
 */
public class ModelScore {

	// model identity
	int	idx;				// just to identify our method as #1
	String target;			// e.g. T0512
	String methodName;		// e.g. BAKER-ROBETTA_TS1
	
	// model scores
	public double gdt;		// GDT_TS
	public double acc;		// accuracy of contact map compared to native
	public double cov;		// coverage of contact map compared to native
	public double acc12;	// acc for contacts with seq.sep. >= 12
	public double cov12; 	// cov for contacts with seq.sep. >= 12
	public double acc24; 	// acc for contacts with seq.sep. >= 24
	public double cov24; 	// cov for contacts with seq.sep. >= 24
	public int cts;			// total number of contacts
	public int cts12;		// number of contacts with seq.sep >= 12
	public int cts24;		// number of contacts with seq.sep >= 24
	//public float accL5; 	// acc for the L/5 highest confidence contacts
	//public float accL10; 	// acc for the L/10 highest confidence contacts
	
	// optionally store ranks
	public int rankGdt;		// rank of this model by gdt_ts score
	public int rankCm;		// rank of this model by (acc+cov)/2
	public int rankCm12;	// rank if this model by acc/cov LR(>=12)
	public int rankCm24;	// rank of this model by acc/cov LR(>=24)
	
	/*----------------------------- constructors ----------------------------*/
	
	public ModelScore(String target, String modelName) {
		this.target = target;
		this.methodName = modelName;
	}

	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Writes this model score to database.
	 * @param conn an active db connection object
	 * @param tableName the target database and table name
	 * @throws SQLException if something went wrong
	 */
	public void writeToDb(MySQLConnection conn, String tableName) throws SQLException {
		String query = String.format("INSERT INTO " + tableName + "(target, model, gdt_ts, acc, cov, acc12, cov12, acc24, cov24) VALUES(\"%s\", \"%s\", %f, %f, %f, %f, %f, %f, %f);", target, methodName, gdt, acc, cov, acc12, cov12, acc24, cov24);
		//	System.out.println(query); // DEBUG
		conn.executeSql(query);
	}
	
	/**
	 * Writes this model score to database including the rank columns.
	 * @param conn an active db connection object
	 * @param tableName the target database and table name
	 * @throws SQLException if something went wrong
	 */
	public void writeToDbWithRanks(MySQLConnection conn, String tableName) throws SQLException {
		String query = String.format("INSERT INTO " + tableName + "(idx, target, method, gdt_ts, acc, cov, acc12, cov12, acc24, cov24, rank_gdt, rank_cm, rank_cm_12, rank_cm_24) VALUES(%d, \"%s\", \"%s\", %f, %f, %f, %f, %f, %f, %f, %d, %d, %d, %d);", idx, target, methodName, gdt, acc, cov, acc12, cov12, acc24, cov24, rankGdt, rankCm, rankCm12, rankCm24);
		//	System.out.println(query); // DEBUG
		conn.executeSql(query);
	}
	
	/**
	 * Updates an existing score table with the number of edges for each model
	 * @param conn an active db connection object
	 * @param tableName the target database and table name
	 * @throws SQLException if something went wrong
	 */	
	public void updateNumContacts(MySQLConnection conn, String tableName) throws SQLException {
		String query = String.format("UPDATE " + tableName + " set cts=%d, cts12=%d, cts24=%d WHERE target='%s' AND method='%s'", cts, cts12, cts24, target, methodName);
		// System.out.println(query); // DEBUG
		conn.executeSql(query);
	}
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Creates the table to hold the scoring results (if not exists).
	 * @param conn an active db connection object
	 * @param tableName the name of the target database and table
	 * @throws SQLException if something went wrong
	 */
	public static void createDbTable(MySQLConnection conn, String tableName) throws SQLException {
		
		//TODO: Not safe against SQL injection!
		String query = "CREATE TABLE IF NOT EXISTS " + tableName + " (" +
				"idx int(5) auto_increment key," +
				"target char(5)," +
				"model varchar(50)," +
				"gdt_ts float," +
				"acc float," +
				"cov float," +
				"acc12 float," +
				"cov12 float," +
				"acc24 float," +
				"cov24 float," +
				//"acc_l5 float," +
				//"acc_l10 float" +
				"rank_gdt int," +
				"rank_cm int," +
				"rank_cm_12 int," +
				"rank_cm_24 int," +
				"cts int," +
				"cts12 int," +
				"cts24 int" + 
				");";
		
		conn.executeSql(query);
	}
	
	/**
	 * Creates the table to hold the scoring results (if not exists).
	 * @param conn an active db connection object
	 * @param tableName the name of the target database and table
	 * @throws SQLException if something went wrong
	 */
	public static void createDbTableWithRanks(MySQLConnection conn, String tableName) throws SQLException {
		
		//TODO: Not safe against SQL injection!
		String query = "CREATE TABLE IF NOT EXISTS " + tableName + " (" +
				"idx int(5)," +
				"target char(5)," +
				"method varchar(50)," +
				"gdt_ts float," +
				"acc float," +
				"cov float," +
				"acc12 float," +
				"cov12 float," +
				"acc24 float," +
				"cov24 float," +
				"rank_gdt int," +
				"rank_cm int," +
				"rank_cm_12 int," +
				"rank_cm_24 int" +
				");";
		
		conn.executeSql(query);
	}
}

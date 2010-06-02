package examples;
import owl.core.util.MySQLConnection;
import java.sql.*;


/* 
 * HelloWorld program for determining the rank 
 *  
 *  */

public class helloRank {
	static int maxRank=21,  VL=2; // Verbosity Level 
	static String backgrndDB="cullpdb_90", host="localhost", port="3306", user="lappe", password="apple";
	static String verzeichnis, datei, pdb, chain="A", c_type, cutoff, generation, individual ;
	static MySQLConnection conn;
	

	private static int getRank( String nbs, String centRes) throws SQLException {
		String sql, res; 
		Statement stmt;  
		ResultSet rsst;
		int counter=0, c=0, lastc=0, rank=0, sRank=0; 
		boolean seenCentRes = false; 

		if (VL>=2) {
			System.out.println("getCountRank for ");
			System.out.println("nbs    : "+nbs); 
			System.out.println("centRes: "+centRes);
		}

		sql = "select res, sum(c) as t from "+backgrndDB+".nb_equals group by res order by t DESC;";
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		if (VL>=2) System.out.println("###\tres\ttotal t");
		sRank = maxRank; 
		while (rsst.next()) {	
			counter++; 
			res = rsst.getString(1); // 1st column -- res
			c = rsst.getInt( 2); // 2nds column : count/residue
			if (VL>=2) System.out.print(counter+"\t"+res+"\t"+c);

			if ((c == lastc) && (lastc>0) && seenCentRes) { // tie 
				if (VL>=2) System.out.print(" <-- TIE!");
				rank = counter; 		
			} // end if 
			if (res.equals(centRes)) { 
				if (VL>=2) System.out.print(" <== " + centRes);
				seenCentRes = true;
				rank = counter; 
			}
			if (VL>=2) System.out.println(".");
		} // end while 
		if (VL>=2) System.out.println("=> rank "+rank); 
		rsst.close(); 
		stmt.close(); 

		if (rank==0) sRank = maxRank;
		else sRank = rank; 
		return sRank; 
	} // end of getCountRank 


	/**
	 * @param args
	 */
	public static void main(String[] args) {
	
		try { // try to load the driver
			Class.forName("com.mysql.jdbc.Driver").newInstance();
			try { // try to connect
				conn = new MySQLConnection("localhost",user,"apple", "GATSE");
				// conn = new MySQLconnection(host+port,user,password,backgrndDB);
				
				getRank( "%P%R%T%x%W%", "L"); 
				
				conn.close(); // closing the Database connection

			} catch (SQLException E) {
				System.out.println("SQLException: " + E.getMessage());
				System.out.println("SQLState:     " + E.getSQLState());
				System.out.println("VendorError:  " + E.getErrorCode());
			} // end try/catch connection

		} // end try load Driver

		catch (Exception E) {
			System.err.println("Unable to load driver.");
			E.printStackTrace();
		} // end catch load driver

	} // end main

} // end class 

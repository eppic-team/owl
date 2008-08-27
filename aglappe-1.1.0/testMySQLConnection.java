import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

public class testMySQLConnection {

	/**
	 * "Hello World" for a MySQLConnection to white
	 * 
	 * @author duarte
	 */
	
	static String user = "duarte"	; // change user name!!
	
	public static void main(String[] args){
				
		MySQLConnection conn = null;
		try {
			conn = new MySQLConnection("white",user,"nieve","newmsdgraph");
		} catch (SQLException e1) {
			e1.printStackTrace();
			System.err.println("Can't connect to the db. Exiting");
			System.exit(1);
		}
		
		try {
			String sql = "SELECT num,res FROM nodes limit 3;";
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			while (rsst.next()) {				
				int num = rsst.getInt(1); // 1st column -- num
				String res = rsst.getString(2); // 2nd column -- res
				System.out.println("serial number: "+num+", residue type: "+res);
			}
		} catch (SQLException e) {
			e.printStackTrace();
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
		}

	}

}

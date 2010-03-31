import tools.MySQLConnection;

import java.sql.SQLException;
import java.sql.Statement;
import java.sql.ResultSet;

/**
 * "Hello World" for a MySQLConnection to our MySQL server (talyn)
 * 
 * @author duarte
 */
public class testMySQLConnection {

	// we throw the SQLException, the program will terminate upon a SQLException and print the stack trace
	public static void main(String[] args) throws SQLException{
				
		// using the empty MySQLConnection constructor the connection settings are read from the ~/.my.cnf file
		MySQLConnection conn = new MySQLConnection();

		String sql = "SELECT num,res FROM cullpdb_20.single_model_node limit 3";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		while (rsst.next()) {				
			int num = rsst.getInt(1); // 1st column -- num
			String res = rsst.getString(2); // 2nd column -- res
			System.out.println("serial number: "+num+", residue type: "+res);
		}
		
		conn.close();
	}

}

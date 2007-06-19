package tools;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

public class MySQLConnectionCheck extends Thread {

	/**
	 * Class to check wheter a MySQL server is alive and contains a certain database, 
	 * to be run in different threads so that it can be used to check all nodes in cluster parallely
	 * Once thread is executed we can call the connTestPassed method which returns a boolean that
	 * tells wheter the check was passed or not.
	 * 
	 * @author duarte
	 */
	
	String dbServer;
	String dbUserName;
	String dbPassword;
	String dbName;
	boolean success;
	
	public MySQLConnectionCheck(String dbServer,String dbUserName, String dbPassword, String dbName){
		this.dbServer=dbServer;
		this.dbUserName=dbUserName;
		this.dbPassword=dbPassword;
		this.dbName=dbName;
		this.success=true;
	}

	@Override
	public void run() {
		MySQLConnection.loadMySQLDriver();		
		String connStr="jdbc:mysql://"+dbServer+"/"+dbName;
		Connection tryConn = null;
		try {
			tryConn = DriverManager.getConnection(connStr, dbUserName, dbPassword);
		}
		catch (SQLException e){
			success = false;
		}
		if (tryConn!=null){
			try {
				tryConn.close();
			} catch (SQLException e) {
				e.printStackTrace();
				System.err.println("Couldn't close the test connection.");
			}
		}
	}
	
	public boolean connTestPassed(){
		return success;		
	}

}

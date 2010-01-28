package tools; 

import java.sql.*;
/**
 * Package:		tools
 * Class: 		utils4DB
 * Author:		Ioannis Filippis, filippis@molgen.mpg.de
 * Date:		23/01/2006
 * 
 * utils4DB contains static methods to facilitate information retrieval from database.
 * For example execute an update Query, get the range of a numeric field from a table,
 * get information for the size of a database
 * 
 * Changelog:
 * 04/04/06 modified by IF - execUpdateQuery method added and code in comments was removed
 * 21/03/06 modified by IF - getRange methods return double array instead of float 
 * 23/01/06 first created by IF
 */
public class utils4DB {

    private utils4DB() {};
    
    public static void execUpdateQuery(Connection C, String query) {

    	Statement S;

    	try {
    	    S = C.createStatement();
    	    S.executeUpdate(query);
    	    S.close();
    	} // end try
    	catch (SQLException E) {
    	    System.out.println("SQLException: " + E.getMessage());
    	    System.out.println("SQLState:     " + E.getSQLState());
    	    System.out.println("VendorError:  " + E.getErrorCode());
    	} // end catch
    	    	
    }
    
    public static double[] getRange(Connection C, String table, String column) {

    	String query = "";
    	Statement S;
    	ResultSet R;

    	double[] range = new double[2];

    	try { 
    	    query = "SELECT MIN("+column+"), MAX("+column+") FROM "+table+";";
    	    S = C.createStatement();
    	    R = S.executeQuery(query);
    	
    	    if (R.next()) {
	    		range[0] = R.getDouble(1);
	    		range[1] = R.getDouble(2);
    	    } 

    	    R.close();
    	    S.close();

    	} // end try
    	catch (SQLException E) {
    	    System.out.println("SQLException: " + E.getMessage());
    	    System.out.println("SQLState:     " + E.getSQLState());
    	    System.out.println("VendorError:  " + E.getErrorCode());
    	} // end catch
    	
    	return range;
    	
    }

    public static double[] getRange(Connection C, String table, String column, String selection) {

    	String query = "";
    	Statement S;
    	ResultSet R;

    	double[] range = new double[2];

    	try { 
    	    query = "SELECT MIN("+column+"), MAX("+column+") FROM "+table+" WHERE ("+selection+");";
    	    S = C.createStatement();
    	    R = S.executeQuery(query);
    	
    	    if (R.next()) {
	    		range[0] = R.getDouble(1);
	    		range[1] = R.getDouble(2);
    	    } 

    	    R.close();
    	    S.close();

    	} // end try
    	catch (SQLException E) {
    	    System.out.println("SQLException: " + E.getMessage());
    	    System.out.println("SQLState:     " + E.getSQLState());
    	    System.out.println("VendorError:  " + E.getErrorCode());
    	} // end catch
    	
    	return range;
    	
    }
    
    public static void getDbSizeInfo(String connFile, String DB) {

		Connection myConnection;
		mySQLConnect SQLC;
	
		double data = 0, index = 0, table_data = 0, table_index = 0, GB = Math.pow(2, 30);
		String Query = null, table = null;
		Statement Stmt = null;
		ResultSet RS = null;
	
		SQLC = new mySQLConnect();
		SQLC.readConnectionFile(connFile);
		myConnection = SQLC.openConnection();
	
		try {
		    Query = "SHOW TABLE STATUS FROM "+DB;
		    Stmt = myConnection.createStatement();
		    RS = Stmt.executeQuery(Query);
		    while (RS.next()) {
				table = RS.getString("Name");
				table_data = RS.getDouble("Data_length");
				table_index = RS.getDouble("Index_length");
				data += RS.getDouble("Data_length");
				index += RS.getDouble("Index_length");
		
				System.out.println("Table "+table+"##data:"+table_data+", index:"+table_index);
		    }
		    RS.close();
		    Stmt.close();
	
		    System.out.println("Database "+DB+" needs "+((data+index)/GB)+ " GB (data:"+(data/GB)+", index:"+(index/GB)+").");
		} 
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end try/catch connection
	
		SQLC.closeConnection(myConnection);

    }
    

}

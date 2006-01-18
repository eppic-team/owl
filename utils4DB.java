package tools.trunk; 

import java.sql.*;
import java.io.*;
import java.util.*;
import java.lang.Math.*;

public class utils4DB {

    private utils4DB() {};

    public static float[] getRange(Connection C, String table, String column) {

	String query = "";
	Statement S;
	ResultSet R;

	float[] range = new float[2];

	try { 
	    query = "SELECT MIN("+column+"), MAX("+column+") FROM "+table+";";
	    S = C.createStatement();
	    R = S.executeQuery(query);
	
	    if (R.next()) {
		range[0] = R.getFloat(1);
		range[1] = R.getFloat(2);
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

    public static float[] getRange(Connection C, String table, String column, String selection) {

	String query = "";
	Statement S;
	ResultSet R;

	float[] range = new float[2];

	try { 
	    query = "SELECT MIN("+column+"), MAX("+column+") FROM "+table+" WHERE ("+selection+");";
	    S = C.createStatement();
	    R = S.executeQuery(query);
	
	    if (R.next()) {
		range[0] = R.getFloat(1);
		range[1] = R.getFloat(2);
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

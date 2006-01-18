package tools.trunk; 

import java.sql.*;
import java.io.*;
import java.util.*;

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

}

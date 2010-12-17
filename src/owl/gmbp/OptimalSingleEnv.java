package owl.gmbp;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Vector;

import owl.core.util.MySQLConnection;

public class OptimalSingleEnv {
	

    private MySQLConnection conn;
	
    // variables for queries
	private String host; // = "talyn";
	private String username; // = "vehlow";
	private String password; // = "nieve";
	private String db; // = "mw";
	
	private String fullnbs;
	private char iRes='A';	

	private Vector<String[]> optNBHStrings;
	
	public OptimalSingleEnv(String nbhString, char ires, String host, String user, String pwd) throws SQLException{
		this.fullnbs = nbhString;
		this.iRes = ires;
		this.host = host;
		this.username = user;
		this.password = pwd;

//		conn = new MySQLConnection(this.host,this.username,this.password,this.db);
	}
	
	public OptimalSingleEnv(String nbhString, char ires, String host, String user, String pwd, String db) throws SQLException{
		this.fullnbs = nbhString;
		this.iRes = ires;
		this.host = host;
		this.username = user;
		this.password = pwd;
		this.db = db;

//		conn = new MySQLConnection(this.host,this.username,this.password,this.db);
	}
	
	public OptimalSingleEnv(String nbhString, char ires) throws SQLException{
		this.fullnbs = nbhString;
		this.iRes = ires;

//		conn = new MySQLConnection(this.host,this.username,this.password,this.db);
	}
	
	public void run() throws SQLException {
		// limit by where locate<nullrank, 2b determined here in advance of the sql query 
	    Statement stmt = conn.createStatement();
	    this.optNBHStrings = new Vector<String[]>();
	    String query = "select nbstring, rvector, support, locate('"+this.iRes+"', rvector) from rvecs10 order by locate('"+this.iRes+"', rvector);"; 
//	    System.out.println(query);
	    ResultSet nbstrings = stmt.executeQuery(query); 
	    String nbstring ="", rvector="", support="", loc=""; 
//	    System.out.println("fullnbs: " + fullnbs); 
	    while (nbstrings.next()) {
	    	nbstring = nbstrings.getString(1); 
	    	rvector = nbstrings.getString(2);
	    	support = nbstrings.getString(3); //.getInt(3)
	    	loc = nbstrings.getString(4); //.getInt(4)
	    	if (isAinB(nbstring, fullnbs)) {
		    	String[] string = {nbstring, rvector, support, loc};
		    	this.optNBHStrings.add(string);
//	    		System.out.println(nbstring+"\t"+rvector+"\t"+support+"\t"+loc);  //( " !!!");
	    	} else {
//	    		System.out.println( " -");
	    	} // end if A is in B 
	    } // end while nbstring 
	    nbstrings.close();
		stmt.close(); // closing the Database connection
		
	}

	private static boolean isAinB( String part, String full) { 
		boolean isAregExp=true; 
		int l = part.length(), lastOcc=0, occ=0;
		String p="";
		for ( int i=0; i<l; i++) {
			p=part.substring(i, i+1); 
			occ =  full.indexOf( p, lastOcc+1);
			if (occ <= 0) {
				isAregExp=false;
			} else lastOcc=occ; 
//			System.out.print(" "+lastOcc+"+"+p); 
		} // next 
		return isAregExp; 
	} // end isAinB 
	
	public Vector<String[]> getOptNBHStrings(){
		return this.optNBHStrings;
	}	

	public void setFullNBHString(String fullnbs) {
		this.fullnbs = fullnbs;
	}

	public void setiRes(char iRes) {
		this.iRes = iRes;
	}
	
	public void setDBaccess(String dbUSER, String dbPWD, String dbHOST, String dbNAME) throws SQLException {
		this.host = dbHOST;
		this.username = dbUSER;
		this.password = dbPWD;
		this.db = dbNAME;
		conn = new MySQLConnection(this.host,this.username,this.password,this.db);
	}
	
}

package tools; 

import java.sql.*;
import java.io.*;
import java.util.*;

public class mySQLConnect 
{
    private String host="", port="", user="", password="", dbname=""; 
    private String driverName = "com.mysql.jdbc.Driver"; 

    public void setDriverName( String xDriver) { driverName = xDriver; }
    public String getDriverName() { return driverName; }

    public void setHost( String xHost) { host = xHost; }
    public String getHost() { return host; }
    
    public void setPort( String xPort) { port = xPort; }
    public String getPort() { return port; }

    public void setUser( String xUser) { user = xUser; }
    public String getUser() { return user; }

    public void setPassword( String xPW) { password = xPW; }
    public String getPassword() { return password; }

    public void setDB( String xDB) { dbname = xDB; }
    public String getDB() { return dbname; }

    public mySQLConnect() {

	loadSQLDriver();
	
    }

    /** this function loads the SQL driver specified by driverName */ 
    public void loadSQLDriver() {
	try {
	    // The newInstance() call is a work around for some broken Java implementations
	    //System.out.print("Loading "+driverName +" ... "); 
	    Class.forName(driverName).newInstance(); 
	    //System.out.println("done."); 
	} 
	catch (Exception E) {
	    System.out.println("\nUnable to load SQL-driver "+driverName);
	    E.printStackTrace();
	} // end catch 
    }// end loadSQLDriver 
    
    /** closes the current connection to the SQL-database */ 
    public void closeConnection( Connection C ) {
	try {
	    C.close(); 
	} catch (SQLException E) {
	    System.out.println("SQLException:\t" + E.getMessage());
	    System.out.println("SQLState:\t" + E.getSQLState());
	    System.out.println("VendorError: \t" + E.getErrorCode());
	} // end try/catch connection 
    } // end closeConnection 

    public Connection openConnection() { 
	// open a connection with the parameters retrieved to the database <dbname> 
	Connection Conn=null;
	String connector = "";
	try {
	    Class.forName( driverName).newInstance(); 	     
	    try {
		// System.out.println("// open Connection()");
		// System.out.println("host "+host);
		// System.out.println("port "+port);
		// System.out.println("dbname "+dbname);
		// System.out.println("user "+user);
		// System.out.println("password "+password);
		
		if (!password.equals("")) {
		    connector = "jdbc:mysql://"+host+port+"/"+dbname+"?user="+user+"&password="+password+"&allowMultiQueries=true";
		} else {
		    connector = "jdbc:mysql://"+host+port+"/"+dbname+"?user="+user+"&allowMultiQueries=true";
		}
		//System.out.println("Opening a connection to "+connector); 
		Conn = DriverManager.getConnection( connector);
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
	return Conn; 
    } // end openConnection()


    public void readConnectionFile( String connFile) {
	// reads the values of the connFile into the static variables; 
	String homedir = System.getProperty("user.home"); 
	// System.out.println("HOME@ "+homedir); // HOME=/home/lappe
	if( connFile.length()==0) { // no file was specified
	    connFile=homedir+"/.my.cnf"; // assume default configuration file 
	} // end if
	// else the location of the connection file was given 

	//System.out.println("Reading settings from file "+connFile);
	// Open the configuration file
	FileReader theFile = null;
	BufferedReader fileIn = null;
	StringTokenizer str;
	String item, oneLine; 
	
	// list the entries in the file and decompose them 
	try {
	    File inputFile = new File(connFile);
	    theFile = new FileReader(inputFile); // open the File
	    fileIn = new BufferedReader( theFile); // open BufferedReader 
	    while ((oneLine = fileIn.readLine() ) != null ) {
		// Write the line at hand to stdout, just for testing purposes 
		// System.out.println("["+oneLine+"]");
		// Construct a stringTokenizer for the line that we read with : delimited
		str = new StringTokenizer( oneLine, " ="); // true sets returnDelimiters flag 
		while ( str.hasMoreTokens()) {
		    item = str.nextToken();
		    // System.out.println("item:"+item);
		    if( item.equals("host")) {
			host=str.nextToken();
			// System.out.println("host:"+host);
			break; 
		    } // end if host 
		    if( item.equals("port")) {
			port=":"+str.nextToken();
			// System.out.println("port:"+port);
			break; 
		    } // end if port
		    if( item.equals("user")) {
			user=str.nextToken();
			// System.out.println("user:"+user);
			break; 
		    } // end if password 
		    if( item.equals("password")) {
			password=str.nextToken();
			// System.out.println("password:"+password);
			break; 
		    } // end if password
		    if( item.equals("database")) {
			dbname=str.nextToken();
			// System.out.println("database:"+dbname);
			break; 
		    } // end if password 
		    
		} // next token in this line 
	    } // next line in the file 
	} // end try opening the file 
	catch ( Exception e ) { System.out.println( e); }  
	
	try { // closing the file
	    if( fileIn != null) fileIn.close();
	    if( theFile != null) theFile.close();
	} catch ( Exception e ) { System.out.println( e); }
	
    } // end class readConnectionFile

} // end class mySQLconnect 


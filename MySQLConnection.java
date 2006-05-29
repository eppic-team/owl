package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.StringTokenizer;
import java.io.IOException;
import java.util.ArrayList;

public class MySQLConnection {

	/*--------------------- constants -----------------------*/
	
    // -- constants for database connection --
    static final String    HOST = 			"white";	
    static final String    USER = 		    "";
    static final String    PASSWORD = 		"nieve";
	
	/*------------------- member variables --------------------*/
	
	public Connection conn; 
	private String host;
	private String user;
	private String password=PASSWORD;
	private String port;
	private String dbname;
	
	/*-------------------- constructors -----------------------*/

	/** 
	 * Connect to database using the given server, user and password
	 * @param dbServer
	 * @param dbUserName
	 * @param dbPassword
	 */
	public MySQLConnection(String dbServer, String dbUserName, String dbPassword) {
		loadMySQLDriver();
		host=dbServer;
		user=dbUserName;
		password=dbPassword;
		port="";
		dbname="";		
		String connStr="jdbc:mysql://"+host+port+"/"+dbname;
		try {
			conn = DriverManager.getConnection(connStr, user, password);
		} catch (SQLException e) {
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
			System.err.println("VendorError:  " + e.getErrorCode());
			e.printStackTrace();
		} // end try/catch connection 		 		 
	}

	/**
	 * Connect to database using the given server, user, password and dbname.
	 * Please always use this constructor in preference rather than constructing without specifying a database
	 * @param dbServer
	 * @param dbUserName
	 * @param dbPassword
	 * @param dbName 
	 */
	public MySQLConnection(String dbServer, String dbUserName, String dbPassword, String dbName) {
		loadMySQLDriver(); 
		host=dbServer;
		user=dbUserName;
		password=dbPassword;
		port="";
		dbname=dbName;		
		String connStr="jdbc:mysql://"+host+port+"/"+dbname;
		try {
			conn = DriverManager.getConnection(connStr, user, password);
		} catch (SQLException e) {
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
			System.err.println("VendorError:  " + e.getErrorCode());
			e.printStackTrace();
		} // end try/catch connection 		 		 		
	}

	/**
	 * Connect to database using the given server, user, password, dbname and port
	 * Only needed if mysql server uses a port different from the standard 3306
	 * @param dbServer
	 * @param dbUserName
	 * @param dbPassword
	 * @param dbName
	 * @param portNum
	 */
	public MySQLConnection(String dbServer, String dbUserName, String dbPassword, String dbName, int portNum) {
		loadMySQLDriver();
		host=dbServer;
		user=dbUserName;
		password=dbPassword;
		port=":"+portNum;
		dbname=dbName;		
		String connStr="jdbc:mysql://"+host+port+"/"+dbname;
		try {
			conn = DriverManager.getConnection(connStr, user, password);
		} catch (SQLException e) {
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
			System.err.println("VendorError:  " + e.getErrorCode());
			e.printStackTrace();
		} // end try/catch connection 		 		 
	}
		
	/**
	 * Connect to database giving a connection file
	 * @param connFile the connection file's name
	 */
	public MySQLConnection(String connFile) {		
		loadMySQLDriver();
		readConnectionFile(connFile);
		String connStr="jdbc:mysql://"+host+port+"/"+dbname;
		try {
			conn = DriverManager.getConnection(connStr, user, password);
		} catch (SQLException e) {
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
			System.err.println("VendorError:  " + e.getErrorCode());
			e.printStackTrace();
		} // end try/catch connection 		 		 
	}
	
	
	/*---------------------- methods -------------------------*/
	
	private void loadMySQLDriver() {
		try {
			Class.forName("com.mysql.jdbc.Driver").newInstance();
		}
		catch(Exception e) {
			e.printStackTrace();
			System.err.println("An exception occurred while loading the mysql jdbc driver, exiting.");
			System.exit(1);
		}
	}

	/**
	 * Used in the constructor that gets a connFile as argument. To read the connection parameters from a connection file.
	 * @param connFile
	 */
	private void readConnectionFile(String connFile) {
		// reads the values of the connFile into the static variables; 
		String homedir = System.getProperty("user.home"); 
		if (connFile.length()==0) { // no file was specified
			connFile=homedir+"/.my.cnf"; // assume default configuration file 
		}		
		// else the location of the connection file was given 		
		// Open the configuration file		
		BufferedReader fileIn = null;
		StringTokenizer str;
		String item, oneLine;
		// to control if the minimum mandatory 4 parameters are given in file
		int paramCount=0; 
		// setting default blank values for port and dbname, they are set to blank unless fields specified in file
		port="";
		dbname="";
		// list the entries in the file and decompose them 
		try {
			fileIn = new BufferedReader(new FileReader(new File(connFile))); // open BufferedReader to file connFile 
			while ((oneLine = fileIn.readLine()) != null ) {
				// Construct a stringTokenizer for the line that we read with : delimited
				str = new StringTokenizer(oneLine, "="); // true sets returnDelimiters flag 
				while ( str.hasMoreTokens()) {
					item = str.nextToken();
					if( item.equals("host")) { // mandatory parameter
						host=str.nextToken();
						paramCount++;
						break; 
					} // end if host 
					if( item.equals("port")) { // optional parameter
						port=":"+str.nextToken();
						break; 
					} // end if port
					if( item.equals("user")) { // mandatory parameter
						user=str.nextToken();
						paramCount++;
						break; 
					} // end if password 
					if( item.equals("password")) { // mandatory parameter
						password=str.nextToken();
						paramCount++;
						break; 
					} // end if password
					if( item.equals("database")) { // mandatory parameter
						dbname=str.nextToken();
						paramCount++;
						break; 
					} // end if password 					
				} // next token in this line 
			} // next line in the file
			if (paramCount<4){
				System.err.println("Not all mandatory parameters are given in connection file "+connFile+". Can't connect to mysql server, exiting.");
				System.exit(1);
			}			
		} 
		catch (IOException e) {
			System.err.println("Couldn't open file "+connFile);			
			e.printStackTrace();
			System.exit(1);
		}  
		
		try { // closing the file
			if (fileIn != null) fileIn.close();			
		} catch (IOException e) {
			System.err.println("Couldn't close file "+connFile);
			e.printStackTrace(); 
		}
		
	}

	public Statement createStatement() throws SQLException {
		return this.conn.createStatement();
	}
	
	public void executeSql(String query) throws SQLException{
		Statement stmt;		
		stmt = conn.createStatement();	
		stmt.execute(query);
		stmt.close();		    
	}
		
	/** 
	 * @param query
	 * @return the first column of the first row of the result of the given query as a string
	 * or null if no results were returned
	 */	
	public String getStringFromDb(String query) {
		Statement    stmt;
		ResultSet    rs;
		String       result = null;
		
		try { 
			
		    stmt = conn.createStatement();
		    rs = stmt.executeQuery(query);
		    if(rs.next()) {
		    	result = rs.getString(1);
		    }
		    rs.close();
		    stmt.close(); 
		    
		} // end try
		catch (SQLException e) {
		    System.err.println("SQLException: " + e.getMessage());
		    System.err.println("SQLState:     " + e.getSQLState());
		    System.err.println("VendorError:  " + e.getErrorCode());
		    e.printStackTrace();
		} // end catch				
		
		return result;
	}
	
	/** 
	 * @param query
	 * @return the first column of the first row of the result of the given query as an integer
	 * or -1 if no results were returned
	 */
	public int getIntFromDb(String query) {
		Statement    stmt;
		ResultSet    rs;
		int          result = -1;

		try { 
			
		    stmt = conn.createStatement();
		    rs = stmt.executeQuery(query);
		    if(rs.next()) {
		    	result = rs.getInt(1);
		    }
		    rs.close();
		    stmt.close(); 
		    
		} // end try
		catch (SQLException e) {
		    System.err.println("SQLException: " + e.getMessage());
		    System.err.println("SQLState:     " + e.getSQLState());
		    System.err.println("VendorError:  " + e.getErrorCode());
		    e.printStackTrace();
		} // end catch			
		
		return result;
	}
	
	public void close() {
		
		try {
			conn.close();
	    } catch (SQLException e) {
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
			System.err.println("VendorError:  " + e.getErrorCode());
			e.printStackTrace();
	    } // end try/catch connection 
	}
	
	/**
	 * Get the actual Connection object within this MySQLConnection. In case that an external library need a Connection object
	 * @return
	 */
	public Connection getConnectionObject() {
		return this.conn;
	}

    /**
     * To print the db size info for the given db of this MySQLConnection.
     * @param dbName
     */
	public void printDbSizeInfo (String dbName) {
		double data = 0, index = 0, table_data = 0, table_index = 0, GB = Math.pow(2, 30);
		String Query = null, table = null;
		Statement Stmt = null;
		ResultSet RS = null;		
		try {
		    Query = "SHOW TABLE STATUS FROM "+dbName;
		    Stmt = this.conn.createStatement();
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
		    System.out.println("Database "+dbName+" needs "+((data+index)/GB)+ " GB (data:"+(data/GB)+", index:"+(index/GB)+").");
		} 
		catch (SQLException e) {
		    System.err.println("SQLException: " + e.getMessage());
		    System.err.println("SQLState:     " + e.getSQLState());
		    System.err.println("VendorError:  " + e.getErrorCode());
		    e.printStackTrace();
		} // end try/catch connection
    }

    /**
     * To get range of a column from a table in this MySQLConnection. Note that if the connection was created without pointing 
     * to a particular database then the argument table must be specified as dbname.tablename.
     * @param table
     * @param column
     * @return
     */
    public double[] getRange(String table, String column) {
    	String query = "";
    	Statement S;
    	ResultSet R;
    	double[] range = new double[2];
    	try { 
    	    query = "SELECT MIN("+column+"), MAX("+column+") FROM "+table+";";
    	    S = this.conn.createStatement();
    	    R = S.executeQuery(query);    	
    	    if (R.next()) {
	    		range[0] = R.getDouble(1);
	    		range[1] = R.getDouble(2);
    	    } 
    	    R.close();
    	    S.close();
    	} // end try
    	catch (SQLException e) {
    	    System.err.println("SQLException: " + e.getMessage());
    	    System.err.println("SQLState:     " + e.getSQLState());
    	    System.err.println("VendorError:  " + e.getErrorCode());
    	    e.printStackTrace();
    	} // end catch    	
    	return range;    	
    }

    /**
     * To get range of a column from a table in this MySQLConnection using a given WHERE condition. Note that if 
     * the connection was created without pointing to a particular database then the argument table must 
     * be specified as dbname.tablename.
     * @param table
     * @param column
     * @param whereStr
     * @return
     */
    public double[] getRange(String table, String column, String whereStr) {
    	String query = "";
    	Statement S;
    	ResultSet R;
    	double[] range = new double[2];
    	try { 
    	    query = "SELECT MIN("+column+"), MAX("+column+") FROM "+table+" WHERE ("+whereStr+");";
    	    S = this.conn.createStatement();
    	    R = S.executeQuery(query);    	
    	    if (R.next()) {
	    		range[0] = R.getDouble(1);
	    		range[1] = R.getDouble(2);
    	    } 
    	    R.close();
    	    S.close();
    	} // end try
    	catch (SQLException e) {
    	    System.err.println("SQLException: " + e.getMessage());
    	    System.err.println("SQLState:     " + e.getSQLState());
    	    System.err.println("VendorError:  " + e.getErrorCode());
    	    e.printStackTrace();
    	} // end catch    	
    	return range;    	
    }
	
    /**
     * To get all indexes names for a certain table. Note the MySQLConnection object must be created with a non-blank database.
     * Using INFORMATION_SCHEMA db, only works from mysql 5.0
     * @param table
     * @return
     */
    public String[] getAllIndexes4Table(String table) {
    	ArrayList<String> indexesAL=new ArrayList<String>();
    	String query;
    	Statement S;
    	ResultSet R;
    	try { 
    	    query = "SELECT DISTINCT INDEX_NAME FROM INFORMATION_SCHEMA.STATISTICS WHERE TABLE_SCHEMA='"+dbname+"' AND TABLE_NAME='"+table+"';";
    	    S = this.conn.createStatement();
    	    R = S.executeQuery(query);    	
    	    while (R.next()) {
	    		indexesAL.add(R.getString(1));	    		
    	    } 
    	    R.close();
    	    S.close();
    	} // end try     		
    	catch (SQLException e) {
    	    System.err.println("SQLException: " + e.getMessage());
    	    System.err.println("SQLState:     " + e.getSQLState());
    	    System.err.println("VendorError:  " + e.getErrorCode());
    	    e.printStackTrace();
    	} // end catch
    	String[] indexes=new String[indexesAL.size()];
    	int i=0;
    	for (String index:indexesAL) {
    		indexes[i]=index;
    		i++;
    	}
    	return indexes;
    }
    
    /**
     * Gets an array of Strings with all queries necessary to create all the indexes for a certain table 
     * @param table
     * @return
     */
    public String[] getCreateIndex4Table(String table){
    	String[] indexes=this.getAllIndexes4Table(table);
    	String[] createIndexQueries=new String[indexes.length];     	
    	for (int i=0;i<indexes.length;i++){
    		String index=indexes[i];
        	try { 
        		Statement S;
        		ResultSet R;
        	    String query = "SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.STATISTICS " +
        	    				"WHERE TABLE_SCHEMA='"+dbname+"' AND TABLE_NAME='"+table+"' AND INDEX_NAME='"+index+"' " +
        	    						"ORDER BY SEQ_IN_INDEX;";
        	    S = this.conn.createStatement();
        	    R = S.executeQuery(query);
        	    String createIndexStr="CREATE INDEX "+index+" ON "+table+" (";
        	    while (R.next()) {
        	    	String colName = R.getString(1);
        	    	createIndexStr+=colName+",";
        	    }
        	    createIndexStr=createIndexStr.substring(0,createIndexStr.lastIndexOf(","));
        	    createIndexStr+=");";
        	    createIndexQueries[i]=createIndexStr;
        	    R.close();
        	    S.close();
        	} // end try     		
        	catch (SQLException e) {
        	    System.err.println("SQLException: " + e.getMessage());
        	    System.err.println("SQLState:     " + e.getSQLState());
        	    System.err.println("VendorError:  " + e.getErrorCode());
        	    e.printStackTrace();
        	} // end catch    		
    	}
    	return createIndexQueries;
    }
    
    /**
     * To get the column type for a certain column and table
     * @param table
     * @param column
     * @return
     */
    public String getColumnType(String table,String column){
    	String query = "SELECT COLUMN_TYPE FROM INFORMATION_SCHEMA.COLUMNS " +
    					"WHERE TABLE_SCHEMA='"+this.dbname+"' AND TABLE_NAME='"+table+"' AND COLUMN_NAME='"+column+"';";
    	String colType = this.getStringFromDb(query);
    	return colType;
    }
    
    /**
     * To findout whether a key (i.e. column) is numeric-based or text-based
     * @param table
     * @param key
     * @return true if is numeric-based, false if is text-based
     */
    public boolean isKeyNumeric(String table, String key){
    	boolean isNumeric = false;
    	String colType = getColumnType(table,key);
    	if (colType.contains("int") || colType.contains("INT")){
    		isNumeric = true;
    	}
    	else if (colType.contains("char") || colType.contains("CHAR")){
    		isNumeric = false;
    	}
    	else {
    		System.err.println("The key '"+key+"' from table '"+table+"' is neither numeric-based (int) nor text-based (char/varchar). Check what's wrong!");
    	}
    	return isNumeric;
    }
    
    /**
     * To get all tables for this MySQLConnection's database.
     * @return an array of String with all table names
     */
	public String[] getTables4Db(){
		String[] tables=null;
		ArrayList<String> tablesAL=new ArrayList<String>();
		String query="SELECT TABLE_NAME FROM INFORMATION_SCHEMA.TABLES WHERE TABLE_SCHEMA='"+dbname+"' ORDER BY TABLE_NAME DESC;";
		try {			
			Statement S = this.conn.createStatement();
			ResultSet R=S.executeQuery(query);
			while (R.next()){
				tablesAL.add(R.getString(1));
			}
			S.close();
			R.close();
			tables=new String[tablesAL.size()];
			for (int i=0;i<tablesAL.size();i++) {
				tables[i]=tablesAL.get(i);
			}
		}
		catch(SQLException e){			
			System.err.println("Couldn't get table names from "+host+" for db="+dbname);
    	    System.err.println("SQLException: " + e.getMessage());
    	    System.err.println("SQLState:     " + e.getSQLState());
    	    System.err.println("VendorError:  " + e.getErrorCode());
			e.printStackTrace();
		}
		return tables;
	}
	
	/** 
	 * To get all distinct ordered ids from a certain key and table from this MySQLConnection
	 * @param key the key name
	 * @param table the table name
	 * @return int array containing all ids 
	 */
	public Integer[] getAllNumIds4KeyAndTable(String key, String table){
		Integer[] allIds=null;		
		try {
			Statement S=conn.createStatement();
			String query="SELECT DISTINCT "+key+" FROM "+table+" ORDER BY "+key+";";
			ResultSet R=S.executeQuery(query);
			ArrayList<Integer> idsAL=new ArrayList<Integer>();
			while (R.next()){
				idsAL.add(R.getInt(1));
			}
			allIds=new Integer[idsAL.size()];
			for (int i=0;i<idsAL.size();i++) {
				allIds[i]=idsAL.get(i);
			}
			R.close();
			S.close();
		}
		catch (SQLException e){
			e.printStackTrace();
		}
		return allIds;
	}

	/** 
	 * To get all distinct ordered text (i.e. char/varchar column) ids from a certain key and table from this MySQLConnection
	 * @param key the key name
	 * @param table the table name
	 * @return int array containing all ids 
	 */
	public String[] getAllTxtIds4KeyAndTable(String key, String table){
		String[] allIds=null;		
		try {
			Statement S=conn.createStatement();
			String query="SELECT DISTINCT "+key+" FROM "+table+" ORDER BY "+key+";";
			ResultSet R=S.executeQuery(query);
			ArrayList<String> idsAL=new ArrayList<String>();
			while (R.next()){
				idsAL.add(R.getString(1));
			}
			allIds=new String[idsAL.size()];
			for (int i=0;i<idsAL.size();i++) {
				allIds[i]=idsAL.get(i);
			}
			R.close();
			S.close();
		}
		catch (SQLException e){
			e.printStackTrace();
		}
		return allIds;
	}

	/** 
	 * To get all distinct ordered ids from a certain key and table from this MySQLConnection
	 * @param key the key name
	 * @param table the table name
	 * @return int array containing all ids 
	 */
	public Object[] getAllIds4KeyAndTable(String key, String table){
		Object[] allIds=null;
		try {
			Statement S=conn.createStatement();
			String query="SELECT DISTINCT "+key+" FROM "+table+" ORDER BY "+key+";";
			ResultSet R=S.executeQuery(query);
			ArrayList<String> idsAL=new ArrayList<String>();
			while (R.next()){
				idsAL.add(R.getString(1));
			}
			if (isKeyNumeric(table,key)){
				allIds=new Integer[idsAL.size()];
				for (int i=0;i<idsAL.size();i++) {
					allIds[i]=Integer.parseInt(idsAL.get(i));
				}
			} else {
				allIds=new String[idsAL.size()];
				for (int i=0;i<idsAL.size();i++){
					allIds[i]=idsAL.get(i);
				}
			}
			R.close();
			S.close();
		}
		catch (SQLException e){
			e.printStackTrace();
		}
		return allIds;
	}

    /**
     * To set the sql_mode of this connection. 
     * @param sqlmode either NO_UNSIGNED_SUBTRACTION or blank
     */
    public void setSqlMode(String sqlmode) {
    	String query="SET SESSION sql_mode='"+sqlmode+"';";
    	try {
    		this.executeSql(query);
    	}
    	catch (SQLException e){
    		System.err.println("Couldn't change the sql mode to "+sqlmode);
    	    System.err.println("SQLException: " + e.getMessage());
    	    System.err.println("SQLState:     " + e.getSQLState());
    	    System.err.println("VendorError:  " + e.getErrorCode());
    	    e.printStackTrace();    		
    	}
    }
    
}

//package clusterconnection;
import java.sql.*;

/**
 * ClusterConnection class to wrap the master/node mysql servers so that is transparent to other programs
 * @author Jose Duarte
 */
//TODO Maybe we should always pass the key for CreateStatement and the constructor and drop the 2 methods without key.
//TODO I think that should make things easier. To have to remember to be switching keys is not ideal, I'd rather pass it always as a parameter in CreateStatement
//TODO The only problem with this approach is when constructing the object the key must be passed. That shouldn't be a major problem anyway. Revise all this

public class ClusterConnection {

	private final String url= "jdbc:mysql://";
	private final String masterDb="key_master"; 
	private final String masterHost="white";
	private Connection nCon;
	private Connection mCon;
	public String keyTable;
	public String key;
	public String idxColumn;
	public String host;
	public String db;
	private String user;
	private String password;
	
	//create a ClusterConnection passing a key. Then we need to call createNConStatement only passing the idx as argument.
	public ClusterConnection (String db,String key, String user,String password) {
		loadMySQLDriver();
		setDb(db);
		setKey(key);
		setIdxColumn();
		setMasterTable(db);
		setUser(user);
		setPassword(password);
		try {
			// For nCon we create a connection to the master too. 
			// This is just a place holder because the actual node connection is not created until we create the statement
			// If we don't do this then when closing the two connections an exception might occurr because we try to close a non existing object
			this.nCon = DriverManager.getConnection(url+masterHost+"/"+masterDb,user,password); 
			this.mCon = DriverManager.getConnection(url+masterHost+"/"+masterDb,user,password);
		}
		catch(SQLException e){
    	    System.out.println("SQLException: " + e.getMessage());
    	    System.out.println("SQLState:     " + e.getSQLState());
    	    System.out.println("VendorError:  " + e.getErrorCode());
			System.out.println("Couldn't get connection to master host "+masterHost+", db="+masterDb+", exiting.");
			System.exit(2);			
		}
	}

	// create a cluster connection with an empty key field. In this case we need to call createNConStatement passing the key as argument 
	public ClusterConnection (String db, String user,String password) { 
		this.loadMySQLDriver();
		setDb(db);
		setUser(user);
		setPassword(password);
		try {
			// For nCon we create a connection to the master too. 
			// This is just a place holder because the actual node connection is not created until we create the statement
			// If we don't do this then when closing the two connections an exception might occurr because we try to close a non existing object
			this.nCon = DriverManager.getConnection(url+masterHost+"/"+masterDb,user,password);
			this.mCon = DriverManager.getConnection(url+masterHost+"/"+masterDb,user,password);
		}
		catch(SQLException e){
    	    System.out.println("SQLException: " + e.getMessage());
    	    System.out.println("SQLState:     " + e.getSQLState());
    	    System.out.println("VendorError:  " + e.getErrorCode());
			System.out.println("Couldn't get connection to master host "+masterHost+", db="+masterDb+", exiting.");
			System.exit(2);			
		}
	}

	
	public void loadMySQLDriver() {
		try {
			Class.forName("com.mysql.jdbc.Driver");
		}
		catch(Exception e) {
			System.out.println(e.getMessage());
			System.out.println("An exception occurred while loading the mysql jdbc driver, exiting.");
			System.exit(1);
		}
	}
	public void close() {
		try {
			this.nCon.close();
			this.mCon.close();
		}
		catch(SQLException e) {
			System.out.println("Couldn't close database connections for master: "+masterHost+" and node: "+this.host+", exiting.");
    	    System.out.println("SQLException: " + e.getMessage());
    	    System.out.println("SQLState:     " + e.getSQLState());
			System.exit(3);						
		}
	}
	
	public String getHost4Idx (int idx) {
		String host="";
		Statement S;
		ResultSet R;
		try {
			S=mCon.createStatement();
			String query="SELECT client_name FROM "+keyTable+" AS m INNER JOIN clients_names AS c "+
						"ON (m.client_id=c.client_id) WHERE "+idxColumn+"="+idx+";";
			R=S.executeQuery(query);
			if (R.next()){
				host=R.getString(1);
			}
			S.close();
			R.close();
		}
		catch(SQLException e) {
			System.out.println("Couldn't get the host name for idx "+idx+", exiting");
    	    System.out.println("SQLException: " + e.getMessage());
    	    System.out.println("SQLState:     " + e.getSQLState());
			System.exit(3);									
		}
		return host;
	}
	
	public void setHostFromIdx(int idx){
		setHost(getHost4Idx(idx));
	}
	
	public void setHost(String host) {
		this.host=host;
		try {
			this.nCon=DriverManager.getConnection(url+host+"/"+db,user,password);
		}
		catch (SQLException e){
    	    System.out.println("SQLException: " + e.getMessage());
    	    System.out.println("SQLState:     " + e.getSQLState());
			System.out.println("Couldn't get connection to node "+host+", exiting.");
			System.exit(2);						
		}
	}
	
	public Statement createStatement(int idx) { // to use when the field "key" is already set
		setMasterTable();
		setIdxColumn();
		Statement S=null;
		this.setHostFromIdx(idx);
		try {
			S=this.nCon.createStatement();
		}
		catch (SQLException e){
			System.out.println("SQLException: " + e.getMessage());
			System.out.println("SQLState:     " + e.getSQLState());
			System.out.println("Couldn't create statement for the node connection, idx= "+idx+", exiting.");
			System.exit(2);									
		}
		return S;
	}
	
	public Statement createStatement(String key,int idx) {
		setKey(key);
		setMasterTable();
		setIdxColumn();
		return createStatement(idx);
	}
	
	public void setMasterTable(String db) { // to set masterTable field in constructor (i.e. first time)
		this.keyTable=db+"_"+this.key+"_list_master";
	}
	
	public void setMasterTable() { // to set masterTable field when db is already set 
		this.keyTable=this.db+"_"+this.key+"_list_master";
	}
	
	public String getMasterTable() {
		return this.keyTable;
	}
	
	public void setIdxColumn() {
		this.idxColumn=this.key+"_id";
	}

	public void setUser(String user) {
		this.user=user;
	}

	public void setPassword(String password) {
		this.password=password;
	}
	
	public void setDb(String db){
		this.db=db;
	}
	
	public void setKey(String key){
		this.key=key;
	}

	public void switchKey(String key){
		this.key=key;
		setMasterTable();
		setIdxColumn();
	}

	public Statement createMasterStatement() { 
		Statement S=null;
		try {
			S=this.mCon.createStatement();
		}
		catch (SQLException e){
			System.out.println("SQLException: " + e.getMessage());
			System.out.println("SQLState:     " + e.getSQLState());
			System.out.println("Couldn't create statement for the master connection, exiting.");
			System.exit(2);									
		}
		return S;
	}

	public ResultSet getAllIdxFromMaster(String key) { 
		this.setKey(key);
		this.setIdxColumn();
		this.setMasterTable();
		String query;
		Statement S;
		ResultSet R=null;
		try {
			query="SELECT "+idxColumn+" FROM "+keyTable+";";
			S=this.mCon.createStatement();
			R=S.executeQuery(query);
			//S.close(); // apparently it doesn't work if we close the Statement!! Don't know why!
		}
		catch (SQLException e){
			System.out.println("SQLException: " + e.getMessage());
			System.out.println("SQLState:     " + e.getSQLState());
			System.out.println("Couldn't get all indices from columnn "+idxColumn+" in table "+keyTable+" from "+masterDb+" database in "+masterHost+", exiting.");
			System.exit(2);									
		}
		return R;
	}

	//	 to get client_id for a certain idx and key
	//TODO will need to change the query. In general this method would return more than 1 client_id if the idx is not unique
	public int getHostId4Idx (String key,int idx) { 
		int hostId=0;
		this.setKey(key);
		this.setIdxColumn();
		this.setMasterTable();
		Statement S;
		ResultSet R;
		String query;
		int countCids=0;
		try {
			S=mCon.createStatement();
			query="SELECT count(client_id) FROM "+keyTable+" WHERE "+idxColumn+"="+idx+";";
			R=S.executeQuery(query);
			if (R.next()){
				countCids=R.getInt(1);
			}
			if (countCids!=1){
				System.out.println("Error! the number of client_id for idx "+idxColumn+"= "+idx+" is bigger than 1. Check what's wrong! Exiting now.");
				System.exit(2);
			}
			else {
				query="SELECT client_id FROM "+keyTable+" WHERE "+idxColumn+"="+idx+";";
				R=S.executeQuery(query);
				if (R.next()){
					hostId=R.getInt(1);
				}
			}
			S.close();
			R.close();
		}
		catch(SQLException e) {
			System.out.println("Couldn't get the host id for idx "+idxColumn+"="+idx+", exiting");
    	    System.out.println("SQLException: " + e.getMessage());
    	    System.out.println("SQLState:     " + e.getSQLState());
			System.exit(3);									
		}
		return hostId;
	}
	
	public void insertIdxInMaster(String key, int clientId) {
		Statement S;
		String query;
		this.setKey(key);
		this.setMasterTable();
		this.setIdxColumn();
		try {
			S=this.mCon.createStatement();
			query="INSERT INTO "+this.keyTable+" (client_id) VALUES ("+clientId+");";
			S.executeUpdate(query);
			S.close();
		} 
		catch (SQLException E) {
			System.out.println("SQLException: " + E.getMessage());
			System.out.println("SQLState:     " + E.getSQLState());
			System.out.println("Couldn't insert new "+this.idxColumn+" in master table "+this.getMasterTable()+". The client_id for it was "+clientId+". Exiting.");
			System.exit(2);
		} 
	}

	public void insertIdxInMaster(String keySrc,String keyDest,int idxSrc) {
		this.setKey(keySrc);
		this.setMasterTable();
		this.setIdxColumn();
		int clientId=0;
		clientId=this.getHostId4Idx(keySrc,idxSrc);
		insertIdxInMaster(keyDest,clientId);
	}

	public int getLastInsertId(String key) {
		int lastIdx=0;
		this.setKey(key);
		this.setMasterTable();
		this.setIdxColumn();
		Statement S;
		ResultSet R;
		String query = "";
		try {
			S = this.mCon.createStatement();
			query = "SELECT LAST_INSERT_ID() FROM "+this.keyTable+" LIMIT 1;";
			R = S.executeQuery(query);
			if (R.next()) {
				lastIdx=R.getInt(1);
			}
			R.close();
			S.close();
		} 
		catch (SQLException E) {
			System.out.println("Couldn't get the last insert id for key type "+this.idxColumn+" from table "+this.keyTable+". Exiting");
			System.out.println("SQLException: " + E.getMessage());
			System.out.println("SQLState:     " + E.getSQLState());
			System.exit(3);
		} // end try/catch connection 
		return lastIdx;
	} // end getGraphId
}

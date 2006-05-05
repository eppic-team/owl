package tools;
import java.sql.*;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * ClusterConnection class to wrap the master/node mysql servers so that is transparent to other programs
 * @author Jose Duarte
 */

public class ClusterConnection {

	private String MASTERDB=DataDistribution.KEYMASTERDB; 
	private final String MASTERHOST=DataDistribution.MASTER;
	private MySQLConnection nCon;
	private MySQLConnection mCon;
	public String keyTable;
	public String key;
	public String host;
	public String db;
	private String user;
	private String password;

	/**
	 * Create a ClusterConnection passing a key. 
	 * @param db the database name
	 * @param key the key name: e.g. asu_id
	 * @param user the user name for connection to both master and nodes
	 * @param password the password for connection to both master and nodes
	 */
	public ClusterConnection (String db,String key, String user,String password) {
		setDb(db);
		setUser(user);
		setPassword(password);
		// For nCon we create a connection to the master too. 
		// This is just a place holder because the actual node connection is not created until we create the statement
		// If we don't do this then when closing the two connections an exception might occurr because we try to close a non existing object
		this.nCon = new MySQLConnection(MASTERHOST,user,password,MASTERDB); 
		this.mCon = new MySQLConnection(MASTERHOST,user,password,MASTERDB);
		this.key=key; //can't use the setKey method here before we've got the db field initialized
		setKeyTable(db);
	}

	/**
	 * Create a ClusterConnection without passing a key. The key will be set later when we call createStatement(key,idx) 
	 * @param db the database name
	 * @param user the user name for connection to both master and nodes
	 * @param password the password for connection to both master and nodes
	 */ 
	public ClusterConnection (String db, String user,String password) { 
		setDb(db);
		setUser(user);
		setPassword(password);		
		// For nCon we create a connection to the master too. 
		// This is just a place holder because the actual node connection is not created until we create the statement
		// If we don't do this then when closing the two connections an exception might occurr because we try to close a non existing object
		this.nCon = new MySQLConnection(MASTERHOST,user,password,MASTERDB); 
		this.mCon = new MySQLConnection(MASTERHOST,user,password,MASTERDB);
	}
	
	public void close() {
		this.nCon.close();
		this.mCon.close();
	}
	
	public String getHost4Idx (int idx) {
		String host="";
		String query = "SELECT client_name FROM "+keyTable+" AS m INNER JOIN clients_names AS c "+
					   "ON (m.client_id=c.client_id) WHERE "+key+"="+idx+";";
		host=this.mCon.getStringFromDb(query);
		return host;
	}
	
	public void setHostFromIdx(int idx){
		setHost(getHost4Idx(idx));
	}
	
	public void setHost(String host) {
		this.host=host;
		//Closing previous connection is essential
		//If we don't close it a lot of connections stay open after using a ClusterConnection object for a while
		this.nCon.close(); 
		this.nCon = new MySQLConnection(host,user,password,db);
	}
	
	/**
	 * This method is strictly private. We shouldn't call this from another class as a key might not be set when we call it
	 * and thus we can't get the client_id from the master key table. Only to be called from createStatement(key,idx)
	 * @param idx the value of the id for a certain key already set
	 * @return a Stament object with a connection to the node that contains idx for key
	 */
	private Statement createStatement(int idx) { // to use when the field "key" is already set
		setKeyTable();
		Statement S=null;
		this.setHostFromIdx(idx);
		try {
			S=this.nCon.createStatement();
		}
		catch (SQLException e){
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
			System.err.println("Couldn't create statement for the node connection, idx= "+idx+", exiting.");
			System.exit(2);									
		}
		return S;
	}
	
	/**
	 * This method is used to create a statement passing the key and idx. It will create a connection the the right node 
	 * and return a Statement for that connection
	 * @param idx the key name
	 * @param idx the id value for that key
	 * @return a Statement object with a connection to the node that contains idx for key
	 */
	public Statement createStatement(String key,int idx) {
		setKey(key);
		return createStatement(idx);
	}
	
	/**
	 * To execute a sql update/insert query in the right node given a query, key and idx. Just a shortcut not to have to do the create statement and execute
	 * @param query the SQL query
	 * @param key the name of the key
	 * @param idx the id value for that key
	 */
	public void executeSql(String query,String key, int idx) throws SQLException {
		Statement stmt;
		stmt = this.createStatement(key,idx);	
		stmt.execute(query);
		stmt.close();		    
	}
	
	/**
	 * To change the MASTERDB String, i.e. the name of the key master database. To be used in testing.
	 * @param db the name of the key master db we want to use instead of the default defined in the MASTERDB field
	 */
	public void setKeyDb(String db) { 
		this.MASTERDB=db;
		//Closing previous connection is essential
		this.mCon.close(); 
		this.mCon = new MySQLConnection(MASTERHOST,user,password,MASTERDB);
	}
	
	/**
	 * To set keyTable field in constructor (i.e. first time). Only to be used in constructor.
	 * @param db the database name
	 */
	public void setKeyTable(String db) { 
		String query="SELECT key_master_table FROM dbs_keys WHERE db=\'"+db+"\' AND key_name=\'"+this.key+"\';";
		this.keyTable=this.mCon.getStringFromDb(query);
	}
	
	/**
	 * To set the keyTable field when db is already set
	 * The value of keyTable is taken from the dbs_keys table in the database given the db and key.
	 */
	public void setKeyTable() {  
		String query="SELECT key_master_table FROM dbs_keys WHERE db=\'"+this.db+"\' AND key_name=\'"+this.key+"\';";
		this.keyTable=this.mCon.getStringFromDb(query);
	}
	
	public String getKeyTable() {
		return this.keyTable;
	}
	
	/**
	 * To get the name of the target table where splitted data is stored in nodes, e.g. for keyTable pdbgraph__asu_list, we get asu_list
	 * @return
	 */
	public String getTableOnNode(){
		String table="";
		if (this.keyTable.contains("__")) {
			String[] tokens=this.keyTable.split("__");
			table=tokens[1];
		}
		else {
			System.err.println("Error! The keyTable field is not set in this ClusterConnection object.");
		}
		return table;
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
		setKeyTable();
	}

	public Statement createMasterStatement() { 
		Statement S=null;
		try {
			S=this.mCon.createStatement();
		}
		catch (SQLException e){
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
			System.err.println("Couldn't create statement for the master connection, exiting.");
			System.exit(2);									
		}
		return S;
	}
	
	public int[] getAllIdxFromMaster(String key) { 
		this.setKey(key);
		int[] ids=null;
		ArrayList<Integer> idsAL=new ArrayList<Integer>();
		try {			
			String query="SELECT "+key+" FROM "+keyTable+";";
			Statement S=this.mCon.createStatement();
			ResultSet R=S.executeQuery(query);
			while (R.next()){
				idsAL.add(R.getInt(1));
			}
			R.close();
			S.close();
		}
		catch (SQLException e){
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
			System.err.println("Couldn't get all indices from columnn "+key+" in table "+keyTable+" from "+MASTERDB+" database in "+MASTERHOST+", exiting.");
			System.exit(2);									
		}
		ids=new int[idsAL.size()];
		for (int i=0;i<idsAL.size();i++){
			ids[i]=idsAL.get(i);
		}
		return ids;
	}
	/**
	 * To get all ids and clients_names pairs for a certain key. Useful when need to submit to all hosts using qsub -q
	 * @param key the name of the key
	 * @return HashMap with keys = indices, and values = node names where the corresponding index is stored
	 */
	public HashMap<Integer,String> getAllIdxAndClients (String key){
		this.setKey(key);
		HashMap<Integer,String> idsAndClients=new HashMap<Integer,String>();
		try {
			String query="SELECT a."+key+",c.client_name FROM "+keyTable+" AS a INNER JOIN clients_names AS c ON (a.client_id=c.client_id);";
			Statement S=this.mCon.createStatement();
			ResultSet R=S.executeQuery(query);
			while (R.next()){
				idsAndClients.put(R.getInt(1),R.getString(2));
			}
			R.close();
			S.close();
		}
		catch (SQLException e){
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
			System.err.println("Couldn't get all indices/client_names pairs from table "+keyTable+" from "+MASTERDB+" database in "+MASTERHOST+", exiting.");
			System.exit(2);									
		}
		return idsAndClients;
	}
	
	/**
	 * To get client_id for a certain idx and key
	 * @param key the key name
	 * @param idx the id value for that key
	 * @return the client_id for node that has the data for the idx for that key
	 */
	//TODO will need to change the query. In general this method would return more than 1 client_id if the idx is not unique
	public int getHostId4Idx (String key,int idx) { 
		int hostId=0;
		this.setKey(key);
		String query;
		int countCids=0;
		query="SELECT count(client_id) FROM "+keyTable+" WHERE "+key+"="+idx+";";
		countCids=this.mCon.getIntFromDb(query);
		if (countCids!=1){
			System.err.println("the query was: "+query);
			System.err.println("Error! the count of client_ids for idx "+key+"= "+idx+" is " +countCids+
					". It must be 1! The values were taken from host: "+MASTERHOST+", database: "+MASTERDB+", table: "+keyTable+". Check what's wrong! Exiting now.");
			System.exit(2);
		}
		else {
			query="SELECT client_id FROM "+keyTable+" WHERE "+key+"="+idx+";";
			hostId=this.mCon.getIntFromDb(query);
		}
		return hostId;
	}
	
	public void insertIdxInMaster(String key, int clientId) {
		String query;
		this.setKey(key);
		try {
			query="INSERT INTO "+this.keyTable+" (client_id) VALUES ("+clientId+");";
			this.mCon.executeSql(query);
		} 
		catch (SQLException E) {
			System.err.println("SQLException: " + E.getMessage());
			System.err.println("SQLState:     " + E.getSQLState());
			System.err.println("Couldn't insert new "+this.key+" in master table "+this.getKeyTable()+". The client_id for it was "+clientId+". Exiting.");
			System.exit(2);
		} 
	}

	public void insertIdxInMaster(String keySrc,String keyDest,int idxSrc) {
		this.setKey(keySrc);
		int clientId=0;
		clientId=this.getHostId4Idx(keySrc,idxSrc);
		insertIdxInMaster(keyDest,clientId);
	}

	public int getLastInsertId(String key) {
		int lastIdx=0;
		this.setKey(key);
		String query = "";
		query = "SELECT LAST_INSERT_ID() FROM "+this.keyTable+" LIMIT 1;";
		lastIdx=this.mCon.getIntFromDb(query);
		return lastIdx;
	}
	
	public int[][] getIdxSet(String key) { 
		int[][] indMatrix=null;
		this.setKey(key);
		String query;
		Statement S;
		ResultSet R;
		try {
			// STEP 1 -- getting set of all client_ids			
			query="SELECT count(distinct client_id) FROM "+keyTable+";";
			int count=0;
			count=this.mCon.getIntFromDb(query);
			S=this.mCon.createStatement();
			query="SELECT DISTINCT client_id FROM "+keyTable+" ORDER BY client_id;";
			R=S.executeQuery(query);
			
			// STEP 2 -- putting sets of indices counts into temp tables c_<client_id> with a serial auto_increment field
			int[] clids=new int[count]; //array to store all client_ids. To be used in loops later
			int i=0;
			while (R.next()){
				Statement Sloop=this.mCon.createStatement();
				int clid=R.getInt(1);
				query="CREATE TEMPORARY TABLE c_"+clid+" (serial int(11) NOT NULL AUTO_INCREMENT,"+key+" int(11),client_id int(11), PRIMARY KEY(serial));";				
				Sloop.executeUpdate(query);
				query="INSERT INTO c_"+clid+" ("+key+",client_id) SELECT "+key+",client_id FROM "+keyTable+" WHERE client_id="+clid+";";
				Sloop.executeUpdate(query);
				clids[i]=clid;
				i++;
				Sloop.close();
			}
			
			// STEP3 -- merging all c_<client_id> tables into a temp table tmp_allcs and selecting the client_id with the maximum count
			//query="SELECT client_id,count(*) as c FROM c_34 GROUP BY client_id UNION SELECT client_id,count(*) as c FROM c_32 GROUP BY client_id;";
			query="DROP TABLE IF EXISTS tmp_allcs;"; 
			S.executeUpdate(query);
			//this table must be permanent! otherwise cannot do the select max(c) later
			query="CREATE TABLE IF NOT EXISTS tmp_allcs (client_id int(11), c int(11)) ENGINE=MEMORY;";
			S.executeUpdate(query);
			String unionStr="SELECT client_id,count(*) AS c FROM c_"+clids[0]+" GROUP BY client_id";
			for (i=1;i<clids.length;i++) {
				unionStr+=" UNION SELECT client_id,count(*) AS c FROM c_"+clids[i]+" GROUP BY client_id";
			}
			query="INSERT INTO tmp_allcs "+unionStr+";";
			S.executeUpdate(query);
			query="SELECT client_id,c FROM tmp_allcs WHERE c=(SELECT max(c) FROM tmp_allcs);";
			R=S.executeQuery(query);
			int clidMaxIdxCount=0;
			int maxIdxCount=0;
			if (R.next()) {
				clidMaxIdxCount=R.getInt(1);
				maxIdxCount=R.getInt(2);
			}
			query="DROP TABLE tmp_allcs;";
			S.executeUpdate(query);
			
			// STEP 4 -- join all c_<client_id> tables into a table with a serial column, and c_<client_id> columns each of them with the indices for each client_id 
			//query="SELECT c_34.serial,c_34.asu_id AS c_34,c_32.asu_id AS c_32,c_36.asu_id AS c_36 FROM c_34 LEFT JOIN c_32 ON (c_34.serial=c_32.serial) LEFT JOIN c_36 ON (c_34.serial=c_36.serial);";
			String selectStr="c_"+clidMaxIdxCount+".serial, c_"+clidMaxIdxCount+"."+key+" AS c_"+clidMaxIdxCount;
			String fromStr="c_"+clidMaxIdxCount;
			for (i=0;i<clids.length;i++) {
				if (clids[i]!=clidMaxIdxCount){
					selectStr+=", c_"+clids[i]+"."+key+" AS c_"+clids[i];
					fromStr+=" LEFT JOIN c_"+clids[i]+" ON (c_"+clidMaxIdxCount+".serial=c_"+clids[i]+".serial)";
				}
			} 
			query="CREATE TEMPORARY TABLE indices_matrix "+"SELECT "+selectStr+" FROM "+fromStr+";";
			S.executeUpdate(query);
			
			// STEP 5 -- put the table into a 2-dimensional array and return it
			indMatrix = new int[maxIdxCount][clids.length];
			query="SELECT * FROM indices_matrix";
			R=S.executeQuery(query);
			i=0;
			while (R.next()) {
				for (int j=0;j<clids.length;j++){
					indMatrix[i][j]=R.getInt(j+2);
				}
				i++;
			}
			R.close();
			S.close();
		}
		catch (SQLException e){
			System.err.println("SQLException: " + e.getMessage());
			System.err.println("SQLState:     " + e.getSQLState());
			System.err.println("Couldn't get the indices set from columnn "+key+" in table "+keyTable+" from "+MASTERDB+" database in "+MASTERHOST+", exiting.");
			System.exit(2);									
		}
		return indMatrix;
	}

	
}


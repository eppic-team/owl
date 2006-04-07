package tools;

import java.io.*;
import java.sql.*;
import java.util.ArrayList;
import java.util.HashMap;

public class DataDistribution {

	final static String GLOBALDIR="/project/snow/global/tmp";
	final static String ADMINDIR="/project/StruPPi/Cluster/admin";
	public final static String MASTER="white";
	final static String HFILE=ADMINDIR+"/hosts_ss.txt";
	final static String KEYMASTERDB="key_master";

	boolean debug=false; // if set to true only mysql commands written, no actual dump or load, also dump directory not removed. Use setDebug method to change it
	String dumpdir;
	private String db;
	private String user;
	private String pwd;
	

	public DataDistribution(String db,String user,String pwd) {
		this.dumpdir=GLOBALDIR+"/dumps_tmp_"+System.currentTimeMillis();
		this.db=db;
		this.user=user;
		this.pwd=pwd;
	}
	
	private MySQLConnection getConnectionToMaster() {
		MySQLConnection conn=this.getConnectionToNode(MASTER);
		return conn;
	}

	private MySQLConnection getConnectionToMasterKeyDb() {
		MySQLConnection conn=new MySQLConnection(MASTER,user,pwd,KEYMASTERDB);
		return conn;
	}
	
	private MySQLConnection getConnectionToNode(String node) {
		MySQLConnection conn=new MySQLConnection(node,user,pwd,db);
		return conn;
	}
	
	private MySQLConnection getConnectionToNode(String node, String dbName){
		MySQLConnection conn=new MySQLConnection(node,user,pwd,dbName);
		return conn;
	}
	
	public void setDebug(boolean debug) {
		this.debug=debug;
	}
	
	public static String[] getNodes() {
		String[] nodes=null;	    
	    try {
			File inputFile = new File(HFILE);
	    	BufferedReader hostsFile = new BufferedReader(new FileReader(inputFile)); // open BufferedReader to the file
	    	String nodesstr=hostsFile.readLine();
	    	nodes=nodesstr.split(" ");
	    	hostsFile.close();
	    }
	    catch (IOException e){
	    	e.printStackTrace();
	    	System.err.println("Couldn't read from file "+HFILE);
	    }
	    return nodes;
	}

	public String[] getTables4Db(){
		String[] tables=null;
		ArrayList<String> tablesAL=new ArrayList<String>();
		String query="SELECT table_name FROM information_schema.tables WHERE table_schema='"+db+"' ORDER BY table_name DESC;";
		try {
			MySQLConnection conn = this.getConnectionToMaster();
			Statement S=conn.createStatement();
			ResultSet R=S.executeQuery(query);
			while (R.next()){
				tablesAL.add(R.getString(1));
			}
			S.close();
			R.close();
			conn.close();
			tables=new String[tablesAL.size()];
			for (int i=0;i<tablesAL.size();i++) {
				tables[i]=tablesAL.get(i);
			}
		}
		catch(SQLException e){
			e.printStackTrace();
			System.err.println("Couldn't get table names from "+MASTER+" for db="+db);
			System.exit(1);
		}
		return tables;
	}
	
	public void initializeDirs() {	
	    if (!((new File(dumpdir)).mkdir())) {
	    	System.err.println("Couldn't create directory "+dumpdir);
	    	System.exit(1);
	    }
		SystemCmd.exec("chmod go+rw "+dumpdir);
	}
	
	public void finalizeDirs() {
	    if (debug) {
	    	System.out.println("Temporary directory "+dumpdir+" was not removed. You must remove it manually");
	    } else {
	    	System.out.println("Removing temporary directory "+dumpdir);
	    	//TODO must capture exit state and print to error if problems deleting dir
	    	SystemCmd.exec("rm -rf "+dumpdir);
	    }
	}
	
	public void dumpData(String[] srchosts, String[] tables) {
		// initialising temporary dump directory
		initializeDirs();
		for (String node: srchosts) {
			dumpData(node,tables);
		}
		System.out.println ("Dump finished.");
	}
	
	private void dumpData(String srchost, String[] tables) {
		//String quickpar="--quick"; // not to buffer to memory, needed for big tables
		if (!((new File(dumpdir+"/"+srchost+"/"+db)).mkdirs())) {
			System.err.println("Couldn't create directory "+dumpdir+"/"+srchost+"/"+db);
			System.exit(1);
		}
		SystemCmd.exec("chmod -R go+rw "+dumpdir+"/"+srchost);
		for ( String tbl: tables) {
			String outfile=dumpdir+"/"+srchost+"/"+db+"/"+tbl+".txt";
			String wherestr="";
			String sqldumpstr="SELECT * FROM `"+tbl+"` "+wherestr+" INTO OUTFILE '"+outfile+"';";
			//String dumpstr="$MYSQLDIR/bin/mysql $srcconnpar $quickpar -e \"$sqldumpstr\" ";
			if (debug) {System.out.println ("HOST="+srchost+", sqldumpstr="+sqldumpstr);} 
			else {
				try {
					MySQLConnection conn = this.getConnectionToNode(srchost);
					Statement S=conn.createStatement();
					S.executeQuery(sqldumpstr);
					System.out.println ("Dumped from host="+srchost+", database="+db+", table="+tbl+" to outfile="+outfile);
					S.close();
					conn.close();
				}
				catch(SQLException e){
					e.printStackTrace();
					System.err.println("Couldn't dump from host="+srchost+", database="+db+", table="+tbl+" to outfile="+outfile);
					System.exit(1);
				}		
			} //end else if debug
		} // end foreach tbl
	}
	
	public void loadData(String[] srchosts, String[] desthosts,String[] tables, String destDb) {
	    for (String desthost: desthosts) {
	    	for (String srchost: srchosts) {
	    		loadData(srchost,desthost,tables,destDb);
	    	} 
	    } 
	    System.out.println ("Load finished.");
	    finalizeDirs();
	}

	private void loadData(String srchost, String desthost,String[] tables, String destDb) {
		for (String tbl:tables) {
			String dumpfile=dumpdir+"/"+srchost+"/"+db+"/"+tbl+".txt";
			String sqlloadstr="LOAD DATA INFILE '"+dumpfile+"' INTO TABLE `"+tbl+"`;";
			if (debug) {System.out.println ("SRCHOST="+srchost+", DESTHOST="+desthost+", sqlloadstr="+sqlloadstr);} 
			else {
				try {
					MySQLConnection conn = this.getConnectionToNode(desthost,destDb);
					conn.executeSql(sqlloadstr);
					System.out.println ("HOST: "+desthost+". Loaded from file="+dumpfile+", into database="+destDb+", table="+tbl);
					conn.close();
				}
				catch(SQLException e){
					e.printStackTrace();
					System.err.println("Errors occurred while loading data from file="+dumpfile+" into host="+desthost+", database="+destDb+", table="+tbl);
					System.exit(1);
				}
			}	    		
		} // end foreach tbl
	}
	
	private void loadSplitData(String srchost, String[] desthosts, String tableName, String destDb) {
		for (String desthost:desthosts) {
			String dumpfile=dumpdir+"/"+srchost+"/"+db+"/"+tableName+"_split_"+desthost+".txt";
			String sqlloadstr="LOAD DATA INFILE '"+dumpfile+"' INTO TABLE `"+tableName+"`;";
			if (debug) {System.out.println ("SRCHOST="+srchost+", DESTHOST="+desthost+", sqlloadstr="+sqlloadstr);} 
			else {
				try {
					MySQLConnection conn = this.getConnectionToNode(desthost,destDb);
					conn.executeSql(sqlloadstr);
					System.out.println ("HOST: "+desthost+". Loaded from file="+dumpfile+", database="+db+", table="+tableName);
					conn.close();
				}
				catch(SQLException e){
					e.printStackTrace();
					System.err.println("Errors occurred while loading data from file="+dumpfile+" into host="+desthost+", database="+destDb+", table="+tableName);
					System.exit(1);
				}
			}	    		
		} // end foreach desthosts
	    System.out.println ("Load finished.");
	    finalizeDirs();
	}
	
	public boolean checkCountsAllTables (){
		boolean checkResult=true;
		String[] nodes = getNodes();
		String[] tables = getTables4Db();
		// getting hashmap of all counts from all tables from nodes
		HashMap<String,HashMap<String,Integer>> countsNodes=new HashMap<String,HashMap<String,Integer>>();
		for (String node:nodes){
			HashMap<String,Integer> tableCounts = new HashMap<String,Integer>();		
			for (String tbl:tables){
				String query="SELECT count(*) FROM "+tbl+";";				
				try {
					MySQLConnection conn = this.getConnectionToNode(node);
					Statement S=conn.createStatement();
					ResultSet R=S.executeQuery(query);
					if (R.next()){
						tableCounts.put(tbl,R.getInt(1));
					}
					S.close();
					R.close();
					conn.close();
				}
				catch(SQLException e){
					e.printStackTrace();
					System.err.println("Couldn't execute query in host="+node+", database="+db);
					System.exit(1);
				}
			}
			countsNodes.put(node,tableCounts);
		}
		// getting hashmap of all counts of all tables from master
		HashMap<String,Integer> countsMaster= new HashMap<String,Integer>();
		for (String tbl:tables){
			String query="SELECT count(*) FROM "+tbl+";";			
			try {
				MySQLConnection conn = this.getConnectionToMaster();
				Statement S=conn.createStatement();
				ResultSet R=S.executeQuery(query);
				if (R.next()){
					countsMaster.put(tbl,R.getInt(1));
				}
				S.close();
				R.close();
				conn.close();
			}
			catch(SQLException e){
				e.printStackTrace();
				System.err.println("Couldn't execute query in host="+MASTER+", database="+db);
				System.exit(1);
			}
		}
		// comparing the nodes counts with the master counts
		for (String tbl:countsMaster.keySet()){
			int masterCount=countsMaster.get(tbl);
			for (String node:countsNodes.keySet()){
				int thisNodeCount=countsNodes.get(node).get(tbl);
				if (masterCount!=thisNodeCount) {					
					System.out.println("Count difers for table "+tbl+" in database "+db+" of node "+node+". MASTER COUNT="+masterCount+", NODE COUNT="+thisNodeCount);
					checkResult=false; // if one count difers then the check fails, we return false for checkResult. If no difers at all we return true
				}
				else {
					System.out.println("Count check passed for node "+node+", table "+tbl);
				}
				
			}
		}		
		return checkResult;
	} // end checkCounts

	/**
	 * To check the key counts in master and nodes for a certain key.  
	 * @param key the name of the key
	 * @return boolean: true if check passed, false if not passed
	 */
	public boolean checkKeyCounts(String key) {
		boolean checkResult=true;
		String keyTable=key+"_list";
		String masterKeyTable=db+"_"+key+"_list_master";
		String keyColumn =key+"_id";
		// getting hashmap of counts of keys from nodes
		String[] nodes = getNodes();
		HashMap<String,int[]> countsNodes=new HashMap<String,int[]>();
		String query="SELECT count("+keyColumn+"),count(DISTINCT "+keyColumn+") FROM "+keyTable+";";
		for (String node:nodes){
			try {				
				MySQLConnection conn = this.getConnectionToNode(node);
				Statement S=conn.createStatement();
				ResultSet R=S.executeQuery(query);
				int[] thisNodeKeyCount=new int[2];
				if (R.next()){
					thisNodeKeyCount[0]=R.getInt(1);
					thisNodeKeyCount[1]=R.getInt(2);
				}
				countsNodes.put(node,thisNodeKeyCount);
				S.close();
				R.close();
				conn.close();
			}
			catch(SQLException e){
				e.printStackTrace();
				System.err.println("Couldn't execute query: "+query+"in host="+node+", database="+db);
				System.exit(1);
			}
		}
		// getting hashmap of counts of keys from master
		HashMap<String,Integer> countsMaster= new HashMap<String,Integer>();
		for (String node:nodes){
			String queryM="SELECT count(*) FROM "+masterKeyTable+" AS a,clients_names as c WHERE a.client_id=c.client_id AND c.client_name='"+node+"';";
			try {
				MySQLConnection conn = this.getConnectionToMasterKeyDb();
				Statement S=conn.createStatement();
				ResultSet R=S.executeQuery(queryM);
				if (R.next()){
					countsMaster.put(node,R.getInt(1));
				}
				S.close();
				R.close();
				conn.close();
			}
			catch(SQLException e){
				e.printStackTrace();
				System.err.println("Couldn't execute query in host="+MASTER+", database="+KEYMASTERDB);
				System.exit(1);
			}
		}
		//compare the two hashmaps of key counts
		for (String node:countsMaster.keySet()){
			int masterCount=countsMaster.get(node);
			int[] thisNodeCount=countsNodes.get(node);
			if (thisNodeCount[0]!=thisNodeCount[1]) {
				System.out.println("Key count and distinct key count do not coincide for key "+keyColumn+" in node "+node+". Key count="+thisNodeCount[0]+", distinct key count="+thisNodeCount[1]);
				checkResult=false;
			}
			else if (thisNodeCount[0]!=masterCount) {
				System.out.println("Key counts do not coincide for key "+keyColumn+" in master and node "+node+". MASTER COUNT="+masterCount+", NODE COUNT="+thisNodeCount[0]);
				System.out.print("Differing "+keyColumn+"'s are: ");
				int[] diffKeys = getDifferingKeys(key,node);
				for (int k:diffKeys){
					System.out.print(k+" ");
				}
				System.out.println();
				checkResult=false;
			}
			else {
				System.out.println("Key counts check passed for key "+keyColumn+" in node "+node+". The count is: "+masterCount);
			}
		}				
		return checkResult;
	}
	
	/**
	 * Method to get the differing keys for a certain key and node. Used by checkKeycounts method. Shouldn't be used out of this class.
	 * @param key the key name
	 * @param node the host name of the cluster node
	 * @return array of ints with all differing keys for this key and node
	 */
	public int[] getDifferingKeys (String key,String node) {
		ArrayList<Integer> diffKeys = new ArrayList<Integer>();
		int[] diffKeysAr;
		String keyTable=key+"_list";
		String masterKeyTable=db+"_"+key+"_list_master";
		String keyColumn=key+"_id";
		String query="SELECT DISTINCT "+keyColumn+" FROM "+keyTable+" ORDER BY "+keyColumn+";";
		MySQLConnection mconn=null;
		try {				
			MySQLConnection nconn = this.getConnectionToNode(node);
			mconn = this.getConnectionToMasterKeyDb();
			Statement S=nconn.createStatement();
			ResultSet R=S.executeQuery(query);
			mconn.executeSql("CREATE TEMPORARY TABLE tmp_keys ("+keyColumn+" int(11) default NULL) ENGINE=MEMORY;");
			int thisKey=0;
			while (R.next()){
				thisKey=R.getInt(1);
				query="INSERT INTO tmp_keys VALUES ("+thisKey+");";
				mconn.executeSql(query);
			}
			S.close();
			R.close();
			nconn.close();
		}
		catch(SQLException e){
			e.printStackTrace();
			System.err.println("Couldn't execute query: "+query+"in host="+node+", database="+db);
			System.exit(1);
		}
		try {	
			query="SELECT c.k " +
					"FROM " +
						"(SELECT u.id AS k,count(u.id) AS cnt " +
							"FROM " +
							"(SELECT "+keyColumn+" AS id FROM tmp_keys UNION ALL SELECT kt."+keyColumn+" AS id FROM "+masterKeyTable+" AS kt LEFT JOIN clients_names AS cn ON kt.client_id=cn.client_id WHERE cn.client_name='"+node+"') AS u GROUP BY u.id) " +
					"AS c " +
					"WHERE c.cnt=1;";
			Statement S=mconn.createStatement();
			ResultSet R=S.executeQuery(query);
			while (R.next()){
				diffKeys.add(R.getInt(1));
			}
			S.close();
			R.close();
		}
		catch(SQLException e){
			e.printStackTrace();
			System.err.println("Couldn't execute query: "+query+"in host="+MASTER+", database="+KEYMASTERDB);
			System.exit(1);
		}
		diffKeysAr= new int[diffKeys.size()];
		for (int i=0;i<diffKeys.size();i++) {
			diffKeysAr[i]=diffKeys.get(i);
		}
		return diffKeysAr;
	}
	
	/**
	 * Method used by splitIdsIntoSet method. To get all ordered ids from a certain key and table from this db in master server
	 * @param key the key name
	 * @param table the table name
	 * @return int array containing all ids 
	 */
	public int[] getAllIds4KeyAndTable(String key, String table){
		int[] allIds=null;
		String keyColumn=key+"_id";		
		try {
			MySQLConnection conn=this.getConnectionToMaster();
			Statement S=conn.createStatement();
			String query="SELECT DISTINCT "+keyColumn+" FROM "+table+" ORDER BY "+keyColumn+";";
			ResultSet R=S.executeQuery(query);
			ArrayList<Integer> idsAL=new ArrayList<Integer>();
			while (R.next()){
				idsAL.add(R.getInt(1));
			}
			allIds=new int[idsAL.size()];
			for (int i=0;i<idsAL.size();i++) {
				allIds[i]=idsAL.get(i);
			}
			R.close();
			S.close();
			conn.close();
		}
		catch (SQLException e){
			e.printStackTrace();
		}
		return allIds;
	}
	
	/**
	 * For a certain key and table returns a HashMap containing an int array per cluster node
	 * @param key
	 * @param table
	 * @return HashMap, keys are node names, values: int array with the ids for each node
	 */
	public HashMap<String,int[]> splitIdsIntoSets(String key, String table){
		HashMap<String,int[]> idSets =new HashMap<String,int[]>();
		String[] nodes=DataDistribution.getNodes();
		int numNodes=nodes.length;
		int[] allIds=this.getAllIds4KeyAndTable(key,table);
		int numIds=allIds.length;
		int setSize=numIds/numNodes;
		int remainder=numIds%numNodes;
		for (int i=0;i<numNodes;i++){
			if (i<remainder){ // for the first "remainder" number of nodes we put setSize+1 ids in the node
				int[] thisnodeidset=new int[setSize+1];
				for (int j=0;j<thisnodeidset.length;j++){
					thisnodeidset[j]=allIds[j+i*(setSize+1)];
				}
				idSets.put(nodes[i],thisnodeidset);
			} else {         // for the rest we put only setSize ids
				int[] thisnodeidset=new int[setSize]; 
				for (int j=0;j<thisnodeidset.length;j++){
					thisnodeidset[j]=allIds[j+remainder*(setSize+1)+(i-remainder)*setSize];
				}
				idSets.put(nodes[i],thisnodeidset); 
			}
		}		
		return idSets;
	}

	/**
	 * To split a given table in chunks based on a key, split tables remain in same database and server
	 * @param key
	 * @param table 
	 */
	public void splitTable (String key,String table){
		String keyColumn=key+"_id";
		String query;
		HashMap<String,int[]> idSets = this.splitIdsIntoSets(key,table);
		String[] splitTables=new String[idSets.size()];
		try {
			MySQLConnection conn=this.getConnectionToMaster();
			int i=0;
			for (String node:idSets.keySet()) {
				String splitTbl=table+"_split_"+node;
				splitTables[i]=splitTbl;
				i++;
				// we create permanent tables
				query="CREATE TABLE "+splitTbl+" LIKE "+table+";";
				conn.executeSql(query);
				// drop the indexes if there was any, indexes will slow down the creation of split tables
				String[] indexes=conn.getAllIndexes4Table(table);
				for (String index:indexes) { 
					conn.executeSql("DROP INDEX "+index+" ON "+splitTbl+";");
				}
				int idmin=idSets.get(node)[0];
				int idmax=idSets.get(node)[idSets.get(node).length-1];
				query="INSERT INTO "+splitTbl+" SELECT * FROM "+table+" WHERE "+keyColumn+">="+idmin+" AND "+keyColumn+"<="+idmax+";";				
				conn.executeSql(query);
				//TODO recreate indexes, use method getCreateIndex4Table from MySQLConnection
			}
			conn.close();
		}
		catch (SQLException e){
			e.printStackTrace();
		}
	}

	/**
	 * To split a given table in chunks based on a key, split tables go to different nodes of cluster 
	 * @param key
	 * @param table
	 * @param destDb name of destination db 
	 */
	public void splitTableToCluster (String key,String table, String destDb){
		String keyColumn=key+"_id";
		String query;
		HashMap<String,int[]> idSets = this.splitIdsIntoSets(key,table);
		String[] splitTables=new String[idSets.size()];
		try {
			MySQLConnection conn=this.getConnectionToMaster();
			int i=0;
			for (String node:idSets.keySet()) {
				String splitTbl=table+"_split_"+node;
				splitTables[i]=splitTbl;
				i++;
				// we create permanent tables, later we drop them. Can't be temporary as we use another connection for dumpData
				query="CREATE TABLE "+splitTbl+" LIKE "+table+";";
				conn.executeSql(query);
				// drop the indexes if there was any, indexes will slow down the creation of split tables
				String[] indexes=conn.getAllIndexes4Table(table);
				for (String index:indexes) { 
					conn.executeSql("DROP INDEX "+index+" ON "+splitTbl+";");
				}
				// make the table a memory table (won't be feasible in general case where tables can be VERY big, even white won't cope) 
				//query="ALTER TABLE "+splitTbl+" TYPE=MEMORY;";
				//conn.executeSql(query);
				int idmin=idSets.get(node)[0];
				int idmax=idSets.get(node)[idSets.get(node).length-1];
				query="INSERT INTO "+splitTbl+" SELECT * FROM "+table+" WHERE "+keyColumn+">="+idmin+" AND "+keyColumn+"<="+idmax+";";				
				conn.executeSql(query);
			}
			// transfering data across
			String[] srchosts={MASTER};
			String[] desthosts=getNodes();
			dumpData(srchosts,splitTables);
			// using here loadSplitData rather than loadData because table names are not the same on source and destination, i.e. source: table_split_tla01, dest: table
			loadSplitData(MASTER,desthosts,table,destDb);
			// droping table, we don't want them anymore after loading data to nodes
			for (String tbl:splitTables){					
				query="DROP TABLE "+tbl+";";
				conn.executeSql(query);
			}
			conn.close();
		}
		catch (SQLException e){
			e.printStackTrace();
		}
	}

	/**
	 * Executes a query in all nodes in cluster.
	 * TODO Right now it is serial, must parallelize this with threads
	 * TODO This can be used in lots of methods in this class (all the loadData and dumpData ones)
	 * @param query
	 */
	public void clusterExecuteQuery(String query){
		String[] nodes = getNodes();
		for (String node: nodes){
			try {
				MySQLConnection conn = this.getConnectionToNode(node);
				conn.executeSql(query);
				conn.close();
			}
			catch(SQLException e){
				e.printStackTrace();
				System.err.println("Couldn't execute query="+query+", in node="+node);
				System.exit(1);
			}
		}

	}
}

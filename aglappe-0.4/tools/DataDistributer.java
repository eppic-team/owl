package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;

public class DataDistributer {
	
	public final static String GLOBALDIR="/project/snow/global/tmp";
	public final static String MASTER = DataDistribution.MASTER;
	public final static String KEYMASTERDB = DataDistribution.KEYMASTERDB;
	
	// Following parameters are optimized after some testing. I found these to be good values
	// They represent the number of nodes reading/writing from/to nfs concurrently when loading/dumping
	// 1 means no concurrency at all. Higher numbers mean more concurrency i.e. more nodes writing/reading to/from nfs concurrently
	final static int NUM_CONCURRENT_READ_QUERIES  = 9;
	final static int NUM_CONCURRENT_WRITE_QUERIES = 6;
	// Following parameter only for concurrent writing queries in same host. This number has to be low
	// as the only limitation is the i/o capacity of the host. 
	final static int NUM_CONCURRENT_SAMEHOST_QUERIES = 3; 
	
	String dumpdir;
	String srcDb;
	String destDb;
	String user;
	String pwd;
	boolean debug = false;
	String rmvtmp = "force"; // 3 possible values: force, noremove and prompt
	
	/**
	 * Construct DataDistributer object passing only source db. The destination db will be same as source.
	 * @param srcDb
	 * @param user
	 * @param pwd
	 */
	public DataDistributer(String srcDb, String user, String pwd) {
		this.srcDb=srcDb;
		this.destDb=srcDb;
		this.user=user;
		this.pwd=pwd;
		initializeDirs();
	}

	/**
	 * Construct DataDistributer object passing both source db and destination db
	 * @param srcDb
	 * @param destDb
	 * @param user
	 * @param pwd
	 */
	public DataDistributer(String srcDb, String destDb, String user, String pwd) {
		this.srcDb=srcDb;
		this.destDb=destDb;
		this.user=user;
		this.pwd=pwd;
		initializeDirs();
	}	
	
	public void setDebug(boolean debug){
		this.debug=debug;
	}
	
	public void setRmv(String rmvtmp){
		this.rmvtmp=rmvtmp;
	}
	
	public void setDumpDir(String dumpdir){
		// before reseting the dumpdir, we get rid of the existent dumpdir created when this DataDistributer object was constructed
		// in the case where we set the dumpdir using setDumpDir clearly the randomly named dumpdir created while constructing is not needed
		SystemCmd.exec("rmdir "+this.dumpdir);
		this.dumpdir=dumpdir;
	}
	
	public void initializeDirs() {
		dumpdir=GLOBALDIR+"/dumps_tmp_"+System.currentTimeMillis();
	    if (!((new File(dumpdir)).mkdir())) {
	    	System.err.println("Couldn't create directory "+dumpdir);
	    	System.exit(1);
	    }
		SystemCmd.exec("chmod go+rw "+dumpdir);
	}
	
	public void initializeDirs(String[] hosts){
		for (String host: hosts) {
			if (!((new File(dumpdir+"/"+host+"/"+srcDb)).mkdirs())) {
				System.err.println("Couldn't create directory "+dumpdir+"/"+host+"/"+srcDb);
				System.exit(1);
			}
			SystemCmd.exec("chmod -R go+rw "+dumpdir+"/"+host);
		}
	}
	
	public void finalizeDirs() {
	    if (debug) {
	    	System.out.println("Temporary directory "+dumpdir+" was not removed. You must remove it manually");
	    } else {
	    	if (rmvtmp.equals("force")) {
	    		System.out.println("Removing temporary directory "+dumpdir);
	    		//TODO must capture exit state and print to error if problems deleting dir
	    		SystemCmd.exec("rm -rf "+dumpdir);
	    	} else if (rmvtmp.equals("prompt")){
	    		System.out.println("Would you like to remove the temporary data directory '"+dumpdir+"' ? (y/n)");
	    		BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
	    		try {
	    			String answer = br.readLine();
	    			if (answer.equals("y") || answer.equals("yes") || answer.equals("Y") || answer.equals("YES") || answer.equals("Yes")) {
	    	    		System.out.println("Removing temporary directory "+dumpdir);
	    	    		SystemCmd.exec("rm -rf "+dumpdir);
	    			} else {
	    		    	System.out.println("Temporary directory "+dumpdir+" was not removed.");
	    			}
	    		} catch (IOException e) {
	    			System.err.println("I/O error while reading user input.");
	    			e.printStackTrace();
	    			System.exit(2);
	    		}
	    	} else if (rmvtmp.equals("noremove")){
	    		System.out.println("Temporary directory "+dumpdir+" was not removed.");
	    	}
	    }
	}
	
	private MySQLConnection getConnectionToNode(String node) {
		MySQLConnection conn=new MySQLConnection(node,user,pwd,srcDb);
		return conn;
	}
	
	private MySQLConnection getConnectionToMaster() {
		MySQLConnection conn=this.getConnectionToNode(MASTER);
		return conn;
	}
	
	private MySQLConnection getConnectionToMasterKeyDb() {
		MySQLConnection conn=new MySQLConnection(MASTER,user,pwd,KEYMASTERDB);
		return conn;
	}

	public void dumpData(String[] srchosts, String[] tables) {
		int concurrentQueries = 1;
		if (srchosts.length>1) {
			concurrentQueries = NUM_CONCURRENT_WRITE_QUERIES;
		}
		int i = 0;
		QueryThread[] qtGroup = new QueryThread[concurrentQueries];
		initializeDirs(srchosts);
		for ( String tbl: tables) {
			for (String srchost: srchosts) {				
				String outfile=dumpdir+"/"+srchost+"/"+srcDb+"/"+tbl+".txt";
				String wherestr="";
				String sqldumpstr="SELECT * FROM `"+tbl+"` "+wherestr+" INTO OUTFILE '"+outfile+"';";
				if (debug) {System.out.println ("HOST="+srchost+", database="+srcDb+", sqldumpstr="+sqldumpstr);} 
				else {
					if (i!=0 && i%(concurrentQueries) == 0) {
						try {
							for (int j = 0;j<qtGroup.length;j++){
								qtGroup[j].join(); // wait until previous thread group is finished
							}
						} catch (InterruptedException e) {
							e.printStackTrace();
						}
						i = 0;
						qtGroup = new QueryThread[concurrentQueries];
					}			
					qtGroup[i] = new QueryThread(sqldumpstr,srchost,user,pwd,srcDb);	    				
				} //end else if debug
				i++;
			} // end foreach srchost
		} // end foreach table
		try {
			for (int j = 0;j<qtGroup.length;j++){
				if (qtGroup[j]!=null){ // some slots of the array may be null, check for those before trying the join
					qtGroup[j].join(); // wait until the last thread group is finished
				}
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}		
		System.out.println ("Dump finished.");
	}

	public <T> void dumpSplitData(String srchost, String[] tables, String key, HashMap<String,T[]> idSets) {
		int concurrentQueries = NUM_CONCURRENT_SAMEHOST_QUERIES; 
		int i = 0;
		QueryThread[] qtGroup = new QueryThread[concurrentQueries];
		String[] srchosts = {srchost};
		initializeDirs(srchosts);
		for ( String tbl: tables) {
			for (String node:idSets.keySet()) {
				String outfile=dumpdir+"/"+srchost+"/"+srcDb+"/"+tbl+"_split_"+node+".txt";
				String wherestr="";
				// if number of ids for the key was less than number of nodes, some of the node slots will be empty. We check for those before trying to access the array
				if (idSets.get(node).length>0) {
					T idmin= idSets.get(node)[0];
					T idmax= idSets.get(node)[idSets.get(node).length-1];				
					wherestr="WHERE "+key+">='"+idmin+"' AND "+key+"<='"+idmax+"'";
					String sqldumpstr="SELECT * FROM `"+tbl+"` "+wherestr+" INTO OUTFILE '"+outfile+"';";
					if (debug) {System.out.println ("HOST="+srchost+", database="+srcDb+", sqldumpstr="+sqldumpstr);} 
					else {
						if (i!=0 && i%(concurrentQueries) == 0) {
							try {
								for (int j = 0;j<qtGroup.length;j++){
									qtGroup[j].join(); // wait until previous thread group is finished
								}
							} catch (InterruptedException e) {
								e.printStackTrace();
							}
							i = 0;
							qtGroup = new QueryThread[concurrentQueries];
						}			
						qtGroup[i] = new QueryThread(sqldumpstr,srchost,user,pwd,srcDb);	    				
					} //end else if debug
					i++;
				}
			} // end foreach node			
		} // end foreach table
		try {
			for (int j = 0;j<qtGroup.length;j++){
				if (qtGroup[j]!=null){ // some slots of the array may be null, check for those before trying the join
					qtGroup[j].join(); // wait until the last thread group is finished
				}
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}		
		System.out.println ("Dump finished.");
	}

	public void loadData(String[] srchosts, String[] desthosts,String[] tables) {
		int concurrentQueries = 1;
		if (srchosts.length==1 && desthosts.length>1) {
			concurrentQueries = NUM_CONCURRENT_READ_QUERIES;
		}
		int i = 0;
		QueryThread[] qtGroup = new QueryThread[concurrentQueries];
		for (String srchost:srchosts){
			for (String tbl:tables) {
				String dumpfile=dumpdir+"/"+srchost+"/"+srcDb+"/"+tbl+".txt";
				String sqlloadstr="LOAD DATA INFILE '"+dumpfile+"' INTO TABLE `"+tbl+"`;";
				for (String desthost: desthosts) {
					if (debug) {System.out.println ("SRCHOST="+srchost+", DESTHOST="+desthost+", DESTDB="+destDb+", sqlloadstr="+sqlloadstr);} 
					else {						
						if (i!=0 && i%(concurrentQueries) == 0) {
							try {
								for (int j = 0;j<qtGroup.length;j++){
									qtGroup[j].join(); // wait until previous thread group is finished
								}
							} catch (InterruptedException e) {
								e.printStackTrace();
							}
							i = 0;
							qtGroup = new QueryThread[concurrentQueries];
						}			
						qtGroup[i] = new QueryThread(sqlloadstr,desthost,user,pwd,destDb);	    				
					}	    						
					i++;
				} // end foreach desthost
			} // end foreach tbl
		} // end foreach srchost
		try {
			for (int j = 0;j<qtGroup.length;j++){
				if (qtGroup[j]!=null){ // some slots of the array may be null, check for those before trying the join
					qtGroup[j].join(); // wait until the last thread group is finished
				}
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}		
	    System.out.println ("Load finished.");
	    finalizeDirs();
	}

	public void loadSplitData(String srchost, String[] desthosts, String tableName) {
		int concurrentQueries = 1;
		if (desthosts.length>1) {
			concurrentQueries = NUM_CONCURRENT_READ_QUERIES;
		}
		int i = 0;
		QueryThread[] qtGroup = new QueryThread[concurrentQueries];
		for (String desthost:desthosts) {
			String dumpfile=dumpdir+"/"+srchost+"/"+srcDb+"/"+tableName+"_split_"+desthost+".txt";
			// we test first if dumpfile exists, if not all the slots (nodes) are filled then some nodes won't have to load anything and the files won't be there because dumpSplitData didn't create them 
			if ((new File(dumpfile)).exists()) { 
				String sqlloadstr="LOAD DATA INFILE '"+dumpfile+"' INTO TABLE `"+tableName+"`;";
				if (debug) {System.out.println ("SRCHOST="+srchost+", DESTHOST="+desthost+", DESTDB="+destDb+", sqlloadstr="+sqlloadstr);} 
				else {
					if (i!=0 && i%(concurrentQueries) == 0) {
						try {
							for (int j = 0;j<qtGroup.length;j++){
								qtGroup[j].join(); // wait until previous thread group is finished
							}
						} catch (InterruptedException e) {
							e.printStackTrace();
						}
						i = 0;
						qtGroup = new QueryThread[concurrentQueries];
					}			
					qtGroup[i] = new QueryThread(sqlloadstr,desthost,user,pwd,destDb);
				} // end else if debug
				i++;
			}
		} // end foreach desthosts
		try {
			for (int j = 0;j<qtGroup.length;j++){
				if (qtGroup[j]!=null){ // some slots of the array may be null, check for those before trying the join
					qtGroup[j].join(); // wait until the last thread group is finished
				}
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
		}		
	    System.out.println ("Load finished.");
	    finalizeDirs();
	}

	/**
	 * For a certain key (text or numeric) and table returns a "data distribution" (kind of evenly distributed) of the data to the nodes
	 * To be used when we have a table that we are going to split to the nodes
	 * TODO eventually the code could be cleverer so that the data is actually evenly distributed, right now is only evenly distributed on key ids
	 * @param key
	 * @param table
	 * @return idSets HashMap, keys are node names, values: Integer/String array with the ids for each node
	 */
	public HashMap<String,Object[]> splitIdsIntoSets(String key, String table){
		HashMap<String,Object[]> idSets =new HashMap<String,Object[]>();		
		String[] nodes=DataDistribution.getMySQLNodes();
		int numNodes=nodes.length;
		MySQLConnection conn = this.getConnectionToMaster();
		Object[] allIds=conn.getAllIds4KeyAndTable(key,table);
		conn.close();
		int numIds=allIds.length;
		int setSize=numIds/numNodes;
		int remainder=numIds%numNodes;
		for (int i=0;i<numNodes;i++){
			if (i<remainder){ // for the first "remainder" number of nodes we put setSize+1 ids in the node
				Object[] thisnodeidset=new Object[setSize+1];
				for (int j=0;j<thisnodeidset.length;j++){
					thisnodeidset[j]=allIds[j+i*(setSize+1)];
				}
				idSets.put(nodes[i],thisnodeidset);
			} else {         // for the rest we put only setSize ids
				Object[] thisnodeidset=new Object[setSize]; 
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
		String[] nodes=DataDistribution.getMySQLNodes();
		MySQLConnection conn=this.getConnectionToMaster();
		String[] splitTables = new String[nodes.length]; // we create an array that will contain the name of all split tables
		String[] indexes=conn.getAllIndexes4Table(table);
		try {
			// we create split tables and drop indexes before inserting
			for (int i=0;i<nodes.length;i++) {
				String splitTbl=table+"_split_"+nodes[i];
				splitTables[i]=splitTbl; // will be used later while looping over splitTables
				// we create permanent tables
				String query="CREATE TABLE "+splitTbl+" LIKE "+table+";";
				conn.executeSql(query);
				// drop the indexes if there was any, indexes will slow down the creation of split tables
				for (String index:indexes) { 
					conn.executeSql("DROP INDEX "+index+" ON "+splitTbl+";");
				}
			}			
			HashMap<String,Object[]> idSets = this.splitIdsIntoSets(key,table);
			for (int i=0;i<nodes.length;i++) {
				Object idmin=idSets.get(nodes[i])[0];
				Object idmax=idSets.get(nodes[i])[idSets.get(nodes[i]).length-1];
				String query="INSERT INTO "+splitTables[i]+" SELECT * FROM "+table+" WHERE "+key+">='"+idmin+"' AND "+key+"<='"+idmax+"';";				
				conn.executeSql(query);
				//TODO recreate indexes, use method getCreateIndex4Table from MySQLConnection
			}				
		}
		catch (SQLException e){
			e.printStackTrace();
		}
		conn.close();
	}

	/**
	 * To split a given table in chunks based on a key, split tables go to different nodes of cluster 
	 * @param key
	 * @param table 
	 */
	public DataDistribution splitTableToCluster (String key,String table){
		System.out.println("Splitting table "+table+" to cluster based on key "+key+"...");
		String[] tables={table};
		String[] desthosts=DataDistribution.getMySQLNodes();
		HashMap<String,Object[]> idSets = this.splitIdsIntoSets(key,table);
		// dumping data with the dumpSplitData method, a modified version of dumpData
		dumpSplitData(MASTER,tables,key,idSets);
		// putting the ids in the key_master database so we keep track of where everything is
		insertIdsToKeyMaster(key,table,idSets);
		// using here loadSplitData rather than loadData because table names are not the same on source and destination, 
		// i.e. source: table_split_tla01, dest: table
		loadSplitData(MASTER,desthosts,table);	
		DataDistribution dataDist = new DataDistribution(destDb,user,pwd);
		System.out.println("Done with splitting.");
		return dataDist;
	}
	
	/**
	 * Insert all ids to the key_master database creating a new table for this destDb/given table combination if not exists
	 * @param key name of key on which distribution of table is based
	 * @param table name of table that we are distributing
	 * @param idSets as returned from splitIdsIntoSets or getIdSetsFromNodes from a DataDistribution object
	 */
	public <T> void insertIdsToKeyMaster(String key,String table,HashMap<String,T[]> idSets) {
		System.out.println("Updating key_master database with ids to nodes mapping...");
		MySQLConnection conn = this.getConnectionToMasterKeyDb();		
		String keyMasterTbl = createNewKeyMasterTbl(key,table);
		removePK(keyMasterTbl,key); // attention removing primary keys, duplicates won't be checked!!!
		// getting first mapping between nodes names and node ids
		HashMap<String,Integer> nodes2nodeids = new HashMap<String,Integer>();
		for (String node:idSets.keySet()){
			String query="SELECT client_id FROM clients_names WHERE client_name='"+node+"';";
			int id = conn.getIntFromDb(query);
			nodes2nodeids.put(node,id);
		}
		for (String node:idSets.keySet()){
				T[] thisNodeIds=idSets.get(node);				
				for (T id:thisNodeIds){
					String query="INSERT INTO "+keyMasterTbl+" ("+key+",client_id) VALUES ('"+id+"',"+nodes2nodeids.get(node)+");";
					try {
						conn.executeSql(query);
					} catch (SQLException e) {
						e.printStackTrace();
					}
				}
		}
		conn.close();
		removeZeros(keyMasterTbl,key); // we only have inserted 0s for records that we didn't want, it's safe now to get rid of them
		addPK(keyMasterTbl,key); // if there were duplicates, this should barf
		System.out.println("Done with updating key_master database.");
	}

	/**
	 * To create a new key master table for destination db in the key_master database given a key and table. Used by insertIdsToKeyMaster
	 * Eventually this method on other key_master related should go into their own class, shouldn't they?
	 * @param key
	 * @param table
	 * @return the name of the key master table created
	 */
	public String createNewKeyMasterTbl(String key,String table) {
		// find out whether key is numeric or text and setting accordingly query strings
		String nodes[] = DataDistribution.getMySQLNodes(); // we need the list of nodes only to get one of them no matter which
		MySQLConnection conn=new MySQLConnection(nodes[0],user,pwd,destDb); // here we connect to destDb in one node, needed to getColumnType 
		String colType = conn.getColumnType(table,key);
		String autoIncr = "";
		if (colType.contains("int")){
			autoIncr = "auto_increment";
		}
		conn.close();
		// key master table name and connection to key master db
		String keyMasterTbl=destDb+"__"+table;
		conn=this.getConnectionToMasterKeyDb();
		try {
			String query="CREATE TABLE IF NOT EXISTS "+keyMasterTbl+" ("+
							key+" "+colType+" NOT NULL "+autoIncr+", " +
							"client_id smallint(6) NOT NULL default '0', " +
							"PRIMARY KEY (`"+key+"`) " +
							") ENGINE=MyISAM DEFAULT CHARSET=ascii COLLATE=ascii_bin;";
			Statement S=conn.createStatement();
			S.executeUpdate(query);
			S.close();
		} catch (SQLException e) {
			System.err.println("Couldn't create table "+keyMasterTbl);
			e.printStackTrace();
		}
		try {
			Statement S=conn.createStatement();
			String query="INSERT INTO dbs_keys (key_name,db,key_master_table) VALUES (\'"+key+"\',\'"+destDb+"\',\'"+keyMasterTbl+"\');";
			S.executeUpdate(query);
			S.close();
		} catch (SQLException e) {
			System.err.println("Didn't insert new record into table dbs_keys of database: "+KEYMASTERDB+". The record for key: "+key+", table: "+table+" existed already. This is usually a harmless error!");
			System.err.println("SQLException: " + e.getMessage());
		}
		conn.close();
		return keyMasterTbl;
	}
	
	public void removePK (String keyMasterTbl,String key){
		MySQLConnection conn=this.getConnectionToMasterKeyDb();
		boolean isNumeric = conn.isKeyNumeric(keyMasterTbl,key);
		String colType = conn.getColumnType(keyMasterTbl,key);
		try {
			if (isNumeric){ // removing the auto_increment, only in numeric keys
				String query="ALTER TABLE "+keyMasterTbl+" MODIFY "+key+" "+colType+" NOT NULL;";
				conn.executeSql(query);
			}
			// removing primary key (same sql code for both numeric or text keys
			String query="ALTER TABLE "+keyMasterTbl+" DROP PRIMARY KEY;";
			conn.executeSql(query);			
		} catch (SQLException e) {
			e.printStackTrace();
		}
		conn.close();
	}
	
	public void addPK (String keyMasterTbl, String key){
		MySQLConnection conn=this.getConnectionToMasterKeyDb();
		boolean isNumeric = conn.isKeyNumeric(keyMasterTbl,key);
		String colType = conn.getColumnType(keyMasterTbl,key);
		try {
			// adding primary key (same sql code for both numeric or text keys
			String query="ALTER TABLE "+keyMasterTbl+" ADD PRIMARY KEY("+key+");";				
			conn.executeSql(query);
			if (isNumeric){ // adding auto_increment, only in numeric keys
				query="ALTER TABLE "+keyMasterTbl+" MODIFY "+key+" "+colType+" NOT NULL auto_increment;";
				conn.executeSql(query);
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
		conn.close();		
	}
	
	public void removeZeros (String keyMasterTbl, String key){
		MySQLConnection conn=this.getConnectionToMasterKeyDb();		
		try {
			// attention! the quotes around 0 are very important, otherwise if key is text-based all records get deleted
			// mysql somehow considers all text records = 0, using '0' is ok as is considered as a text literal for char fields and as a number for int fields
			String query="DELETE FROM "+keyMasterTbl+" WHERE "+key+"='0';";				
			conn.executeSql(query);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		conn.close();				
	}

	/**
	 * Executes a query in all nodes in cluster.
	 * Not in use at the moment
	 * TODO Right now it is serial, must parallelize this with threads
	 * @param query
	 */
	public void clusterExecuteQuery(String query){
		String[] nodes = DataDistribution.getMySQLNodes();
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
	
	/**
	 * Executes a query in all nodes in cluster given a HashMap containing a set of queries (one per node)
	 * Not in use at the moment
	 * TODO Right now it is serial, must parallelize this with threads
	 * TODO This can be used in the load/dump methods in this class where queries are different for each node
	 * @param queries a HashMap containing a query per node
	 */
	public void clusterExecuteQuery(HashMap<String,String> queries){
		String[] nodes = DataDistribution.getMySQLNodes();
		for (String node: nodes){
			String query="";
			try {
				query=queries.get(node);
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

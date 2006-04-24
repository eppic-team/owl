package tools;

import java.io.File;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;

public class DataDistributer {
	
	public final static String GLOBALDIR="/project/snow/global/tmp";
	public final static String MASTER = DataDistribution.MASTER;
	public final static String KEYMASTERDB = DataDistribution.KEYMASTERDB;
	
	String dumpdir;
	String srcDb;
	String destDb;
	String user;
	String pwd;
	boolean debug = false;
	
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
	
	public void initializeDirs() {
		dumpdir=GLOBALDIR+"/dumps_tmp_"+System.currentTimeMillis();
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
	
	private MySQLConnection getConnectionToNode(String node) {
		MySQLConnection conn=new MySQLConnection(node,user,pwd,srcDb);
		return conn;
	}
	
	private MySQLConnection getConnectionToNode(String node, String dbName){
		MySQLConnection conn=new MySQLConnection(node,user,pwd,dbName);
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
		for (String srchost: srchosts) {
			if (!((new File(dumpdir+"/"+srchost+"/"+srcDb)).mkdirs())) {
				System.err.println("Couldn't create directory "+dumpdir+"/"+srchost+"/"+srcDb);
				System.exit(1);
			}
			SystemCmd.exec("chmod -R go+rw "+dumpdir+"/"+srchost);
			for ( String tbl: tables) {
				String outfile=dumpdir+"/"+srchost+"/"+srcDb+"/"+tbl+".txt";
				String wherestr="";
				String sqldumpstr="SELECT * FROM `"+tbl+"` "+wherestr+" INTO OUTFILE '"+outfile+"';";
				if (debug) {System.out.println ("HOST="+srchost+", database="+srcDb+", sqldumpstr="+sqldumpstr);} 
				else {
					try {
						MySQLConnection conn = this.getConnectionToNode(srchost);
						Statement S=conn.createStatement();
						S.executeQuery(sqldumpstr);
						System.out.println ("Dumped from host="+srchost+", database="+srcDb+", table="+tbl+" to outfile="+outfile);
						S.close();
						conn.close();
					}
					catch(SQLException e){
						e.printStackTrace();
						System.err.println("Couldn't dump from host="+srchost+", database="+srcDb+", table="+tbl+" to outfile="+outfile);
						System.exit(1);
					}		
				} //end else if debug
			} // end foreach tbl
		}
		System.out.println ("Dump finished.");
	}
	
	public void loadData(String[] srchosts, String[] desthosts,String[] tables) {
	    for (String desthost: desthosts) {
	    	for (String srchost: srchosts) {
	    		for (String tbl:tables) {
	    			String dumpfile=dumpdir+"/"+srchost+"/"+srcDb+"/"+tbl+".txt";
	    			String sqlloadstr="LOAD DATA INFILE '"+dumpfile+"' INTO TABLE `"+tbl+"`;";
	    			if (debug) {System.out.println ("SRCHOST="+srchost+", DESTHOST="+desthost+", DESTDB="+destDb+", sqlloadstr="+sqlloadstr);} 
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
	    } 
	    System.out.println ("Load finished.");
	    finalizeDirs();
	}

	public void loadSplitData(String srchost, String[] desthosts, String tableName) {
		for (String desthost:desthosts) {
			String dumpfile=dumpdir+"/"+srchost+"/"+srcDb+"/"+tableName+"_split_"+desthost+".txt";
			String sqlloadstr="LOAD DATA INFILE '"+dumpfile+"' INTO TABLE `"+tableName+"`;";
			if (debug) {System.out.println ("SRCHOST="+srchost+", DESTHOST="+desthost+", DESTDB="+destDb+", sqlloadstr="+sqlloadstr);} 
			else {
				try {
					MySQLConnection conn = this.getConnectionToNode(desthost,destDb);
					conn.executeSql(sqlloadstr);
					System.out.println ("HOST: "+desthost+". Loaded from file="+dumpfile+", database="+srcDb+", table="+tableName);
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

	/**
	 * For a certain key and table returns a "data distribution" (kind of evenly distributed) of the data to the nodes
	 * To be used when we have a table that we are going to split to the nodes
	 * TODO eventually the code could be cleverer so that the data is actually evenly distributed, right now is only evenly distributed on key ids
	 * @param key
	 * @param table
	 * @return idSets HashMap, keys are node names, values: int array with the ids for each node
	 */
	public HashMap<String,int[]> splitIdsIntoSets(String key, String table){
		HashMap<String,int[]> idSets =new HashMap<String,int[]>();
		String[] nodes=DataDistribution.getNodes();
		int numNodes=nodes.length;
		MySQLConnection conn = this.getConnectionToMaster();
		int[] allIds=conn.getAllIds4KeyAndTable(key,table);
		conn.close();
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
				query="INSERT INTO "+splitTbl+" SELECT * FROM "+table+" WHERE "+key+">="+idmin+" AND "+key+"<="+idmax+";";				
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
	 */
	public DataDistribution splitTableToCluster (String key,String table){
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
				query="INSERT INTO "+splitTbl+" SELECT * FROM "+table+" WHERE "+key+">="+idmin+" AND "+key+"<="+idmax+";";				
				conn.executeSql(query);
			}
			// transfering data across
			String[] srchosts={MASTER};
			String[] desthosts=DataDistribution.getNodes();
			dumpData(srchosts,splitTables);
			// using here loadSplitData rather than loadData because table names are not the same on source and destination, i.e. source: table_split_tla01, dest: table
			loadSplitData(MASTER,desthosts,table);
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
		// putting the ids in the key_master database so we keep track of where everything is
		insertIdsToKeyMaster(key,table,idSets);
		DataDistribution dataDist = new DataDistribution(destDb,user,pwd);
		return dataDist;
	}
	
	/**
	 * Insert all ids to the key_master database creating a new table for this destDb/given table combination if not exists
	 * @param key name of key on which distribution of table is based
	 * @param table name of table that we are distributing
	 * @param idSets as returned from splitIdsIntoSets or getIdSetsFromNodes from a DataDistribution object
	 */
	public void insertIdsToKeyMaster(String key,String table,HashMap<String,int[]> idSets) {
		MySQLConnection conn = this.getConnectionToMasterKeyDb();		
		String keyMasterTbl = createNewKeyMasterTbl(key,table);
		removePK(keyMasterTbl,key); // attention removing primary keys, duplicates won't be checked!!!
		for (String node:idSets.keySet()){
			int[] thisNodeIds=idSets.get(node);
			for (int id:thisNodeIds){
				String query="INSERT INTO "+keyMasterTbl+" ("+key+",client_id) " +
							"SELECT "+id+",c.client_id FROM clients_names AS c WHERE client_name='"+node+"';";
				try {
					conn.executeSql(query);
				} catch (SQLException e) {
					e.printStackTrace();
				}
			}
		}
		removeZeros(keyMasterTbl,key); // we only have inserted 0s for records that we didn't want, it's safe now to get rid of them
		addPK(keyMasterTbl,key); // if there were duplicates, this should barf
	}

	/**
	 * To create a new key master table for destination db in the key_master database given a key and table. Used by insertIdsToKeyMaster
	 * Eventually this method on other key_master related should go into their own class, shouldn't they?
	 * @param key
	 * @param table
	 * @return the name of the key master table created
	 */
	public String createNewKeyMasterTbl(String key,String table) {
		String keyMasterTbl=destDb+"__"+table;
		MySQLConnection conn=this.getConnectionToMasterKeyDb();
		try {
			String query="CREATE TABLE IF NOT EXISTS "+keyMasterTbl+" ("+
						key+" int(11) NOT NULL auto_increment, " +
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
			String query="INSERT INTO dbs_keys (key_name,db,key_master_table) VALUES (\'"+key+"\',\'"+srcDb+"\',\'"+keyMasterTbl+"\');";
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
		try {
			String query="ALTER TABLE "+keyMasterTbl+" MODIFY "+key+" int(11) NOT NULL default '0';";
			conn.executeSql(query);
			query="ALTER TABLE "+keyMasterTbl+" DROP PRIMARY KEY;";
			conn.executeSql(query);			
		} catch (SQLException e) {
			e.printStackTrace();
		}
		conn.close();
	}
	
	public void addPK (String keyMasterTbl, String key){
		MySQLConnection conn=this.getConnectionToMasterKeyDb();		
		try {
			String query="ALTER TABLE "+keyMasterTbl+" ADD PRIMARY KEY("+key+");";				
			conn.executeSql(query);
			query="ALTER TABLE "+keyMasterTbl+" MODIFY "+key+" int(11) NOT NULL auto_increment;";
			conn.executeSql(query);			
		} catch (SQLException e) {
			e.printStackTrace();
		}
		conn.close();		
	}
	
	public void removeZeros (String keyMasterTbl, String key){
		MySQLConnection conn=this.getConnectionToMasterKeyDb();		
		try {
			String query="DELETE FROM "+keyMasterTbl+" WHERE "+key+"=0;";				
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
		String[] nodes = DataDistribution.getNodes();
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
		String[] nodes = DataDistribution.getNodes();
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

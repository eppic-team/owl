package tools;

import java.io.*;
import java.sql.*;
import java.util.ArrayList;
import java.util.HashMap;

public class DataDistribution {

	final static String ADMINDIR="/project/StruPPi/Cluster/admin";
	final static String HOSTSFILE=ADMINDIR+"/hosts_ss.txt";
	public final static String MASTER="white";
	final static String KEYMASTERDB="key_master";

	public String db;
	private String user;
	private String pwd;
	public boolean isSplit; // for future use 
	

	public DataDistribution(String db,String user,String pwd) {
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
	
	public static String[] getNodes() {
		String[] nodes=null;	    
	    try {
			File inputFile = new File(HOSTSFILE);
	    	BufferedReader hostsFile = new BufferedReader(new FileReader(inputFile)); // open BufferedReader to the file
	    	String nodesstr=hostsFile.readLine();
	    	nodes=nodesstr.split(" ");
	    	hostsFile.close();
	    }
	    catch (IOException e){
	    	e.printStackTrace();
	    	System.err.println("Couldn't read from file "+HOSTSFILE);
	    }
	    return nodes;
	}
	
	public boolean checkCountsAllTables (){
		boolean checkResult=true;
		String[] nodes = getNodes();
		MySQLConnection mconn = this.getConnectionToMaster();
		String[] tables = mconn.getTables4Db();
		mconn.close();
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
		ClusterConnection cconn = new ClusterConnection(db,key,user,pwd);
		String keyTable=cconn.getTableOnNode();
		String masterKeyTable=cconn.getKeyTable();
		cconn.close();
		// getting hashmap of counts of keys from nodes
		String[] nodes = getNodes();
		HashMap<String,int[]> countsNodes=new HashMap<String,int[]>();
		String query="SELECT count("+key+"),count(DISTINCT "+key+") FROM "+keyTable+";";
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
				System.out.println("Key count and distinct key count do not coincide for key "+key+" in node "+node+". Key count="+thisNodeCount[0]+", distinct key count="+thisNodeCount[1]);
				checkResult=false;
			}
			else if (thisNodeCount[0]!=masterCount) {
				System.out.println("Key counts do not coincide for key "+key+" in master and node "+node+". MASTER COUNT="+masterCount+", NODE COUNT="+thisNodeCount[0]);
				System.out.print("Differing "+key+"'s are: ");
				int[] diffKeys = getDifferingKeys(key,node);
				for (int k:diffKeys){
					System.out.print(k+" ");
				}
				System.out.println();
				checkResult=false;
			}
			else {
				System.out.println("Key counts check passed for key "+key+" in node "+node+". The count is: "+masterCount);
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
		ClusterConnection cconn = new ClusterConnection(db,key,user,pwd);
		String keyTable=cconn.getTableOnNode();
		String masterKeyTable=cconn.getKeyTable();
		cconn.close();
		String query="SELECT DISTINCT "+key+" FROM "+keyTable+" ORDER BY "+key+";";
		MySQLConnection mconn=null;
		try {				
			MySQLConnection nconn = this.getConnectionToNode(node);
			mconn = this.getConnectionToMasterKeyDb();
			Statement S=nconn.createStatement();
			ResultSet R=S.executeQuery(query);
			mconn.executeSql("CREATE TEMPORARY TABLE tmp_keys ("+key+" int(11) default NULL) ENGINE=MEMORY;");
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
							"(SELECT "+key+" AS id FROM tmp_keys UNION ALL SELECT kt."+key+" AS id FROM "+masterKeyTable+" AS kt LEFT JOIN clients_names AS cn ON kt.client_id=cn.client_id WHERE cn.client_name='"+node+"') AS u GROUP BY u.id) " +
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
	 * For a certain key and table finds out how the table is splitted in nodes and returns the "data distribution"
	 * Take care when using this method as an argument of insertIdsToKeyMaster: if table is not in chunks (but rather all data in all)
	 * then ids can't be inserted in key_master as there will be duplication, i.e. for key_id=1 as data is in all nodes there 
	 * will be 40 copies of it and thus 40 equal ids will try to be inserted in key_master, which a) makes no sense and 
	 * b) mysql won't do it anyway  
	 * @param key
	 * @param table 
	 * @return idSets HashMap, keys are node names, values: int array with the ids for each node
	 */
	public HashMap<String,int[]> getIdSetsFromNodes(String key, String table){
		HashMap<String,int[]> idSets =new HashMap<String,int[]>();
		String[] nodes = getNodes();
		for (String node:nodes){
			MySQLConnection conn = this.getConnectionToNode(node);
			idSets.put(node,conn.getAllIds4KeyAndTable(key,table));
			conn.close();		
		}
		return idSets;
	}
	
}

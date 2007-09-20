package proteinstructure;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collections;
import java.util.TreeMap;

import tools.MySQLConnection;

/**
 * A residue interaction graph derived from a single chain pdb protein structure loaded from a graph database in aglappe's format 
 * 
 * @author 		Jose Duarte
 * Class:		DbGraph
 * Package:		proteinstructure
 */
public class DbGraph extends Graph {

	private final static String MYSQLSERVER="white";
	private final static String MYSQLUSER=MySQLConnection.getUserName();
	private final static String MYSQLPWD="nieve";
	
	private final static String DEFAULT_CR ="(true)"; // default contact range (CR field in graph db)
	private final static String DEFAULT_CW ="1";      // default contact weight (CW field in graph db)
	
	private int graphid=0;
	//private int sm_id=0; // for future use
	
	private String dbname;
	private MySQLConnection conn;

	/**
	 * Constructs Graph object from graph db, given the dbname, pdbCode, pdbChainCode (classic pdb chain code), ct and cutoff
	 * and passing a MySQLConnection
	 * @param dbname
	 * @param conn
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param cutoff
	 * @param ct
	 * @throws GraphIdNotFoundError 
	 * @throws SQLException 
	 */
	public DbGraph(String dbname, MySQLConnection conn, String pdbCode, String pdbChainCode, double cutoff, String ct) throws GraphIdNotFoundError, SQLException {
		this.dbname=dbname;
		this.conn=conn;
		this.cutoff=cutoff;
		this.pdbCode=pdbCode.toLowerCase();				// our convention: pdb codes are lower case
		this.pdbChainCode=pdbChainCode.toUpperCase();	// our convention: chain codes are upper case
		this.ct=ct;
		this.directed=false;
		// we set the sequence to empty when we read from graph db. We don't have the full sequence in graph db
		// when we pass the sequence in getCM to the ContactMap constructor we want to have either a full sequence (with unobserveds) or a blank in case we don't have the info
		this.sequence=""; 
		//TODO graphs in db are never directed, so this doesn't really apply here. Must solve all this!
		if (ct.contains("/")){
			directed=true;
		}
		
		getgraphid(); // initialises graphid, sm_id and chainCode
		read_graph_from_db(); // gets contacts, nodes and sequence
		
		this.obsLength=nodes.size();
		if (!sequence.equals("")){
			this.fullLength=sequence.length();
		} else {
			// if nodes TreeMap has correct residue numbering then this should get the right full length,
			// we will only miss: gaps (unobserved residues) at the end of the sequence. Those we can't know unless full sequence is given
			this.fullLength=Collections.max(nodes.keySet());
		}
		this.numContacts=contacts.size();
		this.modified=false;
	}
	
	/**
	 * Constructs Graph object from graph db, given the dbname, pdbCode, pdbChainCode (classic pdb chain code), ct and cutoff
	 * MySQLConnection is taken from defaults in DbGraph class: MYSQLSERVER, MYSQLUSER, MYSQLPWD
	 * @param dbname
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param cutoff
	 * @param ct
	 * @throws GraphIdNotFoundError
	 * @throws SQLException
	 */
	public DbGraph(String dbname, String pdbCode, String pdbChainCode, double cutoff, String ct) throws GraphIdNotFoundError, SQLException{ 
		this(dbname,new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD),pdbCode,pdbChainCode,cutoff,ct);
	}
	
	
	/**
	 * Constructs Graph object from graph db, given the graphid and passing a MySQLConnection
	 * @param dbname
	 * @param conn
	 * @param graphid
	 * @throws GraphIdNotFoundError 
	 * @throws SQLException 
	 */
	public DbGraph(String dbname, MySQLConnection conn, int graphid) throws GraphIdNotFoundError, SQLException {
		this.dbname=dbname;
		this.conn=conn;
		this.graphid=graphid;
		this.directed=false;
		// we set the sequence to empty when we read from graph db. We don't have the full sequence in graph db
		// when we pass the sequence in getCM to the ContactMap constructor we want to have either a full sequence (with unobserveds) or a blank in case we don't have the info
		this.sequence="";
		
		read_graph_from_db(); // gets contacts, nodes and sequence
		get_db_graph_info(); // gets pdbCode, pdbChainCode, chainCode, ct and cutoff from db (from graph_id)
		
		//TODO graphs in db are never directed, so this doesn't really apply here. Must solve all this!
		if (ct.contains("/")){
			directed=true;
		}
		this.obsLength=nodes.size();
		if (!sequence.equals("")){
			this.fullLength=sequence.length();
		} else {
			// if nodes TreeMap has correct residue numbering then this should get the right full length,
			// we will only miss: gaps (unobserved residues) at the end of the sequence. Those we can't know unless full sequence is given
			this.fullLength=Collections.max(nodes.keySet());
		}
		this.numContacts=contacts.size();
		this.modified=false;
	}

	/**
	 * Constructs Graph object from graph db, given the graphid
	 * MySQLConnection is taken from defaults in DbGraph class: MYSQLSERVER, MYSQLUSER, MYSQLPWD
	 * @param dbname
	 * @param graphid
	 * @throws GraphIdNotFoundError
	 * @throws SQLException
	 */
	public DbGraph(String dbname, int graphid) throws GraphIdNotFoundError, SQLException{
		this(dbname,new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD), graphid);
	}
	
	/**
	 * Reads contacts and nodes from db.
	 * The db must be a graph db following our standard format, i.e. must have tables: 
	 * chain_graph, single_model_graph, single_model_node, single_model_edge
	 * We don't care here about the origin of the data (msdsd, pdbase, predicted) for the generation of the graph as long as it follows our data format
	 * We read both edges and nodes from single_model_edge and single_model_node.
	 * The sequence is set to blank, as we can't get the full sequence from graph db
	 * @param conn
	 * @throws SQLException 
	 */
	private void read_graph_from_db() throws SQLException{
		contacts = new EdgeSet();
		weights = new TreeMap<Edge, Double>();
		nodes = new TreeMap<Integer, String>();

		// we read only half of the matrix (contacts in one direction only) so that we have the same type of contacts as when creating Graph from Pdb object
		String sql="SELECT i_num,j_num,weight FROM "+dbname+".single_model_edge WHERE graph_id="+graphid+" AND j_num>i_num ORDER BY i_num,j_num ";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		while (rsst.next()) {
			int i=rsst.getInt(1);
			int j=rsst.getInt(2);
			double weight=rsst.getDouble(3);
			Edge cont = new Edge(i,j);
			contacts.add(cont);
			weights.put(cont,weight);
		}
		rsst.close();
		stmt.close();
		sql="SELECT num,res FROM "+dbname+".single_model_node WHERE graph_id="+graphid+" ORDER BY num ";
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		while (rsst.next()){
			int num=rsst.getInt(1);
			String res=rsst.getString(2);
			nodes.put(num, AAinfo.oneletter2threeletter(res));
		}
		rsst.close();
		stmt.close();
	}
	
	private void getgraphid () throws GraphIdNotFoundError, SQLException{
		// input is pdbChainCode i.e. pdb chain code
        // we take chainCode (internal chain identifier, pchain_code for msdsd and asym_id for pdbase) from pchain_code field in chain_graph 
        // (in the chain_graph table the internal chain identifier is called 'pchain_code')
		int pgraphid=0;
		String chainstr="='"+pdbChainCode+"' ";
		if (pdbChainCode.equals("NULL")){
			chainstr=" IS NULL ";
		}

		String sql="SELECT graph_id, pchain_code FROM "+dbname+".chain_graph WHERE accession_code='"+pdbCode+"' AND chain_pdb_code"+chainstr+" AND dist="+cutoff;
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		int check=0;
		while (rsst.next()) {
			check++;
			pgraphid=rsst.getInt(1);
			chainCode=rsst.getString(2);
		}
		if (check!=1){
			System.err.println("No pgraph_id match or more than 1 match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode+", dist="+cutoff);
		}
		rsst.close();
		stmt.close();
		// we set the ctstr to the same as ct except in ALL case, where it is BB+SC+BB/SC
		String ctstr=ct;
		if (ct.equals("ALL")){
			ctstr="BB+SC+BB/SC";
		}
		sql="SELECT graph_id,single_model_id FROM "+dbname+".single_model_graph WHERE pgraph_id="+pgraphid+" AND CT='"+ctstr+"' AND dist="+cutoff+" AND CR='"+DEFAULT_CR+"' AND CW="+DEFAULT_CW;
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		check=0;
		while (rsst.next()){
			check++;
			graphid=rsst.getInt(1);
			//sm_id=rsst.getInt(2); // we might want to use it in the future
		}
		if (check!=1){
			//System.err.println("No graph_id match or more than 1 match for pgraph_id="+pgraphid+", CT="+ctstr+" and cutoff="+cutoff);
			throw new GraphIdNotFoundError("No graph_id match or more than 1 match for pgraph_id="+pgraphid+", CT="+ctstr+" and cutoff="+cutoff);
		}
	}
	
	private void get_db_graph_info() throws GraphIdNotFoundError, SQLException {
			int pgraphid=0;
			String sql="SELECT pgraph_id,CT,dist FROM "+dbname+".single_model_graph WHERE graph_id="+graphid;
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			int check=0;
			while (rsst.next()) {
				check++;
				pgraphid=rsst.getInt(1);
				ct=rsst.getString(2);
				if (ct.equals("BB+SC+BB/SC")) ct="ALL";
				cutoff=rsst.getDouble(3);
			}
			if (check!=1){
				//System.err.println("No pgraph_id match or more than 1 match for graph_id="+graphid);
				throw new GraphIdNotFoundError("No pgraph_id match or more than 1 match for graph_id="+graphid+" in db"+conn.getDbname());
			}
			rsst.close();
			stmt.close();
			sql="SELECT accession_code, chain_pdb_code, pchain_code FROM "+dbname+".chain_graph WHERE graph_id="+pgraphid;
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			check=0;
			while (rsst.next()){
				check++;
				pdbCode=rsst.getString(1);
				pdbChainCode=rsst.getString(2);
				// java returns a null if the field is a database null, we want actually the "NULL" string in that case
				if (pdbChainCode==null) pdbChainCode="NULL";
				chainCode=rsst.getString(3);
			}
			if (check!=1){
				System.err.println("No accession_code+chain_pdb_code+pchain_code match or more than 1 match for graph_id="+pgraphid+" in chain_graph table");
			}
			rsst.close();
			stmt.close();
	}

}

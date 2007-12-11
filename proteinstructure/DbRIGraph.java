package proteinstructure;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collections;
import java.util.TreeMap;

import edu.uci.ics.jung.graph.util.EdgeType;

import tools.MySQLConnection;

/**
 * A residue interaction graph derived from a single chain pdb protein structure loaded from a graph database in aglappe's format 
 * 
 * @author 		Jose Duarte
 * Class:		DbGraph
 * Package:		proteinstructure
 */
public class DbRIGraph extends RIGraph {

	private static final long serialVersionUID = 1L;

	private final static String MYSQLSERVER="white";
	private final static String MYSQLUSER=MySQLConnection.getUserName();
	private final static String MYSQLPWD="nieve";
	
	private final static String DEFAULT_CR ="(true)"; 	// default contact range (CR field in graph db)
	private final static String DEFAULT_CW ="1";      	// default contact weight (CW field in graph db)
	
	private final static int DEFAULT_MODEL = 1;			// default model serial (NMR structures)

	private int graphid;

	private boolean directed;
	
	private String dbname;
	private MySQLConnection conn;

	/**
	 * Constructs RIGraph object from graph db, given the dbname, pdbCode, pdbChainCode (classic pdb chain code), ct and cutoff
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
	public DbRIGraph(String dbname, MySQLConnection conn, String pdbCode, String pdbChainCode, double distCutoff, String contactType, boolean directed, int model) throws GraphIdNotFoundError, SQLException {
		this.dbname=dbname;
		this.conn=conn;
		this.distCutoff=distCutoff;
		this.pdbCode=pdbCode.toLowerCase();				// our convention: pdb codes are lower case
		this.pdbChainCode=pdbChainCode.toUpperCase();	// our convention: chain codes are upper case
		this.contactType=contactType;
		this.model=model;
		// we set the sequence to empty when we read from graph db. We don't have the full sequence in graph db
		this.sequence=""; 
		this.directed = directed;

		getgraphid(); // initialises graphid, sm_id and chainCode
		read_graph_from_db(); // gets contacts and nodes and sets fullLength
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
	public DbRIGraph(String dbname, String pdbCode, String pdbChainCode, double cutoff, String ct, boolean directed, int model) throws GraphIdNotFoundError, SQLException{ 
		this(dbname,new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD),pdbCode,pdbChainCode,cutoff,ct,directed,model);
	}
	
	public DbRIGraph(String dbname, MySQLConnection conn, String pdbCode, String pdbChainCode, double cutoff, String ct, boolean directed) throws GraphIdNotFoundError, SQLException {
		this(dbname,conn,pdbCode,pdbChainCode,cutoff,ct,directed,DEFAULT_MODEL);
	}
	
	public DbRIGraph(String dbname, String pdbCode, String pdbChainCode, double cutoff, String ct, boolean directed) throws GraphIdNotFoundError, SQLException {
		this(dbname,new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD),pdbCode,pdbChainCode,cutoff,ct,directed,DEFAULT_MODEL);
	}
	
	/**
	 * Constructs Graph object from graph db, given the graphid and passing a MySQLConnection
	 * @param dbname
	 * @param conn
	 * @param graphid
	 * @throws GraphIdNotFoundError 
	 * @throws SQLException 
	 */
	public DbRIGraph(String dbname, MySQLConnection conn, int graphid) throws GraphIdNotFoundError, SQLException {
		this.dbname=dbname;
		this.conn=conn;
		this.graphid=graphid;
		// we set the sequence to empty when we read from graph db. We don't have the full sequence in graph db
		this.sequence="";
				
		get_db_graph_info(); // gets pdbCode, pdbChainCode, chainCode, ct, cutoff, directed from db (from graph_id)
		read_graph_from_db(); // gets contacts and nodes and sets fullLength
	}

	/**
	 * Constructs Graph object from graph db, given the graphid
	 * MySQLConnection is taken from defaults in DbGraph class: MYSQLSERVER, MYSQLUSER, MYSQLPWD
	 * @param dbname
	 * @param graphid
	 * @throws GraphIdNotFoundError
	 * @throws SQLException
	 */
	public DbRIGraph(String dbname, int graphid) throws GraphIdNotFoundError, SQLException{
		this(dbname,new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD), graphid);
	}
	
	/**
	 * Reads contacts and nodes from db.
	 * The db must be a graph db following our standard format, i.e. must have tables: 
	 * chain_graph, single_model_graph, single_model_node, single_model_edge
	 * We don't care here about the origin of the data (msdsd, pdbase, predictions) for 
	 * the generation of the graph as long as it follows our data format.
	 * The sequence is set to blank, as we can't get the full sequence from graph db.
	 * @param conn
	 * @throws SQLException 
	 */
	private void read_graph_from_db() throws SQLException{

		serials2nodes = new TreeMap<Integer,RIGNode>();

		// reading secondary structure
		secondaryStructure = new SecondaryStructure();
		String sql = "SELECT sstype, ssid, min(num), max(num) FROM "+dbname+".single_model_node " +
				" WHERE graph_id="+graphid+" AND sstype IS NOT NULL "+
				" GROUP BY ssid";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		while (rsst.next()){
			SecStrucElement sselem = new SecStrucElement(rsst.getString(1).charAt(0), rsst.getInt(3), rsst.getInt(4), rsst.getString(2));
			secondaryStructure.add(sselem);
		}
				
		// reading nodes
		sql="SELECT num,res FROM "+dbname+".single_model_node WHERE graph_id="+graphid;
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		int checkCount = 0;
		while (rsst.next()){
			checkCount++;
			int num=rsst.getInt(1);
			String res=rsst.getString(2);
			RIGNode node = new RIGNode(num,AAinfo.oneletter2threeletter(res),secondaryStructure.getSecStrucElement(num));
			serials2nodes.put(num,node);
			this.addVertex(node);
		}
				
		if (checkCount==0) { // no nodes: empty graph, we return
			this.fullLength = 0;
			return;
		}
		
		// reading edges
		// if undirected we have to prefilter and read only half of the matrix (contacts in one direction only) 
		EdgeType et = EdgeType.DIRECTED;
		String filterStr = "";
		if (!directed) {
			filterStr = " AND j_num>i_num ";
			et = EdgeType.UNDIRECTED;
		}
		sql="SELECT i_num,j_num,weight,distance FROM "+dbname+".single_model_edge WHERE graph_id="+graphid+" "+filterStr;
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		while (rsst.next()) {
			int i=rsst.getInt(1);
			int j=rsst.getInt(2);
			int atomWeight=rsst.getInt(3);
			double distance=rsst.getDouble(4);
			RIGEdge e = new RIGEdge(atomWeight);
			this.addEdge(e, serials2nodes.get(i),serials2nodes.get(j),et);
			e.setDistance(distance);
		}
		rsst.close();
		stmt.close();

		// if db has correct residue numbering then this should get the right full length,
		// we will only miss: gaps (unobserved residues) at the end of the sequence. Those we can't know unless full sequence is given
		this.fullLength=Collections.max(serials2nodes.keySet());

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
		
		String CW = DEFAULT_CW;
		String CR = DEFAULT_CR;
		String EXPBB = "0";
		String ctStr = contactType;
		String weightedStr = "0";
		String directedStr = directed?"1":"0";
		
		if (contactType.contains("_CAGLY")) {
			ctStr = contactType.replaceAll("_CAGLY", "");
		}
		// we set the ctstr to the same as ct except in ALL case, where it is BB+SC+BB/SC
		if (ctStr.equals("ALL")) {
			ctStr = "BB+SC+BB/SC";
		}
		if (AAinfo.isValidMultiAtomContactType(contactType, directed)) {
			CW = ctStr;
			weightedStr = "1";
		}
		if (contactType.contains("_CAGLY") || contactType.contains("Cb")) {
			EXPBB = "-1";
		}		
				
		String sql = "SELECT graph_id, pchain_code FROM "+dbname+".chain_graph " +
					" WHERE accession_code='"+pdbCode+"' AND chain_pdb_code "+chainstr+" " +
					" AND model_serial = "+model+" AND dist = "+distCutoff+" AND expBB = "+EXPBB+ 
					" AND method = 'rc-cutoff';";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		int check=0;
		while (rsst.next()) {
			check++;
			pgraphid=rsst.getInt(1);
			chainCode=rsst.getString(2);
		}
		if (check!=1){
			System.err.println("No pgraph_id match or more than 1 match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode+", dist="+distCutoff);
		}
		rsst.close();
		stmt.close();
		
		sql="SELECT graph_id, single_model_id FROM "+dbname+".single_model_graph "+
			" WHERE pgraph_id="+pgraphid+" AND dist="+distCutoff+" AND expBB="+EXPBB+
			" AND CW='"+CW+"' AND CT='"+ctStr+"' AND CR='"+CR+"' "+
			" AND w = "+weightedStr+" AND d = "+directedStr+";";
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		check=0;
		while (rsst.next()){
			check++;
			graphid=rsst.getInt(1);
			//sm_id=rsst.getInt(2); // we might want to use it in the future
		}
		//System.out.println(graphid);
		if (check!=1){
			//System.err.println("No graph_id match or more than 1 match for pgraph_id="+pgraphid+", CT="+ctstr+" and cutoff="+cutoff);
			throw new GraphIdNotFoundError("No graph_id match or more than 1 match for pgraph_id="+pgraphid+", CT="+ctStr+" and cutoff="+distCutoff);
		}
	}
	
	private void get_db_graph_info() throws GraphIdNotFoundError, SQLException {
			int pgraphid=0;
			String sql="SELECT pgraph_id,dist,expBB,CT,d FROM "+dbname+".single_model_graph WHERE graph_id="+graphid;
			String contactType="";
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			int check=0;
			while (rsst.next()) {
				check++;
				pgraphid=rsst.getInt(1);
				distCutoff=rsst.getDouble(2);
				int expBB=rsst.getInt(3);
				contactType=rsst.getString(4);
				if (contactType.equals("BB+SC+BB/SC")) contactType="ALL";
				if (expBB == -1) contactType.replaceAll("SC","SC_CAGLY");
				directed = (rsst.getInt(5)==1);
			}
			rsst.close();
			stmt.close();
			if (check!=1){
				//System.err.println("No pgraph_id match or more than 1 match for graph_id="+graphid);
				throw new GraphIdNotFoundError("No pgraph_id match or more than 1 match for graph_id="+graphid+" in db"+conn.getDbname());
			}
			
			sql="SELECT accession_code, chain_pdb_code, pchain_code, model_serial FROM "+dbname+".chain_graph WHERE graph_id="+pgraphid;
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
				model=rsst.getInt(4);
			}
			if (check!=1){
				System.err.println("No accession_code+chain_pdb_code+pchain_code match or more than 1 match for graph_id="+pgraphid+" in chain_graph table");
			}
			rsst.close();
			stmt.close();
	}

}

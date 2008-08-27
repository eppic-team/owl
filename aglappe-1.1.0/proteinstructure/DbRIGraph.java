package proteinstructure;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collections;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.uci.ics.jung.graph.util.EdgeType;

import tools.MySQLConnection;

/**
 * A residue interaction graph derived from a single chain pdb or scop protein structure loaded from a graph database in aglappe's format 
 * 
 * NOTE: If someone wants to load graphs with contact range != "(true)", 
 * 		then the constructor with the graph id should be used.
 * 		Another way is to load the whole graph if exists with any of the constructors
 * 		and then apply the corresponding methods to the graph object to get the desired contact range. 
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
	private boolean weighted;
	
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
	public DbRIGraph(String dbname, MySQLConnection conn, String pdbCode, String pdbChainCode, double distCutoff, String contactType, boolean directed, boolean weighted, int model) throws GraphIdNotFoundError, SQLException {
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
		this.weighted = weighted;

		getgraphid(pdbCode, pdbChainCode); // initialises graphid, sm_id and chainCode
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
	public DbRIGraph(String dbname, String pdbCode, String pdbChainCode, double cutoff, String ct, boolean directed, boolean weighted, int model) throws GraphIdNotFoundError, SQLException{ 
		this(dbname,new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD),pdbCode,pdbChainCode,cutoff,ct,directed,weighted,model);
	}
	
	public DbRIGraph(String dbname, MySQLConnection conn, String pdbCode, String pdbChainCode, double cutoff, String ct, boolean directed, boolean weighted) throws GraphIdNotFoundError, SQLException {
		this(dbname,conn,pdbCode,pdbChainCode,cutoff,ct,directed,weighted,DEFAULT_MODEL);
	}
	
	public DbRIGraph(String dbname, String pdbCode, String pdbChainCode, double cutoff, String ct, boolean directed, boolean weighted) throws GraphIdNotFoundError, SQLException {
		this(dbname,new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD),pdbCode,pdbChainCode,cutoff,ct,directed,weighted,DEFAULT_MODEL);
	}
	
	/**
	 * Constructs RIGraph object from graph db, given the dbname, sid, ct and cutoff
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
	public DbRIGraph(String dbname, MySQLConnection conn, String sid, double distCutoff, String contactType, boolean directed, boolean weighted, int model) throws GraphIdNotFoundError, SQLException {
		this.dbname=dbname;
		this.conn=conn;
		this.sid=sid;
		this.distCutoff=distCutoff;
		this.contactType=contactType;
		this.model=model;
		// we set the sequence to empty when we read from graph db. We don't have the full sequence in graph db
		this.sequence=""; 
		this.directed = directed;
		this.weighted = weighted;

		getgraphid(sid); // initialises graphid, sm_id, chainCode, pdbCode and pdbChainCode
		read_graph_from_db(); // gets contacts and nodes and sets fullLength
	}
	
	/**
	 * Constructs Graph object from graph db, given the dbname, sid, ct and cutoff
	 * MySQLConnection is taken from defaults in DbGraph class: MYSQLSERVER, MYSQLUSER, MYSQLPWD
	 * @param dbname
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param cutoff
	 * @param ct
	 * @throws GraphIdNotFoundError
	 * @throws SQLException
	 */
	public DbRIGraph(String dbname, String sid, double cutoff, String ct, boolean directed, boolean weighted, int model) throws GraphIdNotFoundError, SQLException{ 
		this(dbname,new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD),sid,cutoff,ct,directed,weighted, model);
	}
	
	public DbRIGraph(String dbname, MySQLConnection conn, String sid, double cutoff, String ct, boolean directed, boolean weighted) throws GraphIdNotFoundError, SQLException {
		this(dbname,conn,sid,cutoff,ct,directed,weighted,DEFAULT_MODEL);
	}
	
	public DbRIGraph(String dbname, String sid, double cutoff, String ct, boolean directed, boolean weighted) throws GraphIdNotFoundError, SQLException {
		this(dbname,new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD),sid,cutoff,ct,directed,weighted,DEFAULT_MODEL);
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
	 * @throws SQLException 
	 */
	private void read_graph_from_db() throws SQLException{

		// reading secondary structure
		secondaryStructure = new SecondaryStructure(this.sequence);
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
			this.addVertex(node); // this takes care of updating the serials2nodes map
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
			this.addEdge(e, getNodeFromSerial(i), getNodeFromSerial(j),et);
			e.setDistance(distance);
		}
		rsst.close();
		stmt.close();

		// if db has correct residue numbering then this should get the right full length,
		// we will only miss: gaps (unobserved residues) at the end of the sequence. Those we can't know unless full sequence is given
		if (sid == null) {
			this.fullLength=Collections.max(getSerials());
		} else {
			this.fullLength=getVertexCount();
		}

	}
	
	private void getgraphid (String pdbCode, String pdbChainCode) throws GraphIdNotFoundError, SQLException{
		// input is pdbChainCode i.e. pdb chain code
        // we take chainCode (internal chain identifier, pchain_code for msdsd and asym_id for pdbase) from pchain_code field in chain_graph 
        // (in the chain_graph table the internal chain identifier is called 'pchain_code')
		int pgraphid=0;
		String chainstr="='"+pdbChainCode+"' ";
		if (pdbChainCode.equals(Pdb.NULL_CHAIN_CODE)){
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
		if (AAinfo.isValidMultiAtomContactType(contactType, directed) && weighted) {
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
			" WHERE pgraph_id="+pgraphid+" AND graph_type='chain' AND dist="+distCutoff+" AND expBB="+EXPBB+
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
	
	private void getgraphid (String sid) throws GraphIdNotFoundError, SQLException{
		int pgraphid=0;
		
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
		if (AAinfo.isValidMultiAtomContactType(contactType, directed) && weighted) {
			CW = ctStr;
			weightedStr = "1";
		}
		if (contactType.contains("_CAGLY") || contactType.contains("Cb")) {
			EXPBB = "-1";
		}		
				
		String sql = "SELECT graph_id, pchain_code, accession_code, chain_pdb_code FROM "+dbname+".scop_graph " +
					" WHERE scop_id='"+sid+"' " +
					" AND model_serial = "+model+" AND dist = "+distCutoff+" AND expBB = "+EXPBB+ 
					" AND method = 'rc-cutoff';";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		int check=0;
		while (rsst.next()) {
			check++;
			pgraphid=rsst.getInt(1);
			this.chainCode=rsst.getString(2);
			this.pdbCode=rsst.getString(3);
			this.pdbChainCode=rsst.getString(4);
			// java returns a null if the field is a database null, we want actually the Pdb.NULL_CHAIN_CODE string in that case
			if (this.pdbChainCode==null) this.pdbChainCode=Pdb.NULL_CHAIN_CODE;
		}
		if (check!=1){
			System.err.println("No pgraph_id match or more than 1 match for scop id="+sid+", dist="+distCutoff);
		}
		rsst.close();
		stmt.close();
		
		sql="SELECT graph_id, single_model_id FROM "+dbname+".single_model_graph "+
			" WHERE pgraph_id="+pgraphid+" AND graph_type='scop' AND dist="+distCutoff+" AND expBB="+EXPBB+
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
			String sql="SELECT graph_type,pgraph_id,dist,expBB,CT,CR,d FROM "+dbname+".single_model_graph WHERE graph_id="+graphid;
			String contactRange="";
			String graphType="";
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			int check=0;
			Pattern p = Pattern.compile("^\\(\\(i_cid!=j_cid\\)OR\\(abs\\(i_num-j_num\\)>=(\\d+)\\)\\)$");
			while (rsst.next()) {
				check++;
				graphType=rsst.getString(1);
				pgraphid=rsst.getInt(2);
				distCutoff=rsst.getDouble(3);
				int expBB=rsst.getInt(4);
				contactType=rsst.getString(5);
				if (contactType.equals("BB+SC+BB/SC")) contactType="ALL";
				if (expBB == -1) {
					contactType = contactType.replaceAll("SC","SC_CAGLY");
				}
				contactRange=rsst.getString(6);
				Matcher m = p.matcher(contactRange);
				if (contactRange.equals("((i_sstype!=j_sstype)OR(i_ssid!=j_ssid))")) {
					interSSE = true;
				} else if (m.matches()) {
					minSeqSep = Integer.valueOf(m.group(1));
				}
				directed = (rsst.getInt(7)==1);
			}
			rsst.close();
			stmt.close();
			if (check!=1){
				//System.err.println("No pgraph_id match or more than 1 match for graph_id="+graphid);
				throw new GraphIdNotFoundError("No pgraph_id match or more than 1 match for graph_id="+graphid+" in db "+conn.getDbname());
			}
			
			sql="SELECT accession_code, chain_pdb_code, pchain_code, model_serial,"+(graphType.equals("scop")?"scop_id":"NULL")+" FROM "+dbname+"."+graphType+"_graph WHERE graph_id="+pgraphid;
			stmt = conn.createStatement();
			rsst = stmt.executeQuery(sql);
			check=0;
			while (rsst.next()){
				check++;
				pdbCode=rsst.getString(1);
				pdbChainCode=rsst.getString(2);
				// java returns a null if the field is a database null, we want actually the Pdb.NULL_CHAIN_CODE string in that case
				if (pdbChainCode==null) pdbChainCode=Pdb.NULL_CHAIN_CODE;
				chainCode=rsst.getString(3);
				model=rsst.getInt(4);
				sid=rsst.getString(5);
			}
			if (check!=1){
				System.err.println("No accession_code+chain_pdb_code+pchain_code match or more than 1 match for graph_id="+pgraphid+" in chain_graph table");
			}
			rsst.close();
			stmt.close();
	}

}

package proteinstructure;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.IOException;
import java.util.Collections;
import java.util.Iterator;
import java.util.Locale;
import java.util.Set;
import java.util.TreeMap;
import java.util.HashMap;

import java.sql.SQLException;
import java.sql.Statement;

import java.sql.ResultSet;

import sadp.ContactMap;
import tools.MySQLConnection;

/**
 * A residue interaction graph derived from a single chain pdb protein structure
 * 
 */
public class Graph {

	public final static String GRAPHFILEFORMATVERSION = "1.0";
	
	private final static String SINGLEMODELS_DB = "ioannis";	//TODO: Is this being used??
	
	public EdgeSet contacts; // we keep it public to be able to re-reference the object directly (getContacts() copies it)
	
	protected TreeMap<Integer,String> nodes; // nodes is a TreeMap of residue serials to residue types (3 letter code)
	protected SecondaryStructure secondaryStructure; // secondary structure annotation for this protein graph

	protected String sequence; 				// the full sequence (with unobserved residues and non-standard aas ='X')
	protected String pdbCode;
	protected String chainCode;
	protected String pdbChainCode;
	protected int model;
	protected double cutoff;
	protected String ct;					// the contact type
	protected boolean directed;
	protected int minSeqSep = -1;
	protected int maxSeqSep = -1;
	
	// fullLength is length of full sequence or:
	// -if sequence not provided (when reading from db): length of everything except possible unobserved residues at end of chain
	// -if sequence and nodes not provided (when reading from file and sequence field missing): length except possible unobserved residues at end of chain and possible nodes without contacts at end of chain
	protected int fullLength; 
	protected int obsLength;  // length without unobserved, non standard aas 
	
	protected int numContacts;
	
	protected boolean modified;

	/**
	 * Constructs an empty graph object
	 */
	public Graph() {
		this.contacts=new EdgeSet();
		this.cutoff=0;
		this.nodes=new TreeMap<Integer,String>();
		this.sequence=null;
		this.pdbCode=null;
		this.chainCode=null;
		this.pdbChainCode=null;
		this.model=1;
		this.ct=null;
		this.fullLength=0;
		this.obsLength=0;
		this.numContacts=0;
		this.modified=false;
		this.directed=false;
		this.secondaryStructure = new SecondaryStructure();
	}
	
	/**
	 * Constructs a graph with a sequence but no edges
	 */
	public Graph(String sequence) {
		this.contacts=new EdgeSet();
		this.cutoff=0;
		this.nodes=new TreeMap<Integer,String>();
		for(int i=0; i < sequence.length(); i++) {
			nodes.put(i+1, AAinfo.oneletter2threeletter(Character.toString(sequence.charAt(i))));
		}
		this.sequence=sequence;
		this.pdbCode=null;
		this.chainCode=null;
		this.pdbChainCode=null;
		this.model=1;
		this.ct=null;
		this.fullLength=sequence.length();
		this.obsLength=fullLength;
		this.numContacts=0;
		this.modified=false;
		this.directed=false;
		this.secondaryStructure = new SecondaryStructure();
	}
	
	/**
	 * Constructs Graph object by passing EdgeSet with contacts and TreeMap with nodes (res serials and types)
	 * Must also pass contact type, cutoff, pdbCode and chainCode
	 * @param contacts
	 * @param nodes
	 * @param sequence
	 * @param cutoff
	 * @param ct
	 * @param pdbCode
	 * @param chainCode
	 * @param pdbChainCode
	 * @param model
	 * @param secStruct
	 */
	public Graph (EdgeSet contacts, TreeMap<Integer,String> nodes, String sequence, double cutoff,String ct, String pdbCode, String chainCode, String pdbChainCode, int model, SecondaryStructure secStruct) {
		this.contacts=contacts;
		this.cutoff=cutoff;
		this.nodes=nodes;
		this.sequence=sequence;
		this.pdbCode=pdbCode;
		this.chainCode=chainCode;
		this.pdbChainCode=pdbChainCode;
		this.model=model;
		this.ct=ct;
		// in case of pdb was read from file and there was no SEQRES field then fullLength here shouldn't be sequence length but maximum observed residue (see Pdb class)
		this.fullLength=Math.max(sequence.length(),Collections.max(nodes.keySet()));
		this.obsLength=nodes.size();
		this.numContacts=contacts.size();
		this.modified=false;
		this.directed=false;
		if (ct.contains("/")){
			directed=true;
		}
		if(secStruct == null) {
			// we allow null to be passed to simplify graph construction
			this.secondaryStructure = new SecondaryStructure();
		} else {
			this.secondaryStructure = secStruct;
		}
		
		// do some verification checks
		assert(secondaryStructure != null);
		assert(this.pdbCode.equals(this.pdbCode.toLowerCase()));				// pdb codes should be always lower case 
	}
	
	/**
	 * Constructs a Graph from a sadp.ContactMap.
	 * TODO: This constructor might produce erroneous Graph instances.
	 * 
	 *  @param cm            the ContactMap
	 *  @param sequence
	 *  @param cutoff
	 *  @param ct
	 *  @param pdbCode
	 *  @param chainCode
	 *  @param pdbChainCode
	 *  @param model
	 *  @param secStruct
	 */
	public Graph( ContactMap cm, String sequence, double cutoff, String ct, String pdbCode, String chainCode, String pdbChainCode, int model, SecondaryStructure secStruct) {
	    
	    EdgeSet contacts = new EdgeSet();
	    int[][] adjList  = cm.getAdjacencyList();
	    
	    // construct the contacts
	    for( int i=0; i<adjList.length; ++i ) {
	    	for( int j=0; j<adjList[i].length; ++j ) {
			    // store contacts only once for i < j
			    if( i<adjList[i][j] ) {
					contacts.add(new Edge(i+1,adjList[i][j]+1));
			    }
			}
	    }
	    
	    TreeMap<Integer, String> nodes = new TreeMap<Integer, String>();
	    
	    // construct the mapping from node index to 3-lettercode
	    if( sequence != null ) {
			for( int i=0; i<sequence.length(); ++i ) {
				nodes.put(i+1,AAinfo.oneletter2threeletter(sequence.substring(i,i)));
			}
		} else {
			// as there is no sequence information available we assign the
			// non-defined unknown amino acid 'Xaa' to every residue 
			for( int i=0; i<cm.countNodes(); ++i ) {
				nodes.put(i+1, "Xaa");
			}		
	    }
	    
	    // construct dummy sequence if necessary
	    if( sequence  == null ) {
			StringBuffer sb = new StringBuffer(cm.countNodes());
			for( int i=0; i<sb.capacity(); ++i ) {
				sb.append('X');
			}
			sequence = sb.toString();
	    }
	    
	    this.contacts=contacts;
	    
	    this.cutoff=cutoff;
	    this.nodes=nodes;
	    this.sequence=sequence;
	    this.pdbCode      = pdbCode == null ? "" : pdbCode;
	    this.chainCode    = chainCode == null ? "" : chainCode;
	    this.pdbChainCode = pdbChainCode == null ? "" : pdbChainCode;
	    this.model=model;
	    this.ct=ct;
	    // in case of pdb was read from file and there was no SEQRES field then fullLength here shouldn't be sequence length but maximum observed residue (see Pdb class)
	    this.fullLength=Math.max(sequence.length(),Collections.max(nodes.keySet()));
	    this.obsLength=nodes.size();
	    this.numContacts=contacts.size();
	    this.modified=false;
	    this.directed=false;
	    if (ct.contains("/")){
		directed=true;
	    }
	    if(secStruct == null) {
		// we allow null to be passed to simplify graph construction
		this.secondaryStructure = new SecondaryStructure();
	    } else {
		this.secondaryStructure = secStruct;
	    }

	    // do some verification checks
	    assert(secondaryStructure != null);
	    assert(this.pdbCode.equals(this.pdbCode.toLowerCase()));				// pdb codes should be always lower case 
	    assert(this.pdbChainCode.equals(this.pdbChainCode.toUpperCase()));		// pdb chain codes should be always upper case	    
	}
	
	/**
	 * Write graph to given db, using our db graph aglappe format, 
	 * i.e. tables: chain_graph, single_model_graph, single_model_node, single_model_edge
	 * @param conn
	 * @param db
	 * @throws SQLException
	 */
	public void write_graph_to_db(MySQLConnection conn, String db) throws SQLException{
		
		conn.setSqlMode("NO_UNSIGNED_SUBTRACTION,TRADITIONAL");

		// we are fixing these 3 values to what corresponds to our graphs 
		String CW = "1";
		String CR = "(true)";
		String EXPBB = "0";
		String ctStr = ct;
		String weightedStr = "0";
		String directedStr = directed?"1":"0";
		
		if (ct.endsWith("_CAGLY")) {
			ctStr = ct.replace("_CAGLY", "");
		}
		if (ctStr.equals("ALL")) {
			ctStr = "BB+SC+BB/SC";
		}
		if (AAinfo.isValidMultiAtomContactType(ct)) {
			CW = ctStr;
			weightedStr = "1";
		}
		if (ct.endsWith("_CAGLY") || ct.equals("Cb")) {
			EXPBB = "-1";
		}
		if (minSeqSep != -1) {
			CR = "((i_cid!=j_cid)OR(abs(i_num-j_num)>="+minSeqSep+"))";
		}
				
		int pgraphid=0;
		int graphid=0;
		String sql = "SELECT graph_id FROM "+db+".chain_graph " +
					" WHERE accession_code='"+pdbCode+"' AND pchain_code='"+chainCode+"'" +
					" AND model_serial = "+model+" AND dist = "+cutoff+" AND expBB = "+EXPBB+ 
					" AND method = 'rc-cutoff';";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		if (rsst.next()){	// if the pdbCode + chainCode were already in chain_graph then we take the graph_id as the pgraphid
			pgraphid = rsst.getInt(1);
		} else {			// no pdbCode + chainCode found, we insert them in chain_graph, thus assigning a new graph_id (pgraphid)
			// we are inserting same number for num_obs_res and num_nodes (the difference would be the non-standard aas, but we can't get that number from this object at the moment)
			String pdbChainCodeStr = pdbChainCode;
			if (!pdbChainCode.equals("NULL")) {
				pdbChainCodeStr="'"+pdbChainCode+"'";
			}
			sql = "INSERT INTO "+db+".chain_graph (accession_code,chain_pdb_code,pchain_code,model_serial,dist,expBB,method,num_res,num_obs_res,num_nodes,sses,date) " +
					"VALUES ('"+pdbCode+"', "+pdbChainCodeStr+",'"+chainCode+"', "+model+", "+cutoff+", "+EXPBB+", 'rc-cutoff', "+getFullLength()+", "+getObsLength()+", "+getObsLength()+", "+((secondaryStructure!=null)?secondaryStructure.getNumElements():0)+", now())";
			Statement stmt2 = conn.createStatement();
			stmt2.executeUpdate(sql);
			// now we take the newly assigned graph_id as pgraphid
			sql = "SELECT LAST_INSERT_ID() FROM "+db+".chain_graph LIMIT 1";
			ResultSet rsst2 = stmt2.executeQuery(sql);
			if (rsst2.next()){
				pgraphid = rsst2.getInt(1);
			}
			stmt2.close();
			rsst2.close();
		}
		rsst.close();
		// now we insert the graph info into single_model_graph
		// 1st we grab the single_model_id
		int singlemodelid = 0;
		sql = "SELECT single_model_id FROM "+SINGLEMODELS_DB+".single_model WHERE "+
				" dist="+cutoff+" AND expBB="+EXPBB+" AND CW='"+CW+"' AND CT='"+ctStr+"' AND CR='"+CR+"';";
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			singlemodelid = rsst.getInt(1);
		}
		rsst.close();
		// and then insert to single_model_graph
		sql = "INSERT INTO "+db+".single_model_graph (pgraph_id,graph_type,accession_code,single_model_id,dist,expBB,CW,CT,CR,w,d,num_nodes,date) " +
				" VALUES ("+pgraphid+", 'chain', '"+pdbCode+"', "+singlemodelid+", "+cutoff+", "+EXPBB+", '"+CW+"','"+ctStr+"', '"+CR+"', "+weightedStr+", "+directedStr+", "+getObsLength()+", now())";
		stmt.executeUpdate(sql);
		// and we grab the graph_id just assigned in single_model_graph
		sql = "SELECT LAST_INSERT_ID() FROM "+db+".single_model_graph LIMIT 1";
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			graphid = rsst.getInt(1);
		}
		rsst.close();
		stmt.close();
		
		// inserting nodes
		// get the max node in db
		int maxNodeId = 0;
		sql = "SELECT MAX(node_id) FROM "+db+".single_model_node;";
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			maxNodeId = rsst.getInt(1);
		}
		rsst.close();
		stmt.close();
		
		stmt = conn.createStatement();
		for (int resser:nodes.keySet()) {
			String res = AAinfo.threeletter2oneletter(getResType(resser));
			NodeNbh nbh = getNodeNbh(resser);
			String secStructType = null;
			String secStructId = null;
			String sheetSerial = null;
			String turn = null;
			if (secondaryStructure.getSecStrucElement(resser)!=null){
				secStructType = quote(Character.toString(secondaryStructure.getSecStrucElement(resser).getType()));
				secStructId = quote(secondaryStructure.getSecStrucElement(resser).getId());
				char sheetSerialChar = secondaryStructure.getSecStrucElement(resser).getSheetSerial();
				if (sheetSerialChar != 0) {
					sheetSerial = quote(Character.toString(sheetSerialChar));
				}
				turn = secondaryStructure.getSecStrucElement(resser).isTurn()?"1":"0";
			}
			if (directed){  // we insert k(=k_in+k_out), k_in and k_out
				sql = "INSERT INTO "+db+".single_model_node "+
					" (graph_id, node_id, cid, num, res, "+
					" sstype, ssid, sheet_serial, turn, "+
					" k, k_in, k_out, "+
					" n, nwg, n_num) " +
					" VALUES ("+graphid+", "+(maxNodeId+resser)+", '"+chainCode+"', "+resser+", '"+res+"', "+
					" "+secStructType+", "+secStructId+", "+sheetSerial+", "+turn+", "+
					(getInDegree(resser)+getOutDegree(resser))+", "+getInDegree(resser)+", "+getOutDegree(resser)+", "+
					"'"+nbh.getMotifNoGaps()+"', '"+nbh.getMotif()+"', '"+nbh.getCommaSeparatedResSerials()+"')";
			} else {		// we insert k (and no k_in or k_out)
				sql = "INSERT INTO "+db+".single_model_node "+
				" (graph_id, node_id, cid, num, res, "+
				" sstype, ssid, sheet_serial, turn, "+
				" k, n, nwg, n_num) " +
				" VALUES ("+graphid+", "+(maxNodeId+resser)+", '"+chainCode+"', "+resser+", '"+res+"', "+
				" "+secStructType+", "+secStructId+", "+sheetSerial+", "+turn+", "+
				getDegree(resser)+", '"+nbh.getMotifNoGaps()+"', '"+nbh.getMotif()+"', '"+nbh.getCommaSeparatedResSerials()+"')";
			}
			stmt.executeUpdate(sql);
		}
		
		// inserting edges
		// get the max weight
		double maxWeight = 0;
		for (Edge cont:contacts) {
			maxWeight = (maxWeight<cont.weight)?cont.weight:maxWeight;
		}
		for (Edge cont:contacts){
			String i_res = AAinfo.threeletter2oneletter(getResType(cont.i));
			String j_res = AAinfo.threeletter2oneletter(getResType(cont.j));

			String i_secStructType = null;
			String i_secStructId = null;
			String i_sheetSerial = null;
			String i_turn = null;
			if (secondaryStructure.getSecStrucElement(cont.i)!=null){
				i_secStructType = quote(Character.toString(secondaryStructure.getSecStrucElement(cont.i).getType()));
				i_secStructId = quote(secondaryStructure.getSecStrucElement(cont.i).getId());
				char sheetSerialChar = secondaryStructure.getSecStrucElement(cont.i).getSheetSerial();
				if (sheetSerialChar != 0) {
					i_sheetSerial = quote(Character.toString(sheetSerialChar));
				}
				i_turn = secondaryStructure.getSecStrucElement(cont.i).isTurn()?"1":"0";
			}
			
			String j_secStructType = null;
			String j_secStructId = null;
			String j_sheetSerial = null;
			String j_turn = null;
			if (secondaryStructure.getSecStrucElement(cont.j)!=null){
				j_secStructType = quote(Character.toString(secondaryStructure.getSecStrucElement(cont.j).getType()));
				j_secStructId = quote(secondaryStructure.getSecStrucElement(cont.j).getId());
				char sheetSerialChar = secondaryStructure.getSecStrucElement(cont.j).getSheetSerial();
				if (sheetSerialChar != 0) {
					j_sheetSerial = quote(Character.toString(sheetSerialChar));
				}
				j_turn = secondaryStructure.getSecStrucElement(cont.j).isTurn()?"1":"0";
			}
			
			sql = "INSERT INTO "+db+".single_model_edge "+
					" (graph_id, i_node_id, i_cid, i_num, i_res, i_sstype, i_ssid, i_sheet_serial, i_turn, "+
					" j_node_id, j_cid, j_num, j_res, j_sstype, j_ssid, j_sheet_serial, j_turn, weight, norm_weight) " +
					" VALUES ("+graphid+", "+(maxNodeId+cont.i)+", '"+chainCode+"', "+cont.i+", '"+i_res+"', "+i_secStructType+", "+i_secStructId+", "+i_sheetSerial+", "+i_turn+", "+
					(maxNodeId+cont.j)+", '"+chainCode+"', "+cont.j+", '"+j_res+"', "+j_secStructType+", "+j_secStructId+", "+j_sheetSerial+", "+j_turn+", "+
					Math.round(cont.weight)+", "+(cont.weight/maxWeight)+")";
			stmt.executeUpdate(sql);
			if(!directed) {// we want both side of the matrix in the table to follow Ioannis' convention
				// so we insert the reverse contact by swapping i, j in insertion
				sql = "INSERT INTO "+db+".single_model_edge "+
				" (graph_id, i_node_id, i_cid, i_num, i_res, i_sstype, i_ssid, i_sheet_serial, i_turn, "+
				" j_node_id, j_cid, j_num, j_res, j_sstype, j_ssid, j_sheet_serial, j_turn, weight, norm_weight) " +
				" VALUES ("+graphid+", "+(maxNodeId+cont.j)+", '"+chainCode+"', "+cont.j+", '"+j_res+"', "+j_secStructType+", "+j_secStructId+", "+j_sheetSerial+", "+j_turn+", "+
				(maxNodeId+cont.i)+", '"+chainCode+"', "+cont.i+", '"+i_res+"', "+i_secStructType+", "+i_secStructId+", "+i_sheetSerial+", "+i_turn+", "+
				Math.round(cont.weight)+", "+(cont.weight/maxWeight)+")";
				stmt.executeUpdate(sql);
			}
		}
		
		stmt.close();
	}
		
	/**
	 * Write graph to given db, using our db graph aglappe format, 
	 * i.e. tables: chain_graph, single_model_graph, single_model_node, single_model_edge
	 * @param conn
	 * @param db
	 * @throws SQLException
	 */
	public void write_graph_to_db_fast(MySQLConnection conn, String db) throws SQLException, IOException {
		
		conn.setSqlMode("NO_UNSIGNED_SUBTRACTION,TRADITIONAL");

		// we are fixing these 3 values to what corresponds to our graphs 
		String CW = "1";
		String CR = "(true)";
		String EXPBB = "0";
		String ctStr = ct;
		String weightedStr = "0";
		String directedStr = directed?"1":"0";
		
		if (ct.endsWith("_CAGLY")) {
			ctStr = ct.replace("_CAGLY", "");
		}
		if (ctStr.equals("ALL")) {
			ctStr = "BB+SC+BB/SC";
		}
		if (AAinfo.isValidMultiAtomContactType(ct)) {
			CW = ctStr;
			weightedStr = "1";
		}
		if (ct.endsWith("_CAGLY") || ct.equals("Cb")) {
			EXPBB = "-1";
		}
		if (minSeqSep != -1) {
			CR = "((i_cid!=j_cid)OR(abs(i_num-j_num)>="+minSeqSep+"))";
		}
				
		int pgraphid=0;
		int graphid=0;
		String sql = "SELECT graph_id FROM "+db+".chain_graph " +
					" WHERE accession_code='"+pdbCode+"' AND pchain_code='"+chainCode+"'" +
					" AND model_serial = "+model+" AND dist = "+cutoff+" AND expBB = "+EXPBB+
					" AND method = 'rc-cutoff';";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		if (rsst.next()){	// if the pdbCode + chainCode were already in chain_graph then we take the graph_id as the pgraphid
			pgraphid = rsst.getInt(1);
		} else {			// no pdbCode + chainCode found, we insert them in chain_graph, thus assigning a new graph_id (pgraphid)
			// we are inserting same number for num_obs_res and num_nodes (the difference would be the non-standard aas, but we can't get that number from this object at the moment)
			String pdbChainCodeStr = pdbChainCode;
			if (!pdbChainCode.equals("NULL")) {
				pdbChainCodeStr="'"+pdbChainCode+"'";
			}
			sql = "INSERT INTO "+db+".chain_graph (accession_code,chain_pdb_code,pchain_code,model_serial,dist,expBB,method,num_res,num_obs_res,num_nodes,sses,date) " +
					"VALUES ('"+pdbCode+"', "+pdbChainCodeStr+",'"+chainCode+"', "+model+", "+cutoff+", "+EXPBB+", 'rc-cutoff', "+getFullLength()+", "+getObsLength()+", "+getObsLength()+", "+((secondaryStructure!=null)?secondaryStructure.getNumElements():0)+", now())";
			Statement stmt2 = conn.createStatement();
			stmt2.executeUpdate(sql);
			// now we take the newly assigned graph_id as pgraphid
			sql = "SELECT LAST_INSERT_ID() FROM "+db+".chain_graph LIMIT 1";
			ResultSet rsst2 = stmt2.executeQuery(sql);
			if (rsst2.next()){
				pgraphid = rsst2.getInt(1);
			}
			stmt2.close();
			rsst2.close();
		}
		rsst.close();
		// now we insert the graph info into single_model_graph
		// 1st we grab the single_model_id
		int singlemodelid = 0;
		sql = "SELECT single_model_id FROM "+SINGLEMODELS_DB+".single_model WHERE "+
				" dist="+cutoff+" AND expBB="+EXPBB+" AND CW='"+CW+"' AND CT='"+ctStr+"' AND CR='"+CR+"';";
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			singlemodelid = rsst.getInt(1);
		}
		rsst.close();
		// and then insert to single_model_graph
		sql = "INSERT INTO "+db+".single_model_graph (pgraph_id,graph_type,accession_code,single_model_id,dist,expBB,CW,CT,CR,w,d,num_nodes,date) " +
				" VALUES ("+pgraphid+", 'chain', '"+pdbCode+"', "+singlemodelid+", "+cutoff+", "+EXPBB+", '"+CW+"','"+ctStr+"', '"+CR+"', "+weightedStr+", "+directedStr+", "+getObsLength()+", now())";
		stmt.executeUpdate(sql);
		// and we grab the graph_id just assigned in single_model_graph
		sql = "SELECT LAST_INSERT_ID() FROM "+db+".single_model_graph LIMIT 1";
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			graphid = rsst.getInt(1);
		}
		rsst.close();
		stmt.close();
		
		// inserting nodes
		PrintStream nodesOut = new PrintStream(new FileOutputStream(graphid+"_nodes.txt"));
		// get the max node in db
		int maxNodeId = 0;
		sql = "SELECT MAX(node_id) FROM "+db+".single_model_node;";
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			maxNodeId = rsst.getInt(1);
		}
		rsst.close();
		stmt.close();
		
		for (int resser:nodes.keySet()) {
			String res = AAinfo.threeletter2oneletter(getResType(resser));
			NodeNbh nbh = getNodeNbh(resser);
			String secStructType = "\\N";
			String secStructId = "\\N";
			String sheetSerial = "\\N";
			String turn = null;
			if (secondaryStructure.getSecStrucElement(resser)!=null){
				secStructType = Character.toString(secondaryStructure.getSecStrucElement(resser).getType());
				secStructId = secondaryStructure.getSecStrucElement(resser).getId();
				char sheetSerialChar = secondaryStructure.getSecStrucElement(resser).getSheetSerial();
				if (sheetSerialChar != 0) {
					sheetSerial = Character.toString(sheetSerialChar);
				}
				turn = secondaryStructure.getSecStrucElement(resser).isTurn()?"1":"0";
			}
			if (directed){  // we insert k(=k_in+k_out), k_in and k_out
				nodesOut.println(graphid+"\t"+(maxNodeId+resser)+"\t"+chainCode+"\t"+resser+"\t"+res+"\t"+
					secStructType+"\t"+secStructId+"\t"+sheetSerial+"\t"+turn+"\t"+
					(getInDegree(resser)+getOutDegree(resser))+"\t"+getInDegree(resser)+"\t"+getOutDegree(resser)+"\t"+
					nbh.getMotifNoGaps()+"\t"+nbh.getMotif()+"\t"+nbh.getCommaSeparatedResSerials());
			} else {		// we insert k (and no k_in or k_out)
				nodesOut.println(graphid+"\t"+(maxNodeId+resser)+"\t"+chainCode+"\t"+resser+"\t"+res+"\t"+
						secStructType+"\t"+secStructId+"\t"+sheetSerial+"\t"+turn+"\t"+
						getDegree(resser)+"\t"+"\\N"+"\t"+"\\N"+"\t"+
						nbh.getMotifNoGaps()+"\t"+nbh.getMotif()+"\t"+nbh.getCommaSeparatedResSerials());
			}
		}
		nodesOut.close();
		stmt = conn.createStatement();
		sql = "LOAD DATA LOCAL INFILE '"+graphid+"_nodes.txt' INTO TABLE "+db+".single_model_node "+
			" (graph_id, node_id, cid, num, res, "+
			" sstype, ssid, sheet_serial, turn, "+
			" k, k_in, k_out, n, nwg, n_num);";
		stmt.executeUpdate(sql);
		File fileToDelete = new File(graphid+"_nodes.txt");
		if (fileToDelete.exists()) {
			fileToDelete.delete();
		}
		
		// inserting edges
		PrintStream edgesOut = new PrintStream(new FileOutputStream(graphid+"_edges.txt"));
		// get the max weight
		double maxWeight = 0;
		for (Edge cont:contacts) {
			maxWeight = (maxWeight<cont.weight)?cont.weight:maxWeight;
		}
		for (Edge cont:contacts){
			String i_res = AAinfo.threeletter2oneletter(getResType(cont.i));
			String j_res = AAinfo.threeletter2oneletter(getResType(cont.j));

			String i_secStructType = "\\N";
			String i_secStructId = "\\N";
			String i_sheetSerial = "\\N";
			String i_turn = null;
			if (secondaryStructure.getSecStrucElement(cont.i)!=null){
				i_secStructType = Character.toString(secondaryStructure.getSecStrucElement(cont.i).getType());
				i_secStructId = secondaryStructure.getSecStrucElement(cont.i).getId();
				char sheetSerialChar = secondaryStructure.getSecStrucElement(cont.i).getSheetSerial();
				if (sheetSerialChar != 0) {
					i_sheetSerial = Character.toString(sheetSerialChar);
				}
				i_turn = secondaryStructure.getSecStrucElement(cont.i).isTurn()?"1":"0";
			}
			
			String j_secStructType = "\\N";
			String j_secStructId = "\\N";
			String j_sheetSerial = "\\N";
			String j_turn = null;
			if (secondaryStructure.getSecStrucElement(cont.j)!=null){
				j_secStructType = Character.toString(secondaryStructure.getSecStrucElement(cont.j).getType());
				j_secStructId = secondaryStructure.getSecStrucElement(cont.j).getId();
				char sheetSerialChar = secondaryStructure.getSecStrucElement(cont.j).getSheetSerial();
				if (sheetSerialChar != 0) {
					j_sheetSerial = Character.toString(sheetSerialChar);
				}
				j_turn = secondaryStructure.getSecStrucElement(cont.j).isTurn()?"1":"0";
			}
			
			edgesOut.println(graphid+"\t"+(maxNodeId+cont.i)+"\t"+chainCode+"\t"+cont.i+"\t"+i_res+"\t"+i_secStructType+"\t"+i_secStructId+"\t"+i_sheetSerial+"\t"+i_turn+"\t"+
					(maxNodeId+cont.j)+"\t"+chainCode+"\t"+cont.j+"\t"+j_res+"\t"+j_secStructType+"\t"+j_secStructId+"\t"+j_sheetSerial+"\t"+j_turn+"\t"+
					Math.round(cont.weight)+"\t"+(cont.weight/maxWeight));
			if(!directed) {// we want both side of the matrix in the table to follow Ioannis' convention
				// so we insert the reverse contact by swapping i, j in insertion
				edgesOut.println(graphid+"\t"+(maxNodeId+cont.j)+"\t"+chainCode+"\t"+cont.j+"\t"+j_res+"\t"+j_secStructType+"\t"+j_secStructId+"\t"+j_sheetSerial+"\t"+j_turn+"\t"+
						(maxNodeId+cont.i)+"\t"+chainCode+"\t"+cont.i+"\t"+i_res+"\t"+i_secStructType+"\t"+i_secStructId+"\t"+i_sheetSerial+"\t"+i_turn+"\t"+
						Math.round(cont.weight)+"\t"+(cont.weight/maxWeight));
			}			
		}
		edgesOut.close();
		sql = "LOAD DATA LOCAL INFILE '"+graphid+"_edges.txt' INTO TABLE "+db+".single_model_edge "+
			" (graph_id, i_node_id, i_cid, i_num, i_res, i_sstype, i_ssid, i_sheet_serial, i_turn, "+
			" j_node_id, j_cid, j_num, j_res, j_sstype, j_ssid, j_sheet_serial, j_turn, weight, norm_weight);";
		stmt.executeUpdate(sql);
		stmt.close();
		fileToDelete = new File(graphid+"_edges.txt");
		if (fileToDelete.exists()) {
			fileToDelete.delete();
		}
	}
	
	/**
	 * Write graph to given outfile in aglappe format
	 * @param outfile
	 * @throws IOException
	 */
	public void write_graph_to_file (String outfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		Out.println("#AGLAPPE GRAPH FILE ver: "+GRAPHFILEFORMATVERSION);
		Out.println("#SEQUENCE: "+sequence);
		Out.println("#PDB: "+pdbCode);
		Out.println("#PDB CHAIN CODE: "+pdbChainCode);
		Out.println("#CHAIN: "+chainCode);
		Out.println("#CT: "+ct);
		Out.println("#CUTOFF: "+cutoff);
		for (Edge pair:contacts){
			int i_resser=pair.i;
			int j_resser=pair.j;
			double weight=pair.weight;
			Out.printf(Locale.US,i_resser+"\t"+j_resser+"\t%6.3f\n",weight);
		}
		Out.close();		
	}
	
	/**
	 * Write graph to given outfile in network(Ioannis) format
	 * @param graphId
	 * @param dirName
	 * @throws IOException
	 */
	public void writeUndirUnweightGraphToNetworkFiles (int graphId, String dirName) throws IOException {
		if (directed) {
			System.err.println("This method is only for undirected graphs!");
			return;
		}
		
		String filePrefix = dirName + "/" + String.valueOf(graphId)+"_"+pdbCode+"_"+chainCode+"_"+ct.replaceAll("/", ".")+"_"+String.valueOf(cutoff).replaceAll("\\.", "_")+"_";
		PrintStream Out = new PrintStream(new FileOutputStream(filePrefix+"edges.txt"));
		for (Edge pair:contacts){
			int i_resser=pair.i;
			int j_resser=pair.j;
			if (i_resser < j_resser) {
				Out.printf(Locale.US,i_resser+"\t"+j_resser+"\t"+1+"\t"+"1.000"+"\n");
			}
		}
		Out.close();
		Out = new PrintStream(new FileOutputStream(filePrefix+"nodes.txt"));
		for (int resser:nodes.keySet()) {
			String res = AAinfo.threeletter2oneletter(getResType(resser));
			Out.printf(Locale.US,resser+"\t"+chainCode+"\t"+resser+"\t"+res+"\t"+getDegree(resser)+"\t"+getDegree(resser)+"\n");
			
		}
		Out.close();
	}
	
	
	/**
	 * Gets list of contacts as a new EdgeSet (deep copied)
	 * 
	 */ 
	public EdgeSet getContacts(){
		return this.contacts.copy();
	}
	
	/**
	 * Checks if the graphs has the given edge.
	 * 
	 * @return false if there is no such contact, else true
	 * */
	public boolean hasContact( Edge e ) {
	    
	    return contacts.contains(e);
	}
	
	/**
	 * Gets iterator on the set of edges/contacts.
	 * 
	 * @return contact iterator
	 * */
	public Iterator<Edge> getContactIterator() {
	    
	    return contacts.iterator();
	}
	
	/**
	 * Gets TreeMap of nodes, deep copying  
	 * 
	 */
	public TreeMap<Integer,String> getNodes(){
		TreeMap<Integer,String> newNodes = new TreeMap<Integer,String>();
		for (int resser:nodes.keySet()){
			newNodes.put(resser, nodes.get(resser));
		}
		return newNodes;
	}
	
	/**
	 * Deep copies this Graph object returning new one
	 * @return
	 */
	public Graph copy(){
		return new Graph(getContacts(),getNodes(),sequence,cutoff,ct,pdbCode,chainCode,pdbChainCode,model,secondaryStructure.copy());		
	}
	
	/**
	 * Returns an int matrix with 1s for contacts and 0s for non contacts, i.e. the contact map
	 * In non-crossed cases this should give us the upper half matrix (contacts are only j>i)
	 * In crossed cases this gives us a full matrix (contacts are both j>i and i>j since they are directed)
	 * @return
	 */
	public int[][] getIntMatrix(){
		// this initialises the matrix to 0 (i.e. no contact)
		int[][] cm = new int[fullLength][fullLength];
		// we put a 1 for all given contacts
		for (Edge cont:contacts){
			int i_resser = cont.i;
			int j_resser = cont.j;
			cm[i_resser-1][j_resser-1]=1;
		}
		return cm;
	}

	/**
	 * Gets a node's residue type given the residue serial
	 * @param resser
	 * @return
	 */
	public String getResType(int resser){
		return nodes.get(resser);
	}
	
	/**
	 * Gets node neighbourhood given a residue serial
	 * @param resser
	 * @return
	 */
	public NodeNbh getNodeNbh(int resser){
		NodeNbh nbh = new NodeNbh(resser, getResType(resser));
		//this could be implemented using the contact map matrix and scanning through 1 column/row
		//it would be just slightly faster, here we do 2*numContacts iterations, using matrix would be only fullLength iterations
		//however we would then have the overhead of creating the matrix
		for (Edge cont:contacts){
			if (cont.i==resser) nbh.put(cont.j, nodes.get(cont.j));
			if (cont.j==resser) nbh.put(cont.i, nodes.get(cont.i));
		}
		return nbh;
	}
	
	/**
	 * Gets edge neighbourhood (common neighbourhood) given a residue serial pair
	 * @param i_resser
	 * @param j_resser
	 * @return
	 */
	public EdgeNbh getEdgeNbh(int i_resser, int j_resser){
		EdgeNbh nbh = new EdgeNbh(i_resser, getResType(i_resser), j_resser, getResType(j_resser));
		NodeNbh i_nbhd = getNodeNbh(i_resser);
		NodeNbh j_nbhd = getNodeNbh(j_resser);
		if (j_nbhd.size()>=i_nbhd.size()) { //with this we will be slightly faster, always iterating through smallest TreeMap
			for (int resser:i_nbhd.keySet()) {
				if (j_nbhd.containsKey(resser)) nbh.put(resser, i_nbhd.get(resser));
			}
		} else {
			for (int resser:j_nbhd.keySet()) {
				if (i_nbhd.containsKey(resser)) nbh.put(resser, j_nbhd.get(resser));			
			}
		}
		return nbh;
	}
	
	/**
	 * Gets 2nd shell node neighbourhood
	 * @param resser
	 */
	public NodeNbh get2ndshellNodeNbh(int resser){
		// first we create a NodeNbh object for the second shell, central residue is given resser
		NodeNbh nbh2ndshell = new NodeNbh(resser,getResType(resser));
		// we get 1st neighbourhood
		NodeNbh nbh = this.getNodeNbh(resser);
		for (int nb:nbh.keySet()){
			NodeNbh nbh2 = this.getNodeNbh(nb); // for each first neighbour we take its neighbourhood
			for (int nb2:nbh2.keySet()){
				if (nb2!=resser && !nbh.containsKey(nb2)){ // if the 2nd neighbour nb2 is not the given resser or is not a 1st neighbour
					nbh2ndshell.put(nb2, getResType(nb2));
				}
			}
		}
		return nbh2ndshell;
	}
	
	public void addEdge(Edge cont){
		if (!directed && cont.i>cont.j){
			// we invert in case of undirected and i>j because in undirected we have only the half of the matrix j>i
			// if we added an edge i>j it could happen that the edge was already there but inverted and wouldn't be detected as a duplicate
			cont = new Edge(cont.j,cont.i);	
		}
		contacts.add(cont); // contacts is a TreeSet and thus takes care of duplicates
		int oldNumContacts = numContacts;
		numContacts=getNumContacts();
		// if number of contacts changed that means we actually added a new contact and thus we modified the graph
		if (numContacts!=oldNumContacts) modified=true;
		
	}
	
	public void delEdge(Edge cont){
		if (!directed && cont.i>cont.j){
			// we invert in case of undirected and i>j because in undirected we have only the half of the matrix j>i
			// if we try to delete an edge i>j it won't be there, we have to invert it and then try to delete
			cont = new Edge(cont.j,cont.i); 
		}
		contacts.remove(cont);
		int oldNumContacts = numContacts;
		numContacts=getNumContacts();
		// if number of contacts changed that means we actually added a new contact and thus we modified the graph
		if (numContacts!=oldNumContacts) modified=true;
	}
	
	public void restrictContactsToMaxRange(int range){
		EdgeSet edgesToDelete = new EdgeSet();
		for (Edge cont:contacts){
			if (cont.getRange()>range) edgesToDelete.add(cont);
		}
		for (Edge cont:edgesToDelete){
			delEdge(cont);
		}
		maxSeqSep = range;
	}
	
	public void restrictContactsToMinRange(int range){
		EdgeSet edgesToDelete = new EdgeSet();
		for (Edge cont:contacts){
			if (cont.getRange()<range) edgesToDelete.add(cont);
		}
		for (Edge cont:edgesToDelete){
			delEdge(cont);
		}
		minSeqSep = range;
	}

	/**
	 * Returns a HashMap with all edge neighbourhood sizes (if they are >0) for each cell in the contact map
	 * @return
	 */
	public HashMap<Edge,Integer> getAllEdgeNbhSizes() {
		HashMap<Edge,Integer> sizes = new HashMap<Edge, Integer>();
		if (!directed) {
			for (int i=1; i<fullLength;i++){
				for (int j=i+1; j<fullLength;j++){
					int size = getEdgeNbh(i, j).size();
					if (size>0)	sizes.put(new Edge(i,j), size);
				}
			}			
		} else {
			for (int i=1; i<fullLength;i++){
				for (int j=1; j<fullLength;j++){
					if (i!=j){
						int size = getEdgeNbh(i, j).size();
						if (size>0) sizes.put(new Edge(i,j), size);
					}
				}
			}
		}
		return sizes;
	}

	//TODO not sure what kind of return we want, for now is a HashMap with three graph objects 
	public HashMap<String,Graph> compare(Graph other) throws Exception{
		//first check that other has same sequence than this, otherwise throw exception
		if (this.getFullLength()!=other.getFullLength()) {
			//TODO throw specific exception
			throw new Exception("Sequence of 2 graphs to compare differ, can't compare them.");
		} else {
			for (int resser:this.nodes.keySet()){
				String this_res = this.nodes.get(resser);
				String other_res = other.nodes.get(resser);
				if (!this_res.equals("X") && !other_res.equals("X") && !this_res.equals(other_res)) {
					//TODO throw specific exception
					throw new Exception("Sequence of 2 graphs to compare differ, can't compare them.");
				}
			}
		}
		
		EdgeSet common = new EdgeSet();
		EdgeSet onlythis = new EdgeSet();
		EdgeSet onlyother = new EdgeSet();
		for (Edge cont:this.contacts){
			if (other.contacts.contains(cont)) {
				common.add(cont);
			} else {
				onlythis.add(cont);
			}
		}
		for (Edge cont:other.contacts){
			if (!this.contacts.contains(cont)){
				onlyother.add(cont);
			}
		}
		Graph commongraph = new Graph (common,getNodes(),sequence,cutoff,ct,pdbCode,chainCode,pdbChainCode,model,secondaryStructure.copy());
		Graph onlythisgraph = new Graph (onlythis,getNodes(),sequence,cutoff,ct,pdbCode,chainCode,pdbChainCode,model,secondaryStructure.copy());
		Graph onlyothergraph = new Graph (onlyother,getNodes(),sequence,cutoff,ct,other.pdbCode,other.chainCode,other.pdbChainCode,model,secondaryStructure.copy());
		HashMap<String,Graph> result = new HashMap<String,Graph>();
		result.put("common", commongraph);
		result.put("onlythis", onlythisgraph);
		result.put("onlyother",onlyothergraph);
		return result;
	}
	
	public boolean isModified(){
		return modified;
	}
	
	public boolean isDirected(){
		return directed;
	}
	
	public String getPdbCode() {
		return pdbCode;
	}
	
	public String getPdbChainCode(){
		return pdbChainCode;
	}
	
	public String getChainCode(){
		return chainCode;
	}
	
	public int getModel() {
		return model;
	}
	
	public String getSequence(){
		return sequence;
	}
	
	/**
	 * Sets the sequence of this graph to the given one.
	 * It will also set the nodes labels to the new sequence.
	 * At the moment there are no checks on the input, if sequence of
	 * diferrent length than current is passed there will be unexpected
	 * results 
	 * @param sequence
	 */
	public void setSequence(String sequence) {
		if (sequence.length()!=this.fullLength) {
			//TODO throw exception
			System.err.println("Setting sequence to another sequence with different length. This can cause problems!");
		}
		this.sequence = sequence;
		Set<Integer> allnodes = nodes.keySet();
		for (int node:allnodes) {
			nodes.put(node,AAinfo.oneletter2threeletter(Character.toString(sequence.charAt(node-1))));
		}
	}
	
	public int getFullLength(){
		return fullLength;
	}
	
	public int getObsLength(){
		return obsLength;
	}
	
	public int getNumContacts(){
		// in theory we could return just numContacts, because we have taken care of updating it every time contacts changed
		// however we call directly contacts.size() as I feel is safer
		return contacts.size(); 
	}
	
	/**
	 * Returns the contact type of this graph.
	 * @return the contact type
	 */
	public String getContactType() {
		return ct;
	}
	
	/**
	 * Sets the contact type of this graph.
	 * @param ct the contact type
	 */
	public void setContactType(String ct) {
		this.ct=ct;
	}
	
	/**
	 * Returns the distance cutoff for this graph.
	 * @return the distance cutoff
	 */
	public double getCutoff(){
		return cutoff;
	}
	
	/**
	 * Sets the distance cutoff for this graph.
	 * @param distCutoff the distance cutoff
	 */
	public void setCutoff(double distCutoff) {
		this.cutoff = distCutoff;
	}
	
	public boolean containsContact(Edge cont){
		// be careful with order, this checks strictly whether the cont.i, cont.j is given, strictly in that order!
		// in undirected case contacts are stored only in 1 direction (j>i) and thus if wrong order given it won't be found 
		return contacts.contains(cont);
	}
	
	public void resetContacts(){
		this.contacts = new EdgeSet();
	}
	
	public int getDegree(int resser){
		if (directed) {
			System.err.println("Can't get degree for a directed graph, only in or out degree");
			return 0;
		}
		int k = 0;
		for (Edge cont:contacts){
			if (cont.i==resser || cont.j==resser) {
				k++;
			}
		}
		return k;
	}
	
	public int getInDegree(int resser){
		if (!directed){
			System.err.println("Can't get in degree for an undirected graph");
			return 0;
		}
		int k = 0;
		for (Edge cont:contacts){
			if (cont.j==resser) {
				k++;
			}
		}		
		return k;
	}
	
	public int getOutDegree(int resser){
		if (!directed){
			System.err.println("Can't get out degree for an undirected graph");
			return 0;
		}
		int k = 0;
		for (Edge cont:contacts){
			if (cont.i==resser) {
				k++;
			}
		}		
		return k;
	}

	// secondary structure related methods
	
	/** 
	 * Returns true if secondary structure information is available, false otherwise. 
	 */
	public boolean hasSecondaryStructure() {
		return !this.secondaryStructure.isEmpty();
	}
	
	/**
	 * Returns the secondary structure annotation object of this graph.
	 */
	public SecondaryStructure getSecondaryStructure() {
		return this.secondaryStructure;
	}

	/**
	 * Evaluate this graph (assuming it is a prediction) against an original graph
	 * @param originalGraph
	 * @return
	 */
	public PredEval evaluatePrediction(Graph originalGraph) {
		return evaluatePrediction(originalGraph, 1);
	}
	
	/**
	 * Evaluate this graph (assuming it is a prediction) against an original graph,
	 * considering only edges with sequence separation at least minSeqSep.
	 * @param originalGraph
	 * @return
	 */
	public PredEval evaluatePrediction(Graph originalGraph, int minSeqSep) {
		// total predicted contacts
		int predicted = 0;
		for(Edge e:this.getContacts()) {
			if(Math.abs(e.j-e.i) >= minSeqSep) {
				predicted++;
			}
		}
		
		// total native contacts
		int original = 0;
		for(Edge e:originalGraph.getContacts()) {
			if(Math.abs(e.j-e.i) >= minSeqSep) {
				original++;
			}
		}		
				
		// total size of contact map (potential contacts)
		int cmtotal = 0;
		if (originalGraph.isDirected()){
			cmtotal = (originalGraph.getFullLength()-(minSeqSep-1))*(originalGraph.getFullLength()-minSeqSep);
		} else {
			cmtotal = (int)(((originalGraph.getFullLength()-(minSeqSep-1))*(originalGraph.getFullLength()-minSeqSep))/2);
		}
		int TruePos=0, FalsePos=0, TrueNeg=0, FalseNeg=0;
		
		EdgeSet origContacts = originalGraph.getContacts();
		EdgeSet predictedContacts = this.getContacts();

		// directed/ non-directed graphs should be both fine with this code (as long as in the predicted directed graph we store the contacts as j>i)
		// the only thing that changes between directed/non-directed is the count of total cells in contact map (taken care for above)
		for (Edge predictedCont:predictedContacts){
			if(Math.abs(predictedCont.j-predictedCont.i) >= minSeqSep) {
				//System.out.println(predictedCont);
				if (origContacts.contains(predictedCont)) {
					TruePos++;
				}
				else {
					FalsePos++;
				}
			}
		}

		for (Edge origCont:origContacts) {
			if(Math.abs(origCont.j-origCont.i) >= minSeqSep) {
				//System.out.println(origCont);
				if (!predictedContacts.contains(origCont)) {
					//System.out.println(origCont);
					FalseNeg++;
				}
			}
		}
		TrueNeg=cmtotal-TruePos-FalsePos-FalseNeg;
		PredEval eval = new PredEval(TruePos,FalsePos,TrueNeg,FalseNeg,0,predicted,original,cmtotal);
		return eval;
	}
	
	private static String quote(String s) {
		return ("'"+s+"'");
	}
}


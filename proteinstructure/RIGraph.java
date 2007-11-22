package proteinstructure;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collection;
import java.util.HashMap;
import java.util.Locale;
import java.util.TreeMap;

import tools.MySQLConnection;

import edu.uci.ics.jung.graph.util.Pair;

/**
 * A Residue Interaction Graph
 *
 */
public class RIGraph extends ProtStructGraph<RIGNode,RIGEdge> {

	private static final long serialVersionUID = 1L;
	
	private static final String SINGLEMODELS_DB = "ioannis";
	
	// fields
	protected double distCutoff;
	protected String contactType;				// use AAinfo.isValidContactType() to test for validity

	public RIGraph() {
		super();
		this.distCutoff=0;
		this.contactType=null;
	}

	/**
	 * Constructs a RIGraph with a sequence but no edges
	 * @param sequence
	 */
	public RIGraph(String sequence) {
		super();
		this.sequence = sequence;
		this.fullLength = sequence.length();
		this.distCutoff=0;
		this.contactType=null;
		serials2nodes = new TreeMap<Integer,RIGNode>();
		for(int i=0; i < sequence.length(); i++) {
			RIGNode node = new RIGNode(i+1,AAinfo.oneletter2threeletter(Character.toString(sequence.charAt(i))));
			this.addVertex(node);
			serials2nodes.put(i+1, node);
		}
	}
	
	/**
	 * Returns the contact type of this RIGraph
	 * @return
	 */
	public String getContactType() {
		return contactType;
	}

	/**
	 * Sets the contact type of this RIGraph
	 * @param ct the contact type
	 */
	public void setContactType(String contactType) {
		this.contactType=contactType;
	}
	
	/**
	 * Returns the distance cutoff for this RIGraph.
	 * @return the distance cutoff
	 */
	public double getCutoff(){
		return distCutoff;
	}
	
	/**
	 * Sets the distance cutoff for this RIGraph.
	 * @param distCutoff the distance cutoff
	 */
	public void setCutoff(double distCutoff) {
		this.distCutoff = distCutoff;
	}
	
	/**
	 * Returns a RIGNbhood that contains the neighbourhood of given RIGNode
	 * @param node
	 * @return
	 */
	public RIGNbhood getNbhood (RIGNode node) {
		Collection<RIGNode> nbs = this.getNeighbors(node);
		RIGNbhood nbhood = new RIGNbhood(node);
		for (RIGNode nb:nbs) {
			nbhood.put(nb.getResidueSerial(), nb);
		}
		return nbhood;
	}
	
	/**
	 * Returns a RIGNbhood that contains the 2nd shell neighbourhood of given RIGNode
	 * @param node
	 * @return
	 */
	public RIGNbhood getSecondShellNbhood (RIGNode node) {
		Collection<RIGNode> nbs = this.getNeighbors(node);
		RIGNbhood nbhood = new RIGNbhood(node);
		for (RIGNode nb:nbs) {
			for (RIGNode nb2:this.getNeighbors(nb)) {
				if (nb2!=node) {
					// RIGNbhood is a TreeMap that should take care of not inserting duplicates
					nbhood.put(nb2.getResidueSerial(),nb2);
				}
			}
		}
		return nbhood;
	}

	/**
	 * Returns a RIGCommonNbhood that contains common neighbours of given RIGNodes iNode, jNode
	 * @param iNode
	 * @param jNode
	 * @return
	 */
	public RIGCommonNbhood getCommonNbhood(RIGNode iNode, RIGNode jNode) {
		Collection<RIGNode> iNbs = this.getNeighbors(iNode);
		Collection<RIGNode> jNbs = this.getNeighbors(jNode);
		boolean connected = false;
		//NOTE in DIRECTED case this means strictly an edge from iNode to jNode 
		if (this.findEdge(iNode, jNode)!=null) connected = true;
		RIGCommonNbhood comNbhood = new RIGCommonNbhood(iNode, jNode, connected);
		for (RIGNode iNb: iNbs) {
			if (jNbs.contains(iNb)) {
				comNbhood.put(iNb.getResidueSerial(), iNb);
			}
		}
		return comNbhood;
	}
	
	/**
	 * Returns all common neighborhood sizes for each cell of the contact map (contact or non-contact) 
	 * @return
	 */
	public HashMap<Pair<RIGNode>,Integer> getAllCommonNbhSizes() {
		HashMap<Pair<RIGNode>,Integer> comNbhSizes = new HashMap<Pair<RIGNode>, Integer>();
		boolean directed = this.isDirected();
		for (RIGNode n1:this.getVertices()) {
			for (RIGNode n2:this.getVertices()) {
				if (directed) {
					if (n1!=n2) {
						comNbhSizes.put(new Pair<RIGNode>(n1,n2),getCommonNbhood(n1, n2).size());					
					}					
				} else {
					if (n1.getResidueSerial()<n2.getResidueSerial()) {
						comNbhSizes.put(new Pair<RIGNode>(n1,n2),getCommonNbhood(n1, n2).size());					
					}
				}
			}
		}
		return comNbhSizes;
	}

	public int getContactRange(RIGEdge edge) {
		Pair<RIGNode> pair = this.getEndpoints(edge);
		return Math.abs(pair.getFirst().getResidueSerial()-pair.getSecond().getResidueSerial());
	}

	//TODO evaluatePrediction methods should be in ProtStructGraph. 
	//     But to be able to put them there we would need to pass here a Transformer that gets atom or residue serials depending if we are in AI or RI Graph 
	/**
	 * Evaluate this graph (assuming it is a prediction) against an original graph
	 * @param originalGraph
	 * @return
	 */
	public PredEval evaluatePrediction(RIGraph originalGraph) {
		return evaluatePrediction(originalGraph, 1);
	}
	
	/**
	 * Evaluate this graph (assuming it is a prediction) against an original graph,
	 * considering only edges with sequence separation at least minSeqSep.
	 * @param originalGraph
	 * @param minSeqSep
	 * @return
	 */
	public PredEval evaluatePrediction(RIGraph originalGraph, int minSeqSep) {
		
		Collection<RIGEdge> predictedContacts = this.getEdges();
		Collection<RIGEdge> origContacts = originalGraph.getEdges();
		// total predicted contacts
		int predicted = 0;
		for(RIGEdge e:predictedContacts) {
			if(this.getContactRange(e) >= minSeqSep) {
				predicted++;
			}
		}
		
		// total native contacts
		int original = 0;
		for(RIGEdge e:origContacts) {
			if(originalGraph.getContactRange(e) >= minSeqSep) {
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
		
		// directed/ non-directed graphs should be both fine with this code 
		// the only thing that changes between directed/non-directed is the count of total cells in contact map (taken care for above)
		for (RIGEdge predictedCont:predictedContacts){
			if(this.getContactRange(predictedCont) >= minSeqSep) {
				Pair<RIGNode> predNodePair = this.getEndpoints(predictedCont);
				RIGNode node1inOrig = originalGraph.getNodeFromSerial(predNodePair.getFirst().getResidueSerial());
				RIGNode node2inOrig = originalGraph.getNodeFromSerial(predNodePair.getSecond().getResidueSerial());
				//NOTE order of nodes in findEdge doesn't matter if UNDIRECTED.
				//It does matter if DIRECTED. However even in that case we are fine because we use same order in this graph 
				if (originalGraph.findEdge(node1inOrig, node2inOrig)!=null) {
					TruePos++;
				}
				else {
					FalsePos++;
				}
			}
		}

		for (RIGEdge origCont:origContacts) {
			if(originalGraph.getContactRange(origCont) >= minSeqSep) {
				Pair<RIGNode> origNodePair = originalGraph.getEndpoints(origCont);
				RIGNode node1inPred = this.getNodeFromSerial(origNodePair.getFirst().getResidueSerial());
				RIGNode node2inPred = this.getNodeFromSerial(origNodePair.getSecond().getResidueSerial());
				//NOTE order of nodes in findEdge doesn't matter if UNDIRECTED.
				//It does matter if DIRECTED. However even in that case we are fine because we use same order in originalGraph				
				if (this.findEdge(node1inPred,node2inPred)==null) {
					FalseNeg++;
				}
			}
		}
		TrueNeg=cmtotal-TruePos-FalsePos-FalseNeg;
		PredEval eval = new PredEval(TruePos,FalsePos,TrueNeg,FalseNeg,0,predicted,original,cmtotal);
		return eval;
	}
	
	/**
	 * Write graph to given db, using our db graph aglappe format, 
	 * i.e. tables: chain_graph, single_model_graph, single_model_node, single_model_edge
	 * @param conn
	 * @param db
	 * @throws SQLException
	 */
	//TODO we might want to move this to a graph i/o class
	public void write_graph_to_db(MySQLConnection conn, String db) throws SQLException{
		
		// values we fix to constant 
		String CW = "1";
		String CR = "(true)";
		String EXPBB = "0";
		String ctStr = contactType;
		String weightedStr = "0";
		String directedStr = isDirected()?"1":"0";
		
		if (contactType.endsWith("_CAGLY")) {
			ctStr = contactType.replace("_CAGLY", "");
		}
		if (ctStr.equals("ALL")) {
			ctStr = "BB+SC+BB/SC";
		}
		if (AAinfo.isValidMultiAtomContactType(contactType)) {
			CW = ctStr;
			weightedStr = "1";
		}
		if (contactType.endsWith("_CAGLY") || contactType.equals("Cb")) {
			EXPBB = "-1";
		}
		if (minSeqSep != -1) {
			CR = "((i_cid!=j_cid)OR(abs(i_num-j_num)>="+minSeqSep+"))";
		}
				
		int pgraphid=0;
		int graphid=0;
		String sql = "SELECT graph_id FROM "+db+".chain_graph " +
					" WHERE accession_code='"+pdbCode+"' AND pchain_code='"+chainCode+"'" +
					" AND model_serial = "+model+" AND dist = "+distCutoff+" AND expBB = '"+EXPBB+"'" + 
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
					"VALUES ('"+pdbCode+"', "+pdbChainCodeStr+",'"+chainCode+"', "+model+", "+distCutoff+", "+EXPBB+", 'rc-cutoff', "+getFullLength()+", "+getObsLength()+", "+getObsLength()+", "+secondaryStructure.getNumElements()+", now())";
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
				" dist="+distCutoff+" AND expBB="+EXPBB+" AND CW='"+CW+"' AND CT='"+ctStr+"' AND CR='"+CR+"';";
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			singlemodelid = rsst.getInt(1);
		}
		rsst.close();
		// and then insert to single_model_graph
		sql = "INSERT INTO "+db+".single_model_graph (pgraph_id,graph_type,accession_code,single_model_id,dist,expBB,CW,CT,CR,w,d,num_nodes,date) " +
				" VALUES ("+pgraphid+", 'chain', '"+pdbCode+"', "+singlemodelid+", "+distCutoff+", "+EXPBB+", '"+CW+"','"+ctStr+"', '"+CR+"', "+weightedStr+", "+directedStr+", "+getObsLength()+", now())";
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
		for (int resser:serials2nodes.keySet()) {
			RIGNode node = serials2nodes.get(resser);
			String res = AAinfo.threeletter2oneletter(node.getResidueType());
			RIGNbhood nbh = this.getNbhood(node);
			String secStructType = null;
			String secStructId = null;
			String sheetSerial = null;
			String turn = null;
			SecStrucElement sselem = node.getSecStrucElement();
			if (sselem!=null){
				secStructType = quote(Character.toString(sselem.getType()));
				secStructId = quote(sselem.getId());
				char sheetSerialChar = sselem.getSheetSerial();
				if (sheetSerialChar != 0) {
					sheetSerial = quote(Character.toString(sheetSerialChar));
				}
				turn = sselem.isTurn()?"1":"0";
			}
			if (isDirected()){  // we insert k(=k_in+k_out), k_in and k_out
				sql = "INSERT INTO "+db+".single_model_node "+
					" (graph_id, node_id, cid, num, res, "+
					" sstype, ssid, sheet_serial, turn, "+
					" k, k_in, k_out, "+
					" n, nwg, n_num) " +
					" VALUES ("+graphid+", "+(maxNodeId+resser)+", '"+chainCode+"', "+resser+", '"+res+"', "+
					" "+secStructType+", "+secStructId+", "+sheetSerial+", "+turn+", "+
					(getPredecessorCount(node)+getSuccessorCount(node))+", "+getPredecessorCount(node)+", "+getSuccessorCount(node)+", "+
					"'"+nbh.getMotifNoGaps()+"', '"+nbh.getMotif()+"', '"+nbh.getCommaSeparatedResSerials()+"')";
			} else {		// we insert k (and no k_in or k_out)
				sql = "INSERT INTO "+db+".single_model_node "+
				" (graph_id, node_id, cid, num, res, "+
				" sstype, ssid, sheet_serial, turn, "+
				" k, n, nwg, n_num) " +
				" VALUES ("+graphid+", "+(maxNodeId+resser)+", '"+chainCode+"', "+resser+", '"+res+"', "+
				" "+secStructType+", "+secStructId+", "+sheetSerial+", "+turn+", "+
				getNeighborCount(node)+", '"+nbh.getMotifNoGaps()+"', '"+nbh.getMotif()+"', '"+nbh.getCommaSeparatedResSerials()+"')";
			}
			stmt.executeUpdate(sql);
		}
		
		// inserting edges
		// get the max weight
		double maxWeight = 0;
		for (RIGEdge cont:getEdges()) {
			maxWeight = (maxWeight<cont.getAtomWeight())?cont.getAtomWeight():maxWeight;
		}
		for (RIGEdge cont:getEdges()){
			RIGNode i_node = getEndpoints(cont).getFirst();
			RIGNode j_node = getEndpoints(cont).getSecond();
			String i_res = AAinfo.threeletter2oneletter(i_node.getResidueType());
			String j_res = AAinfo.threeletter2oneletter(j_node.getResidueType());

			String i_secStructType = null;
			String i_secStructId = null;
			String i_sheetSerial = null;
			String i_turn = null;
			SecStrucElement i_sselem = i_node.getSecStrucElement();
			if (i_sselem!=null){
				i_secStructType = quote(Character.toString(i_sselem.getType()));
				i_secStructId = quote(i_sselem.getId());
				char sheetSerialChar = i_sselem.getSheetSerial();
				if (sheetSerialChar != 0) {
					i_sheetSerial = quote(Character.toString(sheetSerialChar));
				}
				i_turn = i_sselem.isTurn()?"1":"0";
			}
			
			String j_secStructType = null;
			String j_secStructId = null;
			String j_sheetSerial = null;
			String j_turn = null;
			SecStrucElement j_sselem = j_node.getSecStrucElement();
			if (j_sselem!=null){
				j_secStructType = quote(Character.toString(j_sselem.getType()));
				j_secStructId = quote(j_sselem.getId());
				char sheetSerialChar = j_sselem.getSheetSerial();
				if (sheetSerialChar != 0) {
					j_sheetSerial = quote(Character.toString(sheetSerialChar));
				}
				j_turn = j_sselem.isTurn()?"1":"0";
			}
			
			sql = "INSERT INTO "+db+".single_model_edge "+
					" (graph_id, i_node_id, i_cid, i_num, i_res, i_sstype, i_ssid, i_sheet_serial, i_turn, "+
					" j_node_id, j_cid, j_num, j_res, j_sstype, j_ssid, j_sheet_serial, j_turn, weight, norm_weight) " +
					" VALUES ("+graphid+", "+(maxNodeId+i_node.getResidueSerial())+", '"+chainCode+"', "+i_node.getResidueSerial()+", '"+i_res+"', "+i_secStructType+", "+i_secStructId+", "+i_sheetSerial+", "+i_turn+", "+
					(maxNodeId+j_node.getResidueSerial())+", '"+chainCode+"', "+j_node.getResidueSerial()+", '"+j_res+"', "+j_secStructType+", "+j_secStructId+", "+j_sheetSerial+", "+j_turn+", "+
					cont.getAtomWeight()+", "+(cont.getAtomWeight()/maxWeight)+")";
			stmt.executeUpdate(sql);
			if(!isDirected()) {// we want both side of the matrix in the table to follow Ioannis' convention
				// so we insert the reverse contact by swapping i, j in insertion
				sql = "INSERT INTO "+db+".single_model_edge "+
				" (graph_id, i_node_id, i_cid, i_num, i_res, i_sstype, i_ssid, i_sheet_serial, i_turn, "+
				" j_node_id, j_cid, j_num, j_res, j_sstype, j_ssid, j_sheet_serial, j_turn, weight, norm_weight) " +
				" VALUES ("+graphid+", "+(maxNodeId+j_node.getResidueSerial())+", '"+chainCode+"', "+j_node.getResidueSerial()+", '"+j_res+"', "+j_secStructType+", "+j_secStructId+", "+j_sheetSerial+", "+j_turn+", "+
				(maxNodeId+i_node.getResidueSerial())+", '"+chainCode+"', "+i_node.getResidueSerial()+", '"+i_res+"', "+i_secStructType+", "+i_secStructId+", "+i_sheetSerial+", "+i_turn+", "+
				cont.getAtomWeight()+", "+(cont.getAtomWeight()/maxWeight)+")";
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
	//TODO we might want to move this to a graph i/o class
	public void write_graph_to_db_fast(MySQLConnection conn, String db) throws SQLException, IOException {
		
		// values we fix to constant 
		String CW = "1";
		String CR = "(true)";
		String EXPBB = "0";
		String ctStr = contactType;
		String weightedStr = "0";
		String directedStr = isDirected()?"1":"0";
		
		if (contactType.endsWith("_CAGLY")) {
			ctStr = contactType.replace("_CAGLY", "");
		}
		if (ctStr.equals("ALL")) {
			ctStr = "BB+SC+BB/SC";
		}
		if (AAinfo.isValidMultiAtomContactType(contactType)) {
			CW = ctStr;
			weightedStr = "1";
		}
		if (contactType.endsWith("_CAGLY") || contactType.equals("Cb")) {
			EXPBB = "-1";
		}
		if (minSeqSep != -1) {
			CR = "((i_cid!=j_cid)OR(abs(i_num-j_num)>="+minSeqSep+"))";
		}
				
		int pgraphid=0;
		int graphid=0;
		String sql = "SELECT graph_id FROM "+db+".chain_graph " +
					" WHERE accession_code='"+pdbCode+"' AND pchain_code='"+chainCode+"'" +
					" AND model_serial = "+model+" AND dist = "+distCutoff+" AND expBB = '"+EXPBB+"'" + 
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
					"VALUES ('"+pdbCode+"', "+pdbChainCodeStr+",'"+chainCode+"', "+model+", "+distCutoff+", "+EXPBB+", 'rc-cutoff', "+getFullLength()+", "+getObsLength()+", "+getObsLength()+", "+secondaryStructure.getNumElements()+", now())";
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
				" dist="+distCutoff+" AND expBB="+EXPBB+" AND CW='"+CW+"' AND CT='"+ctStr+"' AND CR='"+CR+"';";
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			singlemodelid = rsst.getInt(1);
		}
		rsst.close();
		// and then insert to single_model_graph
		sql = "INSERT INTO "+db+".single_model_graph (pgraph_id,graph_type,accession_code,single_model_id,dist,expBB,CW,CT,CR,w,d,num_nodes,date) " +
				" VALUES ("+pgraphid+", 'chain', '"+pdbCode+"', "+singlemodelid+", "+distCutoff+", "+EXPBB+", '"+CW+"','"+ctStr+"', '"+CR+"', "+weightedStr+", "+directedStr+", "+getObsLength()+", now())";
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
		
		for (int resser:serials2nodes.keySet()) {
			
			RIGNode node = serials2nodes.get(resser);
			String res = AAinfo.threeletter2oneletter(node.getResidueType());
			RIGNbhood nbh = this.getNbhood(node);			
			String secStructType = "\\N";
			String secStructId = "\\N";
			String sheetSerial = "\\N";
			String turn = null;
			SecStrucElement sselem = node.getSecStrucElement();
			if (sselem!=null){
				secStructType = Character.toString(sselem.getType());
				secStructId = sselem.getId();
				char sheetSerialChar = sselem.getSheetSerial();
				if (sheetSerialChar != 0) {
					sheetSerial = Character.toString(sheetSerialChar);
				}
				turn = sselem.isTurn()?"1":"0";
			}
			if (isDirected()){  // we insert k(=k_in+k_out), k_in and k_out
				nodesOut.println(graphid+"\t"+(maxNodeId+resser)+"\t"+chainCode+"\t"+resser+"\t"+res+"\t"+
					secStructType+"\t"+secStructId+"\t"+sheetSerial+"\t"+turn+"\t"+
					(getPredecessorCount(node)+getSuccessorCount(node))+"\t"+getPredecessorCount(node)+"\t"+getSuccessorCount(node)+"\t"+
					nbh.getMotifNoGaps()+"\t"+nbh.getMotif()+"\t"+nbh.getCommaSeparatedResSerials());
			} else {		// we insert k (and no k_in or k_out)
				nodesOut.println(graphid+"\t"+(maxNodeId+resser)+"\t"+chainCode+"\t"+resser+"\t"+res+"\t"+
						secStructType+"\t"+secStructId+"\t"+sheetSerial+"\t"+turn+"\t"+
						getNeighborCount(node)+"\t"+"\\N"+"\t"+"\\N"+"\t"+
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
		for (RIGEdge cont:getEdges()) {
			maxWeight = (maxWeight<cont.getAtomWeight())?cont.getAtomWeight():maxWeight;
		}
		for (RIGEdge cont:getEdges()){
			RIGNode i_node = getEndpoints(cont).getFirst();
			RIGNode j_node = getEndpoints(cont).getSecond();
			String i_res = AAinfo.threeletter2oneletter(i_node.getResidueType());
			String j_res = AAinfo.threeletter2oneletter(j_node.getResidueType());
			
			String i_secStructType = "\\N";
			String i_secStructId = "\\N";
			String i_sheetSerial = "\\N";
			String i_turn = null;
			SecStrucElement i_sselem = i_node.getSecStrucElement();
			if (i_sselem!=null){
				i_secStructType = Character.toString(i_sselem.getType());
				i_secStructId = i_sselem.getId();
				char sheetSerialChar = i_sselem.getSheetSerial();
				if (sheetSerialChar != 0) {
					i_sheetSerial = Character.toString(sheetSerialChar);
				}
				i_turn = i_sselem.isTurn()?"1":"0";
			}
			
			String j_secStructType = "\\N";
			String j_secStructId = "\\N";
			String j_sheetSerial = "\\N";
			String j_turn = null;
			SecStrucElement j_sselem = j_node.getSecStrucElement();
			if (j_sselem!=null){
				j_secStructType = Character.toString(j_sselem.getType());
				j_secStructId = j_sselem.getId();
				char sheetSerialChar = j_sselem.getSheetSerial();
				if (sheetSerialChar != 0) {
					j_sheetSerial = Character.toString(sheetSerialChar);
				}
				j_turn = j_sselem.isTurn()?"1":"0";
			}
			
			edgesOut.println(graphid+"\t"+(maxNodeId+i_node.getResidueSerial())+"\t"+chainCode+"\t"+i_node.getResidueSerial()+"\t"+i_res+"\t"+i_secStructType+"\t"+i_secStructId+"\t"+i_sheetSerial+"\t"+i_turn+"\t"+
					(maxNodeId+j_node.getResidueSerial())+"\t"+chainCode+"\t"+j_node.getResidueSerial()+"\t"+j_res+"\t"+j_secStructType+"\t"+j_secStructId+"\t"+j_sheetSerial+"\t"+j_turn+"\t"+
					cont.getAtomWeight()+"\t"+(cont.getAtomWeight()/maxWeight));
			if(!isDirected()) {// we want both side of the matrix in the table to follow Ioannis' convention
				// so we insert the reverse contact by swapping i, j in insertion
				edgesOut.println(graphid+"\t"+(maxNodeId+j_node.getResidueSerial())+"\t"+chainCode+"\t"+j_node.getResidueSerial()+"\t"+j_res+"\t"+j_secStructType+"\t"+j_secStructId+"\t"+j_sheetSerial+"\t"+j_turn+"\t"+
						(maxNodeId+i_node.getResidueSerial())+"\t"+chainCode+"\t"+i_node.getResidueSerial()+"\t"+i_res+"\t"+i_secStructType+"\t"+i_secStructId+"\t"+i_sheetSerial+"\t"+i_turn+"\t"+
						cont.getAtomWeight()+"\t"+(cont.getAtomWeight()/maxWeight));
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
	
	/** Single quotes the given string */
	private String quote(String s) {
		return ("'"+s+"'");
	}

	/**
	 * Write graph to given outfile in aglappe format
	 * @param outfile
	 * @throws IOException
	 */
	//TODO we might want to move this to a graph i/o class
	public void writeToFile (String outfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		Out.println("#AGLAPPE GRAPH FILE ver: "+GRAPHFILEFORMATVERSION);
		Out.println("#SEQUENCE: "+sequence);
		Out.println("#PDB: "+pdbCode);
		Out.println("#PDB CHAIN CODE: "+pdbChainCode);
		Out.println("#CHAIN: "+chainCode);
		Out.println("#CT: "+contactType);
		Out.println("#CUTOFF: "+distCutoff);
		for (RIGEdge cont:getEdges()){
			Pair<RIGNode> pair = getEndpoints(cont);
			int i_resser=pair.getFirst().getResidueSerial();
			int j_resser=pair.getSecond().getResidueSerial();
			//BEWARE!! here we write weights while in writeToDb we write atomWeights (consistent with what we do in FileRIGraph) TODO do we want this behaviour?
			double weight=cont.getWeight();
			Out.printf(Locale.US,i_resser+"\t"+j_resser+"\t%6.3f\n",weight);
		}
		Out.close();		
	}
}

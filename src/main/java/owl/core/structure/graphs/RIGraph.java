package owl.core.structure.graphs;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Locale;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.vecmath.GMatrix;

import owl.core.structure.AminoAcid;
import owl.core.structure.ContactType;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;
import owl.core.util.IntPairComparator;
import owl.core.util.IntPairSet;
import owl.core.util.MySQLConnection;


import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * A Residue Interaction Graph
 *
 */
public class RIGraph extends ProtStructGraph<RIGNode,RIGEdge> {

	private static final long serialVersionUID = 1L;
	
	private static final String DEFAULT_SINGLEMODELS_DB = "ioannis";
	
	// fields
	protected double distCutoff;
	protected String contactType;				// use AAinfo.isValidContactType() to test for validity
	protected String singleModelsDb;
	
	public RIGraph() {
		super();
		this.distCutoff=0;
		this.contactType=null;
		this.singleModelsDb=DEFAULT_SINGLEMODELS_DB;
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
		this.pdbCode=PdbAsymUnit.NO_PDB_CODE;
		this.chainCode=PdbChain.NO_CHAIN_CODE;
		this.pdbChainCode=PdbChain.NO_PDB_CHAIN_CODE;
		this.singleModelsDb=DEFAULT_SINGLEMODELS_DB;
		for(int i=0; i < sequence.length(); i++) {
			RIGNode node = new RIGNode(i+1,AminoAcid.one2three(sequence.charAt(i)));
			this.addVertex(node); // this also updates the serials2nodes map
		}
	}
	
	/**
	 * Returns the db with the single model ids
	 * @return
	 */
	public String getSingleModelsDb () {
		return singleModelsDb;
	}

	/**
	 * Sets the db with the single model ids
	 * @param ct the contact type
	 */
	public void setSingleModelsDb(String db) {
		this.singleModelsDb=db;
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
	 * Gets the number of observed residues for this RIGraph
	 * See also: getFullLength()
	 * @return the observed length of this graph
	 */
	public int getObsLength() {
		return this.getVertexCount();
	}
	
	/**
	 * Gets the full length of this graph
	 * See also: getObsLength()
	 * @return the full (sequence) length of this graph
	 */
	public int getFullLength() {
		return this.fullLength;
	}
	
	/**
	 * Returns the neighbors of the given node ordered by residueSerial.
	 * This is sometimes needed because the order in which neighbors are returned by getNeighbors() is undefined.
	 * @param node
	 * @return
	 */
	public Collection<RIGNode> getOrderedNeighbors(RIGNode node) {
		LinkedList<RIGNode> orderedNeighbours = new LinkedList<RIGNode>();
		orderedNeighbours.addAll(this.getNeighbors(node));
		Collections.sort(orderedNeighbours, new Comparator<RIGNode>() {
			public int compare(RIGNode n1, RIGNode n2) {
				return (new Integer(n1.getResidueSerial()).compareTo(new Integer(n2.getResidueSerial())));
			}
		});	
		return orderedNeighbours;
	}
	
	/**
	 * Returns a RIGNbhood that contains the neighbourhood of given RIGNode
	 * @param node
	 * @return
	 */
	public RIGNbhood getNbhood (RIGNode node) {
		return new RIGNbhood(node,this.getNeighbors(node));
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
	 * Returns all common neighborhood sizes (if they are >0) for each cell of the contact map (contact or non-contact) 
	 * @return
	 */
	public HashMap<Pair<Integer>,Integer> getAllCommonNbhSizes() {
		HashMap<Pair<Integer>,Integer> comNbhSizes = new HashMap<Pair<Integer>, Integer>();
		boolean directed = this.isDirected();
		for (RIGNode n1:this.getVertices()) {
			for (RIGNode n2:this.getVertices()) {
				if (directed) {
					if (n1!=n2) {
						int size = getCommonNbhood(n1, n2).size();
						if (size>0) comNbhSizes.put(new Pair<Integer>(n1.getResidueSerial(),n2.getResidueSerial()),size);					
					}					
				} else {
					if (n1.getResidueSerial()<n2.getResidueSerial()) {
						int size = getCommonNbhood(n1, n2).size();
						if (size>0) comNbhSizes.put(new Pair<Integer>(n1.getResidueSerial(),n2.getResidueSerial()),size);					
					}
				}
			}
		}
		return comNbhSizes;
	}

	/**
	 * Calculates the contact order of this RIGraph
	 * i.e. average contact range divided by number of observed residues
	 * See wikipedia http://en.wikipedia.org/wiki/Contact_order
	 * 
	 * NOTE: not clear whether this is a standard calculation of contact order,
	 * for instance it doesn't filter the contiguous residues 
	 * @return
	 */
	public double getContactOrder() {
		
		int contRangeSum = 0;
		
		for (RIGEdge edge:this.getEdges()) {			
			int contactRange = getContactRange(edge);
			contRangeSum+=contactRange;
		}
		
		return (double)contRangeSum/(double)(getObsLength()*getEdgeCount());
	}
	
	public int getContactRange(RIGEdge edge) {
		Pair<RIGNode> pair = this.getEndpoints(edge);
		return Math.abs(pair.getFirst().getResidueSerial()-pair.getSecond().getResidueSerial());
	}

	public int getResidueSerial(RIGNode node) {
		return node.getResidueSerial();
	}
	
	public int getFirstResidueSerial() {
		return serials2nodes.firstKey();
	}

	public int getLastResidueSerial() {
		return serials2nodes.lastKey();
	}
	
	/**
	 * Adds an edge between residue serials i, j with default weight
	 * If i or j don't map to RIGNodes in this graph nothing will be added
	 * @param i
	 * @param j
	 * @return true if an edge was added, false otherwise
	 */
	public boolean addEdgeIJ(int i, int j) {
		EdgeType et = EdgeType.UNDIRECTED;
		if (isDirected()) {
			et = EdgeType.DIRECTED;
		}
		if (getNodeFromSerial(i)!=null && getNodeFromSerial(j)!=null) {
			return this.addEdge(new RIGEdge(), getNodeFromSerial(i), getNodeFromSerial(j), et);
		}
		return false;
	}
	
	/**
	 * Adds an edge between residue serials i, j with the given weight.
	 * Weight has to be between 0 and 1.
	 * If i or j don't map to RIGNodes in this graph nothing will be added.
	 * @param i
	 * @param j
	 * @return true if an edge was added, false otherwise
	 */
	public boolean addEdgeIJ(int i, int j, double weight) {
		EdgeType et = EdgeType.UNDIRECTED;
		if (isDirected()) {
			et = EdgeType.DIRECTED;
		}
		if (getNodeFromSerial(i)!=null && getNodeFromSerial(j)!=null) {
			return this.addEdge(new RIGEdge(weight), getNodeFromSerial(i), getNodeFromSerial(j), et);
		}
		return false;
	}
	
	/**
	 * Filter this RIGraph by removing the edges with weight below the given minWeight value.
	 * @param minWeight
	 */
	public void filterByMinWeight(double minWeight) {
		for (RIGEdge edge:this.getEdges()) {
			if (edge.getWeight() < minWeight) {
				this.removeEdge(edge);
			}
		}
	}
	
	/**
	 * Discretize this graph such that all edges with weight >= weightCutoff are set to weight 1
	 * and all other edges are discarded. The different to filterByMinWeight is that all edges
	 * which are kept are set to weight 1.
	 * @param weightCutoff the minimum weight to keep the contact
	 */
	public void discretizeByWeightCutoff(double weightCutoff) {
		for (RIGEdge edge:this.getEdges()) {
			if (edge.getWeight() < weightCutoff) {
				this.removeEdge(edge);
			} else {
				edge.setWeight(1);
			}
		}		
	}
	
	/**
	 * Discretize this graph such that the 'top' contacts with highest weight are kept and set
	 * to weight 1 and all other edges are discarded. When edges have equal weight, the order
	 * is undefined.
	 * @param num the number of contacts to keep 
	 */
	public void discretizeByNumContacts(int top) {
		int numEdges = this.getEdgeCount();
		ArrayList<RIGEdge> edges = new ArrayList<RIGEdge>();
		edges.addAll(this.getEdges());
		Collections.sort(edges, new Comparator<RIGEdge>() {
			public int compare(RIGEdge e1,RIGEdge e2) {
				return Double.compare(e2.getWeight(), e1.getWeight());
			}
		});
		for(int i = 0; i < top; i++) {
			edges.get(i).setWeight(1);
		}
		for(int i = top; i < numEdges; i++) {
			this.removeEdge(edges.get(i));
		}
	}
	
	/**
	 * Returns true if this graph has at least one edge in ]0;1[.
	 * @return true if at least one edge is weighted
	 */
	public boolean hasWeightedEdges() {
		for(RIGEdge e: this.getEdges()) {
			if(e.getWeight() > 0 && e.getWeight() < 1) return true;
		}
		return false;
	}
	
	//TODO evaluatePrediction methods should be in ProtStructGraph. 
	//     But to be able to put them there we would need to pass here a Transformer that gets atom or residue serials depending if we are in AI or RI Graph 
	/**
	 * Evaluate this graph (assuming it is a prediction) against an original graph
	 * @param originalGraph
	 * @return
	 */
	public GraphComparisonResult evaluatePrediction(RIGraph originalGraph) {
		return evaluatePrediction(originalGraph, 1);
	}
	
	/**
	 * Evaluate this graph (assuming it is a prediction) against an original graph,
	 * considering only edges with sequence separation at least minSeqSep.
	 * @param originalGraph
	 * @param minSeqSep
	 * @return
	 */
	public GraphComparisonResult evaluatePrediction(RIGraph originalGraph, int minSeqSep) {
		
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
		GraphComparisonResult eval = new GraphComparisonResult(TruePos,FalsePos,TrueNeg,FalseNeg,0,predicted,original,cmtotal);
		return eval;
	}
	
	/**
	 * Write graph to given db, using our db graph format, 
	 * i.e. tables: chain_graph, single_model_graph, single_model_node, single_model_edge
	 * @param conn
	 * @param db
	 * @param weighted whether atom weights will be used or not
	 * @throws SQLException
	 */
	public void writeToDb(MySQLConnection conn, String db, boolean weighted) throws SQLException{
		//TODO we might want to move this to a graph i/o class
		//TODO Get rid of this and only keep fast one??
		
		HashMap<Integer,Integer> resser2nodeid = new HashMap<Integer,Integer>();
		// values we fix to constant 
		String CW = "1";
		String CR = "(true)";
		String EXPBB = "0";
		String ctStr = contactType;
		String weightedStr = "0";
		String directedStr = isDirected()?"1":"0";
		
		if (contactType.contains("_CAGLY")) {
			ctStr = contactType.replaceAll("_CAGLY", "");
		}
		if (ctStr.equals("ALL")) {
			ctStr = "BB+SC+BB/SC";
		}
		
		if (ContactType.isValidMultiAtomContactType(contactType, isDirected()) && weighted) {
			CW = ctStr;
			weightedStr = "1";
		}
		if (contactType.contains("_CAGLY") || contactType.contains("Cb")) {
			EXPBB = "-1";
		}
		if (minSeqSep != -1) {
			CR = "((i_cid!=j_cid)OR(abs(i_num-j_num)>="+minSeqSep+"))";
		} else if (interSSE) {
			CR = "((i_sstype!=j_sstype)OR(i_ssid!=j_ssid))";
		}
				
		int pgraphid=0;
		int graphid=0;
		String sql;		
		if (sid==null) {
			sql = "SELECT graph_id FROM "+db+".chain_graph " +
					" WHERE accession_code='"+pdbCode+"' AND pchain_code='"+chainCode+"'" +
					" AND model_serial = "+model+" AND dist = "+distCutoff+" AND expBB = "+EXPBB+ 
					" AND method = 'rc-cutoff';";
		} else {
			sql = "SELECT graph_id FROM "+db+".scop_graph " +
					" WHERE scop_id = '"+sid+"' "+
					" AND model_serial = "+model+" AND dist = "+distCutoff+" AND expBB = "+EXPBB+ 
					" AND method = 'rc-cutoff';";			
		}		
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		if (rsst.next()){	// if the pdbCode + chainCode were already in chain_graph then we take the graph_id as the pgraphid
			pgraphid = rsst.getInt(1);
		} else {			// no pdbCode + chainCode found, we insert them in chain_graph, thus assigning a new graph_id (pgraphid)
			// we are inserting same number for num_obs_res and num_nodes (the difference would be the non-standard aas, but we can't get that number from this object at the moment)
			String pdbChainCodeStr = "'"+pdbChainCode+"'";
			if (pdbChainCode.equals(PdbAsymUnit.NULL_CHAIN_CODE)) {
				pdbChainCodeStr="NULL";
			}
			if (sid==null) {
				sql = "INSERT INTO "+db+".chain_graph (accession_code,chain_pdb_code,pchain_code,model_serial,dist,expBB,method,num_res,num_obs_res,num_nodes,sses,date) " +
						"VALUES ('"+pdbCode+"', "+pdbChainCodeStr+",'"+chainCode+"', "+model+", "+distCutoff+", "+EXPBB+", 'rc-cutoff', "+getFullLength()+", "+getObsLength()+", "+getObsLength()+", "+((secondaryStructure!=null)?secondaryStructure.getNumElements():0)+", now())";
			} else {
				sql = "INSERT INTO "+db+".scop_graph (scop_id,accession_code,chain_pdb_code,pchain_code,model_serial,dist,expBB,method,num_res,num_obs_res,num_nodes,sses,date) " +
						"VALUES ('"+sid+"', '"+pdbCode+"', "+pdbChainCodeStr+",'"+chainCode+"', "+model+", "+distCutoff+", "+EXPBB+", 'rc-cutoff', "+getFullLength()+", "+getObsLength()+", "+getObsLength()+", "+((secondaryStructure!=null)?secondaryStructure.getNumElements():0)+", now())";				
			}
			Statement stmt2 = conn.createStatement();
			stmt2.executeUpdate(sql);
			// now we take the newly assigned graph_id as pgraphid
			sql = "SELECT LAST_INSERT_ID() FROM "+db+"."+((sid==null)?"chain":"scop")+"_graph LIMIT 1";
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
		sql = "SELECT single_model_id FROM "+singleModelsDb+".single_model WHERE "+
				" dist="+distCutoff+" AND expBB="+EXPBB+" AND CW='"+CW+"' AND CT='"+ctStr+"' AND CR='"+CR+ "' "+
				(singleModelsDb.equals(DEFAULT_SINGLEMODELS_DB)?"":("AND d="+directedStr))+";";
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			singlemodelid = rsst.getInt(1);
		}
		rsst.close();
		// and then insert to single_model_graph
		sql = "INSERT INTO "+db+".single_model_graph (pgraph_id,graph_type,accession_code,single_model_id,dist,expBB,CW,CT,CR,w,d,num_nodes,date) " +
				" VALUES ("+pgraphid+", '"+((sid==null)?"chain":"scop")+"', '"+pdbCode+"', "+singlemodelid+", "+distCutoff+", "+EXPBB+", '"+CW+"','"+ctStr+"', '"+CR+"', "+weightedStr+", "+directedStr+", "+getObsLength()+", now())";
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
		stmt = conn.createStatement();
		for (int resser:this.getSerials()) {			
			RIGNode node = this.getNodeFromSerial(resser);
			String res = Character.toString(AminoAcid.three2one(node.getResidueType()));
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
					" (graph_id, cid, num, res, "+
					" sstype, ssid, sheet_serial, turn, "+
					" k, k_in, k_out, "+
					" n, nwg, n_num) " +
					" VALUES ("+graphid+", '"+chainCode+"', "+resser+", '"+res+"', "+
					" "+secStructType+", "+secStructId+", "+sheetSerial+", "+turn+", "+
					(inDegree(node)+outDegree(node))+", "+inDegree(node)+", "+outDegree(node)+", "+
					"'"+nbh.getMotifNoGaps()+"', '"+nbh.getMotif()+"', '"+nbh.getCommaSeparatedResSerials()+"')";
			} else {		// we insert k (and no k_in or k_out)
				sql = "INSERT INTO "+db+".single_model_node "+
				" (graph_id, cid, num, res, "+
				" sstype, ssid, sheet_serial, turn, "+
				" k, n, nwg, n_num) " +
				" VALUES ("+graphid+", '"+chainCode+"', "+resser+", '"+res+"', "+
				" "+secStructType+", "+secStructId+", "+sheetSerial+", "+turn+", "+
				degree(node)+", '"+nbh.getMotifNoGaps()+"', '"+nbh.getMotif()+"', '"+nbh.getCommaSeparatedResSerials()+"')";
			}
			stmt.executeUpdate(sql);
		}
		
		// get the max node in db
		int maxNodeId = 0;
		sql = "SELECT LAST_INSERT_ID() FROM "+db+".single_model_node LIMIT 1;";
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			maxNodeId = rsst.getInt(1);
		}
		rsst.close();
		int nodeId = maxNodeId-this.getVertexCount();
		for (int resser:this.getSerials()) {
			nodeId++;
			resser2nodeid.put(resser, nodeId);			
		}
		
		// inserting edges
		// get the max weight
		double maxWeight = 0;
		if (weighted) {
			for (RIGEdge cont:getEdges()) {
				maxWeight = (maxWeight<cont.getAtomWeight())?cont.getAtomWeight():maxWeight;
			}
		} else {
			maxWeight = 1;
		}
		for (RIGEdge cont:getEdges()){
			RIGNode i_node = getEndpoints(cont).getFirst();
			RIGNode j_node = getEndpoints(cont).getSecond();
			char i_res = AminoAcid.three2one(i_node.getResidueType());
			char j_res = AminoAcid.three2one(j_node.getResidueType());

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
					" j_node_id, j_cid, j_num, j_res, j_sstype, j_ssid, j_sheet_serial, j_turn, weight, norm_weight, distance) " +
					" VALUES ("+graphid+", "+resser2nodeid.get(i_node.getResidueSerial())+", '"+chainCode+"', "+i_node.getResidueSerial()+", '"+i_res+"', "+i_secStructType+", "+i_secStructId+", "+i_sheetSerial+", "+i_turn+", "+
					resser2nodeid.get(j_node.getResidueSerial())+", '"+chainCode+"', "+j_node.getResidueSerial()+", '"+j_res+"', "+j_secStructType+", "+j_secStructId+", "+j_sheetSerial+", "+j_turn+", "+
					(weighted?cont.getAtomWeight():1)+", "+(weighted?(cont.getAtomWeight()/maxWeight):1)+", "+cont.getDistance()+")";
			stmt.executeUpdate(sql);
			if(!isDirected()) {// we want both side of the matrix in the table to follow Ioannis' convention
				// so we insert the reverse contact by swapping i, j in insertion
				sql = "INSERT INTO "+db+".single_model_edge "+
				" (graph_id, i_node_id, i_cid, i_num, i_res, i_sstype, i_ssid, i_sheet_serial, i_turn, "+
				" j_node_id, j_cid, j_num, j_res, j_sstype, j_ssid, j_sheet_serial, j_turn, weight, norm_weight, distance) " +
				" VALUES ("+graphid+", "+resser2nodeid.get(j_node.getResidueSerial())+", '"+chainCode+"', "+j_node.getResidueSerial()+", '"+j_res+"', "+j_secStructType+", "+j_secStructId+", "+j_sheetSerial+", "+j_turn+", "+
				resser2nodeid.get(i_node.getResidueSerial())+", '"+chainCode+"', "+i_node.getResidueSerial()+", '"+i_res+"', "+i_secStructType+", "+i_secStructId+", "+i_sheetSerial+", "+i_turn+", "+
				(weighted?cont.getAtomWeight():1)+", "+(weighted?(cont.getAtomWeight()/maxWeight):1)+", "+cont.getDistance()+")";
				stmt.executeUpdate(sql);
			}
		}
		
		stmt.close();

	}
	
	/**
	 * Writes graph to given db, using our db graph format,
	 * i.e. tables: chain_graph, single_model_graph, single_model_node, single_model_edge
	 * Atom weights are used. If not desired, use {@link #writeToDb(MySQLConnection, String, boolean)}
	 * @param conn
	 * @param db
	 * @throws SQLException
	 */
	public void writeToDb(MySQLConnection conn, String db) throws SQLException {
		writeToDb(conn, db, true);
	}
	
	/**
	 * Write graph to given db, using our db graph format, 
	 * i.e. tables: chain_graph, single_model_graph, single_model_node, single_model_edge
	 * Uses faster writing through files and LOAD DATA INFILE commands
	 * @param conn
	 * @param db
	 * @param weighted whether atom weights will be used or not
	 * @throws SQLException
	 * @throws IOException
	 */
	public void writeToDbFast(MySQLConnection conn, String db, boolean weighted) throws SQLException, IOException {
		//TODO we might want to move this to a graph i/o class
		
		HashMap<Integer,Integer> resser2nodeid = new HashMap<Integer,Integer>();
		
		// values we fix to constant 
		String CW = "1";
		String CR = "(true)";
		String EXPBB = "0";
		String ctStr = contactType;
		String weightedStr = "0";
		String directedStr = isDirected()?"1":"0";
		
		if (contactType.contains("_CAGLY")) {
			ctStr = contactType.replaceAll("_CAGLY", "");
		}
		if (ctStr.equals("ALL")) {
			ctStr = "BB+SC+BB/SC";
		}
		if (ContactType.isValidMultiAtomContactType(contactType, isDirected()) && weighted) {
			CW = ctStr;
			weightedStr = "1";
		}
		if (contactType.contains("_CAGLY") || contactType.contains("Cb")) {
			EXPBB = "-1";
		}
		if (minSeqSep != -1) {
			CR = "((i_cid!=j_cid)OR(abs(i_num-j_num)>="+minSeqSep+"))";
		} else if (interSSE) {
			CR = "((i_sstype!=j_sstype)OR(i_ssid!=j_ssid))";
		}
				
		int pgraphid=0;
		int graphid=0;
		String sql;		
		if (sid==null) {
			sql = "SELECT graph_id FROM "+db+".chain_graph " +
					" WHERE accession_code='"+pdbCode+"' AND pchain_code='"+chainCode+"'" +
					" AND model_serial = "+model+" AND dist = "+distCutoff+" AND expBB = "+EXPBB+ 
					" AND method = 'rc-cutoff';";
		} else {
			sql = "SELECT graph_id FROM "+db+".scop_graph " +
					" WHERE scop_id = '"+sid+"' "+
					" AND model_serial = "+model+" AND dist = "+distCutoff+" AND expBB = "+EXPBB+ 
					" AND method = 'rc-cutoff';";			
		}		
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		if (rsst.next()){	// if the pdbCode + chainCode were already in chain_graph then we take the graph_id as the pgraphid
			pgraphid = rsst.getInt(1);
		} else {			// no pdbCode + chainCode found, we insert them in chain_graph, thus assigning a new graph_id (pgraphid)
			// we are inserting same number for num_obs_res and num_nodes (the difference would be the non-standard aas, but we can't get that number from this object at the moment)
			String pdbChainCodeStr = "'"+pdbChainCode+"'";
			if (pdbChainCode.equals(PdbAsymUnit.NULL_CHAIN_CODE)) {
				pdbChainCodeStr="NULL";
			}
			if (sid==null) {
				sql = "INSERT INTO "+db+".chain_graph (accession_code,chain_pdb_code,pchain_code,model_serial,dist,expBB,method,num_res,num_obs_res,num_nodes,sses,date) " +
						"VALUES ('"+pdbCode+"', "+pdbChainCodeStr+",'"+chainCode+"', "+model+", "+distCutoff+", "+EXPBB+", 'rc-cutoff', "+getFullLength()+", "+getObsLength()+", "+getObsLength()+", "+((secondaryStructure!=null)?secondaryStructure.getNumElements():0)+", now())";
			} else {
				sql = "INSERT INTO "+db+".scop_graph (scop_id,accession_code,chain_pdb_code,pchain_code,model_serial,dist,expBB,method,num_res,num_obs_res,num_nodes,sses,date) " +
						"VALUES ('"+sid+"', '"+pdbCode+"', "+pdbChainCodeStr+",'"+chainCode+"', "+model+", "+distCutoff+", "+EXPBB+", 'rc-cutoff', "+getFullLength()+", "+getObsLength()+", "+getObsLength()+", "+((secondaryStructure!=null)?secondaryStructure.getNumElements():0)+", now())";				
			}
			Statement stmt2 = conn.createStatement();
			stmt2.executeUpdate(sql);
			// now we take the newly assigned graph_id as pgraphid
			sql = "SELECT LAST_INSERT_ID() FROM "+db+"."+((sid==null)?"chain":"scop")+"_graph LIMIT 1";
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
		sql = "SELECT single_model_id FROM "+singleModelsDb+".single_model WHERE "+
				" dist="+distCutoff+" AND expBB="+EXPBB+" AND CW='"+CW+"' AND CT='"+ctStr+"' AND CR='"+CR+ "' "+
				(singleModelsDb.equals(DEFAULT_SINGLEMODELS_DB)?"":("AND d="+directedStr))+";";
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			singlemodelid = rsst.getInt(1);
		}
		rsst.close();
		// and then insert to single_model_graph
		sql = "INSERT INTO "+db+".single_model_graph (pgraph_id,graph_type,accession_code,single_model_id,dist,expBB,CW,CT,CR,w,d,num_nodes,date) " +
				" VALUES ("+pgraphid+", '"+((sid==null)?"chain":"scop")+"', '"+pdbCode+"', "+singlemodelid+", "+distCutoff+", "+EXPBB+", '"+CW+"','"+ctStr+"', '"+CR+"', "+weightedStr+", "+directedStr+", "+getObsLength()+", now())";
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
		for (int resser:this.getSerials()) {			
			RIGNode node = this.getNodeFromSerial(resser);
			String res = Character.toString(AminoAcid.three2one(node.getResidueType()));
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
				nodesOut.println(graphid+"\t"+chainCode+"\t"+resser+"\t"+res+"\t"+
					secStructType+"\t"+secStructId+"\t"+sheetSerial+"\t"+turn+"\t"+
					(inDegree(node)+outDegree(node))+"\t"+inDegree(node)+"\t"+outDegree(node)+"\t"+
					nbh.getMotifNoGaps()+"\t"+nbh.getMotif()+"\t"+nbh.getCommaSeparatedResSerials());
			} else {		// we insert k (and no k_in or k_out)
				nodesOut.println(graphid+"\t"+chainCode+"\t"+resser+"\t"+res+"\t"+
						secStructType+"\t"+secStructId+"\t"+sheetSerial+"\t"+turn+"\t"+
						degree(node)+"\t"+"\\N"+"\t"+"\\N"+"\t"+
						nbh.getMotifNoGaps()+"\t"+nbh.getMotif()+"\t"+nbh.getCommaSeparatedResSerials());
			}
		}
		nodesOut.close();
		stmt = conn.createStatement();
		sql = "LOAD DATA LOCAL INFILE '"+graphid+"_nodes.txt' INTO TABLE "+db+".single_model_node "+
			" (graph_id, cid, num, res, "+
			" sstype, ssid, sheet_serial, turn, "+
			" k, k_in, k_out, n, nwg, n_num);";
		stmt.executeUpdate(sql);
		stmt.close();
		File fileToDelete = new File(graphid+"_nodes.txt");
		if (fileToDelete.exists()) {
			fileToDelete.delete();
		}
		
		// get the first inserted node id
		int firstNodeId = 0;
		sql = "SELECT LAST_INSERT_ID() FROM "+db+".single_model_node LIMIT 1;";
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		if (rsst.next()){
			firstNodeId = rsst.getInt(1);
		}
		rsst.close();
		stmt.close();
		int nodeId = firstNodeId;
		for (int resser:this.getSerials()) {
			resser2nodeid.put(resser, nodeId);	
			nodeId++;
		}
		
		// inserting edges
		PrintStream edgesOut = new PrintStream(new FileOutputStream(graphid+"_edges.txt"));
		// get the max weight
		double maxWeight = 0;
		if (weighted) {
			for (RIGEdge cont:getEdges()) {
				maxWeight = (maxWeight<cont.getAtomWeight())?cont.getAtomWeight():maxWeight;
			}
		} else {
			maxWeight = 1;			
		}
		for (RIGEdge cont:getEdges()){
			RIGNode i_node = getEndpoints(cont).getFirst();
			RIGNode j_node = getEndpoints(cont).getSecond();
			char i_res = AminoAcid.three2one(i_node.getResidueType());
			char j_res = AminoAcid.three2one(j_node.getResidueType());
			
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
			
			edgesOut.println(graphid+"\t"+resser2nodeid.get(i_node.getResidueSerial())+"\t"+chainCode+"\t"+i_node.getResidueSerial()+"\t"+i_res+"\t"+i_secStructType+"\t"+i_secStructId+"\t"+i_sheetSerial+"\t"+i_turn+"\t"+
					resser2nodeid.get(j_node.getResidueSerial())+"\t"+chainCode+"\t"+j_node.getResidueSerial()+"\t"+j_res+"\t"+j_secStructType+"\t"+j_secStructId+"\t"+j_sheetSerial+"\t"+j_turn+"\t"+
					(weighted?cont.getAtomWeight():1)+"\t"+(weighted?(cont.getAtomWeight()/maxWeight):1)+"\t"+cont.getDistance());
			if(!isDirected()) {// we want both side of the matrix in the table to follow Ioannis' convention
				// so we insert the reverse contact by swapping i, j in insertion
				edgesOut.println(graphid+"\t"+resser2nodeid.get(j_node.getResidueSerial())+"\t"+chainCode+"\t"+j_node.getResidueSerial()+"\t"+j_res+"\t"+j_secStructType+"\t"+j_secStructId+"\t"+j_sheetSerial+"\t"+j_turn+"\t"+
						resser2nodeid.get(i_node.getResidueSerial())+"\t"+chainCode+"\t"+i_node.getResidueSerial()+"\t"+i_res+"\t"+i_secStructType+"\t"+i_secStructId+"\t"+i_sheetSerial+"\t"+i_turn+"\t"+
						(weighted?cont.getAtomWeight():1)+"\t"+(weighted?(cont.getAtomWeight()/maxWeight):1)+"\t"+cont.getDistance());
			}			
		}
		edgesOut.close();
		sql = "LOAD DATA LOCAL INFILE '"+graphid+"_edges.txt' INTO TABLE "+db+".single_model_edge "+
			" (graph_id, i_node_id, i_cid, i_num, i_res, i_sstype, i_ssid, i_sheet_serial, i_turn, "+
			" j_node_id, j_cid, j_num, j_res, j_sstype, j_ssid, j_sheet_serial, j_turn, weight, norm_weight, distance);";
		stmt = conn.createStatement();
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
	 * Write graph to given db, using our db graph format, 
	 * i.e. tables: chain_graph, single_model_graph, single_model_node, single_model_edge
	 * Uses faster writing through files and LOAD DATA INFILE commands
	 * Atom weights will be used.
	 * @param conn
	 * @param db
	 * @throws SQLException
	 * @throws IOException
	 */
	public void writeToDbFast(MySQLConnection conn, String db) throws SQLException, IOException {
		writeToDb(conn, db, true);
	}
	
	/**
	 * Write graph to given outfile in CMView format
	 * @param outfile
	 * @throws IOException
	 */
	public void writeToFile (String outfile) throws IOException {
		//TODO we might want to move this to a graph i/o class
		
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		Out.println("#CMVIEW GRAPH FILE ver: "+GRAPHFILEFORMATVERSION);
		Out.println("#SEQUENCE: "+sequence);
		Out.println("#PDB: "+(pdbCode==null?"":pdbCode));
		Out.println("#PDB CHAIN CODE: "+(pdbChainCode==null?"":pdbChainCode));
		Out.println("#CHAIN: "+(chainCode==null?"":chainCode));
		Out.println("#MODEL: "+(model==0?PdbAsymUnit.DEFAULT_MODEL:model));
		Out.println("#CT: "+contactType);
		Out.println("#CUTOFF: "+distCutoff);
		
		// we use a temp TreeMap to be able to order the output
		TreeMap<Pair<Integer>,Double> pairs = new TreeMap<Pair<Integer>,Double>(new IntPairComparator());
		for (RIGEdge cont:getEdges()){
			Pair<RIGNode> pair = getEndpoints(cont);
			int i_resser=pair.getFirst().getResidueSerial();
			int j_resser=pair.getSecond().getResidueSerial();
			//BEWARE!! here we write weights while in writeToDb we write atomWeights (consistent with what we do in FileRIGraph) TODO do we want this behaviour?
			double weight=cont.getWeight();
			pairs.put(new Pair<Integer>(i_resser,j_resser),weight);
			
		}
		for (Pair<Integer> pair:pairs.keySet()) { 
			Out.printf(Locale.US,pair.getFirst()+"\t"+pair.getSecond()+"\t%6.3f\n",pairs.get(pair));
		}
		Out.close();		
	}

	/**
	 * Write graph to given outfile in format appropriate for the mfinder program
	 * @param outfile
	 * @throws IOException
	 */
	public void writeToMotifFile (String outfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		
		HashMap<Integer,Integer> resser2nodeser = new HashMap<Integer,Integer>();
		int nodeSer = 0;
		for (int resser:this.getSerials()) {
			RIGNode node = this.getNodeFromSerial(resser);
			if (outDegree(node) > 0) {
				nodeSer++;
				resser2nodeser.put(resser, nodeSer);
			}
		}
		
		// we use a temp TreeMap to be able to order the output
		TreeSet<Pair<Integer>> pairs = new TreeSet<Pair<Integer>>(new IntPairComparator());
		for (RIGEdge cont:getEdges()){
			Pair<RIGNode> pair = getEndpoints(cont);
			int i_resser=resser2nodeser.get(pair.getFirst().getResidueSerial());
			int j_resser=resser2nodeser.get(pair.getSecond().getResidueSerial());	
			pairs.add(new Pair<Integer>(i_resser,j_resser));
		}
		for (Pair<Integer> pair:pairs) { 
			Out.printf(Locale.US,pair.getFirst()+"\t"+pair.getSecond()+"\t1\n");
		}
		Out.close();		
	}
	
	/**
	 * Writes the graph to given outfile in "paul" format.
	 * Possibly this is compatible with all other LiSA programs. Here we
	 * could only test it with paul and drawCM  
	 * See https://www.mi.fu-berlin.de/w/LiSA/LiSAlibrary
	 * @param outfile
	 * @throws IOException
	 */
	public void writeToPaulFile(String outfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		
		// the order of the first 2 lines is:
		// for drawCM 1st number of contacts, 2nd number of residues
		// for paul 1st number of residues, 2nd number of contacts
		Out.println(getFullLength());// number of residues
		Out.println(getEdgeCount()); // number of contacts		
		
		TreeMap<Pair<Integer>,Double> pairs = new TreeMap<Pair<Integer>,Double>(new IntPairComparator());
		for (RIGEdge cont:getEdges()){
			Pair<RIGNode> pair = getEndpoints(cont);
			int i_resser=pair.getFirst().getResidueSerial();
			int j_resser=pair.getSecond().getResidueSerial();
			double weight=cont.getWeight();
			pairs.put(new Pair<Integer>(i_resser,j_resser),weight);			
		}
		for (Pair<Integer> pair:pairs.keySet()) { 
			// The paul format has 4 columns: the first 2 are the end points of the edge,
			// 3rd column is weights, 4th column is contact class which we set to 1
			// In paul the numbering of residues starts at 0: that's why we subtract 1 here
			Out.printf(Locale.US,(pair.getFirst()-1)+" "+(pair.getSecond()-1)+" %6.3f 1\n",pairs.get(pair));
		}

		Out.close();
	}
	
	/**
	 * Export graph as a Casp contact prediction file.
	 * Note: Writes weighted edges to file.
	 * @param outFile name of the output file
	 * @throws IOException
	 */
	public void writeToCaspRRFile(String outFile) throws IOException {
		PrintStream out = new PrintStream(new FileOutputStream(outFile));
		CaspRRFileData rrData = new CaspRRFileData();
		rrData.setSequence(this.sequence);
		rrData.setTargetNum(this.targetNum);
		rrData.setModelNum(this.caspModelNum);
		rrData.setAuthor(this.authorStr);
		rrData.setMethod(this.methodStr);
		
		// we use a temp TreeMap to be able to order the output
		TreeMap<Pair<Integer>,CaspRRFileData.RRContact> pairs = new TreeMap<Pair<Integer>,CaspRRFileData.RRContact>(new IntPairComparator());
		for (RIGEdge cont:getEdges()){
			Pair<RIGNode> nodePair = getEndpoints(cont);
			int i = nodePair.getFirst().getResidueSerial();
			int j = nodePair.getSecond().getResidueSerial();
			double minDist = CaspRRFileData.DEFAULT_MIN_DIST;
			double maxDist = this.getCutoff();
			double weight = cont.getWeight();
			CaspRRFileData.RRContact rrCont = rrData.new RRContact(i,j, minDist, maxDist, weight);
			pairs.put(new Pair<Integer>(i, j),rrCont);
		}
		
		for (CaspRRFileData.RRContact rrCont:pairs.values()) {
			rrData.addContact(rrCont);
		}
		rrData.writeToStream(out);
		out.close();
	}
	
	/**
	 * Write the graph to a file comptabile with the SADP program by Johannes Jain.
	 * See: Brijnesh J. Jain, Michael Lappe: Joining Softassign and Dynamic Programming for the Contact Map Overlap Problem. BIRD 2007: 410-423
	 * @param outfile
	 * @throws IOException
	 */
	public void writeToSADPFile (String outfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		Out.println(this.getFullLength());
		for (RIGEdge cont:getEdges()){
			Pair<RIGNode> pair = getEndpoints(cont);
			int i_resser=pair.getFirst().getResidueSerial() - 1;
			int j_resser=pair.getSecond().getResidueSerial() - 1 ;
			//BEWARE!! here we write weights while in writeToDb we write atomWeights (consistent with what we do in FileRIGraph) TODO do we want this behaviour?
			double weight=cont.getWeight();
			Out.printf(Locale.US,i_resser+"\t"+j_resser+"\t%1.0f\t1\n",weight);
		}
		Out.close();		
	}

	
	/**
	 * Write graph to given outfile in network(Ioannis) format
	 * @param graphId
	 * @param dirName
	 * @throws IOException
	 */
	//TODO we might want to move this to a graph i/o class
	//TODO refactor 
	public void writeUndirUnweightGraphToNetworkFiles (int graphId, String dirName) throws IOException {
		if (isDirected()) {
			System.err.println("This method is only for undirected graphs!");
			return;
		}
		
		String filePrefix = dirName + "/" + String.valueOf(graphId)+"_"+pdbCode+"_"+chainCode+"_"+contactType.replaceAll("/", ".")+"_"+String.valueOf(distCutoff).replaceAll("\\.", "_")+"_";
		PrintStream Out = new PrintStream(new FileOutputStream(filePrefix+"edges.txt"));
		for (RIGEdge cont:getEdges()){
			Pair<RIGNode> pair = getEndpoints(cont);
			int i_resser=pair.getFirst().getResidueSerial();
			int j_resser=pair.getSecond().getResidueSerial();
			if (i_resser < j_resser) {
				Out.printf(Locale.US,i_resser+"\t"+j_resser+"\t"+1+"\t"+"1.000"+"\n");
			}
		}
		Out.close();
		Out = new PrintStream(new FileOutputStream(filePrefix+"nodes.txt"));
		for (int resser:this.getSerials()) {
			RIGNode node = this.getNodeFromSerial(resser);
			String res = Character.toString(AminoAcid.three2one(node.getResidueType()));
			Out.printf(Locale.US,resser+"\t"+chainCode+"\t"+resser+"\t"+res+"\t"+degree(node)+"\t"+degree(node)+"\n");
			
		}
		Out.close();
	}
	
	/**
	 * Compares this RIGraph to given RIGraph returning 3 graphs: common, onlythis, onlyother in a HashMap
	 * @param other
	 * @return
	 * @throws Exception
	 * See also: {@link #evaluatePrediction(RIGraph)}
	 */
	//TODO not sure what kind of return we want, for now is a HashMap with three graph objects 
	public HashMap<String,RIGraph> compare(RIGraph other) throws Exception{
		//first check that other has same sequence than this, otherwise throw exception
		if (this.getFullLength()!=other.getFullLength()) {
			//TODO throw specific exception
			throw new Exception("Sequence of 2 graphs to compare differ, can't compare them.");
		} //else {
//			for (int resser:this.serials2nodes.keySet()){
//				String this_res = getNodeFromSerial(resser).getResidueType();
//				RIGNode otherNode = other.getNodeFromSerial(resser);
//				String other_res = AAinfo.NONSTANDARD_AA_THREE_LETTER;
//				if (otherNode!=null) {
//					other_res = otherNode.getResidueType();
//				}
//				if (!this_res.equals(AAinfo.NONSTANDARD_AA_THREE_LETTER) && !other_res.equals(AAinfo.NONSTANDARD_AA_THREE_LETTER) && !this_res.equals(other_res)) {
//					//TODO throw specific exception
//					throw new Exception("Sequence of 2 graphs to compare differ, can't compare them.");
//				}
//			}
//		}
		//NOTE: the common graph will have same node/edge properties as this graph, 
		//      which doesn't make a lot of sense, but anyway one has to choose between this or other, 
		//      or otherwise make some kind of merge, e.g. merge the weights by averaging? 
		RIGraph commongraph = this.copy(); 
		RIGraph onlythisgraph = this.copy();
		RIGraph onlyothergraph = other.copy();

		for (RIGEdge cont:this.getEdges()){
			Pair<RIGNode> pair = this.getEndpoints(cont);
			int i_resser = pair.getFirst().getResidueSerial();
			int j_resser = pair.getSecond().getResidueSerial();
			if (other.findEdge(other.getNodeFromSerial(i_resser), other.getNodeFromSerial(j_resser))!=null) {
				onlythisgraph.removeEdge(onlythisgraph.findEdge(onlythisgraph.getNodeFromSerial(i_resser),onlythisgraph.getNodeFromSerial(j_resser)));
				onlyothergraph.removeEdge(onlyothergraph.findEdge(onlyothergraph.getNodeFromSerial(i_resser),onlyothergraph.getNodeFromSerial(j_resser)));
			} else {
				commongraph.removeEdge(commongraph.findEdge(commongraph.getNodeFromSerial(i_resser),commongraph.getNodeFromSerial(j_resser)));
			}
		}

		HashMap<String,RIGraph> result = new HashMap<String,RIGraph>();
		result.put("common", commongraph);
		result.put("onlythis", onlythisgraph);
		result.put("onlyother",onlyothergraph);
		return result;
	}
	
	/**
	 * Compares this RIGraph to given RIGraph and returns the number of common edges 
	 * @param other
	 * @return
	 */
	public int getCommonEdgesCount(RIGraph other) {
		
		int commonEdges = 0;
		
		if ((this.getEdgeCount() == 0) || (other.getEdgeCount() == 0)) {
			return commonEdges;
		}
		
		for (RIGEdge cont:this.getEdges()){
			Pair<RIGNode> pair = this.getEndpoints(cont);
			int i_resser = pair.getFirst().getResidueSerial();
			int j_resser = pair.getSecond().getResidueSerial();
			if (other.findEdge(other.getNodeFromSerial(i_resser), other.getNodeFromSerial(j_resser))!=null) {
				commonEdges++;				
			}
		}

		return commonEdges;
	}
	
	/**
	 * Returns a new RIGraph copy (deep) of this one
	 * @return
	 */
	public RIGraph copy() {
		RIGraph newGraph = new RIGraph();
		newGraph.setPdbCode(pdbCode);
		newGraph.setPdbChainCode(pdbChainCode);
		newGraph.setChainCode(chainCode);
		newGraph.setModel(model);
		newGraph.setTargetNum(targetNum);
		newGraph.setGroupNum(groupNum);
		newGraph.setCaspModelNum(caspModelNum);
		newGraph.setContactType(contactType);
		newGraph.setCutoff(distCutoff);
		newGraph.setSequence(sequence);
		newGraph.setSingleModelsDb(singleModelsDb);
		
		// copying nodes
		for (RIGNode node:this.getVertices()) {
			RIGNode newNode = node.copy();
			newGraph.addVertex(newNode);
		}
		
		// copying edges
		for (RIGEdge edge:this.getEdges()) {
			Pair<RIGNode> pair = this.getEndpoints(edge);
			int i_resser = pair.getFirst().getResidueSerial();
			int j_resser = pair.getSecond().getResidueSerial();
			// EdgeType enum should copy correctly because enums are treated as ints in copying (always deep copied)
			newGraph.addEdge(edge.copy(), newGraph.getNodeFromSerial(i_resser), newGraph.getNodeFromSerial(j_resser), this.getEdgeType(edge));
		}
		
		// copying the SecondaryStructure object by retrieving all references from the new nodes
		if (this.getSecondaryStructure()!=null) {
			SecondaryStructure secStruct = new SecondaryStructure(this.getSecondaryStructure().getSequence());
			for (RIGNode node:newGraph.getVertices()) {
				SecStrucElement sselem = node.getSecStrucElement();
				if (sselem!=null && !secStruct.contains(sselem)) {
					secStruct.add(sselem);
				}
			}
			newGraph.setSecondaryStructure(secStruct);
		}
		
		return newGraph;
	}
	
	/**
	 * Gets the complement graph of this RIGraph returning it
	 * @return the complement graph of this RIGraph
	 */
	public RIGraph getComplement() {
		RIGraph complGraph = this.copy();
		complGraph.removeAllEdges();
		
		if (!isDirected()) {
			for (int i = 1; i < fullLength; i++) {
				for (int j = (i+1); j <= fullLength; j++) {
					if (this.getEdgeFromSerials(i, j) == null) {
						complGraph.addEdgeIJ(i, j);
					}
				}
			}
		} else {
			for (int i = 1; i <= fullLength; i++) {
				for (int j = 1; j <= fullLength; j++) {
					if (i == j) {
						continue;
					}
					if (this.getEdgeFromSerials(i, j) == null) {
						complGraph.addEdgeIJ(i, j);
					}
				}
			}			
		}
		
		return complGraph;
	}
	
	/**
	 * Turns the RIGraph into the GMatrix 
	 * @return the complement graph of this RIGraph
	 */
	public void convert2GMatrix( GMatrix A) { // Attn : Matrix starts from 0,0 as first entry 
		A.setSize( fullLength, fullLength); 
		A.setIdentity(); 	
		if (!isDirected()) { // undirected Graph 
			for (int i = 1; i < fullLength; i++) {
				for (int j = (i+1); j <= fullLength; j++) {
					if (this.getEdgeFromSerials(i, j) != null) {
						A.setElement( i-1,j-1, 1.0);
						A.setElement( j-1,i-1, 1.0);
					} // end if i,j / j,i is an edge  
				} // next j 
			}	// next i 	
		} else { // directed Graph 
			for (int i = 1; i <= fullLength; i++) {
				for (int j = 1; j <= fullLength; j++) {
					if (i == j) {
						continue;
					}
					if (this.getEdgeFromSerials(i, j) != null) {
						A.setElement( i-1,j-1, 1.0);  
					} // end if i,j is an edge  
				} // next j 
			}	// next i 		
		} // end if directed Graph 
	} // end convert2GMatrix

	
	/**
	 * Removes a vertex from this graph.
	 * This overridden function also updates the serials2nodes map.
	 * @return true if vertex is present in this graph and thus can be removed, false if vertex is not present in this graph
	 */
	@Override
	public boolean removeVertex(RIGNode vertex) {
		if (serials2nodes.containsKey(vertex.getResidueSerial())) {
			serials2nodes.remove(vertex.getResidueSerial());
			return super.removeVertex(vertex);
		} else {
			return false;
		}
	}
	
	/**
	 * Adds a vertex to this graph.
	 * This overridden function also updates the serials2nodes map.
	 * @return true if add is successful (node not present in graph so can be added), false otherwise (node already in graph so not added)
	 */
	@Override
	public boolean addVertex(RIGNode vertex) {
		if (!serials2nodes.containsKey(vertex.getResidueSerial())) {
			serials2nodes.put(vertex.getResidueSerial(),vertex);
			return super.addVertex(vertex);
		} else {
			return false;
		}
	}
	
	/**
	 * Sets the serials2nodes Map from the RIGNodes objects present in this RIGraph
	 * This should be used in cases where a RIGraph has been created using some external 
	 * classes (e.g. JUNG algorithms). In those cases a generic addVertex is called instead of our
	 * overridden addVertex to create the nodes, thus serials2nodes will not have been set correctly.
	 * So this method needs to be called to set the mapping of serials to nodes.  
	 */
	protected void setSerials2NodesMap() {
		this.serials2nodes = new TreeMap<Integer, RIGNode>();
		for (RIGNode node:this.getVertices()) {
			serials2nodes.put(node.getResidueSerial(),node);
		}
	}
	
	/**
	 * Returns a new RIGraph that contains a subset of edges of the given size of this RIGraph.  
	 * @param toSample the fraction of edges to sample  
	 * @return
	 */
	public RIGraph getRandomSubset(double toSample) {
		int numSampledContacts = (int) ((double)(this.getEdgeCount())*toSample);
		return getRandomSubset(numSampledContacts);
	}

	/**
	 * Returns a new RIGraph that contains a subset of edges of the given size of this RIGraph.
	 * @param numSampledContacts the number of edges to sample
	 * @return
	 */
	private RIGraph getRandomSubset(int numSampledContacts) {
		Random rand = new Random();
		
		IntPairSet allcmpairs = new IntPairSet();
		for (RIGEdge edge:getEdges()) {
			Pair<RIGNode> nodePair = this.getEndpoints(edge);
			Pair<Integer> pair = new Pair<Integer>(nodePair.getFirst().getResidueSerial(),nodePair.getSecond().getResidueSerial());
			allcmpairs.add(pair);
		}
		
		TreeSet<Integer> sampledIndices = new TreeSet<Integer>();
		for (int i=0;i<numSampledContacts;i++) {
			// the TreeSet takes care of not having duplicates indices
			// add() returns false if the value is a duplicate, so if it returns false 
			// we loop until the return value is true, i.e. we have a new index
			while (!sampledIndices.add(rand.nextInt(allcmpairs.size())));
		}
		
		RIGraph subsetGraph = this.copy();
		subsetGraph.removeAllEdges();
		
		// the iterator returns the indices in order since it is from a TreeSet
		Iterator<Integer> indexIterator = sampledIndices.iterator();
		int index = 0;
		int sampledIndex = indexIterator.next();
		for (Pair<Integer> pair:allcmpairs) {
			if (index==sampledIndex) {
				subsetGraph.addEdgeIJ(pair.getFirst(), pair.getSecond());
				if (indexIterator.hasNext()) {
					sampledIndex = indexIterator.next();
				} else {
					break;
				}
			}
			index++;
		}
		if (subsetGraph.getEdgeCount()!=numSampledContacts) {
			System.err.println("Error: sampled number of edges does not coincide with required number to sample. Please report the bug.");
		}
		return subsetGraph;
	}

	/**
	 * Returns a new RIGraph that contains this RIGraph's edges plus a number of random 
	 * added edges 
	 * @param noiseToAdd the fraction of edges to add (with respect to the total number of 
	 * edges in this RIGraph)
	 * @return
	 */
	public RIGraph getRandomNoiseMap(double noiseToAdd) {
		int numContactsToAdd = (int) ((double)(this.getEdgeCount())*noiseToAdd);
		return getRandomNoiseMap(numContactsToAdd);
	}
	
	/**
	 * Returns a new RIGraph that contains this RIGraph's edges plus a number of random 
	 * added edges
	 * @param numContactsToAdd the number of edges to add
	 * @return
	 */
	private RIGraph getRandomNoiseMap(int numContactsToAdd) {
		Random rand = new Random();
		RIGraph noiseGraph = this.copy();
		int n = 0;
		while (n<numContactsToAdd) {
			int i = rand.nextInt(this.getVertexCount())+1;
			int j = rand.nextInt(this.getVertexCount())+1;
			if (i==j) continue;
			if (i>j) {
				int t = i;
				i = j;
				j = t;
			}
			if (!noiseGraph.containsEdgeIJ(i, j) && 
					noiseGraph.getNodeFromSerial(i)!=null && // we have to check that the sampled i,j doesn't fall on unobserved residues  
					noiseGraph.getNodeFromSerial(j)!=null) {
				noiseGraph.addEdgeIJ(i, j);
				n++;
				continue;
			}
		}
		
		if (noiseGraph.getEdgeCount()!=this.getEdgeCount()+numContactsToAdd) {
			System.err.println("Error: number of added noise edges do not coincide with required number to add. Please report the bug.");
		}
		return noiseGraph;
	}
}

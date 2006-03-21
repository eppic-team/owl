package tools;

import java.sql.*;
import java.io.*;
import java.util.*;

public class Graph2Pml {

    private Connection conn;
    private PrintWriter out = null;
    private PyMol pml = null;
    private boolean msdsd = true, cgoEdge = false, directed = true;
    
    private HashMap<String, String> chainNodes = new HashMap<String, String>();
    private HashMap<String, String> chainEdges = new HashMap<String, String>();
    private HashMap<Long, String> nodeSetIds = new HashMap<Long, String> ();
    private HashMap<String, String> nodeSets = new HashMap<String, String> ();
    private HashMap<String, String> nodeSet = new HashMap<String, String> ();
    private HashMap<Long, String> edgeSetIds = new HashMap<Long, String> ();
    private HashMap<String, String> edgeSets = new HashMap<String, String> ();
    private HashMap<String, String> edgeSet = new HashMap<String, String> ();

    private String[] dbInfo = new String[] {"newmsdgraph", "list", "nodes", "edges", "graph_id"};
    private String[] nodeInfo = new String[] {"cid", "num"};
    private String[] edgeInfo = new String[] {"i_cid", "j_cid", "i_num", "j_num"};
    private String[] graphInfo = new String[6];
    private String molObjName = null, nodeGraphSel = "", edgeGraphSel = "";
    
    private boolean[] draw = new boolean[] {true, true, false, false};
    private String nodeColorMethod = "chain", nodeSizeMethod = "uniform", edgeColorMethod = "uniform", edgeSizeMethod = "uniform", nodeColor = "blue", edgeColor = "orange", specialResColor = "purple";
    private boolean nodeColDiscr = false, edgeColDiscr = false, nodeSizeRev = false, edgeSizeRev = false;
    private double nodeSize = 0.46, edgeSize = 3.65;
    private double[] nodeSizeRange = {0.18, 0.75};
    private double[] edgeSizeRange = {0.3, 7.0};
    private String edgeGapCondition = "", nodeTranspCondition = "", specialResQuery = null;
    private double edgeGap = 0.25, nodeTransp = 0.6, surfTransp = 0.6;
    private String backgroundColor = "white";
			  
    private String[] colors = {"light_grey", "green", "red", "blue", "yellow", "violet", "cyan", "salmon", "lime", "pink", "slate", "magenta", "orange", "marine", "olive", "purple", "teal", "forest", "firebrick", "chocolate", "wheat", "white", "grey"};

    public Graph2Pml(PrintWriter out, String molObjName, int graphId, String edgeWeight, String nodeContactType, String edgeContactType, String contactRange, String userFilter, boolean msdsd, boolean cgoEdge, boolean directed, Connection conn) {
	
		graphInfo[0] = String.valueOf(graphId);
		graphInfo[1] = edgeWeight;
		graphInfo[2] = nodeContactType;
        graphInfo[3] = edgeContactType;
		graphInfo[4] = contactRange;
		graphInfo[5] = userFilter;
		nodeGraphSel = "("+dbInfo[4]+" = "+graphInfo[0]+") AND ("+graphInfo[2]+" > 0) AND ("+graphInfo[4]+")";
		edgeGraphSel = "("+dbInfo[4]+" = "+graphInfo[0]+") AND ("+graphInfo[3]+" > 0) AND ("+graphInfo[4]+")";
		
		this.conn = conn;
		this.out = out;
		this.pml = new PyMol(out);
		this.molObjName = molObjName;
		this.msdsd = msdsd;
		this.cgoEdge = cgoEdge;
		this.directed = directed;

    }
    
    public boolean draw(boolean nodes, boolean edges, boolean specialRes, boolean surface) {

		boolean status = true;
	
		if (!nodes && !edges) {
		    System.out.println("You should at least draw nodes or edges!!!");
		    status = false;
		} else {
		    draw[0] = nodes;
		    draw[1] = edges;
		    draw[2] = specialRes;
		    draw[3] = surface;
		    status = true;
		}
	
		return status;

    }

    public void setSchema(String dbName, String graphTbl, String nodesTbl, String edgesTbl, String graphIdCol, String nodeChainCol, String nodeNumCol, String edgeINodeChainCol, String edgeJNodeChainCol, String edgeINodeNumCol, String edgeJNodeNumCol) {

		dbInfo[0] = dbName;
		dbInfo[1] = graphTbl;
		dbInfo[2] = nodesTbl;
		dbInfo[3] = edgesTbl;
		dbInfo[4] = graphIdCol;
		nodeGraphSel = "("+dbInfo[4]+" = "+graphInfo[0]+") AND ("+graphInfo[2]+" > 0) AND ("+graphInfo[4]+")";
		edgeGraphSel = "("+dbInfo[4]+" = "+graphInfo[0]+") AND ("+graphInfo[3]+" > 0) AND ("+graphInfo[4]+")";
	
		nodeInfo[0] = nodeChainCol;
		nodeInfo[1] = nodeNumCol;
	
		edgeInfo[0] = edgeINodeChainCol;
		edgeInfo[1] = edgeJNodeChainCol;
		edgeInfo[2] = edgeINodeNumCol;
		edgeInfo[3] = edgeJNodeNumCol;

    }

    public void setUserFilter(String filter) { graphInfo[4] = filter; }

    public void seNodeSizeRange(double[] range) { nodeSizeRange = range; }

    public void setEdgeSizeRange(double[] range) { edgeSizeRange = range; }

    public void setBackgroundColor(String color) { backgroundColor = color; }
    
    public boolean setUniformNodeColor(String color) {

		if (draw[0]) {
		    nodeColor = color;
		    nodeColorMethod = "Uniform";
		}
	
		return draw[0];

    }

    public boolean setUniformEdgeColor(String color) {

		if (draw[1]) {
		    edgeColor = color;
		    edgeColorMethod = "Uniform";
		}

		return draw[1];

    }

    public boolean setSpecialRes(String color, String query) {
	
		if (draw[2]) {
		    specialResColor = color;
		    specialResQuery = query;
		}

		return draw[2];

    }

    public boolean setUniformNodeSize(double size) {

		if (draw[0]) {
		    nodeSize = size;
		    nodeSizeMethod = "Uniform";
		}

		return draw[0];

    }

    public boolean setUniformEdgeSize(double size) {

		if (draw[1]) {
		    edgeSize = size;
		    edgeSizeMethod = "Uniform";
		}

		return draw[1];

    }

    public void setNodeColorMethod(String method, boolean discr) { nodeColorMethod = method; nodeColDiscr = discr; }

    public void setEdgeColorMethod(String method, boolean discr) { edgeColorMethod = method; edgeColDiscr = discr; }
    
    public void setNodeSizeMethod(String method, boolean rev) { nodeSizeMethod = method; nodeSizeRev = rev; }

    public void setEdgeSizeMethod(String method, boolean rev) { edgeSizeMethod = method; edgeSizeRev = rev; }

    public void setEdgeGapCondition(String condition) { edgeGapCondition = condition; }
    
    public void setEdgeGap(double gap) { edgeGap = gap; }

    public void setNodeTranspCondition(String condition) { nodeTranspCondition = condition; }
    
    public void setNodeTransp(double transp) { nodeTransp = transp; }

    public void setSurfTransp(double transp) { surfTransp = transp; }

    public boolean outputGraph(String connFile) {

		boolean status = check();
		if (!status) { return status; }
		exportPML();
		return status;

    }

    private boolean check() {

		boolean status = true;
	
		if (draw[0]) {
		    if ((nodeColorMethod == null) || (nodeColorMethod.equals("Uniform") && (nodeColor == null))) {
				System.out.println("Provide method for coloring your nodes or set the uniform node color!");
				status = false;
		    }
		    if ((nodeSizeMethod == null) || (nodeSizeMethod.equals("Uniform") && (nodeSize == 0))) {
				System.out.println("Provide method for defining your nodes' size or set the uniform node size!");
				status = false;
		    }
		    
		}
	
		if (draw[1]) {
		    if ((edgeColorMethod == null) || (edgeColorMethod.equals("Uniform") && (edgeColor == null))) {
				System.out.println("Provide method for coloring your edges or set the uniform edge color!");
				status = false;
		    }
		    if ((edgeSizeMethod == null) || (edgeSizeMethod.equals("Uniform") && (edgeSize == 0))) {
				System.out.println("Provide method for defining your edges' size or set the uniform edge size!");
				status = false;
		    }
		    if ((!edgeGapCondition.equals("")) && (directed)) {
		    	System.out.println("Gap condition will be ignored!");
		    }
		    
		}
	
		if ((draw[2]) && ((specialResColor == null) || (specialResQuery == null))) {
		    System.out.println("Set the color and query used for highlighting special residues!");
		    status = false;
		}
		
		return status;

    }

    public void exportPML() {

		pml.background(backgroundColor);
		pml.init();
		pml.createColor("light_grey", new double[] {0.8, 0.8, 0.8});
	
		if (draw[3]) { 
		    pml.showWhat("surface", molObjName);
		    pml.set("transparency", surfTransp, molObjName);
		}
	
		chains();
	
		if (draw[2]) { colorSpecialRes(); }
	
		pml.zoom("");
	
		if (draw[0]) {
		    createNodes();
	
		    if (nodeSizeMethod.equals("uniform")) {
		    	setNodesSize("nodes");
		    } else {
		    	setNodesSize();
		    }
	
		    if (nodeColorMethod.equals("uniform")) {
		    	setNodesColor("nodes");
		    } else {
		    	setNodesColor();
		    }
	
		    if (nodeTranspCondition.equals("")) {
		    	setNodesTransp("nodes");
		    } else {
		    	setNodesTransp();
		    }
		}
	
		pml.zoom("");
	
		if (draw[1]) {
	
		    if (cgoEdge) {
		    
			createCgoEdges();
	
		    } else {
			
			createEdges();
			
			/*
			if (edgeSizeMethod.equals("uniform")) {
			    setEdgesSize("edges");
			} else {
			    setEdgesSize();
			}
			*/
	
			if (edgeColorMethod.equals("uniform")) {
			    setEdgesColor("edges");
			} else {
			    //setEdgesColor();
			}
	
			/*
			if (directed) {
			    setGapDir();
			} else {
			    if (edgeGapCondition.equals("")) {
				setEdgesGap("edges");
			    } else {
				setEdgesGap();
			    }
			}
			*/
		    }
	
		}
	
		pml.zoom("");
  
    } // end of exportPML

    private void chains() {

		Statement S;
		ResultSet R;
		String cid = "", query = "", graphChainSet = "", chainSel = "";
		int numChains = -1;
		
		try { 
		    query = "SELECT DISTINCT("+nodeInfo[0]+") FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+dbInfo[4]+" = "+graphInfo[0]+") ORDER BY "+nodeInfo[0]+";";
		    System.out.println(query);
		    S = conn.createStatement();
		    R = S.executeQuery(query);
		    while (R.next()) { 
				cid = R.getString(nodeInfo[0]);
				numChains++;
				pml.createColor("color_cid"+cid, colors[numChains%colors.length]);
				chainSel = pml.selectChain(cid, msdsd);
				pml.setColor("color_cid"+cid, "\""+chainSel+"\"");
				graphChainSet = graphChainSet + chainSel+ " or ";
				pml.initList("edges_cid"+cid);
		    }
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch
	
		pml.selPml("restMol", molObjName+" and not ("+graphChainSet.substring(0,graphChainSet.length()-4)+")", false);
		pml.hide("restMol");

    }

    private void colorSpecialRes() {

		String cid = "", resSel = "";
		int num = 0;
		Statement S;
		ResultSet R;
	
		try {
		    S = conn.createStatement();	
		    R = S.executeQuery(specialResQuery);
		    while (R.next()) {
				cid = R.getString(nodeInfo[0]);
				num = R.getInt(nodeInfo[1]);
				resSel = resSel + pml.selectNode(cid, num, msdsd, false)+ " or ";	    
		    }
		    R.close();
		    S.close();
	
		    out.println();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch
	
		pml.selPml("specialRes", resSel.substring(0,resSel.length()-4), false);
		pml.setColor(specialResColor, "specialRes");

    } // end colorSpecialRes

    public String createNodes() {

		Statement S;
		ResultSet R;
		String query = "", nodes = "", node = "", chainNodesStr = "", cid = "";
	
		try { 
	
		    S = conn.createStatement();
	
		    query = "SELECT "+nodeInfo[0]+", "+nodeInfo[1]+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+nodeInfo[0]+", "+nodeInfo[1]+";";
		    System.out.println(query);
		    R = S.executeQuery(query);
		    while (R.next()) { 
				cid = R.getString(nodeInfo[0]);
				node = pml.addNode(cid, R.getInt(nodeInfo[1]), msdsd);
				chainNodes.put(cid, chainNodes.get(cid)+node+"+");
				nodes = nodes+node+"+";
		    } // end while loop through the resultset R
		    
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 
	
		for (String chain : chainNodes.keySet()) {
		    chainNodesStr = chainNodes.get(chain);
		    System.out.println(chain+" "+chainNodesStr);
		    chainNodesStr = "("+chainNodesStr.substring(4,chainNodesStr.length()-1)+")";
		    chainNodes.put(chain, chainNodesStr);
		    pml.selPml("nodes_cid"+chain, chainNodesStr, false);
		}
		nodes = "("+nodes.substring(0,nodes.length()-1)+")";
		pml.selPml("nodes", nodes, false);
		return nodes;

    } // end writeNodes()

    public void createEdges() {
    
		Statement S;
		ResultSet R;
		String query = "", edge = "", chainEdgesStr = "", i_cid = "";
	
		try { 
	
		    pml.initList("edges");
	
		    S = conn.createStatement();
	
		    query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
		    R = S.executeQuery(query);
		    while (R.next()) {
		    	i_cid = R.getString(edgeInfo[0]);	
		       	edge = pml.addEdge(i_cid, R.getInt(edgeInfo[2]), R.getString(edgeInfo[1]), R.getInt(edgeInfo[3]), msdsd);
				chainEdges.put(i_cid, chainEdges.get(i_cid)+edge+"+"); // to be removed
				pml.appendList("edges_cid"+i_cid, edge);
				pml.appendList("edges", edge);
		    } // end while loop through the resultset R
		    
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 
	
		// to be removed
		for (String chain : chainEdges.keySet()) {
		    chainEdgesStr = chainEdges.get(chain);
		    System.out.println(chain+" "+chainEdgesStr);
		    chainEdgesStr = "("+chainEdgesStr.substring(4,chainEdgesStr.length()-1)+")";
		    chainEdges.put(chain, chainEdgesStr);
		}

    }

    public void setNodesSize(String nodes) {

    	pml.set("sphere_scale", nodeSize, nodes);
	
    }

    public void setNodesSize() {

		Statement S;
		ResultSet R;
		String query = "", node = "";
		double curNodeSize = 0;
		double[] nodeCurSizeRange = new double[2];
	
		try { 
	
		    nodeCurSizeRange = utils4DB.getRange(conn, dbInfo[0]+"."+dbInfo[2], nodeSizeMethod, nodeGraphSel+" AND ("+graphInfo[5]+")");
	
		    S = conn.createStatement();
	
		    query = "SELECT "+nodeInfo[0]+", "+nodeInfo[1]+", "+nodeSizeMethod+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+nodeInfo[0]+", "+nodeInfo[1]+";";
		    R = S.executeQuery(query);
		    while (R.next()) { 
	
				node = "n."+R.getString(nodeInfo[0])+"."+R.getInt(nodeInfo[1]);
				curNodeSize = rescale(R.getDouble(nodeSizeMethod), nodeCurSizeRange, nodeSizeRange, nodeSizeRev);
				pml.set("sphere_scale", curNodeSize, node);
	
		    } // end while loop through the resultset R
		    
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 

    }

    public void setEdgesSize(String edges) {

    	pml.set("dash_width", edgeSize, edges);
	
    }

    public void setEdgesSize() {

		Statement S;
		ResultSet R;
		String query = "", edge = "";
		double curEdgeSize = 0;
		double[] edgeCurSizeRange = new double[2];
	
		try { 
	
		    edgeCurSizeRange = utils4DB.getRange(conn, dbInfo[0]+"."+dbInfo[3], edgeSizeMethod, edgeGraphSel+" AND ("+graphInfo[5]+")");
	
		    S = conn.createStatement();
	
		    query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+", "+edgeSizeMethod+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
		    R = S.executeQuery(query);
		    while (R.next()) { 
		
				edge = "e."+R.getString(nodeInfo[0])+"."+R.getInt(nodeInfo[2])+"."+R.getString(nodeInfo[1])+"."+R.getInt(nodeInfo[3]);
				curEdgeSize = rescale(R.getDouble(edgeSizeMethod), edgeCurSizeRange, edgeSizeRange, edgeSizeRev);
				pml.set("dash_width", curEdgeSize, edge);
	
		    } // end while loop through the resultset R
		    
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 

    }

    public void setNodesTransp(String nodes) {

    	pml.set("sphere_transparency", nodeTransp, nodes);

    }

    public void setNodesTransp() {

		Statement S;
		ResultSet R;
		String query = "", nodes = "";
	
		try { 
	
		    pml.set("sphere_transparency", 0, "nodes");
	
		    S = conn.createStatement();
	
		    query = "SELECT "+nodeInfo[0]+", "+nodeInfo[1]+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+")  AND ("+nodeTranspCondition+") ORDER BY "+nodeInfo[0]+", "+nodeInfo[1]+";";
		    R = S.executeQuery(query);
		    while (R.next()) { 
	
		    	nodes = nodes+"n."+R.getString(nodeInfo[0])+"."+R.getInt(nodeInfo[1])+"+";
			
		    } // end while loop through the resultset R
		    
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 
	
		nodes = "("+nodes.substring(0,nodes.length()-1)+")";
		pml.set("sphere_transparency", nodeTransp, nodes);

    }

    public void setEdgesGap(String edges) {

    	pml.set("dash_gap", edgeGap, edges);

    }

    public void setEdgesGap() {

		Statement S;
		ResultSet R;
		String query = "", edges = "";
	
		try { 
	
		    pml.set("dash_gap", 0, "edges");
	
		    S = conn.createStatement();
	
		    query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") AND ("+edgeGapCondition+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
		    R = S.executeQuery(query);
		    while (R.next()) { 
	
		    	edges = edges+"e."+R.getString(nodeInfo[0])+"."+R.getInt(nodeInfo[2])+"."+R.getString(nodeInfo[1])+"."+R.getInt(nodeInfo[3])+"+";
			
		    } // end while loop through the resultset R
		    
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 
	
		edges = "("+edges.substring(0,edges.length()-1)+")";
		pml.set("dash_gap", edgeGap, edges);

    }

    public void setGapDir() {

		Statement S;
		ResultSet R;
		String query = "", subQuery = "", edges = "";
	
		try { 
	
		    pml.set("dash_gap", 0, "edges");
	
		    S = conn.createStatement();
	
		    subQuery = "SELECT "+edgeInfo[0]+" AS i_cid, "+edgeInfo[2]+" AS i_num, "+edgeInfo[1]+" AS j_cid, "+edgeInfo[3]+" AS j_num FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+")";
	
		    query = "SELECT A.* FROM "+
			"("+subQuery+") AS A LEFT JOIN ("+subQuery+") AS B ON "+
			"(A.i_cid = B.j_cid AND A.i_num = B.j_num AND A.j_cid = B.i_cid AND A.j_num = B.i_num) "+
			"WHERE (B.i_cid IS NULL) ORDER BY i_cid, i_num, j_cid, j_num;";
		    R = S.executeQuery(query);
		    while (R.next()) { 
	
		    	edges = edges+"e."+R.getString("i_cid")+"."+R.getInt("i_num")+"."+R.getString("j_cid")+"."+R.getInt("j_num")+"+";
			
		    } // end while loop through the resultset R
		    
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 
	
		edges = "("+edges.substring(0,edges.length()-1)+")";
		pml.set("dash_gap", edgeGap, edges);

    }
    
    public void setNodesColor(String nodes) {

    	pml.setNodeColor(nodeColor, nodes);

    }

    public void setNodesColor() {

		Statement S;
		ResultSet R;
		String query = "", node = "", nodeSetId = "", nodeSetStr = "";
		double colWeight = 0;
		double[] nodeCurColRange = new double[2];
		double[] nodeRGB = new double[3];
		double[] nodeHSV = new double[3];
		int i = 0;
		
	
		try { 
	
		    S = conn.createStatement();
	
		    if (nodeColorMethod.equals("chain")) {
	
				for (String chain : chainNodes.keySet()) {
				    pml.setNodeColor("color_cid"+chain, "nodes_cid"+chain);
				}
	
		    } else if ((!nodeColorMethod.equals("uniform")) && (nodeColDiscr)) {
	
				query = "SELECT DISTINCT "+nodeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+nodeColorMethod+";";
				R = S.executeQuery(query);
				while (R.next()) {
				    i++;
				    pml.createColor("color_nodeSet."+i, colors[i%colors.length]);
				    nodeSetIds.put(new Long(Double.doubleToLongBits(R.getDouble(nodeColorMethod))), "nodeSet."+i);
				}
				R.close();
		
				query = "SELECT "+nodeInfo[0]+", "+nodeInfo[1]+", "+nodeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+nodeInfo[0]+", "+nodeInfo[1]+";";
				while (R.next()) {
				    node = "n."+R.getString(nodeInfo[0])+"."+R.getInt(nodeInfo[1])+"+";
				    colWeight = R.getDouble(nodeColorMethod);
				    nodeSetId = nodeSetIds.get(new Long(Double.doubleToLongBits(colWeight)));
				    nodeSet.put(node, nodeSetId);
				    nodeSets.put(nodeSetId, nodeSets.get(nodeSetId)+node+"+");
				}
				R.close();
		
				for (String curNodeSet : nodeSets.keySet()) {
				    nodeSetStr = nodeSets.get(curNodeSet);
				    nodeSetStr = "("+nodeSetStr.substring(0,nodeSetStr.length()-1)+")";
				    nodeSets.put(curNodeSet, nodeSetStr);
				    pml.selPml(curNodeSet, nodeSetStr, false);
				    pml.setNodeColor("color_"+curNodeSet, curNodeSet);
				}
	
		    } else if ((!nodeColorMethod.equals("uniform")) && (!nodeColDiscr)) {
	
				nodeCurColRange = utils4DB.getRange(conn, dbInfo[0]+"."+dbInfo[2], nodeColorMethod, nodeGraphSel+" AND ("+graphInfo[5]+")");
				
				query = "SELECT "+nodeInfo[0]+", "+nodeInfo[1]+", "+nodeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+nodeInfo[0]+", "+nodeInfo[1]+";";
				R = S.executeQuery(query);
				while (R.next()) {
				    node = "n."+R.getString(nodeInfo[0])+"."+R.getInt(nodeInfo[1]);
				    colWeight = R.getDouble(nodeColorMethod);
				    nodeHSV[0] = rescale(colWeight, nodeCurColRange, new double[] {0, 240}, false);
				    nodeHSV[1] = 1;
				    nodeHSV[2] = 1;
				    nodeRGB = hsvToRgb(nodeHSV);
				    pml.createColor("color_"+node, nodeRGB);
				    pml.setNodeColor("color_"+node, node);
				}
				R.close();
		
		    }   
	
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 

    }

    public void setEdgesColor(String edges) {

		pml.iterateList(edges, "edge");
		pml.setColor(edgeColor, "edge");

    }

    public void setEdgesColor() {

		Statement S;
		ResultSet R;
		String query = "", edge = "", edgeSetId = "", edgeSetStr = "", node = "";
		double colWeight = 0;
		double[] edgeCurColRange = new double[2];
		double[] edgeRGB = new double[3];
		double[] edgeHSV = new double[3];
		int i = 0;
		
	
		try { 
	
		    S = conn.createStatement();
	
		    if (edgeColorMethod.equals("chain")) {
	
				for (String chain : chainEdges.keySet()) {
				    pml.setColor("color_cid"+chain, "edges_cid"+chain);
				}
	
		    } else if (edgeColorMethod.equals("node")) {
	
				if (nodeColorMethod.equals("uniform")) {
				    pml.setColor(nodeColor, "edges");
				} else if (nodeColorMethod.equals("chain")) {
				    for (String chain : chainEdges.keySet()) {
				    	pml.setColor("color_cid"+chain, "edges_cid"+chain);
				    }		    
				} else if (nodeColDiscr) {
				    query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
				    R = S.executeQuery(query);
				    while (R.next()) {
						node = "n."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2]);
						edge = "e."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2])+"."+R.getString(edgeInfo[1])+"."+R.getInt(edgeInfo[3]);
						pml.setColor("color_"+nodeSet.get(node), "\""+edge+"\"");
				    }
				    R.close();
				} else if (!nodeColDiscr) {
				    query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
				    R = S.executeQuery(query);
				    while (R.next()) {
						node = "n."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2]);
						edge = "e."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2])+"."+R.getString(edgeInfo[1])+"."+R.getInt(edgeInfo[3]);
						pml.setColor("color_"+node, "\""+edge+"\"");
				    }
				    R.close();
				}
	
		    } else if ((!edgeColorMethod.equals("uniform")) && (edgeColDiscr)) {
	
				query = "SELECT DISTINCT "+edgeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeColorMethod+";";
				R = S.executeQuery(query);
				while (R.next()) {
				    i++;
				    pml.createColor("color_edgeSet."+i, colors[i%colors.length]);
				    edgeSetIds.put(new Long(Double.doubleToLongBits(R.getDouble(edgeColorMethod))), "edgeSet."+i);
				}
				R.close();
		
				query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+", "+edgeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
				while (R.next()) {
				    edge = "e."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2])+"."+R.getString(edgeInfo[1])+"."+R.getInt(edgeInfo[3]);
				    colWeight = R.getDouble(edgeColorMethod);
				    edgeSetId = edgeSetIds.get(new Long(Double.doubleToLongBits(colWeight)));
				    edgeSet.put(edge, edgeSetId);
				    edgeSets.put(edgeSetId, edgeSets.get(edgeSetId)+edge+"+");
				}
				R.close();
		
				for (String curEdgeSet : edgeSets.keySet()) {
				    edgeSetStr = edgeSets.get(curEdgeSet);
				    edgeSetStr = "("+edgeSetStr.substring(0,edgeSetStr.length()-1)+")";
				    edgeSets.put(curEdgeSet, edgeSetStr);
				    pml.selPml(curEdgeSet, edgeSetStr, false);
				    pml.setColor("color_"+curEdgeSet, curEdgeSet);
				}
	
		    } else if ((!edgeColorMethod.equals("uniform")) && (!edgeColDiscr)) {
	
				edgeCurColRange = utils4DB.getRange(conn, dbInfo[0]+"."+dbInfo[3], edgeColorMethod, edgeGraphSel+" AND ("+graphInfo[5]+")");
				
				query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+", "+edgeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
				R = S.executeQuery(query);
				while (R.next()) {
				    edge = "e."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2])+"."+R.getString(edgeInfo[1])+"."+R.getInt(edgeInfo[3]);
				    colWeight = R.getDouble(edgeColorMethod);
				    edgeHSV[0] = rescale(colWeight, edgeCurColRange, new double[] {0, 240}, false);
				    edgeHSV[1] = 1;
				    edgeHSV[2] = 1;
				    edgeRGB = hsvToRgb(edgeHSV);
				    pml.createColor("color_"+edge, edgeRGB);
				    pml.setColor("color_"+edge, "\""+edge+"\"");
				}
				R.close();
	
		    }   
	
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 

    }

    public void createCgoEdges() {

    }

    private double rescale(double value, double[] curRange, double[] newRange, boolean rev) {

		if (!rev) {
		    return (newRange[0] + ((value-curRange[0])/(curRange[1]-curRange[0]))*(newRange[1]-newRange[0]));
		} else {
		    return (newRange[1] - ((value-curRange[0])/(curRange[1]-curRange[0]))*(newRange[1]-newRange[0]));
		}
		
    }

    private double[] hsvToRgb(double[] hsv) {
	
		double rgb[] = {0, 0, 0};
		double f = 0, p = 0, q = 0, t = 0;
		int i = 0;
	
		if (hsv[1] == 0) {
		    //achromatic (grey)
		    rgb[0] = hsv[2];
		    rgb[1] = hsv[2];
		    rgb[2] = hsv[2];
		} else {
		    hsv[0] = hsv[0]/(double)60;
		    i = (int) Math.floor(hsv[0]);
		    f = hsv[0] - i;
		    
		    p = hsv[2] * ( 1 - hsv[1] );
		    q = hsv[2] * ( 1 - hsv[1] * f );
		    t = hsv[2] * ( 1 - hsv[1] * ( 1 - f ) );
		    
		    switch (i) {
	                case 0:
			    rgb[0] = hsv[2];
			    rgb[1] = t;
			    rgb[2] = p;
			    break;
	                case 1:
			    rgb[0] = q;
			    rgb[1] = hsv[2];
			    rgb[2] = p;
			    break;
	                case 2:
			    rgb[0] = p;
			    rgb[1] = hsv[2];
			    rgb[2] = t;
			    break;
		        case 3:
			    rgb[0] = p;
			    rgb[1] = q;
			    rgb[2] = hsv[2];
			    break;
	                case 4:
			    rgb[0] = t;
			    rgb[1] = p;
			    rgb[2] = hsv[2];
			    break;
		        case 5:
			    rgb[0] = hsv[2];
			    rgb[1] = p;
			    rgb[2] = q;
			    break;
		        default:
			    rgb[0] = hsv[2];
			    rgb[1] = hsv[2];
			    rgb[2] = hsv[2];
		    }
		}

		return(rgb);

    } // end of hsvToRgb

} // end of class Graph2Pml

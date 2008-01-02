package tools;

import java.sql.*;
import java.io.*;
import java.util.*;


/**
 * Package:		tools
 * Class: 		Graph2Pml
 * Author:		Ioannis Filippis, filippis@molgen.mpg.de
 * Date:		21/03/2006
 *  
 * This class enables the visualization of a contact graph. Contact
 * graph is defined based on the contact weight, the contact type, the
 * contact range as well as a user defined filter. User can also define 
 * chain coloring, color/size/transparency of nodes, color/width/gap of
 * edges, representation of edges like lines or cylinders, directedness
 * of graph, etc ... Moreover, the schema of the tables where the graph
 * is extracted from, can be modified.
 * 
 * 
 * Notes:
 * - The program defines some python/pymol string lists that could be used later for
 * 	 easy manipulation of groups. You may encounter the following lists:
 *   - nodes, edges : all nodes/edges
 * 	 - nodes_cidX, edges_cidX (where X the chain_code) contain all nodes/edges per
 * 	  chain
 * 	 - 
 * - The visualization is not "static". Although, the whole class is defined
 *   for a single visualisation instance, all methods could be used for on-the-fly
 *   visualization changes.
 * - The coding is not optimum since the purposes of the visualization are not
 *   well defined at the moment. The main purpose of the class is to provide
 *   insight how the PyMol class can be used, as well as ready code for custom
 *   implementations.
 * - All properties of nodes and edges are handled separately for 2 reasons
 *	 - Even though nodes(edges) could be created and their properties could be
 *	   set at the same time (traversing one the nodes(edges) from the database), 
 *	   it was preferred to be handled in different functions. In this way, these
 *     functions can be called later just for a specific purpose instead of
 *     having to recreate the nodes from the beginning.
 *   - Despite the fact that the code for setting node properties is the same
 *     with the code for the edge properties (except some minor changes), separate
 *     methods are preferred for the same property, one for the node and one for the
 *     color. This is to accommodate future changes to the code, that might not be
 *     common to edges and nodes. 
 * - There are actually 3 implementations in this code, 2 of them commented
 * 	 These implementations concern how to handle in groups nodes, edges. For
 * 	 example, handle simultaneously all the nodes or all the nodes of a chain
 *   or all the nodes of a connected component for coloring purposes.
 *   -Java string method: The first idea was to keep somehow in memory pymol selection
 *    strings and create PyMol selection objects with that. For example, the string
 *    variable nodesConnectedComponent could have value "n.A.1+n.A.3+n.A.4" and a
 *    pymol selection could be easily made using that. However, pymol command length
 *    is limited and especially for edges the command could easily exceed the limit.
 *   -Pymol string method: The second idea was to create the string in pymol stepwise.
 *    So the string would be a python string actually and the command line would be small
 *    since each time that a node was added in the graph, it would be added also in
 *    the string e.g. nodesConnectedComponent += "n.A.1". However, still the command line
 *    was exceeded. For example when the value of the all edges string exceeded the limit,
 *    calling a pymol command with parameter edges would make pymol crash.
 *   -Final method: A pymol list of objects is created and you just iterate the list giving
 *    the command you wish. The nice thing is that if your graph is not large, then you
 *    can still converge to the previous method by concatenating the list. See example in
 *    testGraph2Pml. Of course, the list solution is slower but it will never crash.
 *   -It must be noticed that the edges modelled as distance objects can not be selected
 *    in pymol since they do not contain any atoms.
 *   -The 3 different implementations are limited in the chains, createNodes, createEdges,
 *    setNodesColorUniform(String, String), setEdgesColorUniform(String, String) methods.
 *   -The only point where I don't use list or where I modify node/edge one-by-one, is in 
 *    the colorSpecialRes method. I expect these residues not to be many so the command line
 *    length limit won't be exceeded.
 * - Never use "." for a list's name!!!!
 * - If you put surface even with high transparency, you can not see the complete graph.
 * 
 * Changelog:
 * 28/03/06 modified by IF (adding more flexibility on coloring, popDefaults and setDefaults methods and more comments)
 * 27/03/06 modified by IF (functional at last - added cgoEdges)
 * 22/03/06 modified by IF (modified to compile based on objectNameQuotes variable changes in PyMol class - still non functional)
 * 21/03/06 first created by IF (non functional)
 */

public class Graph2Pml {

	// check constructor
    private Connection conn;
    private PrintWriter out = null;
    private PyMol pml = null;
    private boolean msdsd = true, cgoEdge = false, directed = true;
    private String molObjName = null;
    
    /*
     * chains: HashSet holding all chain codes
     * nodeSetIds: HashMap containing nodes' possible values based on which
     * 			nodes are separated into sets. Values are the set names. Used
     * 			when the coloring of nodes is based on a field and colors are
     * 			discretised.
     * nodeSet:  HashMap containing nodes as keys and the set they belong to
     * 			as values. This is used when the edges are colored as the nodes
     * 			("node" method) and the coloring of nodes is based on a field 
     * 			with colors to be discretised.
     * edgeSetIds, edgeSet:  similar with nodeSetIds, nodeSet. edgeSet is not
     * 			actually used.
     * defaults: HashMap containing all attributes/variables that can be set
     * 			to their default values using setDefaults method. Values are 
     * 			integers so a switch can be utilised.  
     */
    private HashSet<String> chains = new HashSet<String>();
    private HashMap<Long, String> nodeSetIds = new HashMap<Long, String> ();
    private HashMap<String, String> nodeSet = new HashMap<String, String> ();
    private HashMap<Long, String> edgeSetIds = new HashMap<Long, String> ();
    private HashMap<String, String> edgeSet = new HashMap<String, String> ();   
    //private HashMap<String, String> chainNodes = new HashMap<String, String>(); //java String method
    //private HashMap<String, String> chainEdges = new HashMap<String, String>(); //java String method
    private HashMap<String, Integer> defaults = new HashMap<String, Integer>();
    
    /*
     * dbInfo, nodeInfo and edgeInfo hold the db schema info, while the graphInfo
     * 	concenrns the graph model
     * dbInfo: String array holding the db name, the table name with the list of available
     * 		graphs, the nodes table name, the edges table name and the column name for the
     * 		graph id
     * nodeInfo: String array holding the chain code and the residue serial column names
     * edgeInfo: String array holding the chain code and the residue serial column names
     * 		for the i and j edges
     * graphInfo: String array holding the graph id, the edge weight, the node contact type,
     * 		the edge contact type, the contact range and the user-defined filter
     * nodeGraphSel: String containing the mysql where-condition when quering the nodes table
     * edgeGraphSel: String containing the mysql where-condition when quering the edges table 
     */
    private String[] dbInfo = new String[] {"newmsdgraph", "list", "nodes", "edges", "graph_id"};
    private String[] nodeInfo = new String[] {"cid", "num"};
    private String[] edgeInfo = new String[] {"i_cid", "j_cid", "i_num", "j_num"};
    private String[] graphInfo = new String[6];
    private String nodeGraphSel = "", edgeGraphSel = "";
    
    /*
     * draw: boolean array defining whether nodes, edges, special residues, surface will 
     * 		drawn
     * nodeColorMethod: String containing the method for defining the node color. These can be: 
     * 		- "chain" (default): all nodes of the same chain will have the same color, either defined
     * 			by the user or by the color of the chain they belong to
     * 		- "uniform": all nodes will have the same color
     * 		- db field based: in this case either the values of the field are rescaled in the rgb
     * 			scale (from red to blue) or all nodes of the same value will have specific color
     * 			(discretised case)
     * edgeColorMethod: Similar as above plus one more method:
     * 		- "node": edges are colored based on the color of the i_node.
     * nodeSizeMethod, edgeSizeMethod: the same as above except the "chain" and the "node"
     * nodeColDiscr, edgeColDiscr : whether the coloring is discretised in the db field based case
     * nodeSizeRev, edgeSizeRev : whether values should be considered in the descending order in the
     * 		db field based case (high value -> small sphere)
     * nodeColor, edgeColor, specialResColor, backgroundColor: ...
     * nodeSizeRange, edgeSizeRange, cgoEdgeSizeRange: ...
     * edgeGap, cgoEdgeGap, cgoEdgeLength, nodeTransp, surfTransp: ...
     * 		Cylinder is defined by the cgoEdgeGap, cgoEdgeLength and cgoEdgeSize (actually the radius).
     * 		Only the last one can be defined by the user
     * specialResQuery: a mysql query statement that identifies residues(not nodes - spheres) that 
     * 		should be colored with specialResColor. For example, this could be used in a SC_SC graph to 
     * 		mark the SC_SC=0 regions in the chain, or early folding residues
     * edgeGapCondition: a mysql where-condition based on which specific edges will have gaps. If the graph is
     * 		defined as directed, then this will be ignored
     * nodeTranspCondition: a mysql where-condition based on which specific nodes will be transparent
     * 
     */
    private boolean[] draw = new boolean[] {true, true, false, false};
    private String nodeColorMethod = "chain", nodeSizeMethod = "uniform", edgeColorMethod = "uniform", edgeSizeMethod = "uniform", nodeColor = "blue", edgeColor = "orange", specialResColor = "purple";
    private boolean nodeColDiscr = false, edgeColDiscr = false, nodeSizeRev = false, edgeSizeRev = false;
    private double nodeSize = 0.46, edgeSize = 3.65, cgoEdgeSize = 0.3;
    private double[] nodeSizeRange = {0.2, 0.8}; //{0.18, 0.75};
    private double[] edgeSizeRange = {0.3, 7.0};
    private double[] cgoEdgeSizeRange = {0.2, 0.5};
    private String edgeGapCondition = "", nodeTranspCondition = "", specialResQuery = null;
    private double edgeGap = 0.25, nodeTransp = 0.6, surfTransp = 0.6, cgoEdgeGap = 0, cgoEdgeLength = 0.5;
    private String backgroundColor = "white";
    
    /*
     * All the following string arrays holding colors can be changed by the user!
     * chainColors: default chain colors (chain A will be light_grey, second green and so on).
     * chainNodeColors, chainEdgeColors: String arrays holding the colors to be used for
     * 		coloring the nodes(edges) if "chain" is the coloring methods. By default, these are null
     * 		and the chainColors are used instead.
     * setNodeColors, setEdgeColors: String array holding the colors to be used for nodes/edges
     * 		based on the set they belong to (db field based, discretised coloring). By default, are similar
     * 		to the chainColors (only the first light_grey color is omitted)
     */  
    private String[] chainColors = {"light_grey", "green", "red", "blue", "yellow", "violet", "cyan", "salmon", "lime", "pink", "slate", "magenta", "orange", "marine", "olive", "purple", "teal", "forest", "firebrick", "chocolate", "wheat", "white", "grey"};
    private String[] chainNodeColors = null;
    private String[] chainEdgeColors = null;
    private String[] setNodeColors = {"green", "red", "blue", "yellow", "violet", "cyan", "salmon", "lime", "pink", "slate", "magenta", "orange", "marine", "olive", "purple", "teal", "forest", "firebrick", "chocolate", "wheat", "white", "grey"};
    private String[] setEdgeColors = {"green", "red", "blue", "yellow", "violet", "cyan", "salmon", "lime", "pink", "slate", "magenta", "orange", "marine", "olive", "purple", "teal", "forest", "firebrick", "chocolate", "wheat", "white", "grey"};
    
    /*
     * @param out: PrintWriter where pymol commands will be send to
     * @param molObjName: the object name where the molecule is loaded to
     * @param graphId: the graph IDX
     * @param edgeWeight, nodeContactType, edgeContactType, contactRange, userFilter
     * 		All should be given in the format of a sql expression embedded in a where clause
     * @param edgeWeight: the weight to be used for the edge (non functional)
     * 		This should not be confused with contact type. For example, you may choose 
     * 		all SC_SC edges but the edge weight to be uniform.
     * @param nodeContactType, edgeContact Type: both refer to contact type (SC_SC, BB_BB, .. edges)
     * 		Contact type is defined separately due to different db fields that have to be >0. For 
     * 		example for newmsdgraph	db and for a SC_SC graph, nodeContactType should be "SC_SC_in+SC_SC_out",
     * 		while edgeContactType "SC_SC".
     * @param contactRange: inter-secondary-structure/SC-dominated/long-medium-short range edges
     * 		This filter is applied only to edges. So, the nodes will be the same whether you
     * 		define a contactRange filter or not (in this case you should use "(true)").
     * @param userFilter: a user defined filter
     * 		This is applies both to nodes and edges! If not applicable, just give "(true)".
     * @param msdsd: whether the pdb file is extracted from msdsd
     * @param cgoEdge: whether the edge should be represented as a cylinder
     * @param directed: whether the edge should be directed 
     * 		If edge is line, then it will have gaps if not present in both directions else solid
     * 		If edge is cylinder, then directed edge is represented as cylinder-line combination
     * @param conn: mysql connection
     *  
     */
    public Graph2Pml(PrintWriter out, String molObjName, int graphId, String edgeWeight, String nodeContactType, String edgeContactType, String contactRange, String userFilter, boolean msdsd, boolean cgoEdge, boolean directed, Connection conn) {
	
		graphInfo[0] = String.valueOf(graphId);
		graphInfo[1] = edgeWeight;
		graphInfo[2] = nodeContactType;
        graphInfo[3] = edgeContactType;
		graphInfo[4] = contactRange;
		graphInfo[5] = userFilter;
		nodeGraphSel = "("+dbInfo[4]+" = "+graphInfo[0]+") AND ("+graphInfo[2]+" > 0) AND (true)";
		edgeGraphSel = "("+dbInfo[4]+" = "+graphInfo[0]+") AND ("+graphInfo[3]+" > 0) AND ("+graphInfo[4]+")";
		
		this.conn = conn;
		this.out = out;
		this.pml = new PyMol(out);
		this.molObjName = molObjName;
		this.msdsd = msdsd;
		this.cgoEdge = cgoEdge;
		this.directed = directed;
		
		// initialise defaults hashMap
		popDefaults();

    }
    
    /*
     * define whether nodes, edges, special residues, surface will be drawn
     */
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
    
    /*
     * define the schema of the database from where the contact graph will be extracted
     */
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

    public void setNodeSizeRange(double[] range) { nodeSizeRange = range; }

    public void setEdgeSizeRange(double[] range) { edgeSizeRange = range; }
    
    public void setCgoEdgeSizeRange(double[] range) { cgoEdgeSizeRange = range; }
    
    public void setBackgroundColor(String color) { backgroundColor = color; }
    
    public void setChainColors(String[] colors) { chainColors = colors; }
    
    public void setNodeColors(String[] colors) { setNodeColors = colors; }
    
    public void setEdgeColors(String[] colors) { setEdgeColors = colors; }
    
    public boolean setUniformNodeColor(String color) {

		if (draw[0]) {
		    nodeColor = color;
		    nodeColorMethod = "uniform";
		}
	
		return draw[0];

    }

    public boolean setUniformEdgeColor(String color) {

		if (draw[1]) {
		    edgeColor = color;
		    edgeColorMethod = "uniform";
		}

		return draw[1];

    }
    
    public boolean setChainNodeColor() {

		if (draw[0]) {
			chainNodeColors = null;
		    nodeColorMethod = "chain";
		}
	
		return draw[0];

    }
    
    public boolean setChainNodeColor(String[] colors) {

		if (draw[0]) {
		    chainNodeColors = colors;
		    nodeColorMethod = "chain";
		}
	
		return draw[0];

    }
    
    public boolean setChainEdgeColor() {

		if (draw[1]) {
			chainEdgeColors = null;
		    edgeColorMethod = "chain";
		}

		return draw[1];

    }
    
    public boolean setChainEdgeColor(String[] colors) {

		if (draw[1]) {
		    chainEdgeColors = colors;
		    edgeColorMethod = "chain";
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
		    nodeSizeMethod = "uniform";
		}

		return draw[0];

    }

    public boolean setUniformEdgeSize(double size) {

		if (draw[1]) {
		    edgeSize = size;
		    edgeSizeMethod = "uniform";
		}

		return draw[1];

    }
    
    public boolean setUniformCgoEdgeSize(double size) {

		if (draw[1]) {
		    cgoEdgeSize = size;
		    edgeSizeMethod = "uniform";
		}

		return draw[1];

    }

    /*
     * define the db field to be used for coloring the nodes. If discr (discretised is true), 
     * all nodes of the same value will have specific color, else the values of the field 
     * are rescaled in the rgb scale (from red to blue)
     */
    public void setNodeColorMethod(String method, boolean discr) { nodeColorMethod = method; nodeColDiscr = discr; }
    
    /*
     * define the db field to be used for coloring the edges. If discr (discretised is true), 
     * all edges of the same value will have specific color, else the values of the field 
     * are rescaled in the rgb scale (from red to blue)
     */
    public void setEdgeColorMethod(String method, boolean discr) { edgeColorMethod = method; edgeColDiscr = discr; }
    
    /*
     * define the db field to be used for coloring the nodes. If rev (reversed is true), 
     * values will be considered in the descending order (high value -> small sphere).
     */    
    public void setNodeSizeMethod(String method, boolean rev) { nodeSizeMethod = method; nodeSizeRev = rev; }

    /*
     * define the db field to be used for coloring the edges. If rev (reversed is true), 
     * values will be considered in the descending order (high value -> thin edge).
     */       
    public void setEdgeSizeMethod(String method, boolean rev) { edgeSizeMethod = method; edgeSizeRev = rev; }

    public void setEdgeGapCondition(String condition) { edgeGapCondition = condition; }
    
    public void setEdgeGap(double gap) { edgeGap = gap; }

    public void setNodeTranspCondition(String condition) { nodeTranspCondition = condition; }
    
    public void setNodeTransp(double transp) { nodeTransp = transp; }

    public void setSurfTransp(double transp) { surfTransp = transp; }

    /*
     * visualise the graph
     */
    public boolean outputGraph() {

		boolean status = check();
		if (!status) { return status; }
		exportPML();
		return status;

    }

    /*
     * private method to do some simple checks on the compatibility of user defined preferences
     */
    private boolean check() {

		boolean status = true;
	
		if (draw[0]) {
		    if ((nodeColorMethod == null) || (nodeColorMethod.equals("uniform") && (nodeColor == null))) {
				System.out.println("Provide method for coloring your nodes or set the uniform node color!");
				status = false;
		    }
		    if ((nodeSizeMethod == null) || (nodeSizeMethod.equals("uniform") && (nodeSize == 0))) {
				System.out.println("Provide method for defining your nodes' size or set the uniform node size!");
				status = false;
		    }
		    
		}
	
		if (draw[1]) {
		    if ((edgeColorMethod == null) || (edgeColorMethod.equals("uniform") && (edgeColor == null))) {
				System.out.println("Provide method for coloring your edges or set the uniform edge color!");
				status = false;
		    }
		    if ((edgeColorMethod.equals("node")) && (!draw[0])) {
				System.out.println("Since nodes are not drawn, edges can not be colored based on nodes coloring!");
				status = false;
		    }
		    if ((edgeSizeMethod == null) || (edgeSizeMethod.equals("uniform") && (edgeSize == 0))) {
				System.out.println("Provide method for defining your edges' size or set the uniform edge size!");
				status = false;
		    }
		    if ((!edgeGapCondition.equals("")) && (directed)) {
		    	System.out.println("Gap condition will be ignored!");
		    }
		    if ( (cgoEdge) && ((!edgeGapCondition.equals("")) || (edgeGap != 0)) ) {
		    	System.out.println("Gap is ignored for CGO edges!");
		    }
		}
	
		if ((draw[2]) && ((specialResColor == null) || (specialResQuery == null))) {
		    System.out.println("Set the color and query used for highlighting special residues!");
		    status = false;
		}
		
		return status;

    }

    /*
     * private method that exports all pymol commands for visualization
     */
    private void exportPML() {

		pml.background(backgroundColor);
		pml.init();
		pml.createColor("light_grey", new double[] {0.8, 0.8, 0.8});
	
		chains();
		
		//draw special residues
		if (draw[2]) { colorSpecialRes(); } 
		pml.zoom("graphMol", false);
	
		//draw nodes
		if (draw[0]) {
		    createNodes();
		    pml.zoom("graphMol", false);
		    if (nodeSizeMethod.equals("uniform")) {
		    	setNodesSizeUniform("nodes", nodeSize);
		    } else {
		    	setNodesSize();
		    }
		    pml.zoom("graphMol", false);			
		    if (nodeColorMethod.equals("uniform")) {
		    	setNodesColorUniform("nodes", nodeColor);
		    } else {
		    	setNodesColor();
		    }
		    pml.zoom("graphMol", false);
		    if (nodeTranspCondition.equals("")) {
		    	setNodesTranspUniform("nodes", nodeTransp);
		    } else {
		    	setNodesTransp();
		    }
		    pml.zoom("graphMol", false);
		}
	
		//draw edges
		if (draw[1]) {	
		    if (cgoEdge) {		    
		    	createCgoEdges();
		    	pml.zoom("graphMol", false);
		    } else {			
				createEdges();
				pml.zoom("graphMol", false);
				if (edgeSizeMethod.equals("uniform")) {
				    setEdgesSizeUniform("edges", edgeSize);
				} else {
				    setEdgesSize();
				}
				pml.zoom("graphMol", false);
				if (edgeColorMethod.equals("uniform")) {
				    setEdgesColorUniform("edges", edgeColor);
				} else {
				    setEdgesColor();
				}
				pml.zoom("graphMol", false);
				if (directed) {
				    setGapDir();
				} else {
				    if (edgeGapCondition.equals("")) {
				    	setEdgesGapUniform("edges", edgeGap);
				    } else {
				    	setEdgesGap();
				    }
				}
				pml.zoom("graphMol", false);				
		    }	
		}
		
		//draw surface
		if (draw[3]) { 
		    pml.showWhat("surface", "graphMol", false);
		    pml.set("transparency", surfTransp, "graphMol", false);
		}
		pml.zoom("graphMol", false);
		
    } // end of exportPML

	/*
	 * private method to color the chains, initialise chain based lists and hide the
	 * rest of the macromolecule (restMol) not included in the contact graph (graphMol)
	 */
    private void chains() {

		Statement S;
		ResultSet R;
		String cid = "", query = "", graphChainSet = "", chainSel = "";
		int numChains = -1;
		
		try { 
		    query = "SELECT DISTINCT("+nodeInfo[0]+") FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+dbInfo[4]+" = "+graphInfo[0]+") ORDER BY "+nodeInfo[0]+";";
		    S = conn.createStatement();
		    R = S.executeQuery(query);
		    while (R.next()) { 
				cid = R.getString(nodeInfo[0]);
				chains.add(cid);
				numChains++;
				
				pml.createColor("color_cid"+cid, chainColors[numChains%chainColors.length]);
				chainSel = pml.selectChain(cid, msdsd);
				pml.setColor("color_cid"+cid, chainSel, false);
				
				pml.initList("nodes_cid"+cid);
				pml.initList("edges_cid"+cid);
				//pml.initString("nodes_cid"+cid, ""); // pymol string method
				//pml.initString("edges_cid"+cid, ""); // pymol string method	
				
				graphChainSet = graphChainSet + chainSel+ " or ";				
		    }
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch
	
		pml.selPml("graphMol", graphChainSet.substring(0,graphChainSet.length()-4), false, false);
		pml.selPml("restMol", molObjName+" and not graphMol", false, false);
		pml.hide("restMol", false);

    }

    /*
     * private method to color special residues
     */
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
	
		pml.selPml("specialRes", resSel.substring(0,resSel.length()-4), false, false);
		pml.setColor(specialResColor, "specialRes", false);

    } // end colorSpecialRes

    /*
     * creates all nodes with pymol default settings and appends them in the nodes, nodes_cidX lists
     */
    public void createNodes() {

		Statement S;
		ResultSet R;
		String query = "", node = "", cid = "";
		//String nodes = "", chainNodesStr = ""; //java String method
	
		try { 
	
			pml.initList("nodes");
			//pml.initString("nodes", ""); // pymol String method
			
		    S = conn.createStatement();	
		    query = "SELECT "+nodeInfo[0]+", "+nodeInfo[1]+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+nodeInfo[0]+", "+nodeInfo[1]+";";
		    R = S.executeQuery(query);
		    while (R.next()) { 
				cid = R.getString(nodeInfo[0]);
				node = pml.addNode(cid, R.getInt(nodeInfo[1]), msdsd);
				pml.appendList("nodes_cid"+cid, node);
				pml.appendList("nodes", node);
				//pml.appendString("nodes_cid"+cid, node+"+"); // pymol String method
				//pml.appendString("nodes", node+"+"); // pymol String method
				//chainNodes.put(cid, chainNodes.get(cid)+node+"+"); // java String method
				//nodes = nodes+node+"+"; // java String method				
		    } 
		    
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 
	
		/*
		for (String chain : chainNodes.keySet()) {
			//pymol String method
			pml.rstripString("nodes_cid"+chain, "+");
		    pml.selPml("nodes_cid"+chain, "nodes_cid"+chain, false, true);
			//java String method
		    chainNodesStr = chainNodes.get(chain);
		    chainNodesStr = chainNodesStr.substring(4,chainNodesStr.length()-1);
		    chainNodes.put(chain, chainNodesStr);
		    pml.selPml("nodes_cid"+chain, chainNodesStr, false, false);		    	    
		}
		//pymol String method
		pml.rstripString("nodes", "+");
		pml.selPml("nodes", "nodes", false, true);
		//java String method
		nodes = nodes.substring(0,nodes.length()-1);
		pml.selPml("nodes", nodes, false, false);
		*/
		
    } // end writeNodes()
    
    /*
     * creates all edges with pymol default settings and appends them in the edges, edges_cidX list
     */    
    public void createEdges() {

		Statement S;
		ResultSet R;
		String query = "", edge = "", i_cid = "";
		//String edges = "", chainEdgesStr = ""; //java String method
	
		try { 
	
			pml.initList("edges");
			//pml.initString("edges", ""); // pymol String method
			
		    S = conn.createStatement();	
		    query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
		    R = S.executeQuery(query);
		    while (R.next()) { 
				i_cid = R.getString(edgeInfo[0]);
				edge = pml.addEdge(i_cid, R.getInt(edgeInfo[2]), R.getString(edgeInfo[1]), R.getInt(edgeInfo[3]), msdsd);
				pml.appendList("edges_cid"+i_cid, edge);
				pml.appendList("edges", edge);
				//pml.appendString("edges_cid"+i_cid, edge+"+"); // pymol String method
				//pml.appendString("edges", edge+"+"); // pymol String method
				//chainEdges.put(cid, chainEdges.get(i_cid)+edge+"+"); // java String method
				//edges = edges+edge+"+"; // java String method				
		    } 
		    
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 
	
		/*
		for (String chain : chainEdges.keySet()) {
			//pymol String method
			pml.rstripString("edges_cid"+chain, "+");
			//java String method
		    chainEdgesStr = chainEdges.get(chain);
		    chainEdgesStr = chainEdgesStr.substring(4,chainEdgesStr.length()-1);
		    chainEdges.put(chain, chainEdgesStr);   	    
		}
		//pymol String method
		pml.rstripString("edges", "+");
		//java String method
		edges = edges.substring(0,nodes.length()-1);
		*/
		
    } // end writeEdges()
    
    /*
     * sets the node size for all nodes included in the @param nodes list
     */
    public void setNodesSizeUniform(String nodes, double size) {

    	pml.iterateList(nodes, "node");
		pml.set("sphere_scale", size, "node", true);
	
    }

    /*
     * sets the size for all nodes based on a db field
     */
    public void setNodesSize() {

		Statement S;
		ResultSet R;
		String query = "", node = "";
		double curNodeSize = 0;
		double[] nodeCurSizeRange = new double[2];
	
		try { 
	
			// get min, max values of the db field
		    nodeCurSizeRange = utils4DB.getRange(conn, dbInfo[0]+"."+dbInfo[2], nodeSizeMethod, nodeGraphSel+" AND ("+graphInfo[5]+")");
	
		    S = conn.createStatement();
	
		    query = "SELECT "+nodeInfo[0]+", "+nodeInfo[1]+", "+nodeSizeMethod+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+nodeInfo[0]+", "+nodeInfo[1]+";";
		    R = S.executeQuery(query);
		    while (R.next()) { 
	
				node = "n."+R.getString(nodeInfo[0])+"."+R.getInt(nodeInfo[1]);
				curNodeSize = rescale(R.getDouble(nodeSizeMethod), nodeCurSizeRange, nodeSizeRange, nodeSizeRev);
				pml.set("sphere_scale", curNodeSize, node, false);
	
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
    
    /*
     * sets the edges size for all nodes included in the @param edges list
     */
    public void setEdgesSizeUniform(String edges, double size) {

    	pml.iterateList(edges, "edge");
		pml.set("dash_width", size, "edge", true);
	
    }
    
    /*
     * sets the size for all nodes based on a db field
     */
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
		
		    	edge = "e."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2])+"."+R.getString(edgeInfo[1])+"."+R.getInt(edgeInfo[3]);
		    	curEdgeSize = rescale(R.getDouble(edgeSizeMethod), edgeCurSizeRange, edgeSizeRange, edgeSizeRev);
				pml.set("dash_width", curEdgeSize, edge, false);
	
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

    /*
     * sets the node transparency for all nodes included in the @param nodes list
     */   
    public void setNodesTranspUniform(String nodes, double transp) {

    	pml.iterateList(nodes, "node");
		pml.set("sphere_transparency", transp, "node", true);

    }
    
    /*
     * sets the node transparency for all nodes defined by the nodeTranspCondition
     */
    public void setNodesTransp() {

		Statement S;
		ResultSet R;
		String query = "", node = "";
	
		try { 
	
			//initialise transparency to 0
	    	pml.iterateList("nodes", "node");
			pml.set("sphere_transparency", 0, "node", true);
	
		    S = conn.createStatement();
	
		    query = "SELECT "+nodeInfo[0]+", "+nodeInfo[1]+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+")  AND ("+nodeTranspCondition+") ORDER BY "+nodeInfo[0]+", "+nodeInfo[1]+";";
		    R = S.executeQuery(query);
		    while (R.next()) { 
		    	node = "n."+R.getString(nodeInfo[0])+"."+R.getInt(nodeInfo[1]);
		    	pml.set("sphere_transparency", nodeTransp, node, false);
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
    
    /*
     * sets the gap for all edges included in the @param edges list
     */   
    public void setEdgesGapUniform(String edges, double gap) {

    	pml.iterateList(edges, "edge");
		pml.set("dash_gap", gap, "edge", true);

    }
    
    /*
     * sets the gap for all edges defined by the edgeGapCondition
     */
    public void setEdgesGap() {

		Statement S;
		ResultSet R;
		String query = "", edge = "";
	
		try { 
	
			// initialise gap to 0
			pml.iterateList("edges", "edge");
			pml.set("dash_gap", 0, "edge", true);
				
		    S = conn.createStatement();
	
		    query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") AND ("+edgeGapCondition+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
		    R = S.executeQuery(query);
		    while (R.next()) { 
		    	edge = "e."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2])+"."+R.getString(edgeInfo[1])+"."+R.getInt(edgeInfo[3]);
		    	pml.set("dash_gap", edgeGap, edge, false);
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

    /*
     * edges that exist only in one direction are set to have gap, while all bidirectional ones
     * are solid
     */
    public void setGapDir() {

		Statement S;
		ResultSet R;
		String query = "", subQuery = "", edge = "";
	
		try { 
			
			//initialise gap to 0
			pml.iterateList("edges", "edge");
			pml.set("dash_gap", 0, "edge", true);
	
		    S = conn.createStatement();
	
		    subQuery = "SELECT "+edgeInfo[0]+" AS i_cid, "+edgeInfo[2]+" AS i_num, "+edgeInfo[1]+" AS j_cid, "+edgeInfo[3]+" AS j_num FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+")";
		    // left join the edges with themselves and choose the ones with null on the match
		    query = "SELECT A.* FROM "+
			"("+subQuery+") AS A LEFT JOIN ("+subQuery+") AS B ON "+
			"(A.i_cid = B.j_cid AND A.i_num = B.j_num AND A.j_cid = B.i_cid AND A.j_num = B.i_num) "+
			"WHERE (B.i_cid IS NULL) ORDER BY i_cid, i_num, j_cid, j_num;";
		    R = S.executeQuery(query);
		    while (R.next()) { 
		    	edge = "e."+R.getString("i_cid")+"."+R.getInt("i_num")+"."+R.getString("j_cid")+"."+R.getInt("j_num");
		    	pml.set("dash_gap", edgeGap, edge, false);
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
    
    /*
     * sets the color for all nodes included in the @param nodes list
     */       
    public void setNodesColorUniform(String nodes, String color) {

    	pml.iterateList(nodes, "node");
		pml.setNodeColor(color, "node", true);
    	//java String method passing either a nodes selection string or just "nodes"
    	//pymol String method passing "nodes" (the pymol selection nodes is used)
    	//pml.setNodeColor(nodeColor, nodes, false);
    	//pymol String method passing "nodes" (the pymol string nodes is used) 
    	//pml.setNodeColor(nodeColor, nodes, true);

    }
    
    /*
     * sets the color for all nodes based on the nodeColorMethod
     */
    public void setNodesColor() {

		Statement S;
		ResultSet R;
		String query = "", node = "", nodeSetId = "";
		double colWeight = 0;
		double[] nodeCurColRange = new double[2];
		double[] nodeRGB = new double[3];
		double[] nodeHSV = new double[3];
		int i = -1;
		
	
		try { 
	
		    S = conn.createStatement();
	
		    if (nodeColorMethod.equals("chain")) {
	
		    	if (chainNodeColors == null) {
					for (String chain : chains) {
					    setNodesColorUniform("nodes_cid"+chain, "color_cid"+chain);
					}		    		
		    	} else {
		    		for (String chain : chains) {
		    			pml.createColor("color_nodes_cid"+chain, chainNodeColors[i%chainNodeColors.length]);
					    setNodesColorUniform("nodes_cid"+chain, "color_nodes_cid"+chain);
					}	
		    	} 
	
		    } else if ((!nodeColorMethod.equals("uniform")) && (nodeColDiscr)) {
	
		    	/* The solution to select residues with specific value for the 
		    	 * nodeColorMethod field, while traversing resultset R, wouldn't 
		    	 * work because equalities with floating-point numbers are dangerous.
		    	 */
				query = "SELECT DISTINCT "+nodeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+nodeColorMethod+";";
				R = S.executeQuery(query);
				while (R.next()) {
				    i++;
				    pml.createColor("color_nodes_set"+i, setNodeColors[i%setNodeColors.length]);
				    nodeSetIds.put(new Long(Double.doubleToLongBits(R.getDouble(nodeColorMethod))), "nodes_set"+i);
				    pml.initList("nodes_set"+i);
				}
				R.close();
				
				/*
				 * I could skip the nodeSet.X lists and just set the color for the specific residue
				 * but it is easier if someone wants for example to highlight one nodeSet later.
				 * Also, the idea to use the "iterate list" solution here is mainly to avoid many
				 * lines in the script. 
				 */
				query = "SELECT "+nodeInfo[0]+", "+nodeInfo[1]+", "+nodeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+nodeInfo[0]+", "+nodeInfo[1]+";";
				R = S.executeQuery(query);
				while (R.next()) {
				    node = "n."+R.getString(nodeInfo[0])+"."+R.getInt(nodeInfo[1]);
				    nodeSetId = nodeSetIds.get(new Long(Double.doubleToLongBits(R.getDouble(nodeColorMethod))));
				    nodeSet.put(node, nodeSetId);
				    pml.appendList(nodeSetId, node);
				}
				R.close();
		
				for (String curNodeSet : nodeSetIds.values()) {
					setNodesColorUniform(curNodeSet, "color_"+curNodeSet);
				}
	
		    } else if ((!nodeColorMethod.equals("uniform")) && (!nodeColDiscr)) {
	
				nodeCurColRange = utils4DB.getRange(conn, dbInfo[0]+"."+dbInfo[2], nodeColorMethod, nodeGraphSel+" AND ("+graphInfo[5]+")");
				
				query = "SELECT "+nodeInfo[0]+", "+nodeInfo[1]+", "+nodeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[2]+" WHERE ("+nodeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+nodeInfo[0]+", "+nodeInfo[1]+";";
				R = S.executeQuery(query);
				while (R.next()) {
				    node = "n."+R.getString(nodeInfo[0])+"."+R.getInt(nodeInfo[1]);
				    colWeight = R.getDouble(nodeColorMethod);
				    //HSV (Hue, Saturation, Value) model- Hue is the color type defined by the rescaled db field value
				    nodeHSV[0] = rescale(colWeight, nodeCurColRange, new double[] {0, 240}, false);
				    nodeHSV[1] = 1;
				    nodeHSV[2] = 1;
				    //convert the hsv color to rgb
				    nodeRGB = hsvToRgb(nodeHSV);
				    pml.createColor("color_"+node, nodeRGB);
				    pml.setNodeColor("color_"+node, node, false);
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
    
    /*
     * sets the color for all edges included in the @param edges list
     */    
    public void setEdgesColorUniform(String edges, String color) {
    	
    	pml.iterateList(edges, "edge");
		pml.setColor(color, "edge", true);
    	//java String method passing either a edges selection string or just "edges"
    	//pymol String method passing "edges" (the pymol selection edges is used)
    	//pml.setColor(edgeColor, edges, false);
    	//pymol String method passing "edges" (the pymol string edges is used) 
    	//pml.setColor(edgeColor, edges, true);
		
    }

    /*
     * sets the color for all edges based on the edgeColorMethod
     */
    public void setEdgesColor() {

		Statement S;
		ResultSet R;
		String query = "", edge = "", edgeSetId = "", node = "";
		double colWeight = 0;
		double[] edgeCurColRange = new double[2];
		double[] edgeRGB = new double[3];
		double[] edgeHSV = new double[3];
		int i = -1;
		
	
		try { 
	
		    S = conn.createStatement();
	
		    if (edgeColorMethod.equals("chain")) {
		    	
		    	if (chainEdgeColors == null) {		    	
					for (String chain : chains) {
						setEdgesColorUniform("edges_cid"+chain, "color_cid"+chain);
					}
		    	} else {
		    		for (String chain : chains) {
		    			pml.createColor("color_edges_cid"+chain, chainEdgeColors[i%chainEdgeColors.length]);
					    setEdgesColorUniform("edges_cid"+chain, "color_edges_cid"+chain);
					}	
		    	} 
		    	
		    } else if (edgeColorMethod.equals("node")) {
	
				if (nodeColorMethod.equals("uniform")) {
					setEdgesColorUniform("edges", nodeColor);
				} else if (nodeColorMethod.equals("chain")) {
					if (chainNodeColors == null) {
						for (String chain : chains) {
							setEdgesColorUniform("edges_cid"+chain, "color_cid"+chain);
						}
					} else {
						for (String chain : chains) {
							setEdgesColorUniform("edges_cid"+chain, "color_nodes_cid"+chain);
						}
					}
				} else if (nodeColDiscr) {
				    query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
				    R = S.executeQuery(query);
				    while (R.next()) {
						node = "n."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2]);
						edge = "e."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2])+"."+R.getString(edgeInfo[1])+"."+R.getInt(edgeInfo[3]);
						pml.setColor("color_"+nodeSet.get(node), edge, false);
				    }
				    R.close();
				} else if (!nodeColDiscr) {
				    query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
				    R = S.executeQuery(query);
				    while (R.next()) {
						node = "n."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2]);
						edge = "e."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2])+"."+R.getString(edgeInfo[1])+"."+R.getInt(edgeInfo[3]);
						pml.setColor("color_"+node, edge, false);
				    }
				    R.close();
				}
	
		    } else if ((!edgeColorMethod.equals("uniform")) && (edgeColDiscr)) {
	
				query = "SELECT DISTINCT "+edgeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeColorMethod+";";
				R = S.executeQuery(query);
				while (R.next()) {
				    i++;
				    pml.createColor("color_edges_set"+i, setEdgeColors[i%setEdgeColors.length]);
				    edgeSetIds.put(new Long(Double.doubleToLongBits(R.getDouble(edgeColorMethod))), "edges_set"+i);
				    pml.initList("edges_set"+i);
				}
				R.close();
		
				query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+", "+edgeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
				R = S.executeQuery(query);
				while (R.next()) {
				    edge = "e."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2])+"."+R.getString(edgeInfo[1])+"."+R.getInt(edgeInfo[3]);
				    edgeSetId = edgeSetIds.get(new Long(Double.doubleToLongBits(R.getDouble(edgeColorMethod))));
				    edgeSet.put(edge, edgeSetId);
				    pml.appendList(edgeSetId, edge);
				}
				R.close();
		
				for (String curEdgeSet : edgeSetIds.values()) {
					setEdgesColorUniform(curEdgeSet, "color_"+curEdgeSet);
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
				    pml.setColor("color_"+edge, edge, false);
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

    /*
     * create cgo edges (cylinders)
     */
    public void createCgoEdges() {
    	
		Statement S;
		ResultSet R;
		String query = "", edge = "", i_cid = "", curCgoEdgeColor = "", queryFields = "", edgeSetId = "", node = "";
		double colWeight = 0, curCgoEdgeSize = 0;
		double[] edgeCurSizeRange = new double[2];
		double[] edgeCurColRange = new double[2];
		double[] edgeRGB = new double[3];
		double[] edgeHSV = new double[3];
		int i = 0;
		
		try { 
	
			pml.initList("edges");
			
			S = conn.createStatement();	
			
			if (edgeColorMethod.equals("chain") && (chainEdgeColors != null)) {
				for (String chain : chains) {
	    			pml.createColor("color_edges_cid"+chain, chainEdgeColors[i%chainEdgeColors.length]);
				}
			}			
			if (!edgeSizeMethod.equals("uniform")) {
				edgeCurSizeRange = utils4DB.getRange(conn, dbInfo[0]+"."+dbInfo[3], edgeSizeMethod, edgeGraphSel+" AND ("+graphInfo[5]+")");
				queryFields += ", "+edgeSizeMethod;
			}
			
			if ( !(edgeColorMethod.equals("uniform") || edgeColorMethod.equals("chain") | edgeColorMethod.equals("node")) ) {
				if (edgeColDiscr) {	
					query = "SELECT DISTINCT "+edgeColorMethod+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeColorMethod+";";
					R = S.executeQuery(query);
					while (R.next()) {
					    i++;
					    pml.createColor("color_edges_set"+i, setEdgeColors[i%setEdgeColors.length]);
					    edgeSetIds.put(new Long(Double.doubleToLongBits(R.getDouble(edgeColorMethod))), "edges_set"+i);
					    pml.initList("edges_set"+i);
					}
					R.close();
			
				} else if (!edgeColDiscr) {
					edgeCurColRange = utils4DB.getRange(conn, dbInfo[0]+"."+dbInfo[3], edgeColorMethod, edgeGraphSel+" AND ("+graphInfo[5]+")");
				}
				queryFields += ", "+edgeColorMethod;
			}   
			
			query = "SELECT "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+queryFields+" FROM "+dbInfo[0]+"."+dbInfo[3]+" WHERE ("+edgeGraphSel+") AND ("+graphInfo[5]+") ORDER BY "+edgeInfo[0]+", "+edgeInfo[2]+", "+edgeInfo[1]+", "+edgeInfo[3]+";";
			R = S.executeQuery(query);
		    while (R.next()) { 
				i_cid = R.getString(edgeInfo[0]);
				
				if (!edgeSizeMethod.equals("uniform")) {
					curCgoEdgeSize = rescale(R.getDouble(edgeSizeMethod), edgeCurSizeRange, cgoEdgeSizeRange, edgeSizeRev);
				} else {
					curCgoEdgeSize = cgoEdgeSize;
				}
				
				if (edgeColorMethod.equals("uniform")) {
					curCgoEdgeColor = edgeColor;
				} else if (edgeColorMethod.equals("chain")){
					curCgoEdgeColor = (chainEdgeColors == null)?"color_cid"+i_cid:"color_edges_cid"+i_cid;					
				} else if (edgeColorMethod.equals("node")){
					if (nodeColorMethod.equals("uniform")) {
						curCgoEdgeColor = nodeColor;
					} else if (nodeColorMethod.equals("chain")) {
						curCgoEdgeColor = (chainNodeColors == null)?"color_cid"+i_cid:"color_nodes_cid"+i_cid;
					} else if (nodeColDiscr) {
						node = "n."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2]);
						curCgoEdgeColor = "color_"+nodeSet.get(node);
					} else if (!nodeColDiscr) {
						node = "n."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2]);
						curCgoEdgeColor = "color_"+node;
					}
				} else if ((!edgeColorMethod.equals("uniform")) && (edgeColDiscr)) {
					edge = "e."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2])+"."+R.getString(edgeInfo[1])+"."+R.getInt(edgeInfo[3]);
				    edgeSetId = edgeSetIds.get(new Long(Double.doubleToLongBits(R.getDouble(edgeColorMethod))));
				    edgeSet.put(edge, edgeSetId);
				    pml.appendList(edgeSetId, edge);
				    curCgoEdgeColor = "color_"+edgeSetId;
				} else if ((!edgeColorMethod.equals("uniform")) && (!edgeColDiscr)) {
					edge = "e."+R.getString(edgeInfo[0])+"."+R.getInt(edgeInfo[2])+"."+R.getString(edgeInfo[1])+"."+R.getInt(edgeInfo[3]);
				    colWeight = R.getDouble(edgeColorMethod);
				    edgeHSV[0] = rescale(colWeight, edgeCurColRange, new double[] {0, 240}, false);
				    edgeHSV[1] = 1;
				    edgeHSV[2] = 1;
				    edgeRGB = hsvToRgb(edgeHSV);
				    pml.createColor("color_"+edge, edgeRGB);
				    curCgoEdgeColor = "color_"+edge;
				}
				
				edge = pml.addCgoEdge(i_cid, R.getInt(edgeInfo[2]), R.getString(edgeInfo[1]), R.getInt(edgeInfo[3]), curCgoEdgeColor, cgoEdgeGap, cgoEdgeLength, curCgoEdgeSize, directed, curCgoEdgeColor, msdsd);
				pml.appendList("edges_cid"+i_cid, edge);
				pml.appendList("edges", edge);		
		    } 
		    
		    R.close();
		    S.close();
	
		} // end try
		catch (SQLException E) {
		    System.out.println("SQLException: " + E.getMessage());
		    System.out.println("SQLState:     " + E.getSQLState());
		    System.out.println("VendorError:  " + E.getErrorCode());
		} // end catch 
		
    }

    private double rescale(double value, double[] curRange, double[] newRange, boolean rev) {

		if (!rev) {
		    return (newRange[0] + ((value-curRange[0])/(curRange[1]-curRange[0]))*(newRange[1]-newRange[0]));
		} else {
		    return (newRange[1] - ((value-curRange[0])/(curRange[1]-curRange[0]))*(newRange[1]-newRange[0]));
		}
		
    }
    
    /*
     * convert HSV (Hue, Saturation, Value) model to RGB (red-green-blue)
     */
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

    /*
     * private method to initialise defaults HashMap
     */
    private void popDefaults() {
    
    	defaults.put("nodeSize", new Integer(1));
    	defaults.put("edgeSize", new Integer(2));
    	defaults.put("cgoEdgeSize", new Integer(3));
    	defaults.put("nodeSizeRange", new Integer(4));
    	defaults.put("edgeSizeRange", new Integer(5));
    	defaults.put("cgoEdgeSizeRange", new Integer(6));
    	defaults.put("edgeGap", new Integer(7));
    	defaults.put("nodeTransp", new Integer(8));
    	defaults.put("surfTransp", new Integer(9));
    	defaults.put("chainColors", new Integer(10));
    	defaults.put("setNodeColors", new Integer(11));
    	defaults.put("setEdgeColors", new Integer(12));    	
    	
    }
    
    /*
     * define which attributes to be reset to their default values
     */
    public void setDefaults(String[] attrs) {
    	for(String attr : attrs) {
    		switch(defaults.get(attr).intValue()) {
    			case 1: 
    				nodeSize = 0.46;
    				break;
    			case 2:
    				edgeSize = 3.65;
    				break;
    			case 3:
    				cgoEdgeSize = 0.3;
    				break;
    			case 4:
    				nodeSizeRange = new double[] {0.2, 0.8};
    				break;
    			case 5:
    				edgeSizeRange = new double[] {0.3, 7.0};
    				break;
    			case 6:
    				cgoEdgeSizeRange = new double[] {0.2, 0.5};
    				break;
    			case 7:
    				edgeGap = 0.25;
    				break;
    			case 8:
    				nodeTransp = 0.6;
    				break;
    			case 9:
    				surfTransp = 0.6;
    				break;
    			case 10:
    				chainColors = new String[] {"light_grey", "green", "red", "blue", "yellow", "violet", "cyan", "salmon", "lime", "pink", "slate", "magenta", "orange", "marine", "olive", "purple", "teal", "forest", "firebrick", "chocolate", "wheat", "white", "grey"};
    				break;
    			case 11:
    				setNodeColors = new String[] {"green", "red", "blue", "yellow", "violet", "cyan", "salmon", "lime", "pink", "slate", "magenta", "orange", "marine", "olive", "purple", "teal", "forest", "firebrick", "chocolate", "wheat", "white", "grey"};
    				break;
    			case 12:
    				setEdgeColors = new String[] {"green", "red", "blue", "yellow", "violet", "cyan", "salmon", "lime", "pink", "slate", "magenta", "orange", "marine", "olive", "purple", "teal", "forest", "firebrick", "chocolate", "wheat", "white", "grey"};
    				break;
    			default:
    				System.out.println("Attribute "+attr+" doesn't exist!");    		  		
    		}    		
    	}
    }
    
} // end of class Graph2Pml

package proteinstructure;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Locale;

import tools.MySQLConnection;

/**
 * A light weight graph class consisting of an EdgeSet and a NodeSet.
 * Name chosen because Graph is already taken for a ResidueInteractionGraph.
 * @author Henning Stehr
 */
public class NodesAndEdges {

	/*------------------------------ constants ------------------------------*/
	public static final String DEFAULT_GRAPH_TABLE = "graph";
	public static final String DEFAULT_NODE_TABLE = "node";
	public static final String DEFAULT_EDGE_TABLE = "edge";
	public static final String DEFAULT_GRAPH_ID_COL = "graph_id";
	public static final String DEFAULT_GRAPH_COMMENT_COL = "comment";
	public static final String DEFAULT_NODE_NUM_COL = "num";
	public static final String DEFAULT_EDGE_I_COL = "i_num";
	public static final String DEFAULT_EDGE_J_COL = "j_num";
	
	/*--------------------------- member variables --------------------------*/
	protected NodeSet nodes;
	protected EdgeSet edges;
	protected String comment;
	
	/*----------------------------- constructors ----------------------------*/
	/**
	 * Create a new graph given a node set and an edge set
	 */
	public NodesAndEdges(NodeSet n, EdgeSet e) {
		nodes = n;
		edges = e;
	}
	
	/**
	 * Create a new graph with empty node- and edge sets
	 */
	public NodesAndEdges() {
		nodes = new NodeSet();
		edges = new EdgeSet();
		comment = null;
	}
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Creates tables for a simple graph database with default column names. The database has to exist.
	 */
	public static void createGraphDatabase(MySQLConnection conn, String db) {
		try {
		// graph table
		String query = "CREATE TABLE " + db + "." + DEFAULT_GRAPH_TABLE 
		             + "(" + DEFAULT_GRAPH_ID_COL + " int(10) unsigned primary key auto_increment, " 
		             +       DEFAULT_GRAPH_COMMENT_COL + " varchar(1000)" 
		             + ")";
		conn.executeSql(query);

		// node table
		query = "CREATE TABLE " + db + "." + DEFAULT_NODE_TABLE 
        + "(" + DEFAULT_GRAPH_ID_COL + " int(10) unsigned, " 
        +       DEFAULT_NODE_NUM_COL + " int(10) unsigned" 
        + ")";		
		conn.executeSql(query);
		
		// edge table
		query = "CREATE TABLE " + db + "." + DEFAULT_EDGE_TABLE 
        + "(" + DEFAULT_GRAPH_ID_COL + " int(10) unsigned, " 
        +       DEFAULT_EDGE_I_COL + " int(10) unsigned, "
        +       DEFAULT_EDGE_J_COL + " int(10) unsigned"         
        + ")";		
		conn.executeSql(query);
		
		} catch(SQLException e) {
			System.err.println("Error while creating graph database: " + e.getMessage());
		}
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Returns the comment
	 * @return the comment
	 */
	public String getComment() {
		return comment;
	}

	/**
	 * Sets the comment
	 * @param comment the comment to set
	 */
	public void setComment(String comment) {
		this.comment = comment;
	}
		
	/**
	 * Returns the node set of the graph
	 * @return the node set
	 */
	public NodeSet getNodes() { return nodes; }
	
	/**
	 * Returns the edge set of the graph
	 * @return the edge set
	 */
	public EdgeSet getEdges() { return edges; }
	
	/**
	 * Returns the number of nodes
	 * @return the number of nodes
	 */
	public int getNumNodes() { return nodes.size(); }
	
	/**
	 * Returns the number of edges
	 * @return the number of edges
	 */
	public int getNumEdges() { return edges.size(); }

	/**
	 * Returns a deep copy of this graph
	 * @return a deep copy of this graph
	 */
	public NodesAndEdges copy() {
		return new NodesAndEdges(nodes.copy(), edges.copy());
	}
	
	/**
	 * Returns true iff this graph equals the given other graph. Graphs are considered
	 * equals if they contain equal sets of nodes and edges. Nodes are equal if
	 * the node num is equal and edges are equal if they connect the same node nums.
	 * @param other the graph to compare to
	 * @return true if the graphs are equal, false otherwise
	 */
	public boolean equals(NodesAndEdges other) {
		if(!this.nodes.equals(other.nodes)) return false;
		if(!this.edges.equals(other.edges)) return false;
		return true;
	}
	
	/**
	 * Check whether all edges map to nodes from the node set.
	 * @return true if graph is consistent, false otherwise
	 */
	public boolean isConsistent() {
		boolean consistent = true;
		for(Edge e:edges) {
			if(!nodes.contains(e.i) || !nodes.contains(e.j)) {
				consistent = false;
			}
		}
		return consistent;
	}
	
	/**
	 * Output information about the graph to stdout
	 * @param printNodes whether to output all nodes
	 * @param printEdges whether to output all edges
	 */
	public void printGraph(boolean printNodes, boolean printEdges) {
		System.out.println("Nodes: " + nodes.size() + " Edges: " + edges.size());
		if(printNodes) {
			System.out.println("Nodes:");
			for(Node node:nodes) {
				System.out.print(node + " ");
			}
			System.out.println();
		}
		if(printEdges) {
			System.out.println("Edges:");
			for(Edge e:edges) {
				System.out.print("(" + e.i + "," + e.j + "," + e.weight + ") ");
			}
			System.out.println();
		}
	}
	
	/**
	 * Write graph to database.
	 * @param conn Database connection
	 */
	public void writeToDb(MySQLConnection conn, String db, String graphTable, String graphIdxCol, String nodeTable, String nodeIdxCol, String edgeTable, String edgeICol, String edgeJCol, int graphId) {
		try {
			// graph id
			String query = "SELECT COUNT(*) FROM " + db + "." + graphTable + " WHERE " + graphIdxCol + " = " + graphId;
			if(conn.getIntFromDb(query) > 0) {
				System.err.println("Graph id " + graphId + " already exists in table " + db + "." + graphTable);
				return;
			}
			query = "INSERT INTO " + db + "." + graphTable + " (" + graphIdxCol + ") VALUES (?)";
			PreparedStatement p = conn.getConnectionObject().prepareStatement(query);
			p.setInt(1, graphId);
			p.executeUpdate();
			p.close();

			// nodes
			query = "INSERT INTO " + db + "." + nodeTable + " (" + graphIdxCol + "," + nodeIdxCol + ") VALUES (?,?)";
			p = conn.getConnectionObject().prepareStatement(query);
			for(Node n:nodes) {
				p.setInt(1, graphId);
				p.setInt(2, n.num);
				p.addBatch();
			}
			p.executeBatch();
			p.close();

			// edges
			query = "INSERT INTO " + db + "." + edgeTable + " (" + graphIdxCol + "," + edgeICol + "," + edgeJCol + ") VALUES (?,?,?)";		
			p = conn.getConnectionObject().prepareStatement(query);
			for(Edge e:edges) {
				p.setInt(1, graphId);
				p.setInt(2, e.i);
				p.setInt(3, e.j);
				p.addBatch();
			}
			p.executeBatch();
			p.close();

		} catch (SQLException e) {
			System.err.println("Error while writing graph to database: " + e.getMessage());
		}
	}
	
	/**
	 * Write graph to database using default parameters for tables and columns.
	 * @param conn database connection
	 **/
	public void writeToDb(MySQLConnection conn, String db, int graph_id) {
		writeToDb(conn, db, DEFAULT_GRAPH_TABLE, DEFAULT_GRAPH_ID_COL, DEFAULT_NODE_TABLE, DEFAULT_NODE_NUM_COL, DEFAULT_EDGE_TABLE, DEFAULT_EDGE_I_COL, DEFAULT_EDGE_J_COL, graph_id);
	}
	
	/**
	 * Read graph from database. Nodes and edges read from the database are _added_ to the current graph.
	 * @param conn database connection
	 */
	public void loadFromDb(MySQLConnection conn, String db, String graphTable, String graphIdxCol, String nodeTable, String nodeIdxCol, String edgeTable, String edgeICol, String edgeJCol, int graphId) {
		try {
			// check graphId
			String query = "SELECT * FROM " + db + "." + graphTable + " WHERE " + graphIdxCol + " = ?";
			PreparedStatement p = conn.getConnectionObject().prepareStatement(query);
			p.setInt(1, graphId);
			ResultSet rs = p.executeQuery();
			if(!rs.next()) System.err.println("Graph id not found in database");	// TODO: throw exception
			rs.close();
			p.close();
			
			// read nodes
			query = "SELECT " + nodeIdxCol + " FROM " + db + "." + nodeTable + " WHERE " + graphIdxCol + " = ?";
			p = conn.getConnectionObject().prepareStatement(query);
			p.setInt(1, graphId);
			rs = p.executeQuery();
			while(rs.next()) {
				int num = rs.getInt(1);
				nodes.add(new Node(num));
			}
			rs.close();
			p.close();
			
			// read edges
			query = "SELECT " + edgeICol + "," + edgeJCol + " FROM " + db + "." + edgeTable + " WHERE " + graphIdxCol + " = ?";
			p = conn.getConnectionObject().prepareStatement(query);
			p.setInt(1, graphId);
			rs = p.executeQuery();
			while(rs.next()) {
				int i = rs.getInt(1);
				int j = rs.getInt(2);
				edges.add(new Edge(i,j));
			}
			rs.close();
			p.close();			
		} catch (SQLException e) {
			System.err.println("Error while reading graph from database: " + e.getMessage());			
		}		
	}

	/**
	 * Read graph from database using default parameters for tables and columns.
	 * @param conn database connection
	 **/
	public void loadFromDb(MySQLConnection conn, String db, int graph_id) {
		loadFromDb(conn, db, DEFAULT_GRAPH_TABLE, DEFAULT_GRAPH_ID_COL, DEFAULT_NODE_TABLE, DEFAULT_NODE_NUM_COL, DEFAULT_EDGE_TABLE, DEFAULT_EDGE_I_COL, DEFAULT_EDGE_J_COL, graph_id);
	}	

	/**
	 * Output graph as a GDL file for aiSee
	 * @param fileName output file
	 * @param labelProperty the property of the nodes which is to be used as the label in aiSee or null if none
	 * @param colorProperty the property of the nodes which contains a color string for aiSee or null if none
	 */
	public void writeGdlFile(String fileName, String labelProperty, String colorProperty) {
		try {
			PrintWriter fileOut = new PrintWriter(new FileOutputStream(fileName));
			// ---- init graph ----
			fileOut.println("graph: {");
			fileOut.println("layoutalgorithm:forcedir");
			fileOut.println("     attraction:80");
			fileOut.println("      repulsion:100");
			fileOut.println("        gravity:0.5");
			
			fileOut.println();
			
			// Specifying new colours
			//fileOut.println("colorentry 33: 230 230 230");

			//Setting the background colour
			fileOut.println("color:white");			
			fileOut.println();
			
			// ---- write nodes ----
			
			// node parameters
			fileOut.println("node.shape:circle");		// circle, rhomb
			fileOut.println("node.width:40");
			fileOut.println("node.height:40");
		    fileOut.println("node.color:white");			
		    fileOut.println("node.textcolor:black");
			fileOut.println();
			
			// write nodes
			for(Node node:nodes) {
				fileOut.println("node: { ");
				fileOut.println("\ttitle: \""+ node +"\"");
				if(labelProperty != null) {
					if(node.getProperty(labelProperty) != null) {
						String label = node.getProperty(labelProperty);
						fileOut.println("\tlabel: \""+ label +"\"");
					}
				}
				if(colorProperty != null) {
					if(node.getProperty(colorProperty) != null) {
						String color = node.getProperty(colorProperty);
						fileOut.println("\tcolor: "+ color);
					}
				}
				fileOut.println("      }");				
			}

			// ---- write edges ----
		    
			// edge parameters
			fileOut.println();
			fileOut.println("edge.arrowstyle:none");		// none, solid
			fileOut.println("edge.linestyle:continuous");	// invisible, dotted, continuous
		    fileOut.println("edge.arrowsize:20");
			fileOut.println("edge.thickness:2");
			fileOut.println("edge.color:blue");
			fileOut.println();
			
			// write edges
			for(Edge e:edges) {
				fileOut.println("edge: { source: \"" + e.i + "\" target: \"" + e.j + "\"}");
			}
			
			// ---- finish graph ----
			fileOut.println();
			fileOut.println("}   // end Graph");
			fileOut.close();
			System.out.println("GDL output written to " + fileName);
			
		} catch(IOException e) {
			System.err.println("Error writing to file " + fileName);
			e.printStackTrace();
		}
	}
	
	/*--------------------------------- main --------------------------------*/
	
	/* Run some test on the database methods */
	public static void testDatabaseMethods() {
		String db = "henning_test";
		int graphId = 1;
		NodesAndEdges graph = new NodesAndEdges();
		MySQLConnection conn = null;
		try {
			conn = new MySQLConnection("white", "stehr", "nieve");
		} catch(SQLException e) {
			System.err.println("Error while opening database connection: " + e.getMessage());
		}
		// load graph from rig database
		System.out.println("Loading graph #" + graphId + " from pdb_reps_graph...");		
		graph.loadFromDb(conn, "pdb_reps_graph", "single_model_graph", "graph_id", "single_model_node", "num", "single_model_edge", "i_num", "j_num", graphId);		
		
		// create new simple graph db
		System.out.println("Creating graph database " + db + "...");
		NodesAndEdges.createGraphDatabase(conn, db);
		
		// write graph to simple graph db
		System.out.println("Writing graph # " + graphId + " to db " + db + "...");
		graph.writeToDb(conn, db, graphId);

		// relading graph from simple graph db
		System.out.println("Loading graph # " + graphId + " from db " + db + "...");
		NodesAndEdges graph2 = new NodesAndEdges();
		graph2.loadFromDb(conn, db, graphId);
		graph2.printGraph(false, false);
		if(!graph.equals(graph2)) System.err.println("Error: Graph written to db differs from graph loaded from db.");
		conn.close();
		double f = 0.3333333;
		int num = 5;
		System.out.println(String.format(Locale.GERMANY, "n=%d s=%5.2f", num, f));
		System.out.println("done.");		
	}
	
	/* Run test on the gdl output */
	public static void testGdlOutput() {
		String gdlFileName = "/project/StruPPi/henning/test.gdl";
		int graphId = 1;
		NodesAndEdges graph = new NodesAndEdges();
		MySQLConnection conn = null;
		try {
			conn = new MySQLConnection("white", "stehr", "nieve");
		} catch(SQLException e) {
			System.err.println("Error while opening database connection: " + e.getMessage());
		}
		// load graph from rig database
		System.out.println("Loading graph #" + graphId + " from pdb_reps_graph...");		
		graph.loadFromDb(conn, "pdb_reps_graph", "single_model_graph", "graph_id", "single_model_node", "num", "single_model_edge", "i_num", "j_num", graphId);
		graph.printGraph(false, false);
		graph.writeGdlFile(gdlFileName, null, null);
		conn.close();		
		System.out.println("done.");
	}
	
	/* Testing NodesAndEdges */
	public static void main(String[] args) {
		testDatabaseMethods();
		//		if(args.length < 6) {
//			System.out.println("Usage: NodesAndEdges <GraphId> <Database> <GraphTable> <NodeTable> <EdgeTable> <GdlOutfile>");
//			System.out.println("Reads a graph from a database and writes a gdl file for aiSee");
//			System.exit(1);
//		}
//		int graphId = Integer.parseInt(args[0]);
//		String db = args[1];
//		String graphTable = args[2];
//		String nodeTable = args[3];
//		String edgeTable = args[4];
//		String gdlFile = args[5];
//		
//		MySQLConnection conn = null;
//		try {
//			conn = new MySQLConnection("white", System.getProperty("user.name"), "nieve");
//		} catch(SQLException e) {
//			System.err.println("Error while opening database connection: " + e.getMessage());
//			System.exit(1);
//		}
//		
//		NodesAndEdges graph = new NodesAndEdges();
//		System.out.println("Loading graph #" + graphId + "from db " + db + "...");
//		graph.loadFromDb(conn, db, graphTable, DEFAULT_GRAPH_ID_COL, nodeTable, DEFAULT_NODE_NUM_COL, edgeTable, DEFAULT_EDGE_I_COL, DEFAULT_EDGE_J_COL, graphId);
//		graph.writeGdlFile(gdlFile, null, null);
//		conn.close();
//		System.out.println("done.");
	}
}

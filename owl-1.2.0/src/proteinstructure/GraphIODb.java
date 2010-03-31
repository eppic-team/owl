package proteinstructure;

import java.io.IOException;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.TreeMap;

import org.apache.commons.collections15.Factory;
import org.apache.commons.collections15.Transformer;

import tools.MySQLConnection;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

public class GraphIODb<V,E> {

	/*------------------------------ constants ------------------------------*/
	public static final String DEFAULT_GRAPH_TABLE = "graph";
	public static final String DEFAULT_NODE_TABLE = "node";
	public static final String DEFAULT_EDGE_TABLE = "edge";
	public static final String DEFAULT_GRAPH_ID_COL = "graph_id";
	public static final String DEFAULT_GRAPH_COMMENT_COL = "comment";
	public static final String DEFAULT_NODE_NUM_COL = "num";
	public static final String DEFAULT_EDGE_I_COL = "i_num";
	public static final String DEFAULT_EDGE_J_COL = "j_num";

	/*------------------------------- fields --------------------------------*/

	String db;
	MySQLConnection conn;
	String graphTable;
	String graphIdxCol;
	String nodeTable;
	String nodeIdxCol;
	String edgeTable;
	String edgeICol;
	String edgeJCol;

	/*----------------------------- constructors ----------------------------*/

	/**
	 * Constructs a GraphDbIO using default table/column names
	 * @param db
	 * @param conn
	 */
	public GraphIODb(String db, MySQLConnection conn) {
		super();
		this.db = db;
		this.conn = conn;
		this.graphTable = DEFAULT_GRAPH_TABLE;
		this.graphIdxCol = DEFAULT_GRAPH_ID_COL;
		this.nodeTable = DEFAULT_NODE_TABLE;
		this.nodeIdxCol = DEFAULT_NODE_NUM_COL;
		this.edgeTable = DEFAULT_EDGE_TABLE;
		this.edgeICol = DEFAULT_EDGE_I_COL;
		this.edgeJCol = DEFAULT_EDGE_J_COL;
	}

	/**
	 * Constructs a GraphDbIO using the provided table/column names
	 * @param db
	 * @param conn
	 * @param graphTable
	 * @param graphIdxCol
	 * @param nodeTable
	 * @param nodeIdxCol
	 * @param edgeTable
	 * @param edgeICol
	 * @param edgeJCol
	 */
	public GraphIODb(String db, MySQLConnection conn, String graphTable, String graphIdxCol, String nodeTable,
			String nodeIdxCol, String edgeTable, String edgeICol,
			String edgeJCol) {
		super();
		this.db = db;
		this.conn = conn;
		this.graphTable = graphTable;
		this.graphIdxCol = graphIdxCol;
		this.nodeTable = nodeTable;
		this.nodeIdxCol = nodeIdxCol;
		this.edgeTable = edgeTable;
		this.edgeICol = edgeICol;
		this.edgeJCol = edgeJCol;
	}

	/*---------------------------- static methods ---------------------------*/

	/**
	 * Creates tables for a simple graph database with default column names. The
	 * database has to exist.
	 */
	public static void createGraphDatabase(MySQLConnection conn, String db) {
		try {
			// graph table
			String query = "CREATE TABLE " + db + "." + DEFAULT_GRAPH_TABLE
					+ "(" + DEFAULT_GRAPH_ID_COL
					+ " int(10) unsigned primary key auto_increment, "
					+ DEFAULT_GRAPH_COMMENT_COL + " varchar(1000)" + ")";
			conn.executeSql(query);

			// node table
			query = "CREATE TABLE " + db + "." + DEFAULT_NODE_TABLE + "("
					+ DEFAULT_GRAPH_ID_COL + " int(10) unsigned, "
					+ DEFAULT_NODE_NUM_COL + " int(10) unsigned" + ")";
			conn.executeSql(query);

			// edge table
			query = "CREATE TABLE " + db + "." + DEFAULT_EDGE_TABLE + "("
					+ DEFAULT_GRAPH_ID_COL + " int(10) unsigned, "
					+ DEFAULT_EDGE_I_COL + " int(10) unsigned, "
					+ DEFAULT_EDGE_J_COL + " int(10) unsigned" + ")";
			conn.executeSql(query);

		} catch (SQLException e) {
			System.err.println("Error while creating graph database: "
					+ e.getMessage());
		}
	}

	/*---------------------------- public methods ---------------------------*/

	/**
	 * Write graph to database.
	 * 
	 * @param graph the graph to write out to database
	 * @param nodeSerialTransformer a Transformer that extracts a serial from the Vertex object
	 * @param graphId the graph_id that the graph will get in the database
	 */
	public void writeToDb(Graph<V,E> graph, Transformer<V, Integer> nodeSerialTransformer, int graphId) throws SQLException {
		//TODO eventually we could paas also a Transformer<E, Double> edgeWeightTransformer and write weights as well
		
		// graph id
		String query = "SELECT COUNT(*) FROM " + db + "." + graphTable
		+ " WHERE " + graphIdxCol + " = " + graphId;
		if (conn.getIntFromDb(query) > 0) {
			System.err.println("Graph id " + graphId+ " already exists in table " + db + "." + graphTable);
			return;
		}
		query = "INSERT INTO " + db + "." + graphTable + " (" + graphIdxCol
		+ ") VALUES (?)";
		PreparedStatement p = conn.getConnectionObject().prepareStatement(query);
		p.setInt(1, graphId);
		p.executeUpdate();
		p.close();

		// nodes
		query = "INSERT INTO " + db + "." + nodeTable + " (" + graphIdxCol
		+ "," + nodeIdxCol + ") VALUES (?,?)";
		p = conn.getConnectionObject().prepareStatement(query);
		for (V n : graph.getVertices()) {
			p.setInt(1, graphId);
			p.setInt(2, nodeSerialTransformer.transform(n));
			p.addBatch();
		}
		p.executeBatch();
		p.close();

		// edges
		query = "INSERT INTO " + db + "." + edgeTable + " (" + graphIdxCol
		+ "," + edgeICol + "," + edgeJCol + ") VALUES (?,?,?)";
		p = conn.getConnectionObject().prepareStatement(query);
		for (E e : graph.getEdges()) {
			Pair<V> pair = graph.getEndpoints(e);
			p.setInt(1, graphId);
			p.setInt(2,  nodeSerialTransformer.transform(pair.getFirst()));
			p.setInt(3,  nodeSerialTransformer.transform(pair.getSecond()));
			p.addBatch();
		}
		p.executeBatch();
		p.close();

	}

	/**
	 * Reads graph from database. 
	 * 
	 * @param graphId the graph_id to read from database
	 * @param serialNodeTransformer a Transformer from a serial into a Vertex object
	 * @param edgeFactory an Edge factory
	 */
	public Graph<V,E> loadFromDb(int graphId, Transformer<Integer, V> serialNodeTransformer, Factory<E> edgeFactory) throws SQLException{
		//TODO eventually we could read weights too (would also require a Transformer<Double,E> weightEdgeTransformer

		Graph<V,E> graph = new SparseGraph<V, E>();
		TreeMap<Integer,V> serials2vertices = new TreeMap<Integer, V>();
		// check graphId
		String query = "SELECT * FROM " + db + "." + graphTable + " WHERE " + graphIdxCol + " = ?";
		PreparedStatement p = conn.getConnectionObject().prepareStatement(query);
		p.setInt(1, graphId);
		ResultSet rs = p.executeQuery();
		if (!rs.next())	System.err.println("Graph id not found in database"); // TODO: throw exception
		rs.close();
		p.close();

		// read nodes
		query = "SELECT " + nodeIdxCol + " FROM " + db + "." + nodeTable + " WHERE " + graphIdxCol + " = ?";
		p = conn.getConnectionObject().prepareStatement(query);
		p.setInt(1, graphId);
		rs = p.executeQuery();
		while (rs.next()) {
			int num = rs.getInt(1);
			V vertex = serialNodeTransformer.transform(num);
			serials2vertices.put(num, vertex);
			graph.addVertex(vertex);
		}
		rs.close();
		p.close();

		// read edges
		query = "SELECT " + edgeICol + "," + edgeJCol + " FROM " + db + "."	+ edgeTable + " WHERE " + graphIdxCol + " = ?";
		p = conn.getConnectionObject().prepareStatement(query);
		p.setInt(1, graphId);
		rs = p.executeQuery();
		while (rs.next()) {
			int i = rs.getInt(1);
			int j = rs.getInt(2);
			graph.addEdge(edgeFactory.create(), serials2vertices.get(i), serials2vertices.get(j), EdgeType.UNDIRECTED);
		}
		rs.close();
		p.close();
		return graph;
	}

	
	// tester
	public static void main(String[] args) throws SQLException, IOException {
		int graphId = 1;
		String gdlFileName = "test.gdl";

		MySQLConnection conn = null;
		try {
			conn = new MySQLConnection("talyn", "duarte", "nieve");
		} catch(SQLException e) {
			System.err.println("Error while opening database connection: " + e.getMessage());
			System.exit(1);
		}
		// load graph from rig database
		GraphIODb<RIGNode, RIGEdge> dbloader = 
			new GraphIODb<RIGNode, RIGEdge>("pdb_reps_graph", conn, "single_model_graph", "graph_id", "single_model_node", "num", "single_model_edge", "i_num", "j_num");
		System.out.println("Loading graph #" + graphId + " from pdb_reps_graph...");		
		Graph<RIGNode,RIGEdge> graph = dbloader.loadFromDb(graphId,
				new Transformer<Integer,RIGNode>() {
			public RIGNode transform(Integer arg0) {
				return new RIGNode(arg0);
			}
		},
		new Factory<RIGEdge>() {
			public RIGEdge create() {
				return new RIGEdge();
			}
		}
		);
		
		System.out.println("node count: "+graph.getVertexCount());
		System.out.println("edge count: "+graph.getEdgeCount());
		
		
		
		new GraphIOGDLFile<RIGNode,RIGEdge>().writeGdlFile(graph, gdlFileName, 
		new Transformer<RIGNode,Integer>() {
			public Integer transform(RIGNode node) {
				return node.getResidueSerial();
			}
		}, 
		new Transformer<RIGNode,String>() {
			public String transform(RIGNode node) {
				return node.getResidueType();
			}
		},
		new Transformer<RIGNode,String>() {
			public String transform(RIGNode node) {
				return null;
			}
		});
		
		conn.close();		
		System.out.println("done.");
	}
}

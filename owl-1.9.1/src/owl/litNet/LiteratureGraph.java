package owl.litNet;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.*;
import java.util.*;

import owl.core.util.MySQLConnection;

import edu.uci.ics.jung.graph.DirectedSparseGraph;

public class LiteratureGraph extends DirectedSparseGraph<Article, Reference> {	

	/*------------------------------ constants ------------------------------*/
	private static final long serialVersionUID = 1L;
	private static final String litDb = "literature";
	private static final String tblArticle = "article_view";
	private static final String tblReference = "reference";

	/*--------------------------- member variables --------------------------*/
	private Map<Integer, Article> idx2vertex;
	
	/*----------------------------- constructors ----------------------------*/
	
	public LiteratureGraph() {
		super();
		idx2vertex = new HashMap<Integer, Article>();
	}
	
	/*-------------------------- idxGraph interface -------------------------*/
	public boolean addVertex(Article a) {
		idx2vertex.put(a.getIdx(), a);
		return super.addVertex(a);
	}
	
	public boolean removeVertex(Article a) {
		idx2vertex.remove(a.getIdx());
		return super.removeVertex(a);
	}
	
	/**
	 * @return the vertex with the given idx or null if no such vertex exists
	 */
	public Article getVertexByIdx(int idx) {
		Article a = idx2vertex.get(idx);
		return a;
	}
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Loads the literature network from the database.
	 */
	public void readFromDatabase(MySQLConnection conn) throws SQLException {
		
		String query;
		Statement s;
		ResultSet rs;
		
		// get articles
		query = String.format("SELECT idx, author, year, journal, title FROM %s.%s", litDb, tblArticle);
		s = conn.createStatement();
		rs = s.executeQuery(query);
		while(rs.next()) {
			int id = rs.getInt(1);
			String author = rs.getString(2);
			int year = rs.getInt(3);
			String journal = rs.getString(4);
			String title = rs.getString(5);
			this.addVertex(new Article(id, author, year, journal, title));
		}
		rs.close();
		s.close();
		
		// get references
		query = String.format("SELECT idx, ref_idx FROM %s.%s", litDb, tblReference);
		s = conn.createStatement();
		rs = s.executeQuery(query);
		while(rs.next()) {
			int idx = rs.getInt(1);
			int refIdx = rs.getInt(2);
			Article a = this.getVertexByIdx(idx);
			Article r = this.getVertexByIdx(refIdx);
			if(a != null && r != null) {
				this.addEdge(new Reference(), a, r);
			}
		}
		rs.close();
		s.close();
	}
	
	/**
	 * Writes the literature network to a file in WilmaScope .xwg format.
	 * @param outFile the file to write to
	 */
	public void writeToWilmaFile(File outFile) throws IOException {
		PrintWriter out = new PrintWriter(outFile);
		
		// header
		out.println("<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>");
		out.println("<!DOCTYPE WilmaGraph SYSTEM \"WilmaGraph.dtd\">");
		out.println("<WilmaGraph>");
		out.println("<Cluster>");
		
		// layout
		out.println("<LayoutEngineType Name=\"Force Directed\">");
		out.println("\t<Property Key=\"Spring\" Value=\"5.6\"/>");
		out.println("\t<Property Key=\"DirectedField\" Value=\"9.5\"/>");
		out.println("\t<Property Key=\"Origin\" Value=\"1.0\"/>");
		out.println("\t<Property Key=\"Repulsion\" Value=\"0.5\"/>");
		out.println("\t<Property Key=\"VelocityAttenuation\" Value=\"0.02\"/>");
		out.println("</LayoutEngineType>");
		
		// node properties
		double radius = 0.2;
		Color color = null;
		//Map<Article, Color> node2col = new TreeMap<Article,Color>();
		String label;
		
		// nodes
		for(Article n:this.getVertices()) {
			color = assignColorByYear(n.getYear());
			//if(node2col.containsKey(n)) node2col.get(n);
			radius = 0.05 + this.getPredecessorCount(n)/100.0;
			label = "";
			if(this.getInEdges(n).size() > 20) {
				//label = n.toString();
			}
			label = n.toString();
			out.println(getWilmaNodeStr("N" + Integer.toString(n.getIdx()), label, color, radius));
		}
		
		// edges
		for(Reference e:this.getEdges()) {
			int start = this.getEndpoints(e).getFirst().getIdx();
			int end = this.getEndpoints(e).getSecond().getIdx();
			out.println(getWilmaEdgeStr("N" + Integer.toString(start), "N" + Integer.toString(end)));
		}       
		
		// footer
		out.println("</Cluster>");
		out.println("</WilmaGraph>");
		out.close();
	}
	
	public String getWilmaNodeStr(String id, String label, Color color, double radius) {
		String viewType = "DefaultNodeView";
		//String viewType = "LabelOnly";		
		float[] rgb = null;
		if(color != null) {
			rgb = color.getRGBColorComponents(null);
		}
		String nodeStr = String.format("\t<Node ID=\"%s\">", id);
		nodeStr += String.format("<ViewType Name=\"%s\">", viewType);
		if(label != null) {
			nodeStr += String.format("<Property Key=\"Label\" Value=\"%s\"/>", label);
		}
		if(color != null) {
			nodeStr += String.format("<Property Key=\"Colour\" Value=\"%f %f %f\"/>", rgb[0], rgb[1], rgb[2]);
		}
		if(radius > 0) {
			nodeStr += String.format("<Property Key=\"Radius\" Value=\"%f\"/>", radius);			
		}
		nodeStr += String.format("</ViewType>");		
		nodeStr += String.format("</Node>");
		return nodeStr;
	}
	
	/**
	 * Delete all vertices with less than n incoming edges.
	 * @param n the minimum number of incoming for a vertex to remain in the graph.
	 */
	public void deleteLowConnectedVertices(int n) {
		LinkedList<Article> deleteSet = new LinkedList<Article>();
		// mark all low connected nodes for deletion
		for(Article a:this.getVertices()) {
			Collection<Reference> inEdges = this.getInEdges(a);
			if(inEdges.size() < n) deleteSet.add(a);
		}
		// delete
		for(Article a: deleteSet) {
			this.removeVertex(a);
		}
	}
	
	/**
	 * Delete all vertices with no adjacent edges.
	 */
	public void deleteNonConnectedVertices() {
		LinkedList<Article> deleteSet = new LinkedList<Article>();
		// mark all low connected nodes for deletion
		for(Article a:this.getVertices()) {
			Collection<Reference> inEdges = this.getInEdges(a);
			Collection<Reference> outEdges = this.getOutEdges(a);
			if(inEdges.size() == 0 && outEdges.size() == 0) deleteSet.add(a);
		}
		// delete
		for(Article a: deleteSet) {
			this.removeVertex(a);
		}
	}	
	
	public String getWilmaEdgeStr(String start, String end) {
		String edgeStr = String.format("\t<Edge StartID=\"%s\" EndID=\"%s\">", start, end);
		edgeStr += String.format("<ViewType Name=\"Arrow\">");
		//edgeStr += String.format("<ViewType Name=\"LineEdge\">");
		edgeStr += String.format("</ViewType>");		
		edgeStr += String.format("</Edge>\n", start, end);
		return edgeStr;
	}
	
	/**
	 * Prints number of nodes and edges to stdout.
	 */
	public void info() {
		System.out.printf("Nodes=%d, Edges=%d\n", this.getVertexCount(), this.getEdgeCount());
	}
	
	/**
	 * Given a year of a publication assign a node color.
	 * @param year the publication year
	 * @return
	 */
	public Color assignColorByYear(int year) {
//		int min = 1970;
//		int max = 2007;
//		year = Math.max(min, year);	// everything before 1970 looks like 1970
//		year = Math.min(max, year);	// everything newer than 2007 looks like 2007
//		float scale = 1.0f * (year - min) / (max - min);	// 1=newest, 0=oldest
//		float r = 1-scale;
//		float g = scale;
//		float b = 0;
//		return new Color(r,g,b);
		if(year < 1970) {
			return Color.DARK_GRAY;
		}
		if(year < 1980) {
			return Color.RED;
		}
		if(year < 1990) {
			return Color.ORANGE;
		}
		if(year < 2000) {
			return Color.cyan;
		}
		return Color.green;
		
	}
	
	/*--------------------------------- main --------------------------------*/
	
	public static void main(String[] args) throws SQLException, IOException {
		
		String wilmaFileName = "/project/StruPPi/incoming/litnet.xwg";
		
		String dbHost = "talyn";
		String dbUser = "stehr";
		String dbPwd = "nieve";
		MySQLConnection conn = new MySQLConnection(dbHost, dbUser, dbPwd);
		
		LiteratureGraph litNet = new LiteratureGraph();
		System.out.println("Loading literature network from database...");
		litNet.readFromDatabase(conn);
		System.out.println("done.");
		litNet.info();
		System.out.println("Deleting articles with less than 10 incoming citations.");
		litNet.deleteLowConnectedVertices(10);
		litNet.info();
		System.out.println("Deleting non-connected articles");
		litNet.deleteNonConnectedVertices();		
		litNet.info();
		System.out.println("Writing to file " + wilmaFileName);
		litNet.writeToWilmaFile(new File(wilmaFileName));
		
		conn.close();
		
	}
	
}

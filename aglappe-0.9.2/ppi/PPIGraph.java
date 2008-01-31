package ppi;

import java.awt.Color;
import java.io.*;
import java.sql.SQLException;
import java.util.*;
import java.util.regex.*;

import edu.uci.ics.jung.graph.SparseGraph;
import tools.MySQLConnection;

public class PPIGraph extends SparseGraph<PPINode, PPIEdge> {
	
	private static final long serialVersionUID = 1L;

	public PPIGraph() {
		super();
	}
	
	public boolean containsNodeWithIdx(int idx) {
		for(PPINode n:this.getVertices()) {
			if(n.getIdx() == idx) return true;
		}
		return false;
	}
	
	/**
	 * @return the first vertex with the given idx or null if no such vertex exists
	 */
	public PPINode getNodeByIdx(int idx) {
		for(PPINode n:this.getVertices()) {
			if(n.getIdx() == idx) return n;
		}
		return null;
	}
	
	public void loadNodesFromTabSepFile(File inFile) throws FileNotFoundException, IOException {
		BufferedReader in = new BufferedReader(new FileReader(inFile));
		String line;
		int lineCount = 0;
		int lastIdx = 0;
		PPINode currentNode = null;
		while((line = in.readLine()) != null) {
			lineCount++;
			String[] fields = line.split("\t");
			if(fields.length < 2) {
				System.err.println("Expected " + 4 + " but found " + fields.length + " fields in line " + lineCount + ": " + line);
			} else {
				int idx = Integer.parseInt(fields[0]);
				String name = fields[1];
				String desc = "";
				if(fields.length >= 4) desc = fields[3];
				String gos = "";
				if(fields.length >= 3) gos = fields[2];
				if(idx > lastIdx) {
					// new node
					currentNode = new PPINode(idx, name, desc);
					this.addVertex(currentNode);
					lastIdx = idx;
				}
				// parse GO string
				String regex = "GO; GO:(\\d+); ([CFP]):(\\w*)"; // GO; GO:0005471; F:ATP:ADP antiporter activity; IDA
				Pattern p = Pattern.compile(regex);
				Matcher m = p.matcher(gos);
				if(m.find()) {
					int num = Integer.parseInt(m.group(1));
					String domStr = m.group(2);
					String goDesc = m.group(3);
					GOAnnotation.Domain domain = null;
					if(domStr.equals("F")) domain = GOAnnotation.Domain.F; else
					if(domStr.equals("C")) domain = GOAnnotation.Domain.C; else
					if(domStr.equals("P")) domain = GOAnnotation.Domain.P; else
						System.err.println("Wrong GO domain " + domStr + " in line " + lineCount + ": " + line);
					GOAnnotation newGo = new GOAnnotation(num, domain, goDesc);
					currentNode.addGOAnnotation(newGo);
				} else {
					if(gos.length() > 0) {
						System.err.println("Error parsing GO string in line " + lineCount + ": " + gos);
					}
				}
			}
		}
		in.close();
	}
	
	public void loadEdgesFromTabSepFile(File inFile) throws FileNotFoundException, IOException {
		BufferedReader in = new BufferedReader(new FileReader(inFile));
		String line;
		int lineCount = 0;
		while((line = in.readLine()) != null) {
			lineCount++;
			String[] fields = line.split("\t");
			if(fields.length < 2) {
				System.err.println("Error reading edge from line " + lineCount + ": " + line);
			} else {
				int start = Integer.parseInt(fields[0]);
				int end = Integer.parseInt(fields[1]);
				PPINode startNode = this.getNodeByIdx(start);
				PPINode endNode = this.getNodeByIdx(end);
				if(startNode == null || endNode == null) {
					System.err.printf("Could not add edge (%d,%d): Node not found.\n", start, end);
				} else {
					this.addEdge(new PPIEdge(), startNode, endNode);
				}
			}
		}
		in.close();
	}
	
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
		Map<PPINode, Color> node2col = assignRandomColors();
		
		// nodes
		for(PPINode n:this.getVertices()) {
			Color c = node2col.containsKey(n)?node2col.get(n):Color.black;
			out.println(getWilmaNodeStr(Integer.toString(n.getIdx()), n.getProteinName(), c, radius));
		}
		
		// edges
		for(PPIEdge e:this.getEdges()) {
			int start = this.getEndpoints(e).getFirst().getIdx();
			int end = this.getEndpoints(e).getSecond().getIdx();
			out.println(getWilmaEdgeStr(Integer.toString(start), Integer.toString(end)));
		}       
		
		// footer
		out.println("</Cluster>");
		out.println("</WilmaGraph>");
		out.close();
	}
	
	public String getWilmaNodeStr(String id, String label, Color color, double radius) {
		float[] rgb = color.getRGBColorComponents(null);
		String nodeStr = String.format("\t<Node ID=\"%s\">", id);
		nodeStr += String.format("<ViewType Name=\"DefaultNodeView\">");
		if(label != null) {
			nodeStr += String.format("<Property Key=\"Label\" Value=\"%s\"/>", label);
		}
		if(color != null) {
			nodeStr += String.format("<Property Key=\"Colour\" Value=\"%f %f %f\"/>", rgb[0], rgb[1], rgb[2]);
		}
		if(radius > 0) {
			nodeStr += String.format("<Property Key=\"Radius\" Value=\"%s\"/>", "0.2");			
		}
		nodeStr += String.format("</ViewType>");		
		nodeStr += String.format("</Node>");
		return nodeStr;
	}
	
	public String getWilmaEdgeStr(String start, String end) {
		String edgeStr = String.format("\t<Edge StartID=\"%s\" EndID=\"%s\"/>", start, end);
		return edgeStr;
	}
	
	/**
	 * Assign random colors to different cellular component annotations.
	 */
	public Map<PPINode, Color> assignRandomColors() {
		Map<Integer, Color> go2col = new HashMap<Integer, Color>();	  // assigns colors to GO numbers
		Map<PPINode, Color> node2col = new HashMap<PPINode, Color>(); // assigns colors to nodes
		for(PPINode n : this.getVertices()) {
			ArrayList<GOAnnotation> gos = n.getGoAnnotations();
			for(GOAnnotation go : gos) {
				if(go.getDomain() == GOAnnotation.Domain.C) {		// consider only cellular component annotations
					Color currentColor;
					if(!go2col.containsKey(go.getNumber())) {
						float r = (float) Math.random();
						float g = (float) Math.random();
						float b = (float) Math.random();
						currentColor = new Color(r,g,b);
						go2col.put(go.getNumber(), currentColor);
					} else {
						currentColor = go2col.get(go.getNumber());
					}
					if(!node2col.containsKey(n)) {
						node2col.put(n, currentColor);
					}
				}
			}
		}
		System.out.println("Distinct GO cellular localization terms: " + go2col.size());
		System.out.println("Nodes with color assignments: " + node2col.size());
		return node2col;
	}
	
	public void loadFromDb(MySQLConnection conn, String dbName, int graphId) {
		
	}
	
	/**
	 * Writes the graph to a database.
	 * @param conn
	 * @param dbName
	 * @return the last inserted graph id or 0 on error
	 */
	public int writeToDb(MySQLConnection conn, String dbName) {
		return 1;
	}
	
	public void createDatabase(MySQLConnection conn, String dbName) {
		
	}
	
	public void info() {
		System.out.printf("Nodes=%d, Edges=%d\n", this.getVertexCount(), this.getEdgeCount());

	}
	
	public void printNodes() {
		for(PPINode n:this.getVertices()) {
			System.out.println(n);
		}		
	}
	
	public void printEdges() {
		for(PPIEdge e:this.getEdges()) {
			int start = this.getEndpoints(e).getFirst().getIdx();
			int end = this.getEndpoints(e).getSecond().getIdx();
			System.out.printf("Edge:(%d,%d)\n", start, end);
		}
	}
	
	/*--------------------------------- main --------------------------------*/
	
	public static void main(String[] args) throws SQLException, FileNotFoundException, IOException {
		String nodeFileName = "/project/StruPPi/incoming/fyi_key.txt";
		String edgeFileName = "/project/StruPPi/incoming/fyi.txt";
		String wilmaFileName = "/project/StruPPi/incoming/fyi.xwg";
		
		String dbHost = "white";
		String dbUser = "stehr";
		String dbPwd = "nieve";
		MySQLConnection conn = new MySQLConnection(dbHost, dbUser, dbPwd);
		//String dbName = "henning_ppi_fyi";
		
		// load from file
		PPIGraph ppiGraph = new PPIGraph();
		ppiGraph.loadNodesFromTabSepFile(new File(nodeFileName));
		ppiGraph.loadEdgesFromTabSepFile(new File(edgeFileName));
		ppiGraph.info();
		System.out.println("Writing to " + wilmaFileName);
		ppiGraph.writeToWilmaFile(new File(wilmaFileName));
		System.out.println("done.");
		
		// clean up
		conn.close();
	}
	
}

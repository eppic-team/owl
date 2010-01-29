package proteinstructure;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;

public class GraphIOGDLFile<V,E> {

	
	/**
	 * Output graph as a GDL file for aiSee
	 * @param graph the graph to write out to GDL file
	 * @param fileName output GDL file
	 * @param titleTransformer a Transformer of Vertex to an integer titles for aiSee nodes
	 * @param labelTransformer a Transformer of Vertex to a label string
	 * @param colorTransformer a Transformer of Vertex to an aiSee color string
	 * @throws IOException 
	 */
	public void writeGdlFile(Graph<V,E> graph, String fileName, Transformer<V,Integer> titleTransformer, Transformer<V,String> labelTransformer, Transformer<V,String> colorTransformer) throws IOException {

		PrintWriter fileOut = new PrintWriter(new FileOutputStream(fileName));
		// ---- init graph ----
		fileOut.println("graph: {");
		fileOut.println("layoutalgorithm:forcedir");
		fileOut.println("     attraction:80");
		fileOut.println("      repulsion:300");
		fileOut.println("        gravity:5.0");

		fileOut.println();

		// Specifying new colours
		//fileOut.println("colorentry 33: 230 230 230");

		//Setting the background colour
		fileOut.println("color:white");			
		fileOut.println();

		// ---- write nodes ----

		// node parameters
		fileOut.println("node.shape:circle");		// circle, rhomb
		fileOut.println("node.width:60");
		fileOut.println("node.height:60");
		fileOut.println("node.color:white");			
		fileOut.println("node.textcolor:black");
		fileOut.println();

		// write nodes
		for(V node:graph.getVertices()) {
			fileOut.println("node: { ");
			fileOut.println("\ttitle: \""+ titleTransformer.transform(node) +"\"");

			if(labelTransformer.transform(node) != null) {
				String label = labelTransformer.transform(node);
				fileOut.println("\tlabel: \""+ label +"\"");
			}


			if(colorTransformer.transform(node) != null) {
				String color = colorTransformer.transform(node);
				fileOut.println("\tcolor: "+ color);
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
		for(E e:graph.getEdges()) {
			Pair<V> pair = graph.getEndpoints(e);
			fileOut.println("edge: { source: \"" + titleTransformer.transform(pair.getFirst()) + "\" target: \"" + titleTransformer.transform(pair.getSecond()) + "\"}");
		}

		// ---- finish graph ----
		fileOut.println();
		fileOut.println("}   // end Graph");
		fileOut.close();
		System.out.println("GDL output written to " + fileName);

	}
	
	
	// tester
	public static void main(String[] args) throws SQLException, PdbCodeNotFoundError, PdbLoadError, IOException {
		Pdb pdb = new PdbasePdb("1bxy");
		pdb.load("A");
		RIGraph graph = pdb.getRIGraph("Ca", 8);
		
		String gdlfile = "test.gdl";
		
		GraphIOGDLFile<RIGNode, RIGEdge> gdlfileIO = new GraphIOGDLFile<RIGNode, RIGEdge>();
		
		gdlfileIO.writeGdlFile(graph, gdlfile, 
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
		
		
	}
	
}

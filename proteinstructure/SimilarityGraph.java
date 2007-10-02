package proteinstructure;
import java.io.*;
import java.util.*;
import tools.*;
/**
 * This class represents a similarity graph where nodes correspond to objects and edges to similarities (or distances) between them. 
 * It provides provides methods for
 * - reading a (complete) graph from a similarity matrix,
 * - removing distant edges as a very simple way of clustering
 * - visualizing the resulting graph in aiSee (by outputting an aiSee script)
 * The graph knows whether it is a similarity graph (default) or a distance graph which changes the behaviour
 * of some functions.
 * @author Henning Stehr
 */
public class SimilarityGraph extends NodesAndEdges {

	/*------------------------------ constants ------------------------------*/
	public static final boolean FORCE_UNDIRECTED = true; // if true, lower half of matrix is ignored
	public static final int SIMILARITY = 1;				 // interpretation of edge weights
	public static final int DISTANCE = -1;				 // see scoreType below
	
	/*--------------------------- member variables --------------------------*/
	private int scoreType;	// can be DISTANCE or SIMILARITY (default)
	
	/*----------------------------- constructors ----------------------------*/
	/**
	 * Creates an empty similarity graph
	 */
	public SimilarityGraph() {
		super();
		this.scoreType = SIMILARITY;
	}
	
	/**
	 * Creates an empty similarity or distance graph (depending on the given parameter).
	 * Use SimilarityGraph.SIMILARITY or SimilarityGraph.DISTANCE for scoreType.
	 */
	public SimilarityGraph(int scoreType) {
		this.scoreType = scoreType;
	}
	
	/*---------------------------- public methods ---------------------------*/
	/**
	 * Returns true if the graph is a similarity graph, false otherwise (i.e. it is a distance graph)
	 */
	 public boolean isSimilarityGraph() {
	 	return scoreType == SIMILARITY;
	 }
	
	/**
	 * Read a simple distance matrix text file
	 * @param fileName
	 */
	public void readMatrix(String fileName) {
		try {
		BufferedReader in = new BufferedReader(new FileReader(fileName));
		String firstLine = in.readLine();
		String[] nodeLabels = firstLine.split("\t");
		Vector<Integer> nodeVector = new Vector<Integer>();
		for(String nodeStr:nodeLabels) {
			int node = Integer.parseInt(nodeStr.trim());
			nodeVector.add(node);
			nodes.add(node);
		}
		String line;
		int lineNum = 0;
		while((line = in.readLine()) != null) {
			String[] distStrings = line.split("\t");
			int colNum = 0;
			for(String distStr:distStrings) {
				int i = nodeVector.get(colNum);
				int j = nodeVector.get(lineNum);
				if(!FORCE_UNDIRECTED || i > j) {
					double dist = Double.parseDouble(distStr.trim());
					Edge e = new Edge(i,j,dist);
					edges.add(e);
				}
				colNum++;
			}
			lineNum++;
		}
		
		} catch(FileNotFoundException e) {
			e.printStackTrace();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Read a distance matrix output file of the Maxcluster program
	 * @param fileName
	 */
	public void readMaxclusterMatrix(String fileName) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			System.out.println("Reading from file " + fileName);
			String line;
			while((line = in.readLine()) != null) {
				if(line.startsWith("PDB")) {
					String[] t = line.split("\\s+");
					int node = Integer.parseInt(t[2]);
					nodes.add(node);
				} else
				if(line.startsWith("DIST")) {
					String[] t = line.split("\\s+");
					int i = Integer.parseInt(t[2]);
					int j = Integer.parseInt(t[3]);
					double dist = Double.parseDouble(t[4]);
					edges.add(new Edge(i,j,dist));
				}
			}
		} catch(FileNotFoundException e) {
			e.printStackTrace();
		} catch(IOException e) {
			e.printStackTrace();
		}		
	}
	
	/**
	 * Not implemented: Write graph to database
	 * @param conn
	 */
	public void writeToDb(MySQLConnection conn) {
		
	}
	
	/**
	 * Not implemented: Read graph from a database
	 * @param conn
	 */
	public void loadFromDb(MySQLConnection conn) {
		
	}
	
	/**
	 * Output graph as a GDL file for aiSee
	 * @param fileName output file
	 * @param bestModel a workaround before ranking is working to highlight the best model in the graph
	 * TODO: Move to NodesAndEdges
	 */
	public void writeGdlFile(String fileName, int bestModel) {
		try {
			PrintWriter fileOut = new PrintWriter(new FileOutputStream(fileName));
			// ---- init graph ----
			fileOut.println("graph: {");
			fileOut.println("layoutalgorithm:forcedir \t attraction:80 \t repulsion:100 \t gravity:0.5");
			fileOut.println();
			
			// Specifying new colours
			fileOut.println("colorentry 33: 230 230 230");

			//Setting the background colour
			fileOut.println("color:white");			
			fileOut.println();
			
			// ---- write nodes ----
			
			// node parameters
			fileOut.println("node.shape :circle");		// circle, rhomb
			fileOut.println("node.width :40");
			fileOut.println("node.height :40");
		    fileOut.println("node.color:white");			
		    fileOut.println("node.textcolor:black");
			fileOut.println();
			
			// write nodes
			for(int node:nodes) {
				if(node == bestModel) {
				    fileOut.println("node.color:red");					
				} else {
					fileOut.println("node.color:white");
				}
				fileOut.println("node: { title: \""+ node +"\"");
				fileOut.println("\tlabel: \""+ node +"\"\t}");
				
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
				fileOut.println("edge: { source: \"" + e.i + "\" \t target: \"" + e.j + "\"}");
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
	
	/**
	 * Filters the edgeset to remove edges for distant objects.
	 * If graph is a similarity graph, removes all edges with similarity values below the cutoff.
	 * If graph is a distance graph, removes all edges with distance values above the cutoff.
	 * @param cutoff the filter cutoff
	 */
	public void removeDistantEdges(double cutoff) {
		if(this.scoreType == SIMILARITY) {
			this.edges.filterEdges(cutoff, -1);
		} else
			if(this.scoreType == DISTANCE) {
				this.edges.filterEdges(cutoff, +1);
			} else {
				System.err.println("Unknown score type " + this.scoreType);
			}
	}
		
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		if(args.length < 5) {
			System.out.println("Usage: SimilarityGraph <matrixFile> <outputGdlFile> <cutoff> <s[imilarity]/d[istance]> <numberOfBestModel>");
			System.exit(1);
		}
		String matrixFile = args[0];
		String gdlFile = args[1];
		String cutoffStr = args[2];
		String typeStr = args[3];
		String bestModelStr = args[4];
		double filterCutoff = Double.parseDouble(cutoffStr);
		int bestModel = Integer.parseInt(bestModelStr);
		int scoreType = 0;
		if(typeStr.startsWith("s")) {
			scoreType = SIMILARITY;
		} else
		if(typeStr.startsWith("d")) {
			scoreType = DISTANCE;
		} else {
			System.err.println("Parameter 4 can be only 's' (for similarity) or 'd' (for distance");
			System.exit(1);
		}
		
//		scoreType = SIMILARITY;
//		matrixFile = "/project/StruPPi/henning/projects/ensemb_vis/1fzy_NtoN_rmsd.txt";
//		gdlFile = "/project/StruPPi/henning/projects/ensemb_vis/1fzy_100_rmsd_1_4.gdl";
//		filterCutoff = 1.4;
		
		SimilarityGraph myGraph = new SimilarityGraph(scoreType);

		myGraph.readMaxclusterMatrix(matrixFile);
		System.out.println("Before filtering:");
		myGraph.printGraph(false, false);
		myGraph.removeDistantEdges(filterCutoff);
		System.out.println("After filtering:");
		myGraph.printGraph(false, false);
		System.out.println("Writing GDL file...");
		myGraph.writeGdlFile(gdlFile, bestModel);				
		System.out.println("done.");
	}

}

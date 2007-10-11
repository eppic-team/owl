package proteinstructure;
import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
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
			nodes.add(new Node(node));
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
	public void readFromMaxclusterMatrix(String fileName) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			System.out.println("Reading from file " + fileName);
			String line;
			while((line = in.readLine()) != null) {
				if(line.startsWith("PDB")) {
					String[] t = line.split("\\s+");
					int node = Integer.parseInt(t[2]);
					String fName = t[3];
					Node newNode = new Node(node);
					newNode.setProperty("filename", fName);
					nodes.add(newNode);
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
	 * Reads a maxCluster ranking file and sets the node properties 'rank' and 'score'.
	 * The nodes are identified by the property 'filename', so they have to be set,
	 * for example by calling readFromMaxclsuterMatrix with a corresponding matrix file.
	 */
	public void readLabelsFromMaxClusterRanking(String rankingFile) {
		try {
			BufferedReader in = new BufferedReader(new FileReader(rankingFile));
			System.out.println("Reading from file " + rankingFile);
			String line;
			while((line = in.readLine()) != null) {
				Pattern p = Pattern.compile("INFO  : +(\\d+). (.+) vs. (.+)   (GDT|RMSD)= *(\\d+\\.\\d+)");
				Matcher m = p.matcher(line);
				if(m.find()) {
					String rank = m.group(1);
					//String targetName = m.group(2);
					String fileName = m.group(3);
					//String scoreType = m.group(4);
					String score = m.group(5);
					//System.out.printf("target=%s file=%s rank=%s type=%s score=%s\n", targetName, fileName, rank, scoreType, score);
					Node n = nodes.getNodeByProperty("filename",fileName);
					if(n == null) {
						System.err.println("No node found with property filename=" + fileName);
					} else {
						n.setProperty("rank", rank);
						n.setProperty("score", score);
					}
				}
			}
		} catch(FileNotFoundException e) {
			e.printStackTrace();
		} catch(IOException e) {
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
			System.out.println("Usage: SimilarityGraph <matrixFile> <rankingFile> <outputGdlFile> <cutoff> <s[imilarity]/d[istance]>");
			System.exit(1);
		}
		String matrixFile = args[0];
		String rankingFile = args[1];
		String gdlFile = args[2];
		String cutoffStr = args[3];
		String typeStr = args[4];
		double filterCutoff = Double.parseDouble(cutoffStr);
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

		myGraph.readFromMaxclusterMatrix(matrixFile);
		myGraph.readLabelsFromMaxClusterRanking(rankingFile);
		Node bestModel = myGraph.getNodes().getNodeByProperty("rank", "1");
		if(bestModel != null) {
			System.out.println("Best model: " + bestModel.num + " with " + bestModel.getProperty("score"));
			bestModel.setProperty("color", "red");
		}
		System.out.println("Before edge filtering:");
		myGraph.printGraph(false, false);
		myGraph.removeDistantEdges(filterCutoff);
		System.out.println("After edge filtering:");
		myGraph.printGraph(false, false);
		System.out.println("Writing GDL file...");
		myGraph.writeGdlFile(gdlFile, null, "color");			// labelProperty: 'rank', 'score' or null; colorProperty: 'color' or null		
		System.out.println("done.");
	}

}

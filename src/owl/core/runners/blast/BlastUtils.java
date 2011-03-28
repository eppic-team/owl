package owl.core.runners.blast;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.*;

import org.apache.commons.collections15.Transformer;

import owl.core.runners.MaxClusterRunner;
import owl.core.structure.*;
import owl.core.structure.graphs.GraphIOGDLFile;
import owl.core.util.MySQLConnection;

import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.util.Pair;


/**
 * A collection of little tools related to Blast and Processing Blast output.
 * @author stehr
 *
 */
public class BlastUtils {

	private static final File tempDir = new File("/tmp/");
	private static final String maxClusterExecutable = "/project/StruPPi/bin/maxcluster";
	
	private static final double IDENTITY_SCORE_RMSD = 0.0;
	
	private static String MATRIX_VIS_SCRIPT = "/project/StruPPi/CASP8/scripts/plot_simmatrix.sh";
	
	//private static final double similarityGraphRmsdCutoff = 2.0;
	//private static final double similarityGraphGdtCutoff = 50.0;
	
	class DoubleWrapper {
		public double val;	
		DoubleWrapper(double v) {
			this.val = v;
		}
	}
	
	/**
	 * Uses a perl script to render a PNG image of the given blast output file.
	 * @param blastOutputFile
	 * @param imgFile
	 */
	public static void renderBlast(File blastOutputFile, File imgFile) throws IOException {
		String renderScript = "/project/StruPPi/CASP8/scripts/render_blast.pl";
		String cmdLine = String.format("%s %s > %s", renderScript, blastOutputFile.getAbsolutePath(), imgFile.getAbsolutePath());
		Runtime.getRuntime().exec(cmdLine);
	}
	
	/**
	 * Helper function used by writeClusterGraph writing a similarity matrix to a file.
	 * @param sparseMatrix
	 * @param labels
	 * @param outFile
	 */
	private static void writeMatrixToFile(HashMap<Pair<Integer>, Double> sparseMatrix, String[] labels, File outFile) throws IOException {
				
		// initialize temp matrix
		int mSize = labels.length;	// assuming matrix size equals number of labels
		//System.out.println(mSize);
		double[][] tmpMatrix = new double[mSize][mSize];
		for (int i = 0; i < tmpMatrix.length; i++) {
			tmpMatrix[i][i] = IDENTITY_SCORE_RMSD;
		}
		for(Pair<Integer> pair:sparseMatrix.keySet()) {
			int i = pair.getFirst() - 1;
			int j = pair.getSecond() - 1;
			if(i >= mSize || j >= mSize) {
				System.err.printf("Error: Entry (%d,%d) in matrix exceeds number of labels (%d).", i+1,j+1,mSize);
			} else {
				if(sparseMatrix.containsKey(pair)) {
					tmpMatrix[i][j] = sparseMatrix.get(pair);
					tmpMatrix[j][i] = sparseMatrix.get(pair);
				}
			}
		}
		
		// write matrix to file
		PrintWriter out = new PrintWriter(outFile);
		String sep = "\t";
		for(String label:labels) {
			out.print(sep + label);
		}
		out.println();
		for (int i = 0; i < tmpMatrix.length; i++) {
			out.print(labels[i]);
			for (int j = 0; j < tmpMatrix.length; j++) {
				out.print(sep + tmpMatrix[i][j]);
			}
			out.println();
		}
		out.close();
	}
	
	/**
	 * Calculates a similarity matrix (measured by rmsd) for a set of PDB 
	 * templates and outputs a graph overview for visual inspection.
	 * Now also writes the similarity matrix to a file (to be used by R script).
	 * @param templates
	 * @param conn db connection from where the PDB data will be taken
	 * @param pdbaseDb pdbase database from where the PDB data will be taken
	 * @param graphFile
	 * @param matrixFile
	 * @param similarityGraphRmsdCutoff
	 * @throws IOException if problems writing files or running maxcluster
	 * @throws SQLException if problems getting PDB data from database (only tried if PDB data were not loaded yet in given TemplateList)
	 * @throws PdbLoadException if problems getting PDB data from database (only tried if PDB data were not loaded yet in given TemplateList)
	 */
	public static void writeClusterGraph(TemplateList templates, MySQLConnection conn, String pdbaseDb, File graphFile, File matrixFile, double similarityGraphRmsdCutoff) throws IOException, SQLException, PdbLoadException {
		if (templates.size()<=1) return; // if templates empty or only 1 template there's nothing to compare
		
		if (!templates.isPdbDataLoaded()) {
			templates.loadPdbData(conn, pdbaseDb);
		}
		
		String listFileName = "listfile";
		File listFile = new File(tempDir, listFileName);
		listFile.deleteOnExit();
		PrintWriter out = new PrintWriter(listFile);
		
		// create list file
		Iterator<Template> it = templates.iterator();
		while(it.hasNext()) {
			Template template = it.next();
			// extract pdb and chain code
			
			String pdbCode = template.getId().substring(0, 4);
			String chain = template.getId().substring(4);
			File pdbFile = new File(tempDir, pdbCode + chain + ".pdb");
			pdbFile.deleteOnExit();
			
			if(template.hasPdbData()) {
				
				Pdb pdb = template.getPdb();
				// write to file
				pdb.writeToPDBFile(pdbFile.getAbsolutePath());
				
				// add to listfile
				out.println(pdbFile.getAbsolutePath());
			}
		}
		out.close();
		
		// run maxcluster
		MaxClusterRunner mcr = new MaxClusterRunner(maxClusterExecutable);
		HashMap<Pair<Integer>, Double> matrix = mcr.calculateSequenceIndependentMatrix(listFile.getAbsolutePath(), MaxClusterRunner.ScoreType.RMSD);

		// write similarity matrix file
		String[] templateIds = templates.getIds();
		writeMatrixToFile(matrix, templateIds, matrixFile);
		
		// use shell script to create visualization of matrix file
		String cmdLine = MATRIX_VIS_SCRIPT + " " + matrixFile.getAbsolutePath() + " " + matrixFile.getAbsolutePath() + ".ps";
		System.out.println(cmdLine);
		Runtime.getRuntime().exec(cmdLine);
		
		// generate graph from similarity matrix
		SparseGraph<String, DoubleWrapper> simGraph = new SparseGraph<String, DoubleWrapper>();
		// write nodes
		for(String id:templateIds) {
			simGraph.addVertex(id);
		}
		
		// write edges
		for(Pair<Integer> edge:matrix.keySet()) {
			String start = templateIds[edge.getFirst()-1];
			String end = templateIds[edge.getSecond()-1];
			double weight = matrix.get(edge);
			//System.out.println(weight);
			if(weight < similarityGraphRmsdCutoff) {
				simGraph.addEdge(new BlastUtils().new DoubleWrapper(weight), new Pair<String>(start, end));
			}
		}
		
		// write GDL file for aiSee
		GraphIOGDLFile<String, DoubleWrapper> gdlfileIO = new GraphIOGDLFile<String, DoubleWrapper>();
		gdlfileIO.writeGdlFile(simGraph, graphFile.getAbsolutePath(),
				new Transformer<String, Integer>() {public Integer transform(String s) {return s.hashCode();} },
				new Transformer<String, String>()  {public String transform(String s) {return s;} }, 
				new Transformer<String, String>()  {public String transform(String s) {return "white";} }
		);
	}
	
	
	/**
	 * Testing some of the methods in this class.
	 * @param args
	 */
	public static void main(String[] args) throws IOException, SQLException, PdbLoadException {
//		File blastOutput = new File("");
//		File imgFile = new File("");
//		try {
//			renderBlast(blastOutput, imgFile);
//		} catch(IOException e) {
//			System.out.println("RenderBlast failed: " + e.getMessage());
//		}
		
		// testing writeClusterGraph
		File templateFile = new File("/project/StruPPi/CASP8/results/T0387/T0387.templates");
		File graphFile = new File("/project/StruPPi/CASP8/results/T0387/T0387.templates.gdl");
		File matrixFile = new File("/project/StruPPi/CASP8/results/T0387/T0387.templates.matrix");
		TemplateList templateList = new TemplateList(templateFile);
		
		writeClusterGraph(templateList, new MySQLConnection("talyn","duarte","nieve"),"pdbase",graphFile, matrixFile, 6.0);
		
		
	}
	
}

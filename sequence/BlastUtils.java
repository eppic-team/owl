package sequence;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import org.apache.commons.collections15.Transformer;

import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.util.Pair;

import proteinstructure.*;

/**
 * A collection of little tools related to Blast and Processing Blast output.
 * @author stehr
 *
 */
public class BlastUtils {

	private static final File tempDir = new File("/tmp/");
	private static final String maxClusterExecutable = "/project/StruPPi/bin/maxcluster";
	private static final double similarityGraphGdtCutoff = 50.0;
	
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
	 * Calculates a distance matrix for a set of blast hits and outputs a graph overview for visual inspection.
	 */
	public static void writeClusterGraph(TemplateList templates, File graphFile) throws IOException {
		if (templates.size()==0) return;
		
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

			// write to file
			template.getPdb().dump2pdbfile(pdbFile.getAbsolutePath());

			// add to listfile
			out.println(pdbFile.getAbsolutePath());				
		}
		out.close();
		
		// run maxcluster
		MaxClusterRunner mcr = new MaxClusterRunner(maxClusterExecutable);
		HashMap<Pair<Integer>, Double> matrix = mcr.calculateSequenceIndependentMatrix(listFile.getAbsolutePath(), MaxClusterRunner.ScoreType.GDT);
		
		// generate graph from similarity matrix
		SparseGraph<String, DoubleWrapper> simGraph = new SparseGraph<String, DoubleWrapper>();
		// write nodes
		String[] templateIds = templates.getIds();
		for(String id:templateIds) {
			simGraph.addVertex(id);
		}
		
		// write edges
		for(Pair<Integer> edge:matrix.keySet()) {
			String start = templateIds[edge.getFirst()-1];
			String end = templateIds[edge.getSecond()-1];
			double weight = matrix.get(edge);
			//System.out.println(weight);
			if(weight > similarityGraphGdtCutoff) {
				simGraph.addEdge(new BlastUtils().new DoubleWrapper(weight), new Pair<String>(start, end));
			}
		}
		// create aiSee output
		GraphIOGDLFile<String, DoubleWrapper> gdlfileIO = new GraphIOGDLFile<String, DoubleWrapper>();
		gdlfileIO.writeGdlFile(simGraph, graphFile.getAbsolutePath(),
		new Transformer<String, Integer>() {
				public Integer transform(String s) {return s.hashCode();}
		}, new Transformer<String, String>(){
			public String transform(String s) {return s;}
		}, new Transformer<String, String>(){
			public String transform(String s) {return "white";}
		});
		
	}
	
	/**
	 * Testing some of the method in this class.
	 * @param args
	 */
	public static void main(String[] args) {
		File blastOutput = new File("");
		File imgFile = new File("");
		try {
			renderBlast(blastOutput, imgFile);
		} catch(IOException e) {
			System.out.println("RenderBlast failed: " + e.getMessage());
		}
	}
	
}

package owl.cccp;

import java.io.*;
import java.util.*;

import edu.uci.ics.jung.graph.util.Pair;

import owl.core.structure.*;
import owl.core.structure.graphs.*;

/**
 * Takes a list of pdb files assuming that they are predictions of a single structure
 * and ranks them according to different graph-based measures.
 * See also: /project/StruPPi/henning/projects/graph_averaging/matlab_scripts/benchmark_model_selection.m
 * @author stehr
 */
public class rankModels {
	
	//private static int counter = 1; // for debugging
	
	public static class Scores extends Vector<Double> {
		private static final long serialVersionUID = 1L;

		public Scores(int size) {
			super(size);
		}
		
		public Scores() {
			super();
		}
		
	};
	
	/**
	 * Takes a list of pdb files and create graphs for each of them.
	 * @param listFile
	 */
	private static Vector<RIGraph> loadGraphs(File listFile) {
		Vector<RIGraph> graphs = new Vector<RIGraph>();
		Vector<String> templateFileNames = new Vector<String>();
		
		String contactType = "Cb";
		double distCutoff = 8.0;
		
		// read list of predictions
		try {
			BufferedReader in = new BufferedReader(new FileReader(listFile));
			String line;
			File file;
			Pdb pdb;
			RIGraph graph;
			while ((line =  in.readLine()  ) != null) {
				file = new File(line);
				if(!file.canRead()) {
					System.err.println("File " + line + " not found.");
					System.exit(1);
				} else {
					templateFileNames.add(line);
					try {
						pdb = new PdbfilePdb(file.getAbsolutePath());
						pdb.load(pdb.getChains()[0]);
						graph = pdb.getRIGraph(contactType, distCutoff);
						graphs.add(graph);
					} catch(PdbLoadException e) {
						System.err.println("Error loading pdb structure: " + e.getMessage());
						System.exit(1);
					}
				}
				in.close();
			}

		} catch(FileNotFoundException e) {
			System.err.println("File " + listFile.getAbsolutePath() + " not found.");
			System.exit(1);
		} catch(IOException e) {
			System.err.println("Error reading from file " + listFile.getAbsolutePath());
			System.exit(1);
		}
		
		return graphs;
	}
	
	/**
	 * Returns a vector of the given size containing random numbers between 0.0 and 1.0.
	 * @param size the size of the random vector
	 * @return the random vector
	 */
	private static Scores getRandomVector(int size) {
		Scores randomVector = new Scores(size);
		for (int i = 0; i < size; i++) {
			randomVector.add(Math.random());
		}
		return randomVector;
	}
	
	/**
	 * Returns a vector with the number of contacts for each graph
	 * @return the vector with a contact count for each graph or null on error
	 */
	private static Scores scoreModelsByNumContacts(Vector<RIGraph> graphs) {
		if(graphs == null || graphs.size() == 0) {
			System.err.println("Could not count edges. List of graphs is empty.");
			return null;
		}
		
		Scores scores = new Scores(graphs.size());
		
		for(RIGraph g:graphs) {
			scores.add((double) g.getEdgeCount());
		}
		return scores;
	}
	
	/**
	 * Scores the given models by the fraction of contacts which are consensus contacts.
	 * Thresh < 0 has a special meaning. Instead of counting the number of consensus contacts,
	 * the total number of structures that each contact occurs in is summed up for all edges.
	 * @param graphs the models to be ranked.
	 * @param thresh threshold for calling a consensus contact
	 * @param normalize if true, divide number of consensus contacts by number of total contacts
	 * @returna A scores object or null if the model vector was empty.
	 */
	private static Scores scoreModelsByConsensusContacts(Vector<RIGraph> graphs, double thresh, boolean normalize) {
		
		if(thresh > 1) return null;
		
		if(graphs == null || graphs.size() == 0) return null;
		
		int numGraphs = graphs.size();
		
		Scores scores = new Scores(numGraphs);
		
		int graphSize = graphs.get(1).getFullLength();
		
		HashMap<Pair<Integer>, Integer> contactVotes = new HashMap<Pair<Integer>, Integer>();
		
		// we go through all positions in the alignment
		for (int i=0; i<graphSize; i++){
			for (int j=0; j<graphSize; j++) {
 
				int vote = 0; 
				// scanning all templates to see if they have this contact
				for (RIGraph graph:graphs){					
					if (graph.containsEdgeIJ(i,j)) {
						vote++;
					}
				}
				// putting vote in contactVotes TreeMap
				if (vote>0){
					contactVotes.put(new Pair<Integer>(i,j), vote);
				}				
			}
		}	
		
		//System.out.println("Union edges: " + contactVotes.size());
		
		int voteThreshold = (int) Math.ceil((double)numGraphs*thresh); // i.e. round up of 50%, 40% or 30% (depends on threshold given)
		
		// now iterate over all graphs
		for(RIGraph graph:graphs) {
			int numConsEdges = 0;
			for(RIGEdge e:graph.getEdges()) {
				Pair<RIGNode> p = graph.getEndpoints(e);
				Pair<Integer> pi = new Pair<Integer>(p.getFirst().getResidueSerial(), p.getSecond().getResidueSerial());
				if(contactVotes.containsKey(pi)) {
					if(thresh < 0) {
						// special case, sum up edge weights from contactVotes without applying threshold
						numConsEdges += contactVotes.get(e);
					} else 
 					if(contactVotes.get(e) >= voteThreshold) {
						numConsEdges++;
					}
				}
			}
			double score = (double) numConsEdges;
			// TODO: divide by number of nodes to make scores comparable
			// note, that number of nodes is proportional to number of edges irrespective of the edge type
			// in the case of thresh -1, it may or may not make sense to also divide by the number of templates

			if(normalize) {
				score = score / (double) graph.getEdgeCount();
			}
			
			scores.add(score);
		}
		return scores;
	}
	
	private static void printScores(Scores scores, PrintStream out) {
		//System.out.println(counter++);
		for(Double score:scores) {
			out.printf("%5.2f\t", score);
		}
		out.println();
	}
	
	private static Scores loadGdts(File gdtFile) {
		Scores gdts = new Scores();
		try {
			BufferedReader in = new BufferedReader(new FileReader(gdtFile));
			String line;
			while ((line =  in.readLine()  ) != null) {
				double gdt = Double.parseDouble(line);
				gdts.add(gdt);
			}
		} catch(FileNotFoundException e) {
			System.err.println("File " + gdtFile.getAbsolutePath() + " not found.");
			System.exit(1);
		} catch(IOException e) {
			System.err.println("Error reading from file " + gdtFile.getAbsolutePath());
			System.exit(1);
		}
		return gdts;
	}
		
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		if(args.length < 3) {
			System.out.println("Usage: rankModels <listFile> <gdtFile> <outFile>");
			System.out.println("where");
			System.out.println("listFile contains pdb files names of the models to be scored");
			System.out.println("gdtFile contains gdt scores for the models in the same order");
			System.out.println("the results will be written to outFile");
			System.exit(1);
		}
		
		String listFileName = args[0];
		String gdtFileName = args[1];
		String outFileName = args[2];
		File listFile = new File(listFileName);
		File gdtFile = new File(gdtFileName);
		if(!listFile.canRead()) {
			System.err.println("Error: Could not read from file " + listFileName);
			System.exit(1);
		}		
		if(!gdtFile.canRead()) {
			System.err.println("Error: Could not read from file " + gdtFileName);
			System.exit(1);
		}
		
		Vector<RIGraph> graphs = loadGraphs(listFile);
		Scores gdts = loadGdts(gdtFile);
		Scores scores01r = scoreModelsByConsensusContacts(graphs,0.1, false);
		Scores scores01n = scoreModelsByConsensusContacts(graphs,0.1, true);
		Scores scores02r = scoreModelsByConsensusContacts(graphs,0.2, false);
		Scores scores02n = scoreModelsByConsensusContacts(graphs,0.2, true);		
		Scores scores03r = scoreModelsByConsensusContacts(graphs,0.3, false);
		Scores scores03n = scoreModelsByConsensusContacts(graphs,0.3, true);
		Scores scores04r = scoreModelsByConsensusContacts(graphs,0.4, false);
		Scores scores04n = scoreModelsByConsensusContacts(graphs,0.4, true);
		Scores scores05r = scoreModelsByConsensusContacts(graphs,0.5, false);
		Scores scores05n = scoreModelsByConsensusContacts(graphs,0.5, true);
		//Scores scores06r = scoreModelsByConsensusContacts(graphs,0.6, false);
		//Scores scores06n = scoreModelsByConsensusContacts(graphs,0.6, true);
		Scores scoresSumr = scoreModelsByConsensusContacts(graphs,-1, false);
		Scores scoresSumn = scoreModelsByConsensusContacts(graphs,-1, true);
		Scores numContacts = scoreModelsByNumContacts(graphs);
		Scores randomVector = getRandomVector(graphs.size());
		
		File outFile = new File(outFileName);
		if(!outFile.canWrite()) {
			System.err.println("Error: Could not write to file " + outFileName);
		}
		PrintStream fileOut = System.out;
		try {
			fileOut = new PrintStream(new FileOutputStream(outFile));
		} catch(FileNotFoundException e) {
			System.err.println("File " + outFile.getAbsolutePath() + " not found.");
		}
		
		System.out.println("Cols=" + graphs.size());
		System.out.println("Rows=15 (gdt,0.1r,0.2r,0.3r,0.4r,0.5r,sum_r,num_cont,0.1n,0.2n,0.3n,0.4n,0.5n,sum_n,rand)");
		
		printScores(gdts, fileOut);
		printScores(scores01r, fileOut);
		printScores(scores02r, fileOut);		
		printScores(scores03r, fileOut);
		printScores(scores04r, fileOut);
		printScores(scores05r, fileOut);
		//printScores(scores06r, fileOut);
		printScores(scoresSumr, fileOut);
		printScores(numContacts, fileOut);
		
		printScores(scores01n, fileOut);
		printScores(scores02n, fileOut);
		printScores(scores03n, fileOut);
		printScores(scores04n, fileOut);
		printScores(scores05n, fileOut);
		//printScores(scores06n, fileOut);
		printScores(scoresSumn, fileOut);
		printScores(randomVector, fileOut);
		
		System.out.println("Output written to " + outFileName);
		
		// See also: /project/StruPPi/henning/projects/graph_averaging/matlab_scripts/benchmark_model_selection.m
	}

}

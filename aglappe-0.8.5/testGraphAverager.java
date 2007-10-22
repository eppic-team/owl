import java.io.*;
import java.util.*;
import proteinstructure.*;

public class testGraphAverager {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// Plan:
		// use getopt to parse the following parameters
		// - original structure (optional, for benchmarking)
		// - file with list of structures for averaging (required)
		// - contact type & distance cutoff (both required)
		// - consensus edge threshold (defaults to 0.5?)
		// - run tinker or not
		// - parameters for tinker model picking (optional)
		// - number of tinker models to output (only if run tinker is set, defaults to one)
		// - output benchmark result (graph comparison or structure comparison, only if original structure is specified)
		// - output file for resulting graph (optional)
		
		// take pdb files from command line (later: use getopt)
		// create trivial alignment (check alignment)
		// create graphAverager
		// create consensus graph
		// compare with graph for native structure (first argument)
		
		// do tinker reconstruction
		// compare with native structure
		
		// benchmarking:
		// - use score of best model vs. native (either graph based or structure based)
		// - use matlab to optimize input parameters
		// - create matlab executable for optimizing command line parameters:
		//   * specify command line to be run and constant parameters
		//   * specify parameter ranges to be optimized (discrete numerical, continuous numerical with steps, discrete set)
		//   * specify regexp to parse result from output (stdout or file)
		//   * specify optimization procedure (use default matlab ones, e.g. each parameter separate, genetic algorithm, simul. annealing, ...)
		//   * specify convergence threshold or timeout
		//   * run on cluster?
		
		
		// read command line parameters
		if(args.length < 4) {
			System.out.println("Usage: testGraphAverager <targetPdbFile> <ChainCode> <ListOfPredictionFiles> <edgeThreshold>");
			System.exit(1);
		}
		
		File targetFile = new File(args[0]);
		File listFile = new File(args[2]);
		String chainCode = args[1];
		String thresholdStr = args[3];
		Pdb target = null;
		Graph targetGraph = null;
		Vector<String> modelFileNames = new Vector<String>();
		Vector<Graph> models = new Vector<Graph>();
		String contactType = "Cb";
		double distCutoff = 7.0;
		double graphAveragingThreshold = Double.parseDouble(thresholdStr); 
		
		if(!targetFile.canRead()) {
			System.err.println("Can not read from file " + targetFile.getAbsolutePath());
			System.exit(1);
		}
		if(!listFile.canRead()) {
			System.err.println("Can not read from file " + listFile.getAbsolutePath());
			System.exit(1);
		}
		
		// read input files
		try {
			target = new PdbfilePdb(targetFile.getAbsolutePath(), chainCode);
			targetGraph = target.get_graph(contactType, distCutoff);
		} catch(PdbChainCodeNotFoundError e) {
			System.err.println("Chain code " + chainCode + " not found in file " + targetFile.getAbsolutePath());
			System.exit(1);
		} catch(PdbfileFormatError e) {
			System.err.println("Formating error in file " + targetFile.getAbsolutePath());	
			System.exit(1);
		} catch (FileNotFoundException e) {
			System.err.println("File " + targetFile.getAbsolutePath() + " not found");	
			System.exit(1);			
		} catch(IOException e) {
			System.err.println("Error reading from file " + targetFile.getAbsolutePath());	
			System.exit(1);				
		}
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(listFile));
			String line;
			File file;
			Pdb pdb;
			Graph graph;
			while ((line =  in.readLine()  ) != null) {
				file = new File(line);
				if(!file.canRead()) {
					System.err.println("File " + line + " not found.");
					System.exit(1);
				} else {
					modelFileNames.add(line);
					try {
						pdb = new PdbfilePdb(file.getAbsolutePath(),Pdb.NULL_CHAIN_CODE);
						graph = pdb.get_graph(contactType, distCutoff);
						models.add(graph);
					} catch(PdbChainCodeNotFoundError e) {
						System.err.println("Chain code " + chainCode + " not found in file " + file.getAbsolutePath());
						System.exit(1);
					} catch(PdbfileFormatError e) {
						System.err.println("Formating error in file " + file.getAbsolutePath());	
						System.exit(1);
					} catch (FileNotFoundException e) {
						System.err.println("File " + file.getAbsolutePath() + " not found");	
						System.exit(1);			
					} catch(IOException e) {
						System.err.println("Error reading from file " + file.getAbsolutePath());	
						System.exit(1);				
					}
				}
			}
			in.close();

		} catch(FileNotFoundException e) {
			System.err.println("File " + listFile.getAbsolutePath() + " not found.");
			System.exit(1);
		} catch(IOException e) {
			System.err.println("Error reading from file " + listFile.getAbsolutePath());
			System.exit(1);
		}
		
		// create alignment
		TreeMap<String,String> sequences = new TreeMap<String, String>();
		TreeMap<String, Graph> templateGraphs = new TreeMap<String, Graph>();
		sequences.put(targetFile.getAbsolutePath(), targetGraph.getSequence());
		for(int i=0; i < models.size(); i++) {
			sequences.put(modelFileNames.get(i), models.get(i).getSequence());
			templateGraphs.put(modelFileNames.get(i), models.get(i));
		}
		Alignment al = new Alignment(sequences);
		
		// create GraphAverager
		GraphAverager grav = new GraphAverager(targetGraph.getSequence(), al, templateGraphs, targetFile.getAbsolutePath());
		Graph resultGraph = grav.doAveraging(graphAveragingThreshold);
		
		// compare prediction with target
		PredEval eval = resultGraph.evaluatePrediction(targetGraph);
		eval.print();
		
		// write result graph to file
		try {
			resultGraph.write_graph_to_file("predictedGraph.cm");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
}

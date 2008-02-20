package proteinstructure;

import graphAveraging.GraphAverager;

import java.io.*;
import java.util.*;

/**
 * An ensemble of residue interactions graphs (RIGs) for the same protein structure.
 * Examples: NMR ensembles, structure predictions, folding trajectories.
 * Currently it is assumed that all of these structures share the same underlying
 * sequence. This case would be something like a multi-model RIG, corresponding
 * to a multi-model PDB file. Alternatively, the restriction that all graphs are
 * based on the same structure could be lifted by adding an alignment.
 * @author stehr
 * @date 2007-12-19
 */
public class RIGEnsemble extends ArrayList<RIGraph> {
	
	/*------------------------------ constants ------------------------------*/
	private static final long serialVersionUID = 1L;
	public static final String DEFAULT_EDGE_TYPE = "Cb";
	public static final double DEFAULT_DIST_CUTOFF = 8.0;
	
	/*--------------------------- member variables --------------------------*/
	String edgeType;
	double distCutoff;
	
	/*----------------------------- constructors ----------------------------*/
	/**
	 * Generate an empty RIGEnsemble with the default edgeType and distCutoff.
	 */
	public RIGEnsemble() {
		super();
		edgeType = DEFAULT_EDGE_TYPE;
		distCutoff = DEFAULT_DIST_CUTOFF;
	}
	
	/**
	 * Creates an empty RIGEnsemble with the given edgeType and distCutoff.
	 * @param edgeType
	 * @param distCutoff
	 */
	public RIGEnsemble(String edgeType, double distCutoff) {
		super();
		this.edgeType = edgeType;
		this.distCutoff = distCutoff;
	}
	
	/**
	 * Generate a RIGEnsemble from a listfile, i.e. a text file containing names of data files.
	 * The list file may point to files of different types. If the file type can be recognized
	 * the appropriate loading method will be called. Graphs are being generated from PDB files
	 * using the global edgeType and distanceCutoff. For contact map files, the edgeType/cutoff
	 * has to match the global one, otherwise an exception is thrown. Other erros are supressed
	 * so that if single files contain errors, others will be still loaded. If any file contains
	 * multiple chains, the first one will be read.
	 * @param listFile
	 * @return number of files read
	 */
	public int loadFromListFile(File listFile) throws FileNotFoundException, IOException {
		// for each file in list, load or generate the graph (depending on file type)
			BufferedReader in = new BufferedReader(new FileReader(listFile));
			String line;
			File file;
			Pdb pdb;
			RIGraph graph;
			int fr = 0;
			while ((line =  in.readLine()) != null) {
				file = new File(line);
				if(!file.canRead()) {
					System.err.println("Warning: File " + line + " not found. Skipping.");
				} else {
					int fileType = FileTypeGuesser.guessFileType(file);
					switch(fileType) {
					case(FileTypeGuesser.PDB_FILE):
					case(FileTypeGuesser.RAW_PDB_FILE):
					case(FileTypeGuesser.CASP_TS_FILE):
						try {
							pdb = new PdbfilePdb(file.getAbsolutePath());
							String[] chains = pdb.getChains();
							pdb.load(chains[0]);	// load first chain
							graph = pdb.get_graph(this.edgeType, this.distCutoff);
							this.addRIG(graph);
							fr++;
						} catch(PdbLoadError e) {
							System.err.println("Error loading pdb structure: " + e.getMessage());
							//System.exit(1);
						}
						break;
					case(FileTypeGuesser.CIF_FILE):
						try {
							pdb = new CiffilePdb(file.getAbsolutePath());
							String[] chains = pdb.getChains();
							pdb.load(chains[0]);	// load first chain
							graph = pdb.get_graph(this.edgeType, this.distCutoff);
							this.addRIG(graph);
							fr++;
						} catch(PdbLoadError e) {
							System.err.println("Error loading pdb structure: " + e.getMessage());
							//System.exit(1);
						}
						break;
					case(FileTypeGuesser.AGLAPPE_CM_FILE):
						try {
							graph = new FileRIGraph(file.getAbsolutePath());
							this.addRIG(graph);
							fr++;
						} catch (GraphFileFormatError e) {
							System.err.println("Error loading from contact map file: " + e.getMessage());
							//System.exit(1);
						}
						break;
					case(FileTypeGuesser.CASP_RR_FILE):
						try {
								graph = new CaspRRFileRIGraph(file.getAbsolutePath());
								this.addRIG(graph);
								fr++;
							} catch (GraphFileFormatError e) {
								System.err.println("Error loading from RR file: " + e.getMessage());
								//System.exit(1);
							}
					}
					break;
				}
			}
			in.close();
			return fr;
	}
	
	/**
	 * Generate a RIGEnsemble from a multi-model PDB or mmCIF file. If the file contains
	 * multiple chains, only the first one is read.
	 * @param file the input file
	 * @return the number of models read
	 */
	public int loadFromMultiModelFile(File file) throws IOException, PdbLoadError {
		// for each model in file, generate a graph
		Pdb pdb;
		RIGraph graph;
		Integer[] models;
		String[] chains;
		String chain;
		int mr = 0;
		int fileType = FileTypeGuesser.guessFileType(file);
		switch(fileType) {
		case(FileTypeGuesser.PDB_FILE):
		case(FileTypeGuesser.RAW_PDB_FILE):
			pdb = new PdbfilePdb(file.getAbsolutePath());
			models = pdb.getModels();
			chains = pdb.getChains();
			chain = chains[0];
			for(int mod: models) {
				//pdb = new PdbfilePdb(file.getAbsolutePath());
				pdb.load(chain, mod);
				graph = pdb.get_graph(this.edgeType, this.distCutoff);
				this.addRIG(graph);
				mr++;
			}
			break;
		case(FileTypeGuesser.CIF_FILE):
			pdb = new CiffilePdb(file);
			models = pdb.getModels();
			chains = pdb.getChains();
			chain = chains[0];
			for(int mod: models) {
				//pdb = new CiffilePdb(file.getAbsolutePath());
				pdb.load(chain, mod);
				graph = pdb.get_graph(this.edgeType, this.distCutoff);
				this.addRIG(graph);
				mr++;
			}
			break;			
		}
		return mr;
	}

	/*---------------------------- public methods ---------------------------*/
	
	public int getEnsembleSize() {
		return this.size();
	}
	
	public String getEdgeType() {
		return edgeType;
	}
	
	public double getDistCutoff() {
		return distCutoff;
	}
	
	public void addRIG(RIGraph g) {
		if(!g.getContactType().equals(this.edgeType)) System.err.println("Warning: Contact types do not match.");
		if(g.getCutoff() != this.distCutoff) System.err.println("Warning: Distance cutoffs do not match.");
		this.add(g);
	}
	
	public RIGraph getRIG(int i) {
		return this.get(i);
	}
	
	public RIGraph[] getRIGs() {
		RIGraph[] graphs = new RIGraph[this.size()];
		return this.toArray(graphs);
	}
		
	/*--------------------------------- main --------------------------------*/
	
	// for testing
	public static void main(String[] args) throws IOException, PdbLoadError {
		
		boolean loadFromList = false;
		int filesRead = 0;
		String outFileName = "ensemble.cm";
		
		if(args.length < 2) {
			System.out.println("Usage: RIGEnsemble -l/-m fileName");
			System.out.println("-l	load from list file");
			System.out.println("-m  load from multi-model file (cif or pdb)");
			System.exit(1);
		}
		String opt = args[0];
		if(opt.equals("-l")) {
			loadFromList = true;
		} else
		if(opt.equals("-m")) {
			loadFromList = false;
		} else {
			System.err.println("Unknown option " + opt + ". Expected -l or -m");
		}
		
		String fileName = args[1];
		File inFile = new File(fileName);
		
		RIGEnsemble rigs = new RIGEnsemble();
		if(loadFromList) {
			filesRead = rigs.loadFromListFile(inFile);
		} else {
			filesRead = rigs.loadFromMultiModelFile(inFile);
		}
		
		System.out.println("RIGs read from file:\t" + filesRead);
		System.out.println("RIGs in ensemble:\t" + rigs.getEnsembleSize());
		
		GraphAverager ga = new GraphAverager(rigs);
		RIGraph weightedGraph = ga.getAverageGraph();
		System.out.println("Writing average graph to file "+outFileName);
		weightedGraph.write_graph_to_file(outFileName);
		
	}
}

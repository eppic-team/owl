package proteinstructure;

import graphAveraging.GraphAverager;

import java.io.*;
import java.util.*;

import sequence.Sequence;

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
public class RIGEnsemble {
	
	/*------------------------------ constants ------------------------------*/
	private static final long serialVersionUID = 1L;
	public static final String DEFAULT_EDGE_TYPE = "Cb";
	public static final double DEFAULT_DIST_CUTOFF = 8.0;
	
	/*--------------------------- member variables --------------------------*/
	private ArrayList<RIGraph> ensemble;
	private ArrayList<String> fileNames; // file names of loaded structures
	// TODO: Should this information be stored in the RIG?
	private String edgeType;
	private double distCutoff;
	private boolean loadOnlyFirstModels;	// if true, all models != 1 will be ignored
	
	/*----------------------------- constructors ----------------------------*/
	/**
	 * Generate an empty RIGEnsemble with the default edgeType and distCutoff.
	 */
	public RIGEnsemble() {
		ensemble = new ArrayList<RIGraph>();
		fileNames = new ArrayList<String>();
		edgeType = DEFAULT_EDGE_TYPE;
		distCutoff = DEFAULT_DIST_CUTOFF;
		loadOnlyFirstModels = false;
	}
	
	/**
	 * Creates an empty RIGEnsemble with the given edgeType and distCutoff.
	 * @param edgeType
	 * @param distCutoff
	 */
	public RIGEnsemble(String edgeType, double distCutoff) {
		ensemble = new ArrayList<RIGraph>();
		fileNames = new ArrayList<String>();
		this.edgeType = edgeType;
		this.distCutoff = distCutoff;
		loadOnlyFirstModels = false;
	}
	
	/**
	 * Generate a RIGEnsemble from files in a directory.
	 * The directory may contain files of different types. If the file type can be recognized
	 * the appropriate loading method will be called. Graphs are being generated from PDB files
	 * using the global edgeType and distanceCutoff. For contact map files, the edgeType/cutoff
	 * has to match the global one, otherwise an exception is thrown. Other erros are supressed
	 * so that if single files contain errors, others will be still loaded. If any file contains
	 * multiple chains, the first one will be read.
	 * @param directory
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @return number of files read
	 */
	public int loadFromDirectory(File dir) throws FileNotFoundException, IOException {
		if(!dir.isDirectory()) throw new IOException(dir.getName() + " is not a directory.");
		
		return 0;
	}
	
	/**
	 * Generates a RIGEnsemble from a listfile, i.e. a text file containing names of data files
	 * or from a directory containing data files. For each file in the list or directory, the
	 * method tries to determine the file type and if recognized, calls the appropriate loading
	 * method. Graphs are being generated from PDB files using the global edgeType and distanceCutoff.
	 * For contact map files, the edgeType/cutoff has to match the global one, otherwise an exception
	 * is thrown. Other errors are supressed so that if single files contain errors, others will
	 * still be loaded. If any file contains multiple chains, the first one is taken.
	 * The paramtere commonSequence can be used for the case where the loaded structures are Casp
	 * predictions which contain different numbers of predicted residues but the constructed graphs
	 * should all have the same sequence. In this case the sequence can be passed as a parameter and
	 * will be set for all loaded structures. If commonSequence is null or the file is not a PDB file,
	 * this will be ignored.
	 * @param listFile
	 * @param commonSequence if not null, this sequence will be enforced on all loaded graphs
	 * @return number of files read
	 */
	public int loadFromFileList(File list, Sequence commonSequence) throws FileNotFoundException, IOException {

			// load filenames from listFile or directory
			String[] files;
			if(list.isDirectory()) {
				files = list.list();
			    if (files == null) {
			        // Either dir does not exist or is not a directory
			    	throw new FileNotFoundException("Could not open directory " + list);
			    }
			    for (int i = 0; i < files.length; i++) {
					files[i] = new File(list, files[i]).getAbsolutePath();
				}
			} else {
				BufferedReader in = new BufferedReader(new FileReader(list));
				String line;
				ArrayList<String> tempList = new ArrayList<String>();
				while ((line =  in.readLine()) != null) {
					tempList.add(line);
				}
				in.close();
				files = (String[]) tempList.toArray();
			}
		
			// for each file in list, load or generate the graph (depending on file type)
			File file;
			Pdb pdb;
			RIGraph graph;
			int fr = 0;
			for(String filename:files) {
				file = new File(filename);
				if(!file.canRead()) {
					System.err.println("Warning: File " + filename + " not found. Skipping.");
				} else {
					int fileType = FileTypeGuesser.guessFileType(file);
					switch(fileType) {				
					case(FileTypeGuesser.PDB_FILE):
					case(FileTypeGuesser.RAW_PDB_FILE):
					case(FileTypeGuesser.CASP_TS_FILE):
						try {
							pdb = new PdbfilePdb(file.getAbsolutePath());
							String[] chains = pdb.getChains();
							Integer[] models = pdb.getModels();
							if(loadOnlyFirstModels && models[0] != 1) continue;
							//System.out.println(filename + ":" + chains[0]);
							pdb.load(chains[0], models[0]);	// load first chain and first model
							if(commonSequence != null) ((PdbfilePdb) pdb).setSequence(commonSequence.getSeq());
							graph = pdb.get_graph(this.edgeType, this.distCutoff);
							this.addRIG(graph);
							this.addFileName(filename);
							fr++;
						} catch(PdbLoadError e) {
							System.err.println("Error loading pdb structure " + file.getPath() + ":" + e.getMessage());
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
							this.addFileName(filename);
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
							this.addFileName(filename);
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
								this.addFileName(filename);
								fr++;
							} catch (GraphFileFormatError e) {
								System.err.println("Error loading from RR file: " + e.getMessage());
								//System.exit(1);
							}
					break;
					default: System.err.println("Could not determine filetype of " + filename + ". Skipping.");
					}
				}
			}
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
		default: System.err.println("Error: Could not determine filetype of " + file.getName());
		}
		return mr;
	}

	/**
	 * Creates a new rig ensemble from a map of filenames to rigs.
	 * @param rigs
	 * @return the number of graph loaded
	 */
	public int loadFromGraphMap(Map<String, RIGraph> rigs) {
		int gr = 0;
		for(String name:rigs.keySet()) {
			RIGraph rig = rigs.get(name);
			if(loadOnlyFirstModels && rig.model != 1) continue;	// skip non-first model graphs
			if(!rig.contactType.equals(this.edgeType)) {
				System.err.printf("Warning. Contact type of graph %s (%s) and ensemble (%s) do not match.\n",
						name, rig.contactType, this.edgeType);
			}
			if(rig.distCutoff != this.distCutoff) {
				System.err.printf("Warning. Distance cutoff of graph %s (%3.1f) and ensemble (%3.1f) do not match.\n",
						name, rig.distCutoff, this.distCutoff);				
			}
			this.addRIG(rig);
			this.addFileName(name);
			gr++;
		}
		return gr;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	public int getEnsembleSize() {
		return ensemble.size();
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
		this.ensemble.add(g);
	}
	
	private void addFileName(String filename) {
		this.fileNames.add(filename);
	}
	
	public RIGraph getRIG(int i) {
		return ensemble.get(i);
	}
	
	public RIGraph[] getRIGs() {
		RIGraph[] graphs = new RIGraph[this.getEnsembleSize()];
		return this.ensemble.toArray(graphs);
	}

	/*
	 * If the graphs in this ensemble were read from a list file or directory, this method
	 * returns the filename of the i'th RIG in this ensemble.If the ensemble was created otherwise,
	 * or RIGs have been added manually (using addRIG) the returned value is undefined and may be null.
	 */
	public String getFileName(int i) {
		return this.fileNames.get(i);
	}
	
	/*
	 * If the graphs in this ensemble were read from a list file or directory, this method
	 * returns these filenames as an array where the number corresponds to the numbering of
	 * the RIGs. If the ensemble was created otherwise, or RIGs have been added manually
	 * (using addRIG) the returned value is undefined.
	 */
	public String[] getFilenames() {
		String[] filenames = new String[this.fileNames.size()];
		return this.fileNames.toArray(filenames);
	}
	
	/**
	 * After this function is being called, only first models will be loaded.
	 */
	public void loadOnlyFirstModels() {
		this.loadOnlyFirstModels = true;
	}
	
	/*--------------------------------- main --------------------------------*/
	
	// for testing
	public static void main(String[] args) throws Exception {
		
		boolean loadFromList = false;
		int filesRead = 0;
		String outFileName = "ensemble.cm";
		String outFileName2 = "consensus.cm";
		
		if(args.length < 2) {
			System.out.println("Usage: RIGEnsemble <-d|-l|-m> <fileName> [sequenceFile]");
			System.out.println("-d  load from directory");
			System.out.println("-l  load from list file");
			System.out.println("-m  load from multi-model file (cif or pdb)");
			System.out.println("The optional parameter sequenceFile can be used to set the same sequence for all graphs.");
			System.exit(1);
		}
		String opt = args[0];
		if(opt.equals("-l")) {
			loadFromList = true;
		} else
		if(opt.equals("-d")) {
			loadFromList = true;
		} else 
		if(opt.equals("-m")) {
			loadFromList = false;
		} else {
			System.err.println("Unknown option " + opt + ". Expected -l or -m");
		}
		
		Sequence commonSequence = null;
		if(args.length > 2) {
			String seqFileName = args[2];
			commonSequence = new Sequence();
			try {
				commonSequence.readFromFastaFile(new File(seqFileName));
			} catch (FastaFileFormatError e) {
				System.err.println("Failed to read from Fasta file " + seqFileName + ":" + e.getMessage());
				System.exit(1);
			}
		}
		
		String fileName = args[1];
		File inFile = new File(fileName);
		
		RIGEnsemble rigs = new RIGEnsemble();
		if(loadFromList) {
			filesRead = rigs.loadFromFileList(inFile, commonSequence);
		} else {
			filesRead = rigs.loadFromMultiModelFile(inFile);
		}
		
		System.out.println("RIGs read from file:\t" + filesRead);
		System.out.println("RIGs in ensemble:\t" + rigs.getEnsembleSize());
		
		GraphAverager ga = new GraphAverager(rigs);
		System.out.println("Overall consensus:\t" + ga.getEnsembleConsensusScore());
		RIGraph weightedGraph = ga.getAverageGraph();
		RIGraph consGraph = ga.getConsensusGraph(0.1);
		System.out.println("Writing average graph to file "+outFileName);
		weightedGraph.write_graph_to_file(outFileName);
		System.out.println("Writing consensus graph to file "+outFileName2);
		consGraph.write_graph_to_file(outFileName2);
		
	}
}


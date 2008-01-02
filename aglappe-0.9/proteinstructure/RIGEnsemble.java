package proteinstructure;

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
	 * has to match the global one, otherwise an exception is thrown.
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
	 * Generate a RIGEnsemble from a multi-model PDB or mmCIF file. 
	 */
	public void loadFromMultiModelFile(File file) {
		// for each model in file, generate a graph
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
	
}

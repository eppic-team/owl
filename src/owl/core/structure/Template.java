package owl.core.structure;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import owl.core.connections.ScopConnection;
import owl.core.runners.blast.BlastHit;
import owl.core.structure.features.ScopRegion;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatException;


/**
 * A template protein structure to be used in homology modelling
 *
 */
public class Template {

	private static final String SCOP_DIR = "/project/StruPPi/Databases/SCOP";
	private static final String SCOP_VERSION = "1.73";
	
	// 2 types of Template: 
	// - DB gets PDB data from pdbase db (provided the PDB id) 
	// - FILE gets PDB data from file (provided the file name) 
	private static enum Type {DB, FILE};
	
	// 2 cases for the id: 
	// type==Type.DB -> id is pdbCode+pdbChainCode
	// type==Type.FILE -> id is the file name (without path) 
	private String id;
	
	private File pdbFile; // null if type==Type.DB (i.e. only id provided)
	
	private String scopSccsString;
	private String titleString;
	
	private PdbChain pdb;
	private RIGraph graph;

	private Type type;

	private BlastHit blastHit;
	
	/**
	 * Constructs a Template given an id in format pdbCode+chain, e.g. 1abcA
	 * @param id
	 */
	public Template(String id) {
		this.id = id;
		this.type = Type.DB; 
	}
	
	/**
	 * Constructs a Template given a PDB/CASP TS file
	 * @param pdbFile
	 */
	public Template(File pdbFile) {
		this.pdbFile = pdbFile;
		this.id = pdbFile.getName(); // i.e. the file name without the path
		this.type = Type.FILE;
	}
	
	/**
	 * Constructs a new Template given a BlastHit 
	 * @param hit
	 */
	public Template(BlastHit hit) {
		this.id = hit.getTemplateId();
		this.blastHit = hit;
		this.type = Type.DB;
	}

	/**
	 * Gets the pdb structure coordinates into the pdb object 
	 * reading the data either from CIF dir repository or from FILE depending on the type of this Template
	 * If reading from file only the first chain encountered in the file will be read.
	 * If type is not CIF dir then simply pass nulls for cifRepoDir
	 * @param cifRepoDir a directory containing mmCIF files (.gz) 
	 */
	protected void loadPdbData(String cifRepoDir) throws PdbLoadException, IOException, FileFormatException {
		
		if (type==Type.DB) {
			String pdbCode = id.substring(0, 4);
			String chain = id.substring(4);

			File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".cif");
			cifFile.deleteOnExit();
			PdbAsymUnit.grabCifFile(cifRepoDir, null, pdbCode, cifFile, false);				
			PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile);
			pdb = fullpdb.getChain(chain);
		} 
		else if (type==Type.FILE) {
			PdbAsymUnit fullpdb = null;
			try {
				fullpdb = new PdbAsymUnit(pdbFile);
			} catch (FileFormatException e) {
				System.err.println("Unexpected error while reading pdb file "+pdbFile+": "+e.getMessage());
				System.exit(1);
			}
			pdb = fullpdb.getFirstChain();
		}

	}
	
	/**
	 * Gets the graph from the pdb object for given contact type and cutoff
	 * If there's no PDB data for this template then no graph will be calculated
	 * and graph will remain null, thus calling hasGraphData() would return false.
	 * If there was a graph already loaded for this Template then it will be 
	 * simply overwritten.
	 * @param ct
	 * @param cutoff
	 */
	protected void loadRIGraph(String ct, double cutoff) {
		if (this.hasPdbData()) {
			graph = pdb.getRIGraph(ct, cutoff);
		}
	}
	
	/**
	 * Get all SCOP sccs ids in a comma separated string
	 * pdb must be loaded already. If pdb==null a NullPointerException will be thrown
	 */
	private void getScopInfo() {
		try {
			ScopConnection.parseScop(pdb,SCOP_VERSION, SCOP_DIR);
			Iterator<ScopRegion> it = pdb.getScop().getIterator();
			scopSccsString = "";
			while (it.hasNext()) {
				scopSccsString += it.next().getSccs() +", ";
			}
			if (scopSccsString.contains(",")) // choping off the last comma if the string is not empty
				scopSccsString = scopSccsString.substring(0, scopSccsString.length()-2); 
		} catch (IOException e) {
			System.err.println("Couldn't get SCOP annotation for structure "+id+", error: "+e.getMessage());
		}
	}
	
	/**
	 * Prints info about this template: id and scop sccs string
	 */
	protected void print() {
		System.out.println(id+"\t"+scopSccsString);
	}
	
	/**
	 * Returns the id of this template in the format pdbCode+chain, e.g. 1abcA
	 * @return
	 */
	public String getId() {
		return this.id;
	}
	
	/**  
	 * Returns the PdbChain object with the actual structure
	 * Use {@link #hasPdbData()} to check if this method can be called.
	 * @return the PdbChain object or null if PDB data are not loaded/present
	 */
	public PdbChain getPdb() {
		return this.pdb;
	}
	
	/**
	 * Returns the RIGraph object 
	 * Use {@link #hasGraphData()} to check if this method can be called.
	 * @return
	 */
	public RIGraph getRIGraph(){
		return this.graph;
	}
	
	/**
	 * Sets the PdbChain object of this template.
	 * @param pdb
	 */
	public void setPdb(PdbChain pdb) {
		this.pdb=pdb;
	}
	
	/**
	 * Tells whether there is structure data available for this template
	 * @return
	 */
	public boolean hasPdbData() {
		return (pdb!=null);
	}
	
	/**
	 * Tells wheter this template contains a graph
	 * @return
	 */
	public boolean hasGraphData() {
		return (graph!=null);
	}
	
	/**
	 * Returns the scop sccs string for this Template
	 * The first time this method is called the SCOP data is parsed from file.
	 * Subsequente calls take the SCOP data from the cached variable.
	 * @return
	 * @throws IllegalArgumentException if PDB data is not loaded yet
	 */
	public String getScopSccsString() {
		if (scopSccsString==null) {
			if (!this.hasPdbData()) {
				throw new IllegalArgumentException("PDB data not loaded for this Template. Can't get SCOP data.");
			}
			getScopInfo();				
		}
		return this.scopSccsString;
	}
	
	/**
	 * Returns the pdb title string for this Template
	 * The first time this method is called the title is taken from the pdb object, 
	 * subsequent calls take the title from the cached variable.  
	 * @return the title string
	 * @throws IllegalArgumentException if PDB data is not loaded yet
	 */
	public String getTitle() {
		if (titleString==null) {
			if (!this.hasPdbData()) {
				throw new IllegalArgumentException("PDB data not loaded for this Template. Can't get PDB title.");
			}
			this.titleString = pdb.getParent().getTitle();
		}
		return this.titleString;
	}
	
	/**
	 * Returns the BlastHit corresponding to this Template or null if this 
	 * Template didn't originate from a BlastHit
	 * @return
	 */
	public BlastHit getBlastHit() {
		return this.blastHit;
	}
	
}

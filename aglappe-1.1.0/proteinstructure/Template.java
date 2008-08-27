package proteinstructure;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Iterator;

import sequence.BlastHit;
import sequence.GTGHit;
import tools.MySQLConnection;

/**
 * A template protein structure to be used in homology modelling
 *
 */
public class Template {

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
	
	private Pdb pdb;
	private RIGraph graph;

	private Type type;

	private BlastHit blastHit;
	private GTGHit gtgHit;
	
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
	 * Constructs a new Template given a GTGHit 
	 * @param hit
	 */
	public Template(GTGHit hit) {
		this.id = hit.getTemplateId();
		this.gtgHit = hit;
		this.type = Type.DB;
	}
	
	/**
	 * Gets the pdb structure coordinates into the pdb object 
	 * reading the data either from DB or from FILE depending on the type of this Template
	 * If reading from file only the first chain encountered in the file will be read.
	 * If type is not DB then simply pass nulls for conn and pdbaseDb
	 * @param conn
	 * @param pdbaseDb
	 */
	protected void loadPdbData(MySQLConnection conn, String pdbaseDb) throws SQLException, PdbCodeNotFoundError, PdbLoadError {
		
		if (type==Type.DB) {
			String pdbCode = id.substring(0, 4);
			String chain = id.substring(4);

			pdb = new PdbasePdb(pdbCode, pdbaseDb, conn);
			pdb.load(chain);
		} 
		else if (type==Type.FILE) {
			pdb = new PdbfilePdb(pdbFile.getAbsolutePath());
			String[] chains = pdb.getChains();
			pdb.load(chains[0]);
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
			graph = pdb.get_graph(ct, cutoff);
		}
	}
	
	/**
	 * Get all SCOP sccs ids in a comma separated string
	 * pdb must be loaded already. If pdb==null a NullPointerException will be thrown
	 */
	private void getScopInfo() {
		try {
			pdb.checkScop(SCOP_VERSION, false);
			Iterator<ScopRegion> it = pdb.getScop().getIterator();
			scopSccsString = "";
			while (it.hasNext()) {
				scopSccsString += it.next().sccs +", ";
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
	 * Returns the Pdb object with the actual structure
	 * Use {@link #hasPdbData()} to check if this method can be called.
	 * @return the Pdb object or null if PDB data are not loaded/present
	 */
	public Pdb getPdb() {
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
	 * Sets the Pdb object of this template.
	 * @param pdb
	 */
	public void setPdb(Pdb pdb) {
		this.pdb=pdb;
	}
	
	/**
	 * Tells whether there is structure data available for this template
	 * @return
	 */
	public boolean hasPdbData() {
		return (pdb!=null && pdb.isDataLoaded());
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
	 * If PDB data is not loaded yet it will be loaded from database.
	 * @param conn
	 * @param pdbaseDb
	 * @return
	 */
	public String getScopSccsString(MySQLConnection conn, String pdbaseDb) {
		if (scopSccsString==null) {
			try {
				if (pdb==null) {
					loadPdbData(conn, pdbaseDb);
				}
				getScopInfo();
				
			}  catch (SQLException e) {
				System.err.println("Couldn't get the template structure "+id+" because of SQL error: "+e.getMessage());
				pdb = null;
				scopSccsString = "";
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Couldn't get the template structure "+id+" because pdb code was not found");
				pdb = null;
				scopSccsString = "";
			} catch (PdbLoadError e) {
				System.err.println("Couldn't get the template structure "+id+" because of pdb load error: "+e.getMessage());
				pdb = null;
				scopSccsString = "";
			} 

		}
		return this.scopSccsString;
	}
	
	/**
	 * Returns the pdb title string for this Template
	 * The first time this method is called the title is taken from the pdb object.
	 * If PDB data is not loaded yet it will be loaded from database. 
	 * @param conn
	 * @param pdbaseDb
	 * @return the title string
	 */
	public String getTitle(MySQLConnection conn, String pdbaseDb) {
		if (titleString==null) {
			try {
				if (pdb==null) {
					loadPdbData(conn, pdbaseDb);
				}
				this.titleString = pdb.getTitle();
				
			}  catch (SQLException e) {
				System.err.println("Couldn't get the template structure "+id+" because of SQL error: "+e.getMessage());
				pdb = null;
				scopSccsString = "";
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Couldn't get the template structure "+id+" because pdb code was not found");
				pdb = null;
				scopSccsString = "";
			} catch (PdbLoadError e) {
				System.err.println("Couldn't get the template structure "+id+" because of pdb load error: "+e.getMessage());
				pdb = null;
				scopSccsString = "";
			} 			
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
	
	/**
	 * Returns the GTGHit corresponding to this Template or null if this 
	 * Template didn't originate from a GTGHit 
	 * @return
	 */
	public GTGHit getGTGHit() {
		return this.gtgHit;
	}
}

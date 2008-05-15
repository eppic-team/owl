package proteinstructure;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import sequence.BlastHit;
import sequence.GTGHit;
import tools.MySQLConnection;

/**
 * A template protein structure to be used in homology modelling
 *
 */
public class Template {

	private static final String SCOP_VERSION = "1.73";
	
	private static final String PDBASE_DB = "pdbase";
	
	private String id;
	private String scopSccsString;
	private String titleString;
	private Pdb pdb;
	//private MySQLConnection conn;

	private BlastHit blastHit;
	private GTGHit gtgHit;
	
	/**
	 * Constructs a Template given an id in format pdbCode+chain, e.g. 1abcA
	 * @param id
	 */
	public Template(String id) {
		this.id = id;
		checkId();
		//getPdbAndScopString(); 
	}
	
	/**
	 * Constructs a new Template given a BlastHit 
	 * @param hit
	 */
	public Template(BlastHit hit) {
		this.id = hit.getTemplateId();
		this.blastHit = hit;
		//getPdbAndScopString();
	}

	/**
	 * Constructs a new Template given a GTGHit 
	 * @param hit
	 */
	public Template(GTGHit hit) {
		this.id = hit.getTemplateId();
		this.gtgHit = hit;
		//getPdbAndScopString();
	}
	
	/**
	 * Checks that the id complies to our standard pdbCode+chain, e.g. 1abcA
	 */
	private void checkId() {
		Pattern p = Pattern.compile("\\d\\w\\w\\w\\w");
		Matcher m = p.matcher(id);
		if (!m.matches()) {
			throw new IllegalArgumentException("The given template id: "+id+" is not valid. It must be of the form pdbCode+chain, e.g. 1abcA");
		}
	}
	
	/**
	 * Gets the pdb structure coordinates into the pdb object 
	 * @param conn
	 */
	private void getPdbInfo(MySQLConnection conn) throws SQLException, PdbCodeNotFoundError, PdbLoadError {
		String pdbCode = id.substring(0, 4);
		String chain = id.substring(4);

		pdb = new PdbasePdb(pdbCode, PDBASE_DB, conn);
		pdb.load(chain);

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
	public void print() {
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
	 * Returns the pdb object with the actual structure
	 * The first time the method is called the PDB data is loaded from db. 
	 * Subsequent calls take the PDB data from the  cached variable.
	 * @param conn
	 * @return the Pdb object or null if something goes wrong while retrieving the PDB data
	 */
	public Pdb getPdb(MySQLConnection conn) {
		if (pdb==null) {
			try {
				getPdbInfo(conn);
			}  catch (SQLException e) {
				System.err.println("Couldn't get the template structure "+id+" because of SQL error: "+e.getMessage());
				pdb = null;
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Couldn't get the template structure "+id+" because pdb code was not found");
				pdb = null;
			} catch (PdbLoadError e) {
				System.err.println("Couldn't get the template structure "+id+" because of pdb load error: "+e.getMessage());
				pdb = null;
			} 
		}
		return this.pdb;
	}
	
	/**
	 * Tells whether ther is structure data available for this template
	 * @return
	 */
	public boolean hasStructure() {
		return pdb!=null;
	}
	
	/**
	 * Returns the scop sccs string for this Template
	 * The first time this method is called the SCOP data is parsed from file.
	 * Subsequente calls take the SCOP data from the cached variable.
	 * If PDB data is not loaded yet it will be loaded from database.
	 * @param conn
	 * @return
	 */
	public String getScopSccsString(MySQLConnection conn) {
		if (scopSccsString==null) {
			try {
				if (pdb==null) {
					getPdbInfo(conn);
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
	 * @return the title string
	 */
	public String getTitle(MySQLConnection conn) {
		if (titleString==null) {
			try {
				if (pdb==null) {
					getPdbInfo(conn);
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

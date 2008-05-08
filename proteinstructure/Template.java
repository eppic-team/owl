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
	private Pdb pdb;
	private MySQLConnection conn;

	
	/**
	 * Constructs a Template given an id in format pdbCode+chain, e.g. 1abcA
	 * The blast hit data for this template will be missing
	 * @param id
	 */
	public Template(String id) {
		this.id = id;
		checkId();
		getPdbAndScopString(); 
	}
	
	/**
	 * Constructs a new Template given a BlastHit and a MySQLConnection 
	 * from where PDB data will be taken
	 * @param hit
	 * @param conn
	 */
	public Template(BlastHit hit, MySQLConnection conn) {
		this.id = hit.getTemplateId();
		this.conn = conn;
		getPdbAndScopString();
	}

	/**
	 * Constructs a new Template given a GTGHit and a MySQLConnection 
	 * from where PDB data will be taken
	 * @param hit
	 * @param conn
	 */
	public Template(GTGHit hit, MySQLConnection conn) {
		this.id = hit.getTemplateId();
		this.conn = conn;
		getPdbAndScopString();
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
	 * Gets the pdb structure coordinates into the pdb object and all SCOP 
	 * sccs ids in a comma separated string 
	 */
	private void getPdbAndScopString() {
		String pdbCode = id.substring(0, 4);
		String chain = id.substring(4);
		try {
			pdb = new PdbasePdb(pdbCode, PDBASE_DB, conn);
			pdb.load(chain);
			pdb.checkScop(SCOP_VERSION, false);
			Iterator<ScopRegion> it = pdb.getScop().getIterator();
			scopSccsString = "";
			while (it.hasNext()) {
				scopSccsString += it.next().sccs +", ";
			}
			if (scopSccsString.contains(",")) // choping off the last comma if the string is not empty
				scopSccsString = scopSccsString.substring(0, scopSccsString.length()-2); 
		} catch (SQLException e) {
			System.err.println("Couldn't get the template structure "+id+" because of SQL error: "+e.getMessage());
			pdb = null;
		} catch (PdbCodeNotFoundError e) {
			System.err.println("Couldn't get the template structure "+id+" because pdb code was not found");
			pdb = null;
		} catch (PdbLoadError e) {
			System.err.println("Couldn't get the template structure "+id+" because of pdb load error: "+e.getMessage());
			pdb = null;
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
	 * @return
	 */
	public Pdb getPdb() {
		return this.pdb;
	}
	
	/**
	 * Tells whether ther is structure data available for this template
	 * @return
	 */
	public boolean hasStructure() {
		return pdb!=null;
	}
}

package owl.core.structure;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import owl.core.structure.graphs.RIGEnsemble;
import owl.core.util.MySQLConnection;


/**
 * A set of PDB objects providing various input methods.
 * Note that PDB objects are stored in memory. This is not feasible for large sets of
 * structures. In this case the provided static methods to read pdb+chain codes
 * or filenames can be used to process files individually. See also static methods
 * <code>Pdb.readStructureOrExit</code> and <code>Pdb.readStructureOrNull</code> for
 * loading individual structures.
 * 
 * See also: {@link TemplateList}, {@link RIGEnsemble}
 * @author stehr
 * @date 2009-01-09
 *
 */
public class PdbSet {

	/*------------------------------ constants ------------------------------*/
	private static final String CULLPDB_20_FILE = "/project/StruPPi/Databases/cullpdb/cullpdb_pc20_res1.6_R0.25_d090102_chains1498";

	
	/*--------------------------- member variables --------------------------*/
	Collection<Pdb> pdbs;

	/*----------------------------- constructors ----------------------------*/
	/**
	 * Creates an empty PdbSet.
	 */
	public PdbSet() {
		pdbs = new LinkedHashSet<Pdb>();
	}

	/**
	 * Creates a new pdb set filled with the given pdb objects.
	 * @param inPdbs a collection of pdb objects
	 */
	public PdbSet(Collection<Pdb> inPdbs) {
		pdbs = new HashSet<Pdb>();
		pdbs.addAll(inPdbs);
	}

	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Returns a collection view of the pdb objects in this set.
	 */
	public Collection<Pdb> getPdbs() {
		return pdbs;
	}

	/**
	 * Given a collection of pdb+chain codes or file names, loads the respective structures from PDBase.
	 * @throws PdbCodeNotFoundException 
	 * @throws SQLException 
	 * @throws PdbLoadError 
	 */
	public void readFromList(Collection<String> list) throws SQLException, PdbCodeNotFoundException, PdbLoadError {
		for(String pdbChain:list) {
			String pdbCode = pdbChain.substring(0, 4);
			String chainCode = pdbChain.substring(4, 5);
			Pdb pdb = new PdbasePdb(pdbCode);
			pdb.load(chainCode);
			System.out.print(".");
		}
		System.out.println();
	}
	
	/**
	 * Reads all pdb files from the given directory (<i>not implemented yet</i>).
	 * @param dir
	 */
	public void readFromDirectory(File dir) {
		//TODO: not implemented yet
	}
	
	/**
	 * Generates a set of structures from a text file containing pdb files names or pdb+chain codes (<i>not implemented yet</i>).
	 * @param listFile
	 */
	public void readFromListFile(File listFile) {
		//TODO: not implemented yet
	}
	
	/**
	 * Generates a set of structures from a database table of pdb+chain codes (<i>not implemented yet</i>).
	 */
	public void readFromDatabase(MySQLConnection conn, String db, String table, String pdbTable, String chainTable) {
		//TODO: not implemented yet
	}

	/**
	 * Convenience method to load the cullpdb20 database. Currently, a hard-coded version
	 * from January 2009 is being used. Eventually, the latest version could be downloaded
	 * automatically.
	 * @throws IOException 
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundException 
	 * @throws SQLException 
	 */
	public void readCullPdb20() throws IOException, SQLException, PdbCodeNotFoundException, PdbLoadError {
		LinkedList<String> pdbCodes = new LinkedList<String>();
		
		// read pdb+chain codes from file
		BufferedReader in = new BufferedReader(new FileReader(CULLPDB_20_FILE));
		String line = in.readLine(); // skip header line
		while((line = in.readLine()) != null) {
			String pdbCode = line.split("\\s")[0];
			pdbCodes.add(pdbCode);
		}
		in.close();
		
		// load structures from PDBase
		readFromList(pdbCodes);
	}
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Convenience method to load pdb+chain codes of the cullpdb20 database. Currently, a hard-coded version
	 * from January 2009 is being used. Eventually, the latest version could be downloaded
	 * automatically.
	 * @throws IOException 
	 */
	public static Collection<String> readCullPdb20List() throws IOException {
		LinkedList<String> pdbCodes = new LinkedList<String>();
		
		// read pdb+chain codes from file
		BufferedReader in = new BufferedReader(new FileReader(CULLPDB_20_FILE));
		String line = in.readLine(); // skip header line
		while((line = in.readLine()) != null) {
			String pdbCode = line.split("\\s")[0];
			pdbCodes.add(pdbCode);
		}
		in.close();
		
		return pdbCodes;
	}	
	
	/*--------------------------------- main --------------------------------*/
	
	public static void main(String[] args) {
		
		// prints the lengths of the cullpdb20 proteins
		
		Collection<String> pdbCodes = null;
		try {
			pdbCodes = readCullPdb20List();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("Loading " + pdbCodes.size() + " structures...");
		int n = 0;
		for(String pdbCode:pdbCodes) {
			Pdb pdb = Pdb.readStructureOrNull(pdbCode);
			if(pdb != null) {
				System.out.println(pdb.getFullLength());
				n++;
			}
		}
		System.out.println("Done reading " + n + " structures.");
	}
	
}



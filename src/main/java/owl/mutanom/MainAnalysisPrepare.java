package owl.mutanom;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import owl.core.util.MySQLConnection;
import owl.mutanom.core.TargetList;


/**
 * Tools for preparing the data before running MainAnalysis (which then produces the actual results).
 * @author stehr
 */
public class MainAnalysisPrepare {

	/*------------------------------ constants ------------------------------*/
	
	// old (keep here only for reference):
	//public static final File FOLDX_REPAIRED_PDB_DIR = new File("/project/PyMol/pdbs_download/repairPDB/repair");	
	//public static final File FOLDX_RUN_FILE = 		new File("/project/StruPPi/henning/projects/mutanom/analysis/foldx_runs/run_cluster_sc_3.sh");
	//public static final File FOLDX_RESULT_DIR = 		new File("/project/StruPPi/henning/projects/mutanom/analysis/foldx_runs/out_cx_3/");
	//public static final File NACCESS_RESULT_DIR = 	new File("/project/PyMol/pdbs_download/naccess");
	
	// configurables, TODO: load from config file or environment variable (through static method)
	public static final String BASE_DIR =				"/project/StruPPi/henning/projects/mutanom/analysis/";
	public static final String DATABASE = 				"mutanom3"; 				
	
	// paths
	public static final File FOLDX_REPAIRED_PDB_DIR = 	new File(BASE_DIR, "data/pdb/foldx_repair");
	public static final File FOLDX_RUN_FILE = 			new File(BASE_DIR, "data/pdb/foldx_eval/run_cluster_cx.sh");
	public static final File FOLDX_RESULT_DIR = 		new File(BASE_DIR, "data/pdb/foldx_eval/out_cx");
	public static final File NACCESS_RESULT_DIR = 		new File(BASE_DIR, "data/pdb/naccess");
	public static final File PREDICTED_PDB_DIR = 		new File(BASE_DIR, "data/pdb/predicted");
	public static final File NCBI_SNP_XML_DIR =	 		new File(BASE_DIR, "data/snp/ncbi/xml");
	
	// database
	public static final String RES_TABLE = 				DATABASE + ".res_basic";	// table with basic info about all residues in all substructures
	public static final String FOLDX_RESULT_TABLE = 	DATABASE + ".res_foldx";	// per residue (x19) results of FoldX mutant stability calculations
	public static final String NACCESS_RESULT_TABLE = 	DATABASE + ".res_naccess";	// per residue results of NACCESS surface accessibility
	public static final String NCBI_SNP_TABLE = 		DATABASE + ".snps_ncbi";	// snps parsed from NCBI dbSNP xml files
	
	// parameters
	public static boolean 	WHOLE_COMPLEX = false;	// if false, use whole complex (reverse logic)
	public static boolean 	ONLY_PREDICTED = false;	// if true, load only predicted structures	
	
	// aliases
	public static final boolean ONLINE = true;		// whether to use online databases for loading Uniprot and COSMIC sequence
	public static final boolean OFFLINE = false;	// whether to use online databases for loading Uniprot and COSMIC sequence
	/*---------------------------- static methods ---------------------------*/

	/**
	 * Creates a table with all residues in all substructures. The table is needed
	 * by collectFoldXResults and collectNaccessResults.
	 * @throws SQLException 
	 */
	public static void createResidueTable(MySQLConnection conn, TargetList targets) throws SQLException {
		
		targets.createAllPositionsDbTable(conn, RES_TABLE);
		
	}
	
	/**
	 * Writes the scripts for running FoldX energy calculations on the cluster.
	 * The FoldX repair step has to be run manually before running this method.
	 * See: /project/StruPPi/henning/projects/mutanom/analysis/data/pdb/README
	 */
	public static void prepareFoldXClusterRun(TargetList targets) {
				
		targets.writeFoldXMutScanMasterShellScript(FOLDX_RUN_FILE, FOLDX_REPAIRED_PDB_DIR, FOLDX_RESULT_DIR, WHOLE_COMPLEX);
		System.out.println("FoldXMutScanMasterScriptFile written to: " + FOLDX_RUN_FILE.getPath());	
	}

	/**
	 * Collects the results of a FoldX cluster run and writes them to the database.
	 * Uses default values for paths and database table. FoldX has to be run manually
	 * using the runscript created with prepareFoldXClusterRun().
	 * @param targets
	 * @param conn
	 * @throws SQLException 
	 */
	public static void collectFoldXResults(MySQLConnection conn, TargetList targets) throws SQLException {

		System.out.println("Loading FoldX results from: " + FOLDX_RESULT_DIR.getPath() + " to " + FOLDX_RESULT_TABLE);
		
		targets.collectFoldXResults(conn, RES_TABLE, FOLDX_RESULT_TABLE, FOLDX_RESULT_DIR, WHOLE_COMPLEX);

	}
	
	/**
	 * Collects the results of the NACCESS calculations and writes them to the database.
	 * Creates the target table if not exists. Uses default values for paths and database table.
	 * NACCESS has to be run manually  in /project/StruPPi/henning/projects/mutanom/analysis/data/pdb/
	 * @throws SQLException on errors when reading from or writing to database
	 * @throws IOException on errors when reading from NACCESS output files
	 */
	public static void collectNaccessResults(MySQLConnection conn, TargetList targets) throws SQLException, IOException {
				
		 // parse the exposure for all positions and load to database
		targets.loadExposureToDb(conn, RES_TABLE, NACCESS_RESULT_TABLE, NACCESS_RESULT_DIR);			
	}
	
	public static void collectNcbiSnps(MySQLConnection conn, TargetList targets) throws SQLException {
		
		targets.collectNcbiSnpsFromXmlFiles(conn, NCBI_SNP_TABLE, NCBI_SNP_XML_DIR);
		
	}
	
	// TODO: Is this still useful?
	// write proximity to functional sites to db (to get background distribution)
//	try {
//		targets.writeProximitiesToDb(conn);
//	} catch (SQLException e) {
//		e.printStackTrace();
//	}
	
	/*--------------------------------- main --------------------------------*/
	
	/**
	 * Runs the tools for preparing the data for mutanom analysis.
	 * TODO: Add command line parameters to execute the different steps:
	 * - r create residue table
	 * - n load results of NACCESS calculations to database
	 * - f create FoldX run script for energy calculations on cluster 
	 * - e load results of energy calculations to database
	 * @throws SQLException
	 * @throws IOException 
	 */
	public static void main(String[] args) throws SQLException, IOException {
		
		// suppress logger output (to get rid of JAligner garbage)
		File trashLogFile = File.createTempFile("logger", "trash");
		trashLogFile.deleteOnExit();
		System.setProperty("java.util.logging.config.file",trashLogFile.getAbsolutePath());
		
		MySQLConnection conn = new MySQLConnection();
		
		System.out.println("Loading targets...");
		//TargetList targets = TargetList.loadSingletonList(conn, OFFLINE, "NTRK3");
		//TargetList targets = TargetList.loadTargets(conn, OFFLINE);
		TargetList targets = TargetList.loadPredictedTargets(conn, OFFLINE);
		System.out.println("Loading substructures...");
		targets.loadSubstructures(conn, ONLY_PREDICTED);
		System.out.println("Loading alignments...");
		targets.loadPdbsAndAlignments(PREDICTED_PDB_DIR);
		
		// createResidueTable(conn, targets);
	
		// collectNaccessResults(conn, targets);
		
		// prepareFoldXClusterRun(targets);
		
		collectFoldXResults(conn, targets);
		
		// collectNcbiSnps(conn, targets);
		
		System.out.println("done.");
	}	
}

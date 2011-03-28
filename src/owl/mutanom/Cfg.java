package owl.mutanom;

import java.io.File;

import owl.core.features.StructuralDomainType;

/**
 * Main configuration object for defining static constants and parameters. Use this object for
 * things which are shared by the classes MainAnalysis, MainAnalysisPrepare, GeneInfo.
 * TODO: Currently only used by GeneInfo, merge with MainAnalysis, MainAnalysisPrepare
 * @author stehr
 */
public class Cfg {
	
	/*------------------------------ constants ------------------------------*/
	
	// configurables
	public static final String BASE_DIR =				"/project/StruPPi/henning/projects/mutanom/analysis/";
	public static final String DATABASE = 				"mutanom2"; 
	
	// paths
	public static final File PHOSITE_HTML_PATH = 		new File(BASE_DIR, "data/phosite/html");
	public static final File PDBS_PATH = 				new File(BASE_DIR, "data/pdb/download");	
	public static final File PRED_PDBS_PATH = 			new File(BASE_DIR, "data/pdb/predicted");
	public static final File RESULT_DIR	= 				new File(BASE_DIR, "results");
	public static final File PSE_OUT_DIR = 				new File(RESULT_DIR, "pse");
	
	// database
	public static final String RES_TABLE = 				DATABASE + ".res_basic";	// table with basic info about all residues in all substructures
	public static final String FOLDX_RESULT_TABLE = 	DATABASE + ".res_foldx";	// per residue (x19) results of FoldX mutant stability calculations
	public static final String NACCESS_RESULT_TABLE = 	DATABASE + ".res_naccess";	// per residue results of NACCESS surface accessibility
	
	public static final String QUERY_RES_IDX = 			"SELECT idx FROM " + RES_TABLE + " WHERE pdb_code='%s' AND pdb_chain='%s' AND cif_pos=%d";
	public static final String QUERY_NACCESS_RESULT =	"SELECT rsa_cx FROM " + NACCESS_RESULT_TABLE + " WHERE res_idx=%s";
	public static final String QUERY_FOLDX_RESULT = 	"SELECT ddg_cx FROM " + FOLDX_RESULT_TABLE + " WHERE res_idx=%s AND aa_mut='%s'";
	public static final String QUERY_COUNT_DESTAB = 	"SELECT COUNT(*) FROM " + RES_TABLE + " AS r INNER JOIN " + FOLDX_RESULT_TABLE + " AS f ON r.idx=f.res_idx " +
															"WHERE pdb_code=\"%s\" AND pdb_chain=\"%s\" AND ddg_cx > %f";
	public static final String QUERY_COUNT_NEUTRAL = 	"SELECT COUNT(*) FROM " + RES_TABLE + " AS r INNER JOIN " + FOLDX_RESULT_TABLE + " AS f ON r.idx=f.res_idx " +
															"WHERE pdb_code=\"%s\" AND pdb_chain=\"%s\" AND ddg_cx <= %f";	
	
	// aliases
	public static final boolean ONLINE = true;			// for select online vs. offline access of databases (currently only Uniprot/Cosmic sequence)
	public static final boolean OFFLINE = false;		
	
	// parameters
	public static final String 	RIG_TYPE = 					"Ca";		// for func site proximity
	public static final double 	RIG_CUTOFF = 				8.0;		// for func site proximity
	public static final String 	CLUSTERING_ATOM = 			"Centroid";	// for dist calculation, can be "CA", "CB", or "Centroid"
	public static final boolean CLUST_WITH_REPLACEMENT = 	true; 		// randomizeMutations for clustering
	public static final boolean FUNSITE_WITH_REPLACEMENT = 	true;		// randomizeMutations for random feature count
	public static final int 	NSAMPLES = 					10000;		// sample size for random population
	public static final boolean SEQ_SOURCE = 				OFFLINE;	// whether to use online databases for loading Uniprot and COSMIC sequence
	public static final boolean	WHOLE_COMPLEX = 			false;		// if false, use whole complex for stability and surface (reverse logic)
	public static final double 	EXPOSURE_THRESHOLD = 		15.0;		// cutoff for counting a residue as exposed
	public static final double 	DESTAB_THRESHOLD = 		 	5.0;		// ddg cutoff for counting a mutation as destabilizing
	public static final StructuralDomainType DOMAIN_TYPE = 	StructuralDomainType.NCBI;	// method for domain definition: NCBI (before: DomainParser)

}

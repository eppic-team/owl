package owl.mutanom;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
//import java.sql.SQLException;
//import java.sql.ResultSet;
//import java.sql.SQLException;
//import java.sql.Statement;
import java.sql.SQLException;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.TreeSet;

import javax.vecmath.Point3d;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math.MathException;
import org.apache.commons.math.stat.StatUtils;

import owl.core.features.Feature;
import owl.core.features.FeatureType;
//import owl.core.features.GeneralFeature;
//import owl.core.features.PhosphoSitePlusFeature;
//import owl.core.features.ProteinModificationType;
import owl.core.features.StructuralDomainType;
//import owl.core.features.UniprotFeature;
import owl.core.structure.AminoAcid;
import owl.core.structure.Residue;
import owl.core.structure.graphs.RIGNode;
import owl.core.util.MySQLConnection;
import owl.mutanom.core.Gene;
import owl.mutanom.core.Mutation;
import owl.mutanom.core.Substructure;
import owl.mutanom.core.TargetList;
import owl.mutanom.core.Substructure.SubstructureType;

import edu.uci.ics.jung.graph.util.Pair;



/**
 * Created this class for its main method as TargetList was getting too cluttered
 * @author stehr
 * See also: TargetList.main()
 */
public class MainAnalysis {

	/*------------------------------ constants ------------------------------*/
	
	// configurables, TODO: use config object to load from config file or environment variable (through static method)
	public static final String BASE_DIR =				"/project/StruPPi/henning/projects/mutanom/analysis/";
	public static final String DATABASE = 				"mutanom3"; 
	
	// paths
	public static final File PHOSITE_HTML_PATH = 		new File(BASE_DIR, "data/phosite/html");
	public static final File PRED_PDBS_PATH = 			new File(BASE_DIR, "data/pdb/predicted");
	public static final File RESULT_DIR	= 				new File(BASE_DIR, "results");

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
	public static final String 	RIG_TYPE = 					"Ca";		// for func site proximity (why C-alpha? shall we use all atom?)
	public static final double 	RIG_CUTOFF = 				8.0;		// for func site proximity (4.1 would be 1.0 + vdw radius)
	public static final String 	CLUSTERING_ATOM = 			"Centroid";	// for dist calculation, can be "CA", "CB", or "Centroid"
	//public static final String 	CLUSTERING_ATOM = 			"CA";	// for dist calculation, can be "CA", "CB", or "Centroid"
	public static final boolean CLUST_WITH_REPLACEMENT = 	true; 		// randomizeMutations for clustering
	public static final boolean FUNSITE_WITH_REPLACEMENT = 	true;		// randomizeMutations for random feature count
	public static final int 	NSAMPLES = 					100;		// sample size for random population
	public static final int 	NSAMPLES_TP53 =				100;		// sample size for random population for TP53 (has >800 mut and would be too slow otherwise)
	public static final boolean SEQ_SOURCE = 				OFFLINE;	// whether to use online databases for loading Uniprot and COSMIC sequence
	public static final boolean	WHOLE_COMPLEX = 			false;		// if false, use whole complex for stability and surface (reverse logic)
	public static final double 	EXPOSURE_THRESHOLD = 		15.0;		// cutoff for counting a residue as exposed
	public static final double 	DESTAB_THRESHOLD = 		 	5.0;		// ddg cutoff for counting a mutation as destabilizing
	public static final StructuralDomainType DOMAIN_TYPE = 	StructuralDomainType.NCBI;	// method for domain definition: NCBI (before: DomainParser)
	public static final	boolean USE_ODD_RATIO =				true;		// use observed/expected for scatterplots
	public static final boolean USE_MANUAL_FUNC_ANNO =	    false;		// if true, also load manual annotations for functional sites from file
	
	/*--------------------------- member variables --------------------------*/
	static long timerStart;
	static long timerLastCheckpoint;
	
	/*---------------------------- static methods ---------------------------*/
	
	public static void startTimer() {
		timerStart = System.currentTimeMillis();
		timerLastCheckpoint = timerStart;
	}
	
	public static void takeTime() {
		long current = System.currentTimeMillis();
		System.out.printf("*** TIMER: Since start: %.0f sec\t Since last checkpoint: %.0f sec\n", (current-timerStart)/1000.0, (current-timerLastCheckpoint)/1000.0 );
		timerLastCheckpoint = current;
	}
	
	/*--------------------------------- main --------------------------------*/
	
	/**
	 * Load data for single structure. Used to debug unusual values > 1 for NRAS and BRAF
	 */
	public static void debugSingleStructure() throws SQLException {
		MySQLConnection conn = TargetList.connectToDb();
		TargetList alltargets = TargetList.loadTop29TargetList(conn, SEQ_SOURCE);
		TargetList singleton = new TargetList();
		Gene nras = null;
		for(Gene g:alltargets.getTargets()) {
			if(g.getGeneName().equals("BRAF")) {
				singleton.addTarget(g);
				nras = g;
			}
		}
		singleton.loadSubstructures(conn, false);	// not only predicted
		singleton.loadPdbsAndAlignments(null);
		singleton.loadFeaturesAndGraphsAndDomains(PHOSITE_HTML_PATH, RIG_TYPE, RIG_CUTOFF, DOMAIN_TYPE, USE_MANUAL_FUNC_ANNO);
		singleton.loadMutationsFromDb(conn);
		singleton.loadSnpsFromDb(conn);
		System.out.println(nras.getMutations());
		singleton.printStats(conn);
		System.out.println(reportFeatureCount(singleton, false, Double.NaN, null, null, null));
		conn.close();
	}
	
	/**
	 * Analyze predicted structures (ERBB2 and ...) to test prediction of
	 * oncogene/tumor suppressor
	 * @throws FileNotFoundException 
	 * @throws SQLException 
	 */
	public static void analyzePredictedStructures() throws FileNotFoundException, SQLException {

		MySQLConnection conn = TargetList.connectToDb();

		System.out.println("Loading predicted targets...");
		TargetList targets = TargetList.loadPredictedTargets(conn, SEQ_SOURCE);
		System.out.println(targets.getTargets().size() + " targets loaded.");
		
		System.out.println("Loading substructures...");
		boolean onlyPredicted = false;
		targets.loadSubstructures(conn, onlyPredicted);
		
		// remove non-predicted substructures (not necessary anymore!)
//		int numSS = 0;
//		for(Gene g:targets.getTargets()) {
//			LinkedList<Substructure> toDelete = new LinkedList<Substructure>();
//			for(Substructure ss:g.getSubstructures()) {
//				if(ss.getType() != SubstructureType.PREDICTION) {
//					toDelete.add(ss);
//				} else {
//					numSS++;
//				}
//			}
//			for(Substructure d:toDelete) {
//				g.getSubstructures().remove(d);
//			}
//		}
//		System.out.println(numSS + " substructures loaded.");
		
		//File pdbDir = new File("/home/web/lappe/stehr/mutanom/pdb/");
		//File pdbDir = new File("/project/PyMol/pdbs_download/pred_not_renum/");	// for predicted structures
		targets.loadPdbsAndAlignments(PRED_PDBS_PATH);
		targets.loadFeaturesAndGraphsAndDomains(PHOSITE_HTML_PATH, RIG_TYPE, RIG_CUTOFF, DOMAIN_TYPE, USE_MANUAL_FUNC_ANNO);
		
		targets.loadMutationsFromDb(conn);
		targets.loadSnpsFromDb(conn);		
		targets.printStats(conn);
		//targets.writeOriginalDataFiles(RESULT_DIR, true);
		
		// set up once
		
//		// create database table for results (exposure, functional sites)
//		try {
//			targets.createAllPositionsDbTable(conn, "mutanom.all_observed");
//		} catch (SQLException e) {
//			System.err.println("Could not write positions to database: " + e.getMessage());
//		}
		
		// write surface accessibilities to DB:
		//before: run NACCESS manually on all pdb structure (complex and single chain)
//		try {
//			targets.loadExposureToDb(conn, "mutanom.all_observed", new File("/project/PyMol/pdbs_download/naccess")); // parse the exposure for all positions and load to database
//		} catch(IOException e) {
//			e.printStackTrace();
//		} catch(SQLException e) {
//			e.printStackTrace();
//		}
				
		// analysis
	
		System.out.println("5d. Testgenes");
		boolean singleChain = WHOLE_COMPLEX;	// if true, use single structures from pdbase otherwise, use original pdb files (whole complex)
		double destabThresh = DESTAB_THRESHOLD;	// stability changes above this threshold are considered destabilizing, others neutral
		double exposureThreshold = EXPOSURE_THRESHOLD;		
		
		File resultFile = new File(RESULT_DIR, "5d_testgenes.txt");
		PrintStream out = new PrintStream(resultFile);
		
//		System.out.println("Gene\t#Mut\tExp\tSites\tClust");
//		out.println("Gene\t#Mut\tExp\tSites\tClust");
//		for(Gene g:targets.getTargets()) {
//			System.out.printf("%s\t", g.getGeneName());
//			out.printf("%s\t", g.getGeneName());
//			TargetList singleton = new TargetList();
//			singleton.addTarget(g);
//			// number of mutations
//			int numMut = g.getMutationsInStructures(true, false).size();
//			System.out.printf("%4d\t", numMut);
//			out.printf("%4d\t", numMut);			
//			// exposure
//			Pair<Integer> buriedExposed = countBuriedExposed(conn, singleton, exposureThreshold, singleChain);
//			int numBuried = buriedExposed.getFirst();
//			int numExposed = buriedExposed.getSecond();		
//			double fractionExposed = 1.0 * numExposed / (numExposed+numBuried);
//			System.out.printf("%4.2f\t", fractionExposed);
//			out.printf("%4.2f\t", fractionExposed);
//			// functional sites
//			double fracProxFunctionalSites = reportFeatureCount(singleton, true, Double.NaN, null);
//			System.out.printf("%4.2f\t", fracProxFunctionalSites);
//			out.printf("%4.2f\t", fracProxFunctionalSites);
//			// clustering
//			double clust = calculateClustering(singleton);
//			System.out.printf("%4.2f\t", clust);
//			out.printf("%4.2f\t", clust);
//			
//			System.out.println();
//			out.println();
//		}
		reportAllPerGene(targets, out, conn, singleChain, exposureThreshold, destabThresh);
		
		// ------- generate image for publiation -------
//		File outDir = new File("/home/web/lappe/stehr/mutanom/pub/pse");
//		try {
//			targets.generatePseWithAllMutationsAndSNPs(outDir, pdbDir);
//		} catch (IOException e) {
//			System.err.println("Error creating PyMol session");
//		}
		
		conn.close();
	}

	/**
	 * Keep this only for historic purposes. This was the old data loading part in
	 * analyseMutations(). We may still want to look up why certain genes were
	 * previously exluded from analysis. Otherwise this method does nothing.
	 */
	@SuppressWarnings("null")
	public static void loadDataOld() {
		System.out.println("Loading top20 mutanom targets...");
		//TargetList targets = TargetList.getTop20TargetList(SEQ_SOURCE);

		TargetList targets = null;
		
		// skip SRC because there are no missense mutations in COSMIC and many sequence mismatches
		if(targets.removeTarget("SRC") == null) {
			System.err.println("Could not remove target SRC from list");
		}
		// skip MSH6 because there are no missense mutations in COSMIC
		if(targets.removeTarget("MSH6") == null) {
			System.err.println("Could not remove target MSH6 from list");
		}
		// skip SMO because the only "structured" region is prediction garbage
		if(targets.removeTarget("SMO") == null) {
			System.err.println("Could not remove target SMO from list");
		}
		// skip ERBB2 because there are no missense mutations with known structure in COSMIC, Mut:12 Snp:7
		if(targets.removeTarget("ERBB2") == null) {
			System.err.println("Could not remove target ERBB2 from list");
		}		
		// skip MLH1 because no known structures are available (only 1 predicted) Mut: 5 Snp:100
		if(targets.removeTarget("MLH1") == null) {
			System.err.println("Could not remove target MLH1 from list");
		}
		// skip MSH2 because the 146 SNPs would dominate the whole statistics
		if(targets.removeTarget("MSH2") == null) {
			System.err.println("Could not remove target MSH2 from list");
		}
//		// skip EGFR because the 168 mutations would dominate the whole statistics (for 5 tissues)
//		if(targets.removeTarget("EGFR") == null) {
//			System.err.println("Could not remove target EGFR from list");
//		}
		// skip KIT because the PDB files are not prepared yet
		if(targets.removeTarget("KIT") == null) {
			System.err.println("Could not remove target KIT from list");
		}
		// skip VHL because the PDB files are not prepared yet
		if(targets.removeTarget("VHL") == null) {
			System.err.println("Could not remove target VHL from list");
		}
		// skip PRPS1 because the PDB files are not prepared yet
		if(targets.removeTarget("PRPS1") == null) {
			System.err.println("Could not remove target PRPS1 from list");
		}
		// skip NUP133 because the PDB files are not prepared yet
		if(targets.removeTarget("NUP133") == null) {
			System.err.println("Could not remove target NUP133 from list");
		}
		// skip GRK6 because the PDB files are not prepared yet
		if(targets.removeTarget("GRK6") == null) {
			System.err.println("Could not remove target GRK6 from list");
		}		
		// skip IDH1 because the PDB files are not prepared yet
		if(targets.removeTarget("IDH1") == null) {
			System.err.println("Could not remove target IDH1 from list");
		}
		// skip FGFR3 because the PDB files are not prepared yet
		if(targets.removeTarget("FGFR3") == null) {
			System.err.println("Could not remove target FGFR3 from list");
		}
		// skip MMP2 because the PDB files are not prepared yet
		if(targets.removeTarget("MMP2") == null) {
			System.err.println("Could not remove target MMP2 from list");
		}
		
		// skipping all structures with <= 5 mutations so that APC and NF1 don't mess up classification
		
		// skip APC because it has only 6 mutations
		if(targets.removeTarget("APC") == null) {
			System.err.println("Could not remove target APC from list");
		}	
		// skip NF1 because it has only 3 mutations 
		if(targets.removeTarget("NF1") == null) {
			System.err.println("Could not remove target NF1 from list");
		}
		// skip CTNNB1 because it has only 5 mutations
		if(targets.removeTarget("CTNNB1") == null) {
			System.err.println("Could not remove target CTNNB1 from list");
		}
		// skip NRAS because it has only 5 mutations
		if(targets.removeTarget("NRAS") == null) {
			System.err.println("Could not remove target NRAS from list");
		}
	}
	
	/**
	 * Main method which does the structural impact analysis and outputs the results to stdout and text files.
	 * Source data:
	 * - Genes, Substructures, Mutations, SNPs from database 'mutanom2'
	 * - Uniprot and COSMIC Sequences from text files
	 * - PDB structures from local files
	 * - Sifts mapping from Substructure.LOCAL_SIFTS_FILE
	 * - Feature annotations from
	 *   * Uniprot (online)
	 *   * CSA (local)
	 *   * PhosphoSite (local HTML files)
	 *   * Structural domains from pDomains (online)
	 * Side effects:
	 * - write text files to /project/StruPPi/henning/projects/mutanom/analysis/results
	 * 
	 * @throws FileNotFoundException
	 * @throws SQLException
	 */
	public static void analyzeMutations() throws FileNotFoundException, SQLException {
		
		// TODO: nothing (just to bookmark this code in eclipse)
		
		MySQLConnection conn = TargetList.connectToDb();
		
		startTimer();
		
		// loading targets
		System.out.println("Loading targets...");
		TargetList targets = TargetList.loadTargets(conn, SEQ_SOURCE);
		System.out.println(targets.getTargets().size() + " targets loaded.");
		//TargetList targets = TargetList.getSingletonList(conn, SEQ_SOURCE, "PIK3CA");
		
//		targets.removeTarget("EGFR"); // causing problems with clustering
//		targets.removeTarget("KIT");  // causing problems with clustering
		
		// loading substructures
		System.out.println("Loading substructures...");
		boolean onlyPredicted = false;	// here we want (only) xray structures
		targets.loadSubstructures(conn, onlyPredicted);
		
		// loading PDBs, alingments, features and graphs
		targets.loadPdbsAndAlignments(PRED_PDBS_PATH); // predicted structures not currently used
		takeTime();
		targets.loadFeaturesAndGraphsAndDomains(PHOSITE_HTML_PATH, RIG_TYPE, RIG_CUTOFF, DOMAIN_TYPE, USE_MANUAL_FUNC_ANNO);

		// loading mutations and SNPs
		System.out.println();
		System.out.println("Loading mutations...");
		targets.loadMutationsFromDb(conn);
		targets.loadSnpsFromDb(conn);
		
		// print statistics
		targets.printStats(conn);
		targets.printFunctionalSiteStats();
		
		// System.exit(1);
		
		// filter out unwanted mutations (this should mostly not be necessary because we filter our source data)
		System.out.println();
		System.out.println("Non missense: " + targets.filterMutationsMissense());
		System.out.println("No structure: " + targets.filterMutationsKnownStructure());
		System.out.println("Not observed: " + targets.filterMutationsObserved());
		System.out.println("Pred Struct : " + targets.filterMutationsPredicted());
		System.out.println();
		System.out.println("Mutations left: " + targets.getNumberOfMutations());
		System.out.println("Unique positions: " + targets.getNumberOfMutatedPositions());
		System.out.println();
		
		// now reload everything?
		targets.loadMutationsFromDb(conn);
		
		TargetList oncoGenes = targets.getOncogenes(conn);				// oncogenes (according to db table genes_onc_sup)
		TargetList tumSupGenes = targets.getTumorSuppressors(conn);		// tum sup genes (according to db table genes_onc_sup)
		TargetList otherGenes = targets.getUnAnnotated(conn);			// all others
		System.out.println("Oncogenes:\t\t" + oncoGenes.getTargets().size());
		System.out.println("Tumor Suppressors:\t" + tumSupGenes.getTargets().size());
		System.out.println("Others:\t\t\t" + otherGenes.getTargets().size());
		System.out.println();
		
		// ------ Generate images for publication ------------
		
		// TODO: Move these to other class (or: always write them)
		
		// 1. sequence overview for web
//		File outDir = new File("/home/web/lappe/stehr/mutanom/pub/gene");
//		SubstructureType restrictToType = SubstructureType.XRAY;
//		System.out.println("Writing data files for sequence overview to " + outDir);
//		writeSeqOverviewForWeb(outDir, targets, restrictToType);
		
		// 2. sequence overview for print (SVG)
		// File outFile = new File("/project/StruPPi/henning/projects/mutanom/analysis/results/saved_results/test/figure7_seqoverview.svg");
		File outFile = new File(RESULT_DIR, "svg/seqoverview_onc.svg");
		System.out.println("Writing " + outFile);
		oncoGenes.writeSeqOverviewForPrint(outFile);
		File outFile2 = new File(RESULT_DIR, "svg/seqoverview_sup.svg");
		System.out.println("Writing " + outFile2);
		tumSupGenes.writeSeqOverviewForPrint(outFile2);
		System.out.println();
		
		// 3. write data files for FTP download
		File outDir = new File(RESULT_DIR,"supplementaries");
		targets.writeOriginalDataFiles(outDir, false);
		
		System.exit(1);
		
		PrintStream out;
		PrintStream out2;
		
		File resultFile = null;
		
		takeTime();
		
		// 0. ------ Per gene statistics ------------
		targets.loadMutationsFromDb(conn);
		System.out.println();
		System.out.println("5a. Oncogenes");
		resultFile = new File(RESULT_DIR, "5a_oncogenes.txt");
		out = new PrintStream(resultFile);
		reportAllPerGene(oncoGenes, out, conn, WHOLE_COMPLEX, EXPOSURE_THRESHOLD, DESTAB_THRESHOLD);
		out.close();
		
		System.out.println("5b. Tumour suppressors");
		resultFile = new File(RESULT_DIR, "5b_tumoursuppressors.txt");
		out = new PrintStream(resultFile);		
		reportAllPerGene(tumSupGenes, out, conn, WHOLE_COMPLEX, EXPOSURE_THRESHOLD, DESTAB_THRESHOLD);
		out.close();
		
		System.out.println("5c. Others (unknown or not annotated)");
		resultFile = new File(RESULT_DIR, "5c_others.txt");
		out = new PrintStream(resultFile);		
		reportAllPerGene(otherGenes, out, conn, WHOLE_COMPLEX, EXPOSURE_THRESHOLD, DESTAB_THRESHOLD);
		out.close();
		
		takeTime();
		System.out.println();
		
		//System.exit(1);
		
		// -- Dropping targets from average analysis -- 
		
//		System.out.println("Dropping targets TP53, MSH2...");
//		targets.removeTarget("TP53"); tumSupGenes.removeTarget("TP53");
//		targets.removeTarget("MSH2"); tumSupGenes.removeTarget("MSH2");
		
//		System.out.println("Dropping targets NRAS, KRAS...");
//		targets.removeTarget("NRAS"); tumSupGenes.removeTarget("NRAS");
//		targets.removeTarget("KRAS"); tumSupGenes.removeTarget("KRAS");		
		
		// 1. ------- Exposure -------
		resultFile = new File(RESULT_DIR, "1_exposure.txt");
		out = new PrintStream(resultFile);
		System.out.println("1. Exposure");
		double exposureThreshold = EXPOSURE_THRESHOLD;	// 15.0
		boolean singleChain = false;
		// Rnd
		System.out.println("RND");
		System.out.print("r_a: ");
		targets.loadMutationsFromDb(conn);	// just to be safe
		double expAll = countFractionExposed(conn, targets, exposureThreshold, singleChain, out);		
		System.out.print("r_o: ");
		targets.loadMutationsFromDb(conn);
		double expOnc = countFractionExposed(conn, oncoGenes, exposureThreshold, singleChain, out);
		System.out.print("r_s: ");
		targets.loadMutationsFromDb(conn);
		double expSup = countFractionExposed(conn, tumSupGenes, exposureThreshold, singleChain, out);
		// Snp
		System.out.println("SNP");
		System.out.print("s_a: ");
		targets.loadSNPsAsMutations(conn);
		reportExposure(conn, targets, exposureThreshold, singleChain, expAll, out);
		System.out.print("s_o: ");
		targets.loadSNPsAsMutations(conn);
		reportExposure(conn, oncoGenes, exposureThreshold, singleChain, expOnc, out);
		System.out.print("s_s: ");
		targets.loadSNPsAsMutations(conn);
		reportExposure(conn, tumSupGenes, exposureThreshold, singleChain, expSup, out);		
		// Mut
		targets.loadMutationsFromDb(conn);
		System.out.println("MUT");
		System.out.print("m_a: ");
		reportExposure(conn, targets, exposureThreshold, singleChain, expAll, out);
		System.out.print("m_o: ");
		reportExposure(conn, oncoGenes, exposureThreshold, singleChain, expOnc, out);		
		System.out.print("m_s: ");
		reportExposure(conn, tumSupGenes, exposureThreshold, singleChain, expSup, out);

//		System.out.print("Nuc: ");
//		reportExposure(conn, nuclearGenes, exposureThreshold, singleChain, expAll, out);
//		System.out.print("Cyt: ");
//		reportExposure(conn, cytoplasmicGenes, exposureThreshold, singleChain, expAll, out);
//		System.out.print("Mem: ");
//		reportExposure(conn, membraneGenes, exposureThreshold, singleChain, expAll, out);
//		System.out.print("Trm: ");
//		reportExposure(conn, transMembraneGenes, exposureThreshold, singleChain, expAll, out);
		out.close();
		System.out.println();
		
		takeTime();
		//System.exit(1);
		
		// 2. ------- Clustering -------
		resultFile = new File(RESULT_DIR, "2_clustering.txt");
		out = new PrintStream(resultFile);
		File clustPopFile = new File(RESULT_DIR, "2a_clustering_random_population.txt");
		out2 = new PrintStream(clustPopFile);
		System.out.println("2. Clustering");
		targets.loadMutationsFromDb(conn);
		// per gene
//		System.out.println("All genes:");
//		reportClusteringPerGene(targets);
//		System.out.println("Oncogenes:");
//		reportClusteringPerGene(oncoGenes);
//		System.out.println("Tumor suppressors:");
//		reportClusteringPerGene(tumSupGenes);
		int nSamples = NSAMPLES;
		System.out.println("nSamples = " + nSamples);
		// RND
		System.out.println("RND");
		System.out.print("r_a: ");
		targets.loadMutationsFromDb(conn);
		Double[] rndPop = getClusteringRandomPopulation(targets, nSamples, out, out2); // scrambles mutations
		System.out.print("r_o: ");
		targets.loadMutationsFromDb(conn);
		Double[] rndPopOnc = getClusteringRandomPopulation(oncoGenes, nSamples, out, null); // scrambles mutations
		System.out.print("r_s: ");
		targets.loadMutationsFromDb(conn);
		Double[] rndPopSup = getClusteringRandomPopulation(tumSupGenes, nSamples, out, null); // scrambles mutations
		out2.close();
		// SNP
		System.out.println("SNP");
		System.out.print("s_a: ");
		targets.loadSNPsAsMutations(conn);	
		reportClustering(targets, rndPop, out);	
		System.out.print("s_o: ");
		targets.loadSNPsAsMutations(conn);
		reportClustering(oncoGenes, rndPopOnc, out);	
		System.out.print("s_s: ");
		targets.loadSNPsAsMutations(conn);
		reportClustering(tumSupGenes, rndPopSup, out);	
		// MUT
		System.out.println("MUT");
		System.out.print("m_a: ");		
		targets.loadMutationsFromDb(conn); // is this necessary every time?
		reportClustering(targets, rndPop, out);
		System.out.print("m_o: ");	
		targets.loadMutationsFromDb(conn);
		reportClustering(oncoGenes, rndPopOnc, out);
		System.out.print("m_s: ");	
		targets.loadMutationsFromDb(conn);
		reportClustering(tumSupGenes, rndPopSup, out);
//		targets.loadMutationsFromDb(conn);
//		System.out.print("Nuc: ");
//		reportClustering(nuclearGenes, rndPop, out);
//		targets.loadMutationsFromDb(conn);
//		System.out.print("Cyt: ");
//		reportClustering(cytoplasmicGenes, rndPop, out);
//		targets.loadMutationsFromDb(conn);
//		System.out.print("Mem: ");
//		reportClustering(membraneGenes, rndPop, out);
//		targets.loadMutationsFromDb(conn);
//		System.out.print("Trm: ");
//		reportClustering(transMembraneGenes, rndPop, out);
		out.close();
		System.out.println();
		
		takeTime();
		// System.exit(1);
		
		// 3. ------- Functional sites -------
		resultFile = new File(RESULT_DIR, "3_functionalsites.txt");
		File siteDistFile = new File(RESULT_DIR, "3a_functionalsite_distribution.txt");
		out = new PrintStream(resultFile);
		out2 = new PrintStream(siteDistFile);
		System.out.println("3. Functional sites");
		int nSamples2 = NSAMPLES_TP53;
		System.out.println("nSamples = " + nSamples2);
		System.out.println("RND");
		System.out.println("r_a: ");
		targets.loadMutationsFromDb(conn);
		double fracHits = reportRandomFeatureCount(targets, nSamples2, out, out2, "r_a"); // after this call, mutations will be scrambled!
		System.out.println("r_o: ");
		targets.loadMutationsFromDb(conn);
		double fracHitsOnc = reportRandomFeatureCount(oncoGenes, nSamples2, out, out2, "r_o"); // after this call, mutations will be scrambled!
		System.out.println("r_s: ");
		targets.loadMutationsFromDb(conn);
		double fracHitsSup = reportRandomFeatureCount(tumSupGenes, nSamples2, out, out2, "r_s"); // after this call, mutations will be scrambled!
//		System.out.print("All: ");
//		double fracHits = targets.countFeaturesAllResidues(false, out);
//		System.out.println("a_o:");
//		oncoGenes.countFeaturesAllResidues(false, null);
//		System.out.println("a_s:");
//		tumSupGenes.countFeaturesAllResidues(false, null);
		System.out.println("SNP");
		System.out.println("s_a: ");		
		targets.loadSNPsAsMutations(conn);
		reportFeatureCount(targets, false, fracHits, out, out2, "s_a");
		System.out.println("s_o: ");		
		targets.loadSNPsAsMutations(conn);
		reportFeatureCount(oncoGenes, false, fracHitsOnc, out, out2, "s_o");
		System.out.println("s_s: ");		
		targets.loadSNPsAsMutations(conn);
		reportFeatureCount(tumSupGenes, false, fracHitsSup, out, out2, "s_s");		
		System.out.println("MUT");
		System.out.println("m_a: ");
		targets.loadMutationsFromDb(conn);
		reportFeatureCount(targets, false, fracHits, out, out2, "m_a");
		System.out.println("m_o: ");
		targets.loadMutationsFromDb(conn);
		reportFeatureCount(oncoGenes, false, fracHitsOnc, out, out2, "m_o");
		System.out.println("m_s: ");
		targets.loadMutationsFromDb(conn);
		reportFeatureCount(tumSupGenes, false, fracHitsSup, out, out2, "m_s");
//		System.out.println("Nuc: ");
//		reportFeatureCount(nuclearGenes, false, fracHits, out);
//		System.out.println("Cyt: ");
//		reportFeatureCount(cytoplasmicGenes, false, fracHits, out);
//		System.out.println("Mem: ");
//		reportFeatureCount(membraneGenes, false, fracHits, out);
//		System.out.println("Trm: ");
//		reportFeatureCount(transMembraneGenes, false, fracHits, out);
//		targets.randomizeMutations(WITH_REPLACEMENT);
//		reportFeatureCount(targets, false, fracHits, null);
		out.close();
		out2.close();
		
		takeTime();
		//System.exit(1);
		
		// 4. ------- Energy calculation -------
		resultFile = new File(RESULT_DIR, "4_stability.txt");
		out = new PrintStream(resultFile);
		System.out.println("4. Stability");
		//targets.filterSubstructuresXRay();  // currently only X-ray structures anyways
		//targets.printStats();

		boolean useSingleChain = false;	// if true, use single structures from pdbase otherwise, use original pdb files (whole complex)
		double destabThresh = DESTAB_THRESHOLD;	// stability changes above this threshold are considered destabilizing, others neutral
		// RND
		System.out.println("RND");
		System.out.print("r_a: ");
		targets.loadMutationsFromDb(conn);
		double fracDestabAll = countFractionDestab(conn, targets, destabThresh, useSingleChain, out);
		System.out.print("r_o: ");
		targets.loadMutationsFromDb(conn);
		double fracDestabAllOnc = countFractionDestab(conn, oncoGenes, destabThresh, useSingleChain, out);
		System.out.print("r_s: ");
		targets.loadMutationsFromDb(conn);
		double fracDestabAllSup = countFractionDestab(conn, tumSupGenes, destabThresh, useSingleChain, out);
		// SNP
		System.out.println("SNP");
		System.out.print("s_a: ");
		targets.loadSNPsAsMutations(conn);
		reportStability(conn, targets, destabThresh, useSingleChain, fracDestabAll, out);
		System.out.print("s_o: ");
		targets.loadSNPsAsMutations(conn);
		reportStability(conn, oncoGenes, destabThresh, useSingleChain, fracDestabAllOnc, out);
		System.out.print("s_s: ");
		targets.loadSNPsAsMutations(conn);
		reportStability(conn, tumSupGenes, destabThresh, useSingleChain, fracDestabAllSup, out);		
		//writeStabilitiesToFile(conn, targets, new File("/project/StruPPi/henning/projects/mutanom/analysis/results/stabilities_snp.txt"));
		// MUT
		System.out.println("MUT");
		System.out.print("m_a: ");
		targets.loadMutationsFromDb(conn);
		reportStability(conn, targets, destabThresh, useSingleChain, fracDestabAll, out);
		//writeStabilitiesToFile(conn, targets, new File("/project/StruPPi/henning/projects/mutanom/analysis/results/stabilities_mut.txt"));
		System.out.print("m_o: ");
		targets.loadMutationsFromDb(conn);
		reportStability(conn, oncoGenes, destabThresh, useSingleChain, fracDestabAllOnc, out);
		//writeStabilitiesToFile(conn, oncoGenes, new File("/project/StruPPi/henning/projects/mutanom/analysis/results/stabilities_onc.txt"));
		System.out.print("m_s: ");
		targets.loadMutationsFromDb(conn);
		reportStability(conn, tumSupGenes, destabThresh, useSingleChain, fracDestabAllSup, out);
		//writeStabilitiesToFile(conn, tumSupGenes, new File("/project/StruPPi/henning/projects/mutanom/analysis/results/stabilities_sup.txt"));

//		System.out.print("Nuc: ");
//		reportStability(conn, nuclearGenes, destabThresh, useSingleChain, fracDestabAll, out);
//		System.out.print("Cyt: ");
//		reportStability(conn, cytoplasmicGenes, destabThresh, useSingleChain, fracDestabAll, out);
//		System.out.print("Mem: ");
//		reportStability(conn, membraneGenes, destabThresh, useSingleChain, fracDestabAll, out);
//		System.out.print("Trm: ");
//		reportStability(conn, transMembraneGenes, destabThresh, useSingleChain, fracDestabAll, out);		
		out.close();
		
		takeTime();
		//System.out.println("done.");
		//System.exit(1);
		
		// TODO: nothing (just to bookmark this code in eclipse)
		conn.close();
		
	}

	public static void reportAllPerGene(TargetList targets, PrintStream out, MySQLConnection conn, boolean singleChain, double exposureCutoff, double destabThresh) {
		System.out.println("Gene\t#Mut\tExp\tStab\tSites\tClust\t#Snp\tExp\tStab\tSites\tClust");
		out.println("Gene\t#Mut\tExp\tStab\tSites\tClust\t#Snp\tExp\tStab\tSites\tClust");
		for(Gene g:targets.getTargets()) {
			int nSamples;
			if(g.getGeneName().equals("TP53")) {
				nSamples = NSAMPLES_TP53;
			} else {
				nSamples = NSAMPLES;
			}
			System.out.printf("%s\t", g.getGeneName());
			out.printf("%s\t", g.getGeneName());
			TargetList singleton = new TargetList();
			singleton.addTarget(g);
			
			// MUT - Warning: Below comes the same code for SNPs (code duplication!), changes here must be also applied below!
			{
				singleton.loadMutationsFromDb(conn);
				// number of mutations
				//int numMut = g.getMutationsInStructures(true, true).size();
				int numMut = g.getMutationsInStructures(true, false).size();	// also allow predicted structures
				System.out.printf("%4d\t", numMut);
				out.printf("%4d\t", numMut);			
				
				// exposure
				double expExpected = countFractionExposed(conn, singleton, EXPOSURE_THRESHOLD, WHOLE_COMPLEX, null);
				Pair<Integer> buriedExposed = countBuriedExposed(conn, singleton, exposureCutoff, singleChain);
				int numBuried = buriedExposed.getFirst();
				int numExposed = buriedExposed.getSecond();
				double expObserved = 1.0 * numExposed / (numExposed+numBuried);
				double expValue;
				if(USE_ODD_RATIO) {
					if(expObserved == 0) {
						expValue = 0;
					} else {
						if(expExpected == 0) {
							expValue = -10;
						} else {
							expValue = expObserved / expExpected;
						}
					}		
				} else {
					expValue = expObserved;
				}
				System.out.printf("%6.4f\t", expValue);
				out.printf("%6.4f\t", expValue);
				
				// stability
				double stabExpected = countFractionDestab(conn, singleton, DESTAB_THRESHOLD, WHOLE_COMPLEX, null);
				Pair<Integer> neutralDestab = countNeutralDestab(conn, singleton, destabThresh, singleChain);
				int numNeutral = neutralDestab.getFirst();
				int numDestab = neutralDestab.getSecond();		
				double stabObserved = 1.0 * numDestab / (numDestab+numNeutral);
				double stabValue;
				if(USE_ODD_RATIO) {
					if(stabObserved == 0) {
						stabValue = 0;
					} else {
						if(stabExpected == 0) {
							stabValue = -10;
						} else {
							stabValue = stabObserved / stabExpected;
						}
					}		
				} else {
					stabValue = stabObserved;
				}
				System.out.printf("%6.4f\t", stabValue);
				out.printf("%6.4f\t", stabValue);
				
				// functional sites
				double funcObserved = reportFeatureCount(singleton, true, Double.NaN, null, null, null);
				double funcExpected = reportRandomFeatureCount(singleton, nSamples, null, null, null); // after this call, mutations will be scrambled!
				singleton.loadMutationsFromDb(conn);
				double funcValue;
				if(USE_ODD_RATIO) {
					if(funcObserved == 0) {
						funcValue = 0;
					} else {
						if(funcExpected == 0) {
							funcValue = -10;
						} else {
							funcValue = funcObserved / funcExpected;
						}
					}		
				} else {
					funcValue = funcObserved;
				}
				System.out.printf("%6.4f\t", funcValue);
				out.printf("%6.4f\t", funcValue);
				
				// clustering
				if(g.getGeneName().equals("KIT")) {
					System.out.print("");	// for breakpoint
				}
				double clustObserved = calculateClustering(singleton);
				Double[] rndPop = getClusteringRandomPopulation(singleton, nSamples, null, null);
				singleton.loadMutationsFromDb(conn);
				double clustExpected = StatUtils.mean(ArrayUtils.toPrimitive(rndPop));
				double clustValue;
				if(USE_ODD_RATIO) {
					if(clustObserved == 0) {
						clustValue = 0;
					} else {
						if(clustExpected == 0) {
							clustValue = -1;
						} else {
							clustValue = clustObserved / clustExpected;
						}
					}		
				} else {
					clustValue = clustObserved;
				}
				System.out.printf("%6.4f\t", clustValue);
				out.printf("%6.4f\t", clustValue);
				System.out.println();
			}
			
			// SNP
			{
				singleton.loadSNPsAsMutations(conn);
				// number of SNPs
				//int numMut = g.getMutationsInStructures(true, true).size();
				int numMut = g.getMutationsInStructures(true, false).size();	// also allow predicted structures
				System.out.printf("%4d\t", numMut);
				out.printf("%4d\t", numMut);			
				
				// exposure
				double expExpected = countFractionExposed(conn, singleton, EXPOSURE_THRESHOLD, WHOLE_COMPLEX, null);
				Pair<Integer> buriedExposed = countBuriedExposed(conn, singleton, exposureCutoff, singleChain);
				int numBuried = buriedExposed.getFirst();
				int numExposed = buriedExposed.getSecond();
				double expObserved = 1.0 * numExposed / (numExposed+numBuried);
				double expValue;
				if(USE_ODD_RATIO) {
					if(expObserved == 0) {
						expValue = 0;
					} else {
						if(expExpected == 0) {
							expValue = -10;
						} else {
							expValue = expObserved / expExpected;
						}
					}		
				} else {
					expValue = expObserved;
				}
				System.out.printf("%6.4f\t", expValue);
				out.printf("%6.4f\t", expValue);
				
				// stability
				double stabExpected = countFractionDestab(conn, singleton, DESTAB_THRESHOLD, WHOLE_COMPLEX, null);
				Pair<Integer> neutralDestab = countNeutralDestab(conn, singleton, destabThresh, singleChain);
				int numNeutral = neutralDestab.getFirst();
				int numDestab = neutralDestab.getSecond();		
				double stabObserved = 1.0 * numDestab / (numDestab+numNeutral);
				double stabValue;
				if(USE_ODD_RATIO) {
					if(stabObserved == 0) {
						stabValue = 0;
					} else {
						if(stabExpected == 0) {
							stabValue = -10;
						} else {
							stabValue = stabObserved / stabExpected;
						}
					}		
				} else {
					stabValue = stabObserved;
				}
				System.out.printf("%6.4f\t", stabValue);
				out.printf("%6.4f\t", stabValue);
				
				// functional sites
				double funcObserved = reportFeatureCount(singleton, true, Double.NaN, null, null, null);
				double funcExpected = reportRandomFeatureCount(singleton, NSAMPLES, null, null, null); // after this call, mutations will be scrambled!
				singleton.loadSNPsAsMutations(conn);
				double funcValue;
				if(USE_ODD_RATIO) {
					if(funcObserved == 0) {
						funcValue = 0;
					} else {
						if(funcExpected == 0) {
							funcValue = -10;
						} else {
							funcValue = funcObserved / funcExpected;
						}
					}		
				} else {
					funcValue = funcObserved;
				}
				System.out.printf("%6.4f\t", funcValue);
				out.printf("%6.4f\t", funcValue);
				
				// clustering
				if(g.getGeneName().equals("KIT")) {
					System.out.print("");	// for breakpoint
				}
				double clustObserved = calculateClustering(singleton);
				Double[] rndPop = getClusteringRandomPopulation(singleton, NSAMPLES, null, null);
				singleton.loadSNPsAsMutations(conn);
				double clustExpected = StatUtils.mean(ArrayUtils.toPrimitive(rndPop));
				double clustValue;
				if(USE_ODD_RATIO) {
					if(clustObserved == 0) {
						clustValue = 0;
					} else {
						if(clustExpected == 0) {
							clustValue = -1;
						} else {
							clustValue = clustObserved / clustExpected;
						}
					}		
				} else {
					clustValue = clustObserved;
				}
				// in case that there are not enough SNPs, value can be NaN
				if(Double.isNaN(clustValue)) {
					clustValue=1.0;
				}
				System.out.printf("%6.4f\t", clustValue);
				out.printf("%6.4f\t", clustValue);
				System.out.println();			
			}
			
			// END OF THIS GENE
			out.println();
		}
	}
	
	public static void reportClusteringPerGene(TargetList targets) {
		int numClust = 0;
		double sumNumDist = 0.0;
		double sumClust = 0.0;
		for(Gene g:targets.getTargets()) {
			TargetList singleton = new TargetList();
			singleton.addTarget(g);
			Pair<Double> clustNumDist = calculateClustering2(singleton);
			double clust = clustNumDist.getFirst();
			double numDist =clustNumDist.getSecond();
			numClust++;
			sumClust += clust;
			sumNumDist += numDist;
			System.out.printf("%s\tClustering: %4.2f\tNumDist: %.0f\n", g.getGeneName(), clust, numDist);
		}
		double meanClust = 1.0 * sumClust / numClust;
		double meanNumDist = 1.0 * sumNumDist / numClust;
		System.out.printf("%s\tClustering: %4.2f\tNumDist: %.0f\n", "Mean", meanClust, meanNumDist);
	}
	
	/**
	 * Calculates a population of target lists with randomized mutations and returns for the clustering
	 * values for all individuals as an array. Writes the average clustering to stdout and PrintStream out
	 * and the population to PrintStream out2.
	 * Warning: scrambles mutations!
	 * @param targets
	 * @param out
	 */
	public static Double[] getClusteringRandomPopulation(TargetList targets, int numSamples, PrintStream out, PrintStream out2) {
		double sumClust = 0;
		double clust;
		LinkedList<Double> clustAll = new LinkedList<Double>();
		for (int i = 0; i < numSamples; i++) {
			targets.randomizeMutations(CLUST_WITH_REPLACEMENT); // with replacement
			clust = calculateClustering(targets);
			sumClust += clust;
			if(out2 != null) out2.println(clust);
			clustAll.add(clust);
		}
		double mean = 1.0 * sumClust / numSamples;
		//System.out.printf("Clustering:%5.2f\n", mean);
		if(out != null) out.printf("%f\t%f\t%f\t%e\n", mean, mean, mean, 0.0); // fraction, jackknife, p-value
		Double[] clustArr = new Double[clustAll.size()];
		return clustAll.toArray(clustArr);
	}
	
	/**
	 * Calculates clustering and returns the clustering value. This method is provided for compatibility.
	 * Otherwise, use calculateClustering2.
	 * @param targets
	 * @return
	 */
	public static double calculateClustering(TargetList targets) {
		Pair<Double> clustNumDist = calculateClustering2(targets);
		double clust = clustNumDist.getFirst();
		return clust;
	}

	/**
	 * Calculates clustering and returns the clustering value and the number of distances evaluated.
	 * Distances are measured only within a structural domain. 
	 * @param targets
	 * @return
	 */
	public static Pair<Double> calculateClustering2(TargetList targets) {		
		double sumDistList = 0;
		int numDistList = 0;
		for(Gene g:targets.getTargets()) {
			//System.out.print(g.geneName + " | ");
			double sumDistGene = 0;
			int numDistGene = 0;
			if(g.getFeatures().size()==0) System.err.println("Warning. No features found when calculating clustering for " + g.getGeneName());
			for(Feature f:g.getFeatures()) {
				if(f.getType() == FeatureType.SDOMAIN) {
					double sumDist = 0;
					int numDist = 0;
					TreeSet<Integer> intSet = f.getIntervalSet().getIntegerSet();	// set of positions in current domain
					for(Mutation m:g.getMutations()) {
						int iPos = m.getPos();
						for(Mutation m2:g.getMutations()) {
							int jPos = m2.getPos();
							if(iPos < jPos && intSet.contains(iPos) && intSet.contains(jPos)) {
								if(g.getSubstructure(iPos) != g.getSubstructure(jPos)) System.err.printf("Assertion failed: ss(%d) != ss(%d) in %s", iPos, jPos, f.toString());
								Substructure ss = g.getSubstructure(iPos);
								int iPdbPos = ss.mapUniprotResser2Cif(iPos);
								int jPdbPos = ss.mapUniprotResser2Cif(jPos);
								if(ss.getPdb().hasCoordinates(iPdbPos) && ss.getPdb().hasCoordinates(jPdbPos)) {
									if(CLUSTERING_ATOM.equalsIgnoreCase("Centroid")) {									
										Residue iRes = ss.getPdb().getResidue(iPdbPos);
										Residue jRes = ss.getPdb().getResidue(jPdbPos);
										Point3d px = iRes.getScCentroid();
										Point3d py = jRes.getScCentroid();
										if(px != null && py != null) {
											double dist = 1.0 / px.distance(py);
											sumDist += dist;
											numDist ++;										
										} else {
											//System.out.print("No coordinates for " + iPdbPos + " and/or " + jPdbPos);
										}
									} else {
										if(ss.getPdb().hasCoordinates(iPdbPos, CLUSTERING_ATOM) && ss.getPdb().hasCoordinates(jPdbPos, CLUSTERING_ATOM)) {
												Point3d px = ss.getPdb().getAtomCoord(iPdbPos, CLUSTERING_ATOM);
												Point3d py = ss.getPdb().getAtomCoord(jPdbPos, CLUSTERING_ATOM);
												double dist = 1.0 / px.distance(py);
												sumDist += dist;
												numDist ++;
										} else {
											//System.out.print("No " + CLUSTERING_ATOM + " coordinates for " + iPdbPos + " and/or " + jPdbPos);
										}
									}
								} else {
									//System.out.print("No coordinates for " + iPdbPos + " and/or " + jPdbPos);
								}
								
							}
						}
					}
					if(numDist > 0) {
						//System.out.printf("%s : %6.2f (%d) | ", ss.getRange(), 100.0 * sumDist/numDist, numDist);
						sumDistGene += sumDist;
						numDistGene += numDist;
					}
				} else {
					//System.out.print(ss.getRange() + " | ");					
				}
			}

			
// This is the old way of calculating clustering within whole substructures. After testing the new way, this can be deleted:			
//			for(Substructure ss: g.getSubstructuresWithMutations()) {
//				double sumDist = 0;
//				int numDist = 0;
//				for(Mutation m:g.getMutations()) {
//					int iPos = m.position;
//					for(Mutation m2:g.getMutations()) {
//						int jPos = m2.position;
//						if(iPos < jPos && g.getSubstructure(iPos) == ss && g.getSubstructure(jPos) == ss) {
//							int iPdbPos = ss.mapUniprotResser2Pdb(iPos);
//							int jPdbPos = ss.mapUniprotResser2Pdb(jPos);
//							if(ss.pdb.hasCoordinates(iPdbPos) && ss.pdb.hasCoordinates(jPdbPos)) {
//								if(CLUSTERING_ATOM.equalsIgnoreCase("Centroid")) {									
//									Residue iRes = ss.pdb.getResidue(iPdbPos);
//									Residue jRes = ss.pdb.getResidue(jPdbPos);
//									Point3d px = iRes.getScCentroid();
//									Point3d py = jRes.getScCentroid();
//									if(px != null && py != null) {
//										double dist = 1.0 / px.distance(py);
//										sumDist += dist;
//										numDist ++;										
//									} else {
//										//System.out.print("No coordinates for " + iPdbPos + " and/or " + jPdbPos);
//									}
//								} else {
//									if(ss.pdb.hasCoordinates(iPdbPos, CLUSTERING_ATOM) && ss.pdb.hasCoordinates(jPdbPos, CLUSTERING_ATOM)) {
//											Point3d px = ss.pdb.getAtomCoord(iPdbPos, CLUSTERING_ATOM);
//											Point3d py = ss.pdb.getAtomCoord(jPdbPos, CLUSTERING_ATOM);
//											double dist = 1.0 / px.distance(py);
//											sumDist += dist;
//											numDist ++;
//									} else {
//										//System.out.print("No " + CLUSTERING_ATOM + " coordinates for " + iPdbPos + " and/or " + jPdbPos);
//									}
//								}
//							} else {
//								//System.out.print("No coordinates for " + iPdbPos + " and/or " + jPdbPos);
//							}
//						}
//					}
//				}
//				if(numDist > 0) {
//					//System.out.printf("%s : %6.2f (%d) | ", ss.getRange(), 100.0 * sumDist/numDist, numDist);
//					sumDistGene += sumDist;
//					numDistGene += numDist;
//				} else {
//					//System.out.print(ss.getRange() + " | ");
//				}
//			}
						
			//System.out.println();
			//System.out.printf("%6.2f (%d)\n", 100.0 * sumDistGene / numDistGene, numDistGene);
			sumDistList += sumDistGene;
			numDistList += numDistGene;
		}
		double clust = 100.0 * sumDistList / numDistList;
		Pair<Double> clustNumDist = new Pair<Double>(clust, (double) numDistList);
		return clustNumDist;
	}
	
	public static double reportClustering(TargetList targets, Double[] randPop, PrintStream out) {
		// overall clustering
		double clust = calculateClustering(targets);
		
		// jackknife test
		double maxClust = 0;
		double minClust = 100;
		for(Gene g:targets.getTargets()) {
			// exclude g from list
			TargetList oneLess = targets.getOneLess(g);
			double newClust = calculateClustering(oneLess); 
			maxClust = Math.max(maxClust, newClust);
			minClust = Math.min(minClust, newClust);
		}
		
		// estimated p-value, old version with individual population per set
		//Pair<Double> pair = getClusteringPValue(targets, clust);	// Warning: Scrambles mutations!
		//double mean = pair.getFirst();	// mean over 10000 simulated permutations
		//double pVal = pair.getSecond();
		//System.out.printf("Clustering: %4.2f [%4.2f;%4.2f] Rand: %4.2f P-value: %e\n", clust, minClust, maxClust, mean, pVal);
		
		// estimated p-value, new version based on single population derived from set of cancer mutations
		int numSamples = randPop.length;
		double sumClust = 0;
		int numPositive = 0;
		int numNegative = 0;
		for(double val:randPop) {
			sumClust += val;
			if(val >= clust) numPositive++;
			if(val <= clust) numNegative++;
		}
		double mean = 1.0 * sumClust / numSamples;
		double pVal;
		if(clust >= mean) {
			pVal = (double) numPositive / (double) numSamples;
		} else {
			pVal = -1.0 * (double) numNegative / (double) numSamples;
		}
		System.out.printf("Clustering: %4.2f [%4.2f;%4.2f] P-value: %e\n", clust, minClust, maxClust, pVal);
		out.printf("%f\t%f\t%f\t%e\n", clust, minClust, maxClust, pVal);
		return clust;
	}
	
//	/**
//	 * Estimates a P-value by generating a population of randomly permuted mutations
//	 * and counting how often the given observed value is met or exceeded in the population.
//	 * Side effect: Mutations will be randomly permuted in the given target list!
//	 * @param targets
//	 * @return the mean and the estimated p-value
//	 */
//	public static Pair<Double> getClusteringPValue(TargetList targets, double observedClusteringValue) {
//		int numSamples = 10000;
//		double sumClust = 0;
//		int numEqualOrMore = 0;
//		for (int i = 0; i < numSamples; i++) {
//			targets.randomizeMutations(WITH_REPLACEMENT); // with replacement	
//			double clust = calculateClustering(targets);
//			sumClust += clust;
//			if(clust >= observedClusteringValue) numEqualOrMore++;
//		}
//		double mean = 1.0 * sumClust / numSamples;
//		double pVal = 1.0 * numEqualOrMore / numSamples;
//		return new Pair<Double>(mean,pVal);
//	}
	
//	public static void calculateExposure(TargetList targets) {
//		// calculate exposure
//		int numExposed = 0;
//		int numBuried = 0;
//		for(Gene g:targets.getTargets()) {
//			for(Mutation m:g.getMutations()) {
//				int pos = m.position;
//				Substructure ss = g.getSubstructure(pos);
//				int pdbPos = ss.mapUniprotResser2Pdb(pos);
//				double rsa = ss.getSurfaceAccessibility(pdbPos);
//				if(rsa >= Substructure.EXPOSURE_CUTOFF) numExposed++;
//				if(rsa < Substructure.EXPOSURE_CUTOFF) numBuried++;
//			}
//		}
//		System.out.println("Number of exposed missense mutations: " + numExposed);
//		System.out.println("Number of buried missense mutations: " + numBuried);
//		System.out.printf("Percent exposed: %4.2f\n", 1.0 * numExposed / (numExposed+numBuried));
//		
//		countFractionExposed(targets);
//	}

	/**
	 * Counts the number of buried and exposed mutations in the given target list.
	 * Relies on correct results in the database.
	 */
	public static Pair<Integer> countBuriedExposed(MySQLConnection conn, TargetList targets, double exposureCutoff, boolean singleChain) {
		// calculate exposure
		int numExposed = 0;
		int numBuried = 0;
		for(Gene g:targets.getTargets()) {
			for(Mutation m:g.getMutations()) {
				int pos = m.getPos();
				Substructure ss = g.getSubstructure(pos);
				//if(ss != null && ss.getType() != SubstructureType.PREDICTION) {
				if(ss != null) {
					int cifPos = ss.mapUniprotResser2Cif(pos);
					double rsa = getExposureFromDb(conn, ss.getPdbCode(), ss.getChainCode(), cifPos, singleChain);
					if(rsa > exposureCutoff) numExposed++;
					if(rsa <= exposureCutoff) numBuried++;
				}
			}
		}
		return new Pair<Integer>(numBuried,numExposed);
	}
	
	/**
	 * Counts the number of neutral and destabilizing mutations in the given target list.
	 * Precalculated DeltaDeltaG values are read from the given database connection.
	 */
	public static Pair<Integer> countNeutralDestab(MySQLConnection conn, TargetList targets, double destabThresh, boolean singleChain) {
		// calculate exposure
		int numDestab = 0;
		int numNeutral = 0;
		for(Gene g:targets.getTargets()) {
			for(Mutation m:g.getMutations()) {
				int pos = m.getPos();
				Substructure ss = g.getSubstructure(pos);
				//if(ss != null && ss.getType() != SubstructureType.PREDICTION) {
				if(ss != null) {
					int cifPos = ss.mapUniprotResser2Cif(pos);
					AminoAcid mutAa = m.getMutAA();
					double ddg = getDeltaDeltaGFromDb(conn, ss.getPdbCode(), ss.getChainCode(), cifPos, mutAa, singleChain);
					if(ddg > destabThresh) numDestab++;
					if(ddg <= destabThresh) numNeutral++;
				}
			}
		}
		return new Pair<Integer>(numNeutral,numDestab);
	}
	
	/**
	 * Writes the DeltaDeltaG values for the currently loaded mutations to a file.
	 * Precalculated DeltaDeltaG values are read from the given database connection.
	 */
	public static void writeStabilitiesToFile(MySQLConnection conn, TargetList targets, File outFile) {
		try {
			PrintWriter out = new PrintWriter(new FileWriter(outFile));
			boolean singleChain = false;
			for(Gene g:targets.getTargets()) {
				for(Mutation m:g.getMutations()) {
					int pos = m.getPos();
					Substructure ss = g.getSubstructure(pos);
					if(ss != null && ss.getType() != SubstructureType.PREDICTION) {
						int cifPos = ss.mapUniprotResser2Cif(pos);
						AminoAcid mutAa = m.getMutAA();
						double ddg = getDeltaDeltaGFromDb(conn, ss.getPdbCode(), ss.getChainCode(), cifPos, mutAa, singleChain);
						out.println(ddg);
					}
				}
			}
			out.close();
		} catch(IOException e) {
			System.err.println("Error writing to file " + outFile);
		}
	}
	
	/**
	 * Reports the fraction of exposed mutations, the bootstrap interval and the p-value of the observation for the given target list.
	 * @param conn
	 * @param targets
	 * @param exposureCutoff
	 * @param singleChain
	 * @param fracAll
	 * @return
	 */
	public static double reportExposure(MySQLConnection conn, TargetList targets, double exposureCutoff, boolean singleChain, double fracAll, PrintStream out) {
		// original data
		Pair<Integer> buriedExposed = countBuriedExposed(conn, targets, exposureCutoff, singleChain);
		int numBuried = buriedExposed.getFirst();
		int numExposed = buriedExposed.getSecond();		
		double fractionExposed = 1.0 * numExposed / (numExposed+numBuried);
		double p = getPValue(fracAll,numExposed,numExposed+numBuried);		
		// jackknife test
		double maxFracExposed = 0;
		double minFracExposed = 100;
		for(Gene g:targets.getTargets()) {
			// exclude g from list
			TargetList oneLess = targets.getOneLess(g);
			buriedExposed = countBuriedExposed(conn, oneLess, exposureCutoff, singleChain);
			int newNumBuried = buriedExposed.getFirst();
			int newNumExposed = buriedExposed.getSecond();				
			double newFracExposed = 1.0 * newNumExposed / (newNumExposed+newNumBuried);
			maxFracExposed = Math.max(maxFracExposed, newFracExposed);
			minFracExposed = Math.min(minFracExposed, newFracExposed);
		}		
		System.out.printf("Exposed: %4d Buried: %4d Fraction exposed: %4.2f [%4.2f;%4.2f] P-value: %e\n", numExposed, numBuried, fractionExposed, minFracExposed, maxFracExposed, p);
		out.printf("%f\t%f\t%f\t%e\n", fractionExposed, minFracExposed, maxFracExposed, p); // actual value, error margins, p-value
		return fractionExposed;
	}

	/**
	 * Reports the fraction of destabilizing mutations, the bootstrap interval and the p-value of the observation for the given target list.
	 * @param conn an active connection to the database holding the calculated stability changes
	 * @param targets the list of target genes (with previously loaded substructures and mutations)
	 * @param destabThresh the threshold above which a mutation is considered destabilizing (in terms of DeltaDeltaG)
	 * @param singleChain if true, uses the values calculated for isolated chains, otherwise the ones for pdb complexes
	 * @param fracAll the background probability of observing a destabilizing mutation (used for p-value calculation)
	 * @return the fraction of destabilizing mutations in the set of loaded mutations
	 */
	public static double reportStability(MySQLConnection conn, TargetList targets, double destabThresh, boolean singleChain, double fracAll, PrintStream out) {
		// original data
		Pair<Integer> neutralDestab = countNeutralDestab(conn, targets, destabThresh, singleChain);
		int numNeutral = neutralDestab.getFirst();
		int numDestab = neutralDestab.getSecond();		
		double fractionDestab = 1.0 * numDestab / (numDestab+numNeutral);
		double p = getPValue(fracAll,numDestab,numDestab+numNeutral);		
		// jackknife test
		double maxFracDestab = 0;
		double minFracDestab = 100;
		for(Gene g:targets.getTargets()) {
			// exclude g from list
			TargetList oneLess = targets.getOneLess(g);
			neutralDestab = countNeutralDestab(conn, oneLess, destabThresh, singleChain);
			int newNumNeutral = neutralDestab.getFirst();
			int newNumDestab = neutralDestab.getSecond();				
			double newFracDestab = 1.0 * newNumDestab / (newNumDestab+newNumNeutral);
			maxFracDestab = Math.max(maxFracDestab, newFracDestab);
			minFracDestab = Math.min(minFracDestab, newFracDestab);
		}		
		System.out.printf("Destabilizing: %4d Neutral: %5d Fraction destab: %4.2f [%4.2f;%4.2f] P-value: %e\n", numDestab, numNeutral, fractionDestab, minFracDestab, maxFracDestab, p);
		out.printf("%f\t%f\t%f\t%e\n", fractionDestab, minFracDestab, maxFracDestab, p);
		return fractionDestab;
	}
	
	/**
	 * Gets the naccess result for a single residue from the database.
	 * @param conn
	 * @param pdbCode
	 * @param pdbChain
	 * @param cifPos
	 * @param singleChain this parameter is currently ignored!
	 * @return
	 */
	public static double getExposureFromDb(MySQLConnection conn, String pdbCode, String pdbChain, int cifPos, boolean singleChain) {
		int resIdx = conn.getIntFromDb(String.format(QUERY_RES_IDX, pdbCode, pdbChain, cifPos)); 
		double rsa = conn.getDoubleFromDb(String.format(QUERY_NACCESS_RESULT, resIdx));		
		if(rsa == Double.NaN) {
			System.err.println("Warning: No exposure value found for " + pdbCode + pdbChain + cifPos);
		}
		return rsa;
	}
	
	/**
	 * Gets the foldx mutant stability result for a single mutation from the database.
	 * @param conn
	 * @param pdbCode
	 * @param pdbChain
	 * @param cifPos
	 * @param mutAa
	 * @param singleChain this parameter is currently ignored!
	 * @return
	 */
	public static double getDeltaDeltaGFromDb(MySQLConnection conn, String pdbCode, String pdbChain, int cifPos, AminoAcid mutAa, boolean singleChain) {
		int resIdx = conn.getIntFromDb(String.format(QUERY_RES_IDX, pdbCode, pdbChain, cifPos)); 
		double ddg = conn.getDoubleFromDb(String.format(QUERY_FOLDX_RESULT, resIdx, mutAa.getOneLetterCode()));	
		if(ddg == Double.NaN) {
			System.err.println("Warning: No DeltaDeltaG value found for " + pdbCode + pdbChain + " " + mutAa.getOneLetterCode() + cifPos);
		}
		return ddg;
	}
	
//	/**
//	 * Count the fraction of all residues in the set of substructures which are exposed (rsa above cutoff)
//	 * @param pdbs
//	 */
//	public static double countFractionExposed(TargetList targets) {
//		int numExposed = 0;
//		int numBuried = 0;
//		for(Gene g:targets.getTargets())
//			for(Substructure ss:g.getSubstructuresWithMutations()) {
//				for(int pos:ss.pdb.getAllSortedResSerials()) {
//					double rsa = ss.getSurfaceAccessibility(pos);
//					if(rsa > Substructure.EXPOSURE_CUTOFF) numExposed++;
//					if(rsa <= Substructure.EXPOSURE_CUTOFF) numBuried++;
//				}			
//			}
//		double fractionExposed = 1.0 * numExposed / (numExposed+numBuried);
//		System.out.printf("Exposed: %d Buried: %d Fraction exposed: %4.2f\n", numExposed, numBuried, fractionExposed);
//		return fractionExposed;
//	}
	
	/**
	 * Count the fraction of all residues in the set of substructures which are exposed (rsa above cutoff)
	 * where exposure values are read from the database.
	 * @param pdbs
	 */
	public static double countFractionExposed(MySQLConnection conn, TargetList targets, double exposureCutoff, boolean singleChain, PrintStream out) {
		int numExposed = 0;
		int numBuried = 0;
		for(Gene g:targets.getTargets())
			for(Substructure ss:g.getSubstructuresWithMutations()) {
			//for(Substructure ss:g.getSubstructures()) {
				for(int pos:ss.getPdb().getAllSortedResSerials()) {
					double rsa = getExposureFromDb(conn, ss.getPdbCode(), ss.getChainCode(), pos, singleChain);
					if(rsa > exposureCutoff) numExposed++;
					if(rsa <= exposureCutoff) numBuried++;
				}			
			}		double fractionExposed = 1.0 * numExposed / (numExposed+numBuried);
		//double var = 1.0 * fractionExposed * (1-fractionExposed) * (numExposed+numBuried);
		System.out.printf("Exposed: %d Buried: %d Fraction exposed: %4.2f\n", numExposed, numBuried, fractionExposed);
		if(out != null) out.printf("%f\t%f\t%f\t%e\n",fractionExposed,fractionExposed,fractionExposed, 0.0); // actual value, error margins, p-value
		return fractionExposed;
	}
	
	/**
	 * Count the fraction of all mutations which are destabilizing (delta delta G above threshold)
	 * where ddg values are read from the database.
	 * @param pdbs
	 */
	public static double countFractionDestab_old(MySQLConnection conn, TargetList targets, double destabThresh, boolean singleChain, PrintStream out) {
		int numDestab = 0;
		int numNeutral = 0;
		for(Gene g:targets.getTargets())
			for(Substructure ss:g.getSubstructuresWithMutations()) {
			//for(Substructure ss:g.getSubstructures()) {
				for(int pos:ss.getPdb().getAllSortedResSerials()) {
					AminoAcid wtAa = ss.getPdb().getResidue(pos).getAaType();
					for(AminoAcid mutAa:AminoAcid.values()) {
						if(mutAa != wtAa) {
							double ddg = getDeltaDeltaGFromDb(conn, ss.getPdbCode(), ss.getChainCode(), pos, mutAa, singleChain);
							if(ddg > destabThresh) numDestab++;
							if(ddg <= destabThresh) numNeutral++;
						}
					}
				}			
			}
		double fractionDestab = 1.0 * numDestab / (numDestab+numNeutral);
		System.out.printf("Destabilizing: %d Neutral: %d Fraction destab: %4.2f\n", numDestab, numNeutral, fractionDestab);
		if(out != null) out.printf("%f\t%f\t%f\t%e\n",fractionDestab,fractionDestab,fractionDestab,0.0); // value, jackknife, p-value
		return fractionDestab;
	}
	
	
	
	/**
	 * Count the fraction of all mutations which are destabilizing (delta delta G above threshold)
	 * where ddg values are read from the database.
	 * @param pdbs
	 */
	public static double countFractionDestab(MySQLConnection conn, TargetList targets, double destabThresh, boolean singleChain, PrintStream out) {
		int numDestab = 0;
		int numNeutral = 0;
		for(Gene g:targets.getTargets())
			for(Substructure ss:g.getSubstructuresWithMutations()) {
			//for(Substructure ss:g.getSubstructures()) {
				int d = conn.getIntFromDb(String.format(QUERY_COUNT_DESTAB, ss.getPdbCode(), ss.getChainCode(), destabThresh));
				int n = conn.getIntFromDb(String.format(QUERY_COUNT_NEUTRAL, ss.getPdbCode(), ss.getChainCode(), destabThresh));
				numDestab += d;
				numNeutral += n;
				if(n+d == 0) System.err.println("Warning: No FOLDX results found for " + ss.getPdbCode()+ss.getChainCode());
			}
		double fractionDestab = 1.0 * numDestab / (numDestab+numNeutral);
		System.out.printf("Destabilizing: %d Neutral: %d Fraction destab: %4.2f\n", numDestab, numNeutral, fractionDestab);
		if(out != null) out.printf("%f\t%f\t%f\t%e\n",fractionDestab,fractionDestab,fractionDestab,0.0); // value, jackknife, p-value
		return fractionDestab;
	}
	
	/**
	 * Helper function for reportFeatureCount which does the actual counting and wraps up results in an array
	 * @param targets the target list for which features are to be counted
	 * @return
	 */
	private static int[] getFeatureCount(TargetList targets) {
		
		int actHits = 0, actProx = 0;	// active sites
		int phoHits = 0, phoProx = 0;	// phosphorilation sites
		int ubqHits = 0, ubqProx = 0;	// ubiquitination sites
		int modHits = 0, modProx = 0;	// other modification sites (except phospho + ubiq + glyco)
		int atpHits = 0, atpProx = 0;	// atp binding sites
		int gtpHits = 0, gtpProx = 0;	// gtp binding sites
		int glycHits = 0, glycProx = 0;	// glycosylation sites
		int dnaHits = 0, dnaProx = 0;	// dna binding sites
		int numMut = 0;
				
		//System.out.println("Checking proximity to functional sites");
		for(Gene g:targets.getTargets()) {
			for(Mutation m:g.getMutations()) {
				Substructure ss = g.getSubstructure(m.getPos());
				if(ss != null && ss.getPdb().hasCoordinates(ss.mapUniprotResser2Cif(m.getPos()))) {
//					if(g.getGeneName().equals("PTEN") && m.position == 45) {
//						System.out.println("!");
//					}
					numMut++;
					//System.out.print(".");
					// Remember observed features to prevent overcounting of the same feature
					HashSet<Feature> alreadyObserved = new HashSet<Feature>();
					String fsHit = null;	// flag whether a site has been found already
					String fsNear = null;	// flag whether a proximal site has been found already
					// check whether this position matches a feature
					fsHit = g.getFirstFeatureAndFunctionalClass(m.getPos(), alreadyObserved).getSecond();
					if(fsHit.equals("ACT_SITE")) actHits++;
					if(fsHit.equals("PHO_RES"))  phoHits++;
					if(fsHit.equals("UBQ_RES"))  ubqHits++;
					if(fsHit.equals("MOD_RES"))  modHits++;
					if(fsHit.equals("ATP_BIND")) atpHits++;
					if(fsHit.equals("GTP_BIND")) gtpHits++;
					if(fsHit.equals("GLYC_RES")) glycHits++;
					if(fsHit.equals("DNA_BIND")) dnaHits++;					
//					for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.UNIPROT, m.position)) {
//						alreadyObserved.add(f);
//						UniprotFeature uf = (UniprotFeature) f;
//						if(uf.getUniprotTypeName().equals("ACT_SITE")) {actHits++; fsHit="ACT_SITE"; break;}
//						if(uf.getUniprotTypeName().equals("NP_BIND")) {
//							if(uf.getDescription().startsWith("ATP")) {atpHits++; fsHit="NP_BIND"; break;}
//							else if(uf.getDescription().startsWith("GTP")) {gtpHits++; fsHit="NP_BIND"; break;}
//							else System.err.println("Unknown NP_BIND feature: " + uf.getDescription()); break;	
//						}
//						if(uf.getUniprotTypeName().equals("MOD_RES")) {
//							if(uf.getDescription().startsWith("Phospho")) phoHits++; else
//							if(uf.getDescription().indexOf("nitroso") >= 0) System.err.println(uf);
//							else modHits++; 
//							fsHit="MOD_RES"; break;
//						}
//						if(uf.getUniprotTypeName().equals("CARBOHYD")) {glycHits++;	fsHit="CARBOHYD"; break;}
//						if(uf.getUniprotTypeName().equals("DNA_BIND")) {dnaHits++;	fsHit="DNA_BIND"; break;}
//					}
//					if(fsHit == null) {
//						for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.CSA, m.position)) {
//							alreadyObserved.add(f);
//							actHits++;
//							fsHit="ACT_SITE";
//							break;
//						}					
//					}
//					if(fsHit == null) {
//						for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.GENERAL, m.position)) {
//							alreadyObserved.add(f);
//							GeneralFeature gf = (GeneralFeature) f;
//							if(gf.getDescription().equals("ACT_SITE")) {actHits++; fsHit="ACT_SITE"; break;}
//							if(gf.getDescription().equals("ATP_BIND")) {atpHits++; fsHit="ATP_BIND"; break;}
//							if(gf.getDescription().equals("GTP_BIND")) {gtpHits++; fsHit="GTP_BIND"; break;}
//							if(gf.getDescription().equals("ATY_RES")) {modHits++; fsHit="MOD_RES"; break;}
//							if(gf.getDescription().equals("MOD_RES")) {modHits++; fsHit="MOD_RES"; break;}
//							if(gf.getDescription().equals("PHO_RES")) {phoHits++; fsHit="PHO_RES"; break;}
//							if(gf.getDescription().equals("UBQ_RES")) {ubqHits++; fsHit="UBQ_RES"; break;}
//							if(gf.getDescription().equals("CARBOHYD")) {glycHits++;	fsHit="CARBOHYD"; break;}
//							if(gf.getDescription().equals("DNA_BIND")) {dnaHits++;	fsHit="DNA_BIND"; break;}
//						}
//					}
//					if(fsHit == null) {
//						for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.PHOSPHOSITE, m.position)) {
//							alreadyObserved.add(f);
//							fsHit="MOD_RES";
//							if(((PhosphoSitePlusFeature)f).getModType() == ProteinModificationType.PHOSPHORYLATION) phoHits++;
//							else if(((PhosphoSitePlusFeature)f).getModType() == ProteinModificationType.UBIQUITINATION) ubqHits++;
//							else modHits++;
//							break;
//						}
//					}
					// iterate over all neighbouring positions, check for motif
					int pdbMutPos = ss.mapUniprotResser2Cif(m.getPos());
					if(ss.getGraph() == null) System.err.println("Error: graph not loaded for " + g.getGeneName() + " " + ss.getRange());
					RIGNode mutNode = ss.getGraph().getNodeFromSerial(pdbMutPos);
					if(mutNode == null) {
						// Skip output because this happens too often, this simply means that C-beta (or other contact type) was not found
						//System.err.println("Error: " + ss.getPdbCode()+ss.getChainCode() + " " + pdbMutPos + " has coordinates but graph node is null.");
					} else {
						//if(m.before != AminoAcid.getByThreeLetterCode(mutNode.getResidueType())) 
						//System.err.println("Residue mismatch: " + g.geneName + " " + m.position + " " + m.before.getThreeLetterCode() 
						//		+ " != " + ss.pdbCode+ss.chainCode + " " + pdbMutPos + " " + mutNode.getResidueType());
						if(fsHit == null || fsHit.equals("NONE")) {	// only check proximal hit if no direct hit found
							int pdbNbPos = -1;
							int uniNbPos = -1;
							for(RIGNode n:ss.getGraph().getOrderedNeighbors(mutNode)) {
								// the fact that this call return the neighbours in a nondeterministic order(!) caused a lot of trouble. So now we are ordering them by residueSerial.
								pdbNbPos = n.getResidueSerial();
								uniNbPos = ss.mapCifResser2Uniprot(pdbNbPos);
								if(fsNear == null || fsNear.equals("NONE")) { // only check if no prox hit found yet
									// check whether this position matches a feature
									fsNear = g.getFirstFeatureAndFunctionalClass(uniNbPos, alreadyObserved).getSecond();
									if(fsNear.equals("ACT_SITE")) actProx++;
									if(fsNear.equals("PHO_RES"))  phoProx++;
									if(fsNear.equals("UBQ_RES"))  ubqProx++;
									if(fsNear.equals("MOD_RES"))  modProx++;
									if(fsNear.equals("ATP_BIND")) atpProx++;
									if(fsNear.equals("GTP_BIND")) gtpProx++;
									if(fsNear.equals("GLYC_RES")) glycProx++;
									if(fsNear.equals("DNA_BIND")) dnaProx++;		
									
//									for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.UNIPROT, uniNbPos)) {
//										if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
//											alreadyObserved.add(f);
//											UniprotFeature uf = (UniprotFeature) f;
//											if(uf.getUniprotTypeName().equals("ACT_SITE")) {
//												actProx++; fsNear="ACT_SITE"; break;}
//											if(uf.getUniprotTypeName().equals("NP_BIND")) {
//												if(uf.getDescription().startsWith("ATP")) {atpProx++; fsNear="NP_BIND"; break;}
//												else if(uf.getDescription().startsWith("GTP")) {gtpProx++; fsNear="NP_BIND"; break;}
//												else System.err.println("Unknown NP_BIND feature: " + uf.getDescription()); break;
//											}
//											if(uf.getUniprotTypeName().equals("MOD_RES")) {
//												if(uf.getDescription().startsWith("Phospho")) {phoProx++;}
//												else modProx++; 
//												fsNear="MOD_RES"; break;
//											}
//											if(uf.getUniprotTypeName().equals("CARBOHYD")) {glycProx++; fsNear="CARBOHYD"; break;}
//											if(uf.getUniprotTypeName().equals("DNA_BIND")) {dnaProx++; fsNear="DNA_BIND"; break;}
//										}
//									}
//								}
//								if(fsNear == null) {
//									for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.CSA, uniNbPos)) {
//										if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
//											alreadyObserved.add(f);
//											actProx++;
//											//System.out.print("ac");
//											fsNear="ACT_SITE"; 
//											break;
//										}
//									}
//								}
//								if(fsNear == null) {
//									for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.GENERAL, uniNbPos)) {
//										if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
//											alreadyObserved.add(f);
//											GeneralFeature gf = (GeneralFeature) f;
//											if(gf.getDescription().equals("ACT_SITE")) {actProx++; fsNear="ACT_SITE";break;}
//											if(gf.getDescription().equals("ATP_BIND")) {atpProx++; fsNear="ATP_BIND"; break;}
//											if(gf.getDescription().equals("GTP_BIND")) {gtpProx++; fsNear="GTP_BIND"; break;}
//											if(gf.getDescription().equals("ATY_RES")) {modProx++; fsNear="MOD_RES"; break; }
//											if(gf.getDescription().equals("MOD_RES")) {modProx++; fsNear="MOD_RES"; break; }
//											if(gf.getDescription().equals("UBQ_RES")) {ubqProx++; fsNear="UBQ_RES"; break; }
//											if(gf.getDescription().equals("PHO_RES")) {phoProx++; fsNear="PHO_RES"; break; }
//											if(gf.getDescription().equals("CARBOHYD")) {glycProx++; fsNear="CARBOHYD"; break;}
//											if(gf.getDescription().equals("DNA_BIND")) {dnaProx++; fsNear="DNA_BIND"; break;}
//										}
//									}
//								}
//								if(fsNear == null) {
//									for(Feature f:g.getFeaturesOfTypeForPosition(FeatureType.PHOSPHOSITE, uniNbPos)) {
//										if(!alreadyObserved.contains(f)) {	// avoid double counting of the same feature
//											alreadyObserved.add(f);
//											fsNear="MOD_RES";
//											if(((PhosphoSitePlusFeature)f).getModType() == ProteinModificationType.PHOSPHORYLATION) {phoProx++;}
//											else if(((PhosphoSitePlusFeature)f).getModType() == ProteinModificationType.UBIQUITINATION) ubqProx++;
//											else modProx++;
//											break;
//										}
//									}
								}								
							}
////							if(fsNear != null) {
////								System.err.println(g.getGeneName() + " " + m.position + "->" + uniNbPos + " " + fsNear);
////							}
						}
					}
				} //else System.out.print("-"); // = no structure
			}
		}
		// wrap counts in array
		int[] featureCount = {actHits, actProx, phoHits, phoProx, ubqHits, ubqProx, modHits,modProx ,atpHits,atpProx,gtpHits,gtpProx,glycHits,glycProx,dnaHits,dnaProx, numMut};
		return featureCount;
	}
	
	/**
	 * Generates many sets of randomly scrambled mutations, counts the features for each of them and reports
	 * the sum of the features found. This is used as a test to check whether the background distribution of
	 * functional residues (counted by TargetList.countFeaturesAllResidues) is equivalent to random sets of
	 * mutations or whether there are unforeseen biases for feature size or location.
	 * Warning: As a side-effect, the loaded mutations will be scrambled after this function returns!
	 * @param targets
	 * @return
	 */
	public static double reportRandomFeatureCount(TargetList targets, int numSamples, PrintStream out, PrintStream out2, String name) {
		int[] sumFeatureCount = new int[17];
		for(int i=0; i<numSamples;i++) {
			targets.randomizeMutations(FUNSITE_WITH_REPLACEMENT);	// with replacement
			int[] featureCount = getFeatureCount(targets);
			// sum up results
			for(int j=0;j<featureCount.length;j++) {
				sumFeatureCount[j] += featureCount[j];
			}
		}
		// now sumFeatures contains the summed counts and sumFeatures[16] the summed number of mutations which should be numMutations * numSamples

		int actHits =  sumFeatureCount[0],  actProx =  sumFeatureCount[1];	// active sites
		int phoHits =  sumFeatureCount[2],  phoProx =  sumFeatureCount[3];	// phosphorilation sites
		int ubqHits =  sumFeatureCount[4],  ubqProx =  sumFeatureCount[5];	// ubiquitination sites
		int modHits =  sumFeatureCount[6],  modProx =  sumFeatureCount[7];	// other modification sites (except phospho + ubiq + glyco)
		int atpHits =  sumFeatureCount[8],  atpProx =  sumFeatureCount[9];	// atp binding sites
		int gtpHits =  sumFeatureCount[10], gtpProx =  sumFeatureCount[11];	// gtp binding sites
		int glycHits = sumFeatureCount[12], glycProx = sumFeatureCount[13];	// glycosylation sites
		int dnaHits =  sumFeatureCount[14], dnaProx =  sumFeatureCount[15];	// dna binding sites
		int numMut =   sumFeatureCount[16];
		
		int totalHits = actHits+phoHits+ubqHits+modHits+atpHits+gtpHits+glycHits+dnaHits;
		int totalProx = actProx+phoProx+ubqProx+modProx+atpProx+gtpProx+glycProx+dnaProx;
//		double fracHits = 1.0 * totalHits / numMut;
//		double fracProx = 1.0 * totalProx / numMut;
		double fracTotal = 1.0 * (totalHits+totalProx) / numMut;
		
		// report results
		System.out.printf("Number of mutations:   %6d\n", numMut);
		System.out.printf("Fraction of hits/prox: %6.2f\n", fracTotal);
		System.out.printf("sum=[%3d %3d %3d %3d %3d %3d %3d];\n", actHits+actProx, phoHits+phoProx, ubqHits+ubqProx, glycHits+glycProx+modHits+modProx, atpHits+atpProx, gtpHits+gtpProx, numMut);	
		if(out2 != null) out2.printf("%s\t%3d\t%3d\t%3d\t%3d\t%3d\t%3d\t%3d\n", name, actHits+actProx, phoHits+phoProx, ubqHits+ubqProx, glycHits+glycProx+modHits+modProx, atpHits+atpProx, gtpHits+gtpProx, numMut);
		if(out != null) out.printf("%f\t%f\t%f\t%e\n",fracTotal, fracTotal, fracTotal, 0.0);

		return fracTotal;
	}
	
	/**
	 * Iterates over loaded mutations and counts how many are proximal to or coinciding with functional sites
	 * @param silent if silent=false, prints stats per functional site type, for jackknife test and p-value
	 * @param fracAll background probability that a mutation is coindicing with or proximal to a functional sites (for p-value calculation)
	 * @return the fraction of mutations coinciding with or in proximity of at least one functional site
	 */
	public static double reportFeatureCount(TargetList targets, boolean silent, double fracAll, PrintStream out, PrintStream out2, String name) {
		
		int[] featureCount = getFeatureCount(targets);
		
		int actHits =  featureCount[0],  actProx =  featureCount[1];	// active sites
		int phoHits =  featureCount[2],  phoProx =  featureCount[3];	// phosphorilation sites
		int ubqHits =  featureCount[4],  ubqProx =  featureCount[5];	// ubiquitination sites
		int modHits =  featureCount[6],  modProx =  featureCount[7];	// other modification sites (except phospho + ubiq + glyco)
		int atpHits =  featureCount[8],  atpProx =  featureCount[9];	// atp binding sites
		int gtpHits =  featureCount[10], gtpProx =  featureCount[11];	// gtp binding sites
		int glycHits = featureCount[12], glycProx = featureCount[13];	// glycosylation sites
		int dnaHits =  featureCount[14], dnaProx =  featureCount[15];	// dna binding sites
		int numMut =   featureCount[16];

		int totalHits = actHits+phoHits+ubqHits+modHits+atpHits+gtpHits+glycHits+dnaHits;
		int totalProx = actProx+phoProx+ubqProx+modProx+atpProx+gtpProx+glycProx+dnaProx;
		if(!silent) System.out.printf("Active sites:          %3d hits %3d proximal %3d total\n", actHits, actProx, actHits+actProx);
		if(!silent) System.out.printf("Phosphorilation sites: %3d hits %3d proximal %3d total\n", phoHits, phoProx, phoHits+phoProx);
		if(!silent) System.out.printf("Ubiquitination sites:  %3d hits %3d proximal %3d total\n", ubqHits, ubqProx, ubqHits+ubqProx);
		if(!silent) System.out.printf("Glycosilation sites:   %3d hits %3d proximal %3d total\n", glycHits, glycProx, glycHits+glycProx);		
		if(!silent) System.out.printf("Other modifications:   %3d hits %3d proximal %3d total\n", modHits, modProx, modHits+modProx);
		if(!silent) System.out.printf("ATP binding sites:     %3d hits %3d proximal %3d total\n", atpHits, atpProx, atpHits+atpProx);
		if(!silent) System.out.printf("GTP binding sites:     %3d hits %3d proximal %3d total\n", gtpHits, gtpProx, gtpHits+gtpProx);		
		if(!silent) System.out.printf("DNA binding sites:     %3d hits %3d proximal %3d total\n", dnaHits, dnaProx, dnaHits+dnaProx);		
		if(!silent) System.out.printf("Total:                 %3d hits %3d proximal %3d total\n", totalHits, totalProx, totalHits+totalProx);
		if(!silent) System.out.printf("Number of mutations:   %3d\n", numMut);
		if(!silent) System.out.printf("matlab=[%3d %3d %3d %3d %3d %3d %3d];\n", actHits+actProx, phoHits+phoProx, ubqHits+ubqProx, glycHits+glycProx+modHits+modProx, atpHits+atpProx, gtpHits+gtpProx, numMut);	
		if(out2 != null) out2.printf("%s\t%3d\t%3d\t%3d\t%3d\t%3d\t%3d\t%3d\n", name, actHits+actProx, phoHits+phoProx, ubqHits+ubqProx, glycHits+glycProx+modHits+modProx, atpHits+atpProx, gtpHits+gtpProx, numMut);
		double fracHits = 1.0 * totalHits / numMut;
		double fracProx = 1.0 * totalProx / numMut;
		double fracTotal = 1.0 * (totalHits+totalProx) / numMut;
		//double fracTotal = fracHits;
		if(!silent) {
			// p-value
			double p = getPValue(fracAll,totalHits+totalProx,numMut);
			// jackknife test
			double maxFrac = 0;
			double minFrac = 100;
			for(Gene g:targets.getTargets()) {
				// exclude g from list
				TargetList oneLess = targets.getOneLess(g);
				double newFrac = reportFeatureCount(oneLess, true, Double.NaN, null, null, null);
				maxFrac = Math.max(maxFrac, newFrac);
				minFrac = Math.min(minFrac, newFrac);
			}
			// write output			
			System.out.printf("Fraction of hits:       %6.4f Prox: %6.4f Total: %6.4f [%6.4f;%6.4f] P-value: %e\n", fracHits, fracProx, fracTotal, minFrac, maxFrac, p);
			if(out != null) out.printf("%f\t%f\t%f\t%e\n",fracTotal, minFrac, maxFrac, p);
		}
		if(!silent) System.out.println();
		return fracTotal;
	}
	
	/**
	 * Create pymol sessions with visualization of mutations and other features. These can be
	 * used for generating hypotheses or as examples in the final publication.
	 */
	public static void createPyMolSessions(File outDir) {
		
	}
	
	/**
	 * Based on the analysis of all residues, assess statistical significances for the observations
	 * @returns the p-value
	 */
	public static double getPValue(double probabilityOfPositive, int numberPositiveObserved, int numberTotal) {
		double p = Double.NaN;
		try {
			double fractionPositive = 1.0 * numberPositiveObserved / numberTotal;
			if(fractionPositive > probabilityOfPositive) {
				p = MutationTools.binomialTest2(numberTotal, probabilityOfPositive, numberPositiveObserved);
			} else {
				p = MutationTools.binomialTestUnderrep(numberTotal, probabilityOfPositive, numberPositiveObserved);
			}
		} catch (MathException e) {
			e.printStackTrace();
		}
		return p;
	}
	
	/**
	 * Initial test of the input data required for this analysis
	 * @throws SQLException 
	 */
	public static void testInputData() throws SQLException {

		MySQLConnection conn = TargetList.connectToDb();
		
		System.out.println("Loading top20 mutanom targets...");
		TargetList targets = TargetList.loadTop29TargetList(conn, SEQ_SOURCE);
		
//		// skip PTEN because 1d5rA can not be handled properly with offset (alignment necessary)
//		if(targets.removeTarget("PTEN") == null) {
//			System.err.println("Could not remove target PTEN from list");
//		}
//		// skip NF1 because 1nf1A can not be handled properly with offset (alignment necessary)
//		if(targets.removeTarget("NF1") == null) {
//			System.err.println("Could not remove target NF1 from list");
//		}
//		// skip SMO because the only "structured" region is prediction garbage
//		if(targets.removeTarget("SMO") == null) {
//			System.err.println("Could not remove target SMO from list");
//		}
		
		System.out.println("Loading substructures...");
		targets.loadSubstructures(conn, false);	// not only predicted
		//File pdbDir = new File("/home/web/lappe/stehr/mutanom/pdb/");
		File pdbDir = new File("/project/PyMol/pdbs_download/pred_not_renum/");
		targets.loadPdbsAndAlignments(pdbDir);

		System.out.println();
		System.out.println("Loading mutations...");
		targets.loadMutationsFromDb(conn);
		targets.loadSnpsFromDb(conn);
		//targets.loadSnpsFromDb2(conn);		
		System.out.println("Checking mutations...");
		targets.filterMutationsObserved();
		targets.checkMutations();
		targets.printStats(conn);
		
//		for(Gene g:targets.getTargets()) {
//			System.out.print(g.getGeneName() + "\t" + g.getUniprotId() + "\t");
//			String sql = "SELECT ensg FROM uniprot_15_8.sp2ens where sp_id = \"" + g.getUniprotId() + "\"";
//			String ret = conn.getStringFromDb(sql);
//			String[] ensgs = ret.split("; ");
//			for(String ensg:ensgs) {
//				String query = "SELECT count(*) FROM mutanom.ensvar56_top21 WHERE ensg=\"" + ensg + "\" AND mut RLIKE \"^[A-Z]/[A-Z]$\" AND beg_aa = end_aa";
//				System.out.print(ensg + "\t");
//				int res = conn.getIntFromDb(query);
//				System.out.print(res + " ");
//				query = "SELECT beg_aa, mut FROM mutanom.ensvar56_top21 WHERE ensg=\"" + ensg + "\" AND mut RLIKE \"^[A-Z]/[A-Z]$\" AND beg_aa = end_aa";
//				int mismatches = 0;
//				try {
//					Statement s = conn.createStatement();
//					ResultSet rs = s.executeQuery(query);
//					while(rs.next()) {
//						int pos = rs.getInt(1);
//						String mutStr = rs.getString(2);
//						AminoAcid wt = AminoAcid.getByOneLetterCode(mutStr.charAt(0));
//						AminoAcid mut = AminoAcid.getByOneLetterCode(mutStr.charAt(2));
//						if(g.getUniprotSeq().length() < pos) {
//							//System.out.print(pos + ">" + g.getUniprotSeq().length() + " ");
//							mismatches++;
//						} else if (g.getUniprotSeq().charAt(pos-1) != mutStr.charAt(0)) {
//							//System.out.print(g.getUniprotSeq().charAt(pos-1) +"!="+ mutStr.charAt(0)+ " ");
//							mismatches++;
//						}
//					}
//				} catch(SQLException e) {
//					e.printStackTrace();
//				}
//				System.out.print(mismatches + " ");
//			}
//			System.out.println();
//		}
		
//		System.out.println();
//		System.out.println("Loading cosmic mutations");
//		targets.loadMutationsFromDb(conn);
//		targets.checkMutations();
//		targets.printStats();
		
//		System.out.println();
//		System.out.println("Loading SNPs");
//		targets.loadSnpsFromDb(conn);
//		targets.printStats();
		
		conn.close();
	}
	
	/**
	 * Create pymol session with visualization of mutations and SNPs. These can be
	 * used for generating hypotheses or as examples in the final publication.
	 * Note that currently, all structures are taken from the renumbered PDB files (so no mapping of residue numbers is neccessary).
	 * @throws SQLException 
	 */
	public static void visualizeMutationsInPymol() throws SQLException {
		
		MySQLConnection conn = TargetList.connectToDb();
		
		System.out.println("Loading top20 mutanom targets...");
		TargetList targets = TargetList.loadTop29TargetList(conn, SEQ_SOURCE);

		// skip PTEN because 1d5rA can not be handled properly with offset (alignment necessary)
		if(targets.removeTarget("PTEN") == null) {
			System.err.println("Could not remove target PTEN from list");
		}
		// skip NF1 because 1nf1A can not be handled properly with offset (alignment necessary)
		if(targets.removeTarget("NF1") == null) {
			System.err.println("Could not remove target NF1 from list");
		}
		// skip SMO because the only "structured" region is prediction garbage
		if(targets.removeTarget("SMO") == null) {
			System.err.println("Could not remove target SMO from list");
		}
		
		System.out.println("Loading substructures...");
		targets.loadSubstructures(conn, false);	// not only predicted
		
		File pdbDir = new File("/home/web/lappe/stehr/mutanom/pdb/");
		targets.loadPdbsAndAlignments(pdbDir);

		System.out.println();
		System.out.println("Loading mutations and SNPs...");
		//targets.loadMutations();
		targets.loadMutationsFromDb(conn);
		targets.loadSnpsFromDb(conn);
		targets.printStats(conn);
		
		// filter out unwanted mutations
		System.out.println("Filter mutations:");
		System.out.println("Non missense: " + targets.filterMutationsMissense());
		System.out.println("No structure: " + targets.filterMutationsKnownStructure());
		System.out.println();
		System.out.println("Filter SNPs:");
		System.out.println("Non missense: " + targets.filterSNPsMissense());
		System.out.println("No structure: " + targets.filterSNPsKnownStructure());
		System.out.println();

		File sessionFile = new File("/project/StruPPi/henning/projects/mutanom/analysis/top16.pse");
		System.out.println("Writing session file " + sessionFile);
		
		try {
			targets.generatePseWithAllMutationsAndSNPs(sessionFile, pdbDir);
		} catch (IOException e) {
			System.err.println("Unable to write PyMol session with mutation and SNPs: " + e.getMessage());
		}
	}
	
	/**
	 * Prints a list of the most frequently found prosite motifs in this set
	 */
	public static void findFrequentPrositeMotifs() {
		MySQLConnection conn = TargetList.connectToDb();
		System.out.println("Loading top29 mutanom targets...");
		TargetList targets = TargetList.loadTop29TargetList(conn, SEQ_SOURCE);
		targets.findFrequentPrositeMotifs();
	}
	
	public static void findFrequentUniprotMotifs() {
		MySQLConnection conn = TargetList.connectToDb();
		System.out.println("Loading top29 mutanom targets...");
		TargetList targets = TargetList.loadTop29TargetList(conn, SEQ_SOURCE);
		targets.findFrequentUniprotMotifs();		
	}
	
	/**
	 * @param args
	 * @throws FileNotFoundException 
	 * @throws SQLException 
	 */
	public static void main(String[] args) throws FileNotFoundException, SQLException {
				
		// supress logger output of JAligner
		try {
			File trashLogFile = null;
			trashLogFile = File.createTempFile("JAligner", ".log");
			trashLogFile.deleteOnExit();
			System.setProperty("java.util.logging.config.file",trashLogFile.getAbsolutePath());
		} catch (IOException e) {
			System.err.println("Could not create temporary logfile in default temp directory.");
		}
		
		// perform analysis
		
		//testInputData();
		analyzePredictedStructures();
		analyzeMutations();
		
		//debugSingleStructure();
		//visualizeMutationsInPymol();
		//findFrequentPrositeMotifs();
		//findFrequentUniprotMotifs();
		
		System.out.println("done.");

	}

	/**
	 * Writes the sequence overview with mutations for web, see http://www.molgen.mpg.de/~lappe/stehr/mutanom/
	 */
	public static void writeSeqOverviewForWeb(File outDir, TargetList tl, SubstructureType restrictToType) {
		tl.renderRulers(outDir);
		tl.writeGeneDataFiles(outDir, restrictToType);
		//tl.generatePseWithAllMutations(pseDir, pdbDir);	// for paper: only needed for MLH1 and ERBB2
		
	}

}

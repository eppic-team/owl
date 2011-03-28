package owl.mutanom;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;

import owl.core.features.Feature;
import owl.core.features.FeatureType;
import owl.core.runners.MaxClusterRunner;
import owl.core.runners.MaxClusterRunner.ScoreType;
import owl.core.sequence.alignment.PairwiseSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
//import owl.core.structure.Pdb;
import owl.core.util.MySQLConnection;
import owl.mutanom.core.Gene;
import owl.mutanom.core.Substructure;
import owl.mutanom.core.TargetList;

/**
 * Helper class to print information about candidate genes. We use this for the final gene selection
 * and to detect possible inconsistencies.
 * @author stehr
 */
public class GeneInfo {

	/*------------------------------ constants ------------------------------*/
	
	// see: Cfg.java
	
	// local constants (only used here)
	
	public static final String MAXCLUSTER_EXECUTABLE = "/project/StruPPi/bin/maxcluster";
	public static final String PDB_DIR = "/project/StruPPi/henning/projects/mutanom/analysis/data/pdb/download/singlechains";
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Print pairwise sequence similarities of the targets.
	 */
	public static void printSeqSimilarities(MySQLConnection conn, TargetList targets) throws SQLException {
		
		ArrayList<Gene> orderedGenes = new ArrayList<Gene>();
		orderedGenes.addAll(targets.getTargets());
		
		System.out.println("-- Overall Sequence Similarity --");
		
		// write header
		System.out.print("Gene\t");
		for (int i = 0; i < orderedGenes.size(); i++) {
			Gene g = orderedGenes.get(i);
			System.out.print(g.getGeneName() + "\t");
		}
		System.out.println();
		
		// sequence similarity among cosmic seqs
		for (int i = 1; i < orderedGenes.size(); i++) {
			Gene g1 = orderedGenes.get(i);
			System.out.print(g1.getGeneName() + "\t");
			for (int j = 0; j < i; j++) {
				Gene g2 = orderedGenes.get(j);
				String s1 = g1.getCosmicSeq();
				String s2 = g2.getCosmicSeq();
					try {
						PairwiseSequenceAlignment al = new PairwiseSequenceAlignment(s1,s2,g1.getGeneName(), g2.getGeneName());
						float score = al.getPercentSimilarity();
						if(score < 50) {
							System.out.printf(" %2.0f\t", score);
						} else {
							System.out.printf("*%2.0f*\t", score);
						}
					} catch (PairwiseSequenceAlignmentException e) {
						System.err.println("Error calculating alignment between " + g1.getGeneName() + " and " + g2.getGeneName());
					}
			}
			System.out.println();
		}		
	}
	
	/**
	 * Print pairwise structural similarities of the substructures.
	 */
	public static void printStrucSimilarities(MySQLConnection conn, TargetList targets) throws SQLException {
		
		// load structures
		System.out.println("Loading structures...");
		targets.loadSubstructures(conn, false);	// not only predicted
		targets.loadPdbsAndAlignments(null);
		
		// gather substructures
		ArrayList<Substructure> pdbs = new ArrayList<Substructure>();
		for(Gene g: targets.getTargets()) {
			pdbs.addAll(g.getSubstructures());
		}
		
		System.out.println("-- Structure Similarity --");
		
		// write header
		System.out.print("Pdb\t");
		for (int i = 0; i < pdbs.size(); i++) {
			Substructure ss = pdbs.get(i);
			System.out.print(ss.getPdbCode()+ss.getChainCode() + "\t");
		}
		System.out.println();
		
		// sequence similarity among cosmic seqs
		for (int i = 1; i < pdbs.size(); i++) {
			Substructure ss1 = pdbs.get(i);
			System.out.print(ss1.getPdbCode()+ss1.getChainCode() + "\t");
			for (int j = 0; j < i; j++) {
				Substructure ss2 = pdbs.get(j);
				String pdb1 = new File(PDB_DIR, ss1.getPdbCode()+ss1.getChainCode()+".pdb").getPath();
				String pdb2 = new File(PDB_DIR, ss2.getPdbCode()+ss2.getChainCode()+".pdb").getPath();
				try {
					MaxClusterRunner mc = new MaxClusterRunner("/project/StruPPi/bin/maxcluster");
					double score = mc.calculatePairwiseScore(pdb1, pdb2, ScoreType.GDT);
					if(score < 40) {
						System.out.printf(" %2.0f\t", score);
					} else {
						System.out.printf("*%2.0f*\t", score);
					}
				} catch (IOException e) {
					System.err.println("Error calculating alignment between " + pdb1 + " and " + pdb2);
				}
			}
			System.out.println();
		}		
	}
		
	/**
	 * Convenience method to load substructures, pdbs & alignments, features & graphs & domains, mutations and SNPs.
	 */
	public static void loadEverything(MySQLConnection conn, TargetList targets) throws SQLException {
		
		// loading substructures
		System.out.println("Loading substructures...");
		boolean onlyPredicted = false;	// this is mainly for xray structures
		targets.loadSubstructures(conn, onlyPredicted);
		
		// loading PDBs, alingments, features and graphs
		targets.loadPdbsAndAlignments(Cfg.PRED_PDBS_PATH); // predicted structures not currently used
		targets.loadFeaturesAndGraphsAndDomains(Cfg.PHOSITE_HTML_PATH, Cfg.RIG_TYPE, Cfg.RIG_CUTOFF, Cfg.DOMAIN_TYPE, false);

		// loading mutations and SNPs
		System.out.println();
		System.out.println("Loading mutations...");
		targets.loadMutationsFromDb(conn);
		targets.loadSnpsFromDb(conn);
	}
	
	public static void writeSeqOverviews(MySQLConnection conn, TargetList targets) throws SQLException {
		
		loadEverything(conn, targets);
		
		TargetList oncoGenes = targets.getOncogenes(conn);				// oncogenes (according to db table genes_onc_sup)
		TargetList tumSupGenes = targets.getTumorSuppressors(conn);		// tum sup genes (according to db table genes_onc_sup)
		TargetList otherGenes = targets.getUnAnnotated(conn);			// all others
		
		System.out.println("Oncogenes:\t\t" + oncoGenes.getTargets().size());
		System.out.println("Tumor Suppressors:\t" + tumSupGenes.getTargets().size());
		System.out.println("Others:\t\t\t" + otherGenes.getTargets().size());
		System.out.println();
		
		File outFile;
		// onc
		outFile = new File(Cfg.RESULT_DIR, "svg/seqoverview_onc.svg");
		System.out.println("Writing " + outFile);
		oncoGenes.writeSeqOverviewForPrint(outFile);
		// sup
		outFile = new File(Cfg.RESULT_DIR, "svg/seqoverview_sup.svg");
		System.out.println("Writing " + outFile);
		tumSupGenes.writeSeqOverviewForPrint(outFile);
		// other
		outFile = new File(Cfg.RESULT_DIR, "svg/seqoverview_oth.svg");
		System.out.println("Writing " + outFile);
		otherGenes.writeSeqOverviewForPrint(outFile);

	}
	
	public static void writePymolSessions(MySQLConnection conn, TargetList targets) throws SQLException {
		
		loadEverything(conn, targets);
		
		try {
			targets.generatePseWithMutSnpFuncAndDomains(Cfg.PSE_OUT_DIR, new File(PDB_DIR));	// use single chain files
		} catch (IOException e) {
			System.err.println("Unable to write PyMol session with mutation and SNPs: " + e.getMessage());
		}
		
	}
	
	/**
	 * Prints statistics about the given target list.
	 * TODO: Perform a number of checks (e.g. matching sequence)
	 * @throws SQLException 
	 */
	public static void printGeneInfo(MySQLConnection conn, TargetList targets) throws SQLException {
		
		loadEverything(conn, targets);
		
		// print statistics
		targets.printStats(conn);
		targets.printFunctionalSiteStats();
		System.out.println();
		
	}
	
	/**
	 * Prints the Cosmic2Uniprot Alignments for the genes in the targets list.
	 * @param targets
	 */
	public static void printCosmic2UniprotAlignments(TargetList targets) {
		for(Gene g:targets.getTargets()) {
			System.out.println(g.getGeneName());
			PairwiseSequenceAlignment al;
			try {
				al = new PairwiseSequenceAlignment(g.getUniprotSeq(), g.getCosmicSeq(), "Uniprot", "Cosmic");
				al.printAlignment();
			} catch (PairwiseSequenceAlignmentException e) {
				System.err.println("Error calculating sequence alignment: " + e.getMessage());
			}
			System.out.println();
		}
		
	}
	
	/**
	 * Prints out all stored features.
	 * @throws SQLException 
	 */
	public static void printAllFeatures(MySQLConnection conn, TargetList targets) throws SQLException {
		
		loadEverything(conn, targets);
		
		for(Gene g: targets.getTargets()) {
			System.out.println("--- " + g.getGeneName() + " ---");
			for(FeatureType t:g.getFeatureTypes()) {
				System.out.println(t.toString());
				for(Feature f:g.getFeaturesOfType(t)) {
					System.out.print(f.toString());
				}
				System.out.println();
			}			
			System.out.println();
		}
	}
		
	/*--------------------------------- main --------------------------------*/
	
	public static void main(String[] args) throws SQLException {
		
		// supress logger output of JAligner
		try {
			File trashLogFile = null;
			trashLogFile = File.createTempFile("JAligner", ".log");
			trashLogFile.deleteOnExit();
			System.setProperty("java.util.logging.config.file",trashLogFile.getAbsolutePath());
		} catch (IOException e) {
			System.err.println("Could not create temporary logfile in default temp directory.");
		}
		
		MySQLConnection conn = new MySQLConnection();
		
		System.out.println("Loading targets...");
		//TargetList targets = TargetList.loadTargets(conn, Cfg.OFFLINE);
		TargetList targets = TargetList.loadSingletonList(conn, Cfg.OFFLINE, "ERBB2");		
		//TargetList targets = TargetList.loadPredictedTargets(conn, Cfg.OFFLINE);
		System.out.println(targets.getTargets().size() + " targets loaded.");
		
		// printGeneInfo(conn, targets);	// number of (real) structural mutations, snps, domains, coverage, etc.
		
		// writeSeqOverviews(conn, targets);	// write 3 SVG files for Onc, Sup and others
		
		writePymolSessions(conn, targets);	// write PSE files for each gene with all mutations, SNPs, functional sites and domains
		
		// printAllFeatures(conn, TargetList.getSingletonList(conn, Cfg.OFFLINE, "CTNNB1"));
		
		// printSeqSimilarities(conn, targets);
		
		// printStrucSimilarities(conn, targets);
		
		// printCosmic2UniprotAlignments(targets);
		
		System.out.println("done.");
		
	}
	
}

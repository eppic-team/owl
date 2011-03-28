package owl.mutanom;
import java.io.*;
import java.util.*;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BinomialDistribution;
import org.apache.commons.math.distribution.BinomialDistributionImpl;

import owl.core.runners.NaccessRunner;
import owl.core.structure.Pdb;
import owl.core.structure.PdbSet;
import owl.core.util.Pair;


/**
 * Some (statistical and other) tools for analyzing sets of mutated structures.
 * TODO: Merge with owl.core.util.Statistics
 * @author stehr
 */
public class MutationTools {

	/*------------------------------ constants ------------------------------*/
	private static final String NACCESS_EXECUTABLE = "/project/StruPPi/bin/naccess";
	private static final String NACCESS_PARAMETERS = "";
	private static final double EXPOSURE_CUTOFF = 5.0; // everything above this cutoff is considered exposed
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Performs a binomial test with the given parameters and returns the value P(X >= observations).
	 * @param trials the total number of trials
	 * @param p probability of a positive observation
	 * @param observations the observed number of positive observations
	 * @return the pobability P(X >= observations)
	 * @throws MathException
	 */
	public static double binomialTest(int trials, double p, int observations) throws MathException {
		BinomialDistribution binom = new BinomialDistributionImpl(trials, p);
		double probSum = 0;
		for (int i = observations; i <= trials; i++) {
			probSum += binom.probability(i);
		}
		
		return probSum;
	}

	/**
	 * Performs a binomial test with the given parameters and returns the value P(X >= observations)
	 * using an alternative implementation which differs slightly in the results  due to rounding issues.
	 * @param trials the total number of trials
	 * @param p probability of a positive observation
	 * @param observations the observed number of positive observations
	 * @return the pobability P(X >= observations)
	 * @throws MathException
	 */
	public static double binomialTest2(int trials, double p, int observations) throws MathException {
		BinomialDistribution binom = new BinomialDistributionImpl(trials, p);
		double prob = 1 - binom.cumulativeProbability(observations-1);
		return prob;
	}
	
	/**
	 * Performs a binomial test with the given parameters and returns the value P(X <= observations).
	 * @param trials the total number of trials
	 * @param p probability of a positive observation
	 * @param observations the observed number of positive observations
	 * @return the pobability P(X <= observations)
	 * @throws MathException
	 */
	public static double binomialTestUnderrep(int trials, double p, int observations) throws MathException {
		BinomialDistribution binom = new BinomialDistributionImpl(trials, p);
		double prob = binom.cumulativeProbability(observations-1);
		return -prob;
	}
	
	/**
	 * Performs a binomial test and prints both input parameters and outcome.
	 * @param sampleSize the total number of trials
	 * @param p probability of a positive observation
	 * @param observed the observed number of positive observations
	 * @throws MathException
	 */
	public static void printBinomialTestResult(int sampleSize, double p, int observed) throws MathException {
		double sigma = binomialTest(sampleSize, p, observed);
		System.out.println("Sample size[n]: " + sampleSize);
		System.out.println("Probability of positive outcome[p]: " + p);
		System.out.println("Number of positive outcomes[x]: " + observed);
		System.out.println("Probability to observe x or more positive outcomes: " + sigma);
	}
	
	/**
	 * For a set of pdbs, count the fraction of residues which are exposed.
	 * @param pdbs
	 */
	public static double countFractionExposed(PdbSet pdbs) {
		int exposed = 0;
		int buried = 0;
		for(Pdb pdb:pdbs.getPdbs()) {
			System.out.print(".");
			try {
				NaccessRunner nar = new NaccessRunner(new File(NACCESS_EXECUTABLE), NACCESS_PARAMETERS);
				nar.runNaccess(pdb);
				for(int pos:pdb.getAllSortedResSerials()) {
					if(pdb.getAllRsaFromResSerial(pos) > EXPOSURE_CUTOFF) {
						exposed++;
					} else {
						buried++;
					}
				}
			} catch (IOException e) {
				System.err.println("Error running NACCESS: " + e.getMessage());
				//System.exit(1);
			}
		}
		System.out.println();
		return 1.0 * exposed / (exposed+buried);
	}
	
	public static double countFractionExposed(Collection<String> list) {
		int buried = 0;
		int exposed = 0;
		
		for(String pdbCode:list) {
			Pdb pdb = Pdb.readStructureOrExit(pdbCode);
			Pair<Integer,Integer> pair = countBuriedAndExposed(pdb);
			buried += pair.getFirst();
			exposed += pair.getSecond();
		}		
		return 1.0 * buried/exposed;
		
	}
	
	public static Pair<Integer,Integer> countBuriedAndExposed(Pdb pdb) {
		int exposed = 0;
		int buried = 0;	
		System.out.print(".");
		try {
			NaccessRunner nar = new NaccessRunner(new File(NACCESS_EXECUTABLE), NACCESS_PARAMETERS);
			nar.runNaccess(pdb);
			//pdb.runNaccess(NACCESS_EXECUTABLE, NACCESS_PARAMETERS);
			for(int pos:pdb.getAllSortedResSerials()) {
				if(pdb.getAllRsaFromResSerial(pos) > EXPOSURE_CUTOFF) {
					exposed++;
				} else {
					buried++;
				}
			}
		} catch (IOException e) {
			System.err.println("Error running NACCESS: " + e.getMessage());
			//System.exit(1);
		}
		return new Pair<Integer, Integer>(new Integer(buried), new Integer(exposed));
	}
	
	/**
	 * @param args
	 * @throws MathException 
	 */
	public static void main(String[] args) throws MathException {	
		double fractionExposed = 0.35654072179644847; // taken from cullPdb20 and exposureCutoff 5.0
		//int sampleSize = 20;
		int sampleSize = 111;
		int numExposed = 93;
		
		printBinomialTestResult(sampleSize, fractionExposed, numExposed);
		
		// --- Print significance for a range of outcomes to estimate a reasonable sample size
		// What's the probability that to see i or more exposed residues under the random model?
//		System.out.println("Binomial_1");
//		long millis = System.currentTimeMillis();
//		for(int i=0; i <= sampleSize; i++) {
//			System.out.printf("%2d %4.2f\n", i, 100 * binomialTest(sampleSize, fractionExposed, i));			
//		}
//		System.out.println("Time: " + ((System.currentTimeMillis() - millis)));
//		System.out.println();
//		System.out.println("Binomial_2");
//		millis = System.currentTimeMillis();
//		for(int i=0; i <= sampleSize; i++) {
//			System.out.printf("%2d %4.2f\n", i, 100 * binomialTest2(sampleSize, fractionExposed, i));			
//		}		
//		System.out.println("Time: " + ((System.currentTimeMillis() - millis)));
//		// Result: BinomialTest2 performs better, slightly different results, probably from rounding

		// --- Count fraction exposed over Cullpdb20 set
//		System.out.println("Reading cullpdb_20...");
//		Collection<String> pdbCodes = null;
//		try {
//			pdbCodes = PdbSet.readCullPdb20List();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		System.out.println("Calculating exposed fraction...");
//		System.out.println("Fraction exposed: " + countFractionExposed(pdbCodes));
	}

}

package proteinstructure;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.HashMap;

/**
 * Class to handle statistics of a single decoy set (a native structure plus 
 * a number of decoys of that structure). Contains the statistics for the set:
 * z-score, correlation, etc.
 * Static methods to write sets of statistics to file are also provided.
 * 
 * @author duarte
 *
 */
public class DecoySetStats {

	public String decoyName;
	public int numScoredDecoys;
	public double zscore;
	public double corr;
	public boolean isNativeRank1;

	public DecoySetStats(DecoyScoreSet set, String nativeFileName) {
		this.decoyName = set.getDecoyName();
		this.numScoredDecoys = set.size();
		this.zscore = set.getZscore(nativeFileName);
		this.corr = set.getSpearman();
		this.isNativeRank1 = set.isRank1(nativeFileName);
	}

	/**
	 * Writes the scoring statistics of a group of decoy sets to text file with 5 columns:
	 * decoy name, number of decoys, is native ranked 1, z-score of native, correlation    
	 * @param file
	 * @param stats
	 */
	public static void writeStats(File file, HashMap<String,DecoySetStats> stats) {
		try {
			PrintWriter pw = new PrintWriter(file);
			pw.printf("#%9s\t%6s\t%6s\t%6s\t%6s\n",
					"decoy","nod","rank1","z","corr");
			double sumz = 0,sumcorr = 0;
			int countrank1 = 0;
			for (String decoy:stats.keySet()) {
				int numScDecoys = stats.get(decoy).numScoredDecoys;
				double z = stats.get(decoy).zscore;
				double corr = stats.get(decoy).corr;
				boolean isRank1 = stats.get(decoy).isNativeRank1;
				sumz+=z;
				sumcorr+=corr;
				if (isRank1) countrank1++;
				pw.printf("%10s\t%6d\t%6s\t%6.1f\t%5.2f\n",
						decoy,numScDecoys,isRank1,z,corr);
			}
			int N = stats.size();
			pw.println();
			pw.printf("#%9s\t%6s\t%6s\t%6.1f\t%5.2f\n",
					"means","",countrank1+"/"+N,sumz/N,sumcorr/N);
			pw.close();
		} catch (FileNotFoundException e) {
			System.err.println("Couldn't write stats file "+file);
		}
	}
	
	/**
	 * Writes the scoring statistics (2 scorings: residue-based and atom-based) of a group 
	 * of decoy sets to text file with 8 columns:
	 * decoy name, number of decoys, is native ranked 1 (res-based scoring), z-score of 
	 * native (res-based scoring), correlation (res-based scoring), is native ranked 1 
	 * (atom-based scoring), z-score of native (atom-based scoring), correlation (atom-based 
	 * scoring).     
	 * @param file
	 * @param resStats
	 * @param atomStats
	 */
	public static void writeStats(File file, HashMap<String,DecoySetStats> resStats, HashMap<String,DecoySetStats> atomStats) {
		try {
			PrintWriter pw = new PrintWriter(file);
			pw.printf("#%9s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\n",
					"decoy","nod","res_r1","resz","rescor","atom_r1","atomz","atomcor");
			double sumresz = 0,sumrescor = 0, sumatomz = 0, sumatomcor = 0;
			int countresr1 = 0, countatomr1 = 0;
			for (String decoy:resStats.keySet()) {
				int resNumScDec = resStats.get(decoy).numScoredDecoys;
				double resz = resStats.get(decoy).zscore;
				double rescor = resStats.get(decoy).corr;
				boolean resr1 = resStats.get(decoy).isNativeRank1;
				int atomNumScDec = atomStats.get(decoy).numScoredDecoys;
				double atomz = atomStats.get(decoy).zscore;
				double atomcor = atomStats.get(decoy).corr;
				boolean atomr1 = atomStats.get(decoy).isNativeRank1;
				if (resNumScDec!=atomNumScDec) 
					System.err.println("Warning: number of residue-scored decoys doesn't coincide with number of atom-scored decoys");
				sumresz+=resz;
				sumrescor+=rescor;
				sumatomz+=atomz;
				sumatomcor+=atomcor;
				if (resr1) countresr1++;
				if (atomr1) countatomr1++;
				pw.printf("%10s\t%6d\t%6s\t%6.1f\t%5.2f\t%6s\t%6.1f\t%5.2f\n",
						decoy,resNumScDec,resr1,resz,rescor,atomr1,atomz,atomcor);
			}
			int N = resStats.size();
			pw.println();
			pw.printf("#%9s\t%6s\t%6s\t%6.1f\t%5.2f\t%6s\t%6.1f\t%5.2f\n",
					"means","",countresr1+"/"+N,sumresz/N,sumrescor/N,countatomr1+"/"+N,sumatomz/N,sumatomcor/N);
			pw.close();
		} catch (FileNotFoundException e) {
			System.err.println("Couldn't write stats file "+file);
		}
	}


}

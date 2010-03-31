package owl.decoyScoring;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;

import owl.core.structure.FileFormatError;
import owl.core.util.MySQLConnection;
import owl.core.util.RegexFileFilter;



/**
 * Class to represent a set of decoy sets, i.e. a group of decoys generated 
 * with the same method.
 *  
 * @author duarte
 *
 */
public class DecoyScoreSetsGroup {
	
	ArrayList<DecoyScoreSet> setsGroup;

	public DecoyScoreSetsGroup() {
		setsGroup = new ArrayList<DecoyScoreSet>();
	}
	
	public DecoyScoreSetsGroup(File dir, String decoySetsGroup) {
		setsGroup = new ArrayList<DecoyScoreSet>();
		readGroupOfSets(dir, decoySetsGroup);
	}
	
	public void addDecoyScoreSet(DecoyScoreSet set) {
		setsGroup.add(set);
	}
	
	public DecoyScoreSet getDecoyScoreSet(int i) {
		return setsGroup.get(i);
	}
	
	public void readGroupOfSets(File dir, String decoySetsGroup) {

		File[] files = dir.listFiles(new RegexFileFilter("^"+decoySetsGroup+"_.*\\.scores"));
		for (File file:files) {
			try {
				setsGroup.add(new DecoyScoreSet(file));
			} catch (FileFormatError e) {
				System.err.println(e.getMessage());
				continue;
			} catch (IOException e) {
				System.err.println("Couldn't read scores file "+file+": "+e.getMessage());
				continue;
			}
		}
	}
	
	public int size() {
		return setsGroup.size();
	}
	
	public double getMeanNativeZscore() {
		double sumzscore = 0;
		for (DecoyScoreSet set:setsGroup) {
			sumzscore+=set.getNativeZscore();
		}
		return sumzscore/this.size();
	}
	
	public double getMeanSpearman() {
		double sumcorr = 0;
		for (DecoyScoreSet set:setsGroup) {
			sumcorr+=set.getSpearman();
		}
		return sumcorr/this.size();
	}
	
	/**
	 * Gets the average number of decoys with rmsd below given rmsd for this group.
	 * @param rmsd
	 * @return
	 */
	public double getMeanNumBelowRmsd(double rmsd) {
		int sumnum = 0;
		for (DecoyScoreSet set:setsGroup) {
			sumnum+=set.getNumDecoysBelowRmsd(rmsd);
		}
		return (double)sumnum/(double)this.size();
	}
	
	/**
	 * Gets the average rank of the native for this group.
	 * @return
	 */
	public double getMeanNativeRank() {
		double sumNativeRank = 0;
		for (DecoyScoreSet set:setsGroup) {
			sumNativeRank+=set.getNativeRank();
		}
		return sumNativeRank/(double)this.size();		
	}
	
	/**
	 * Gets the number of times that the native has rank 1 in this group.
	 * @return
	 */
	public int getNumberOfRank1() {
		int numRank1 = 0;
		for (DecoyScoreSet set:setsGroup) {
			if (set.isNativeRank1()) numRank1++;
		}
		return numRank1;
	}
	
	/**
	 * Gets the mean of the delta RMSD measures. See DecoyScoreSet.getDeltaRmsd().
	 * @return
	 */
	public double getMeanDeltaRmsd() {
		double sumDeltaRMSD = 0;
		for (DecoyScoreSet set:setsGroup) {
			sumDeltaRMSD+=set.getDeltaRmsd();
		}
		return sumDeltaRMSD/(double)this.size();				
	}
	
	/**
	 * Gets the mean n% enrichment value. See DecoyScoreSet.getNpercentEnrichment().
	 * @param n
	 * @return
	 */
	public double getMeanNpercentEnrichment(int n) {
		double sumEnrich = 0;
		for (DecoyScoreSet set:setsGroup) {
			sumEnrich+=set.getNpercentEnrichment(n);
		}
		return sumEnrich/(double)this.size();						
	}
	
	public void writeToDb(MySQLConnection conn, String db, String table) throws SQLException {
		for (DecoyScoreSet set:setsGroup) {
			set.writeToDb(conn, db, table);
		}
	}

	/**
	 * Writes the scoring statistics of a group of decoy sets to text file with 5 columns:
	 * decoy name, number of decoys, is native ranked 1, z-score of native, correlation    
	 * @param file
	 * @param stats
	 */
	public void writeStatsToFile(File file) {
		try {
			PrintWriter pw = new PrintWriter(file);
			pw.printf("#%9s\t%6s\t%13s\t%8s\t%8s\t%8s\t%6s\t%6s\t%6s\t%6s\t%9s\t%9s\t%9s\n",
					"decoy","nod","rmsdRange","rmsd<1.0","rmsd<1.5","rmsd<2.0","rank","rank1","z","corr","deltaRMSD","10%enrich","20%enrich");

			for (int i=0;i<setsGroup.size();i++) {
				String decoy = setsGroup.get(i).getDecoyName();
				int numScDecoys = setsGroup.get(i).size();
				String rmsdRange = String.format("%6.3f-%6.3f",setsGroup.get(i).getRmsdRange()[0],setsGroup.get(i).getRmsdRange()[1]);
				int numBelow1 = setsGroup.get(i).getNumDecoysBelowRmsd(1.0);
				int numBelow1_5 = setsGroup.get(i).getNumDecoysBelowRmsd(1.5);
				int numBelow2 = setsGroup.get(i).getNumDecoysBelowRmsd(2.0);
				int rank = setsGroup.get(i).getNativeRank();
				boolean isRank1 = setsGroup.get(i).isNativeRank1();
				double z = setsGroup.get(i).getNativeZscore();
				double corr = setsGroup.get(i).getSpearman();
				double deltaRMSD = setsGroup.get(i).getDeltaRmsd();
				double tenPercentEnrich = setsGroup.get(i).getNpercentEnrichment(10);
				double twentyPercentEnrich = setsGroup.get(i).getNpercentEnrichment(20);

				pw.printf("%10s\t%6d\t%13s\t%8d\t%8d\t%8d\t%6d\t%6s\t%6.1f\t%6.2f\t%9.3f\t%9.2f\t%9.2f\n",
						decoy,numScDecoys,rmsdRange,numBelow1,numBelow1_5,numBelow2,rank,isRank1,z,corr,deltaRMSD,tenPercentEnrich,twentyPercentEnrich);
			}
			
			pw.println();
			pw.printf("#%9s\t%6s\t%13s\t%8.1f\t%8.1f\t%8.1f\t%4.1f\t%6s\t%6.1f\t%6.2f\t%9.3f\t%9.2f\t%9.2f\n",
					"means","","",
					this.getMeanNumBelowRmsd(1.0),
					this.getMeanNumBelowRmsd(1.5),
					this.getMeanNumBelowRmsd(2.0),
					this.getMeanNativeRank(),
					this.getNumberOfRank1()+"/"+this.size(),
					this.getMeanNativeZscore(),
					this.getMeanSpearman(),
					this.getMeanDeltaRmsd(),
					this.getMeanNpercentEnrichment(10),
					this.getMeanNpercentEnrichment(20));
			pw.close();
		} catch (FileNotFoundException e) {
			System.err.println("Couldn't write stats file "+file);
		}
	}

}

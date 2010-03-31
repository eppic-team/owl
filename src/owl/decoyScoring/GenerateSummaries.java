package owl.decoyScoring;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeMap;

import owl.core.util.FileFormatError;
import owl.core.util.RegexFileFilter;



public class GenerateSummaries {

	private static final File outDir = new File("/project/StruPPi/jose/emp_potential");
	private static final File resultsDir = new File("/project/StruPPi/jose/emp_potential");
	public static final String[] DECOYSETS = 
	{"4state_reduced", "fisa", "fisa_casp3", "hg_structal", "ig_structal", "ig_structal_hires", "lattice_ssfit", "lmds", "vhp_mcmd"};

	
	public static void main(String[] args) throws IOException, SQLException {
		
		
		File[] scoresDirs = resultsDir.listFiles(new RegexFileFilter("(res|atom)(type|count|comb)_([^_]+).*"));
		File[] scoresDirsBestCutoffs = resultsDir.listFiles(new RegexFileFilter("(res|atom)(type|count|comb)_(4|8).*"));
		
		File meansFile = new File(outDir,"means.summary");
		PrintWriter mf = new PrintWriter(meansFile);
		
		mf.printf("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
						"dir","type","cutoff","seq separation","size","mean proportion rank1","mean z","mean corr");
		
		// means summary
		for (File scoresDir:scoresDirs) {
		
			double sumz = 0, sumcorr = 0, sumProportionOfRank1 = 0;
			String type = "";
			double cutoff = 0;
			int seqSep = -1;

			for (String decoySet:DECOYSETS) {
				DecoyScoreSetsGroup group = new DecoyScoreSetsGroup(scoresDir,decoySet);
				DecoyScoreSet firstSet = group.getDecoyScoreSet(0);
				type = firstSet.getScoringMethod().getId();
				cutoff = firstSet.getCutoff();
				seqSep = firstSet.getMinSeqSep();
				// calculating means for means.summary
				sumz+=group.getMeanNativeZscore();
				sumcorr+=group.getMeanSpearman();
				sumProportionOfRank1+=(double)group.getNumberOfRank1()/(double)group.size();
			
			}

			int n = DECOYSETS.length;
			mf.printf("%12s\t%10s\t%3.1f\t%2d\t%3d\t%6.2f\t%6.2f\t%6.2f\n",
					scoresDir.getName(),
					type,
					cutoff,
					seqSep,
					n,
					sumProportionOfRank1/n,
					sumz/n,
					sumcorr/n);
			
			
		}
		mf.close();

		// summaries per decoySet
		for (String decoySet:DECOYSETS) {
			
			PrintWriter s = new PrintWriter(new File(outDir,decoySet+".summary"));
			s.printf("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					"dir","type","cutoff","seq separation","size","proportion rank1","z-score","corr");

			for (File scoresDir:scoresDirs) {
				DecoyScoreSetsGroup group = new DecoyScoreSetsGroup(scoresDir,decoySet);
				s.printf("%12s\t%10s\t%3.1f\t%2d\t%3d\t%6.2f\t%6.2f\t%6.2f\n",
						scoresDir.getName(),
						group.getDecoyScoreSet(0).getScoringMethod().getId(),
						group.getDecoyScoreSet(0).getCutoff(),
						group.getDecoyScoreSet(0).getMinSeqSep(),
						group.size(),
						(double)group.getNumberOfRank1()/(double)group.size(),
						group.getMeanNativeZscore(),
						group.getMeanSpearman());
			}
			s.close();
			
			// different methods for individual members per decoy set
			HashMap<String, TreeMap<String,DecoyScoreSet>> decoy2scoresetspermethod = new HashMap<String, TreeMap<String,DecoyScoreSet>>();
			for (File scoresDir:scoresDirsBestCutoffs) {
				File[] files = scoresDir.listFiles(new RegexFileFilter(decoySet+"_.*\\.scores"));
				for (File file:files) {
					try {
						DecoyScoreSet set = new DecoyScoreSet(file);
						if (!decoy2scoresetspermethod.containsKey(set.getDecoyName())) {
							TreeMap<String,DecoyScoreSet> map = new TreeMap<String, DecoyScoreSet>();
							decoy2scoresetspermethod.put(set.getDecoyName(),map);
						}
						decoy2scoresetspermethod.get(set.getDecoyName()).put(set.getScoringMethod().getId(), set);
					} catch (FileFormatError e) {
						System.err.println(e.getMessage());
						continue;
					}
				}
				
			}

			PrintWriter i = new PrintWriter(new File(outDir,decoySet+".indiv.summary"));
			i.printf("#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
					"decoy","rank1","z-score","corr","rank1","z-score","corr","rank1","z-score","corr","rank1","z-score","corr","rank1","z-score","corr","rank1","z-score","corr");
			int recordCount = 0;
			for (String decoy:decoy2scoresetspermethod.keySet()) {
				
				Set<String> scoreTypes = decoy2scoresetspermethod.get(decoy).keySet();
				if (recordCount==0) {
					i.print("#\t");
					for (String scoreType:scoreTypes) {
						i.printf("%s\t",scoreType);
					}
					i.println();
				}
				i.print(decoy);
				for (String scoreType:scoreTypes) {
					DecoyScoreSet set = decoy2scoresetspermethod.get(decoy).get(scoreType);
					i.printf("\t%s\t%6.2f\t%6.2f",set.isNativeRank1(),set.getNativeZscore(),set.getSpearman());
				}
				i.println();
				recordCount++;	
			}
			i.close();
			
		}
 
		
		
	
		

	}



}

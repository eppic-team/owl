import java.io.*;
import java.util.LinkedHashMap;
import java.util.logging.Handler;
import java.util.logging.Logger;

import proteinstructure.PairwiseSequenceAlignment;
import proteinstructure.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;


/**
 * Executable to check the pairwise sequence identity for a list of sequences (multiple fasta file).
 * @author stehr
 *
 */
public class pairwiseSeqId {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		if(args.length < 1) {
			System.out.println("Usage: pairwiseSeqId <fasta_file>");
			System.exit(1);
		}
		
		String fastaFileName = args[0];
		File fastaFile = new File(fastaFileName);
		
		if(!fastaFile.canRead()) {
			System.err.println("Error. Can not read file " + fastaFileName);
			System.exit(1);
		}
		
		// parse sequences
		LinkedHashMap<String, String> seqs = new LinkedHashMap<String,String>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(fastaFile));
			String line;
			String currentTag = null;
			String currentSeq = "";
			while((line=in.readLine()) != null) {
				if(line.startsWith(">")) {
					if(currentTag != null) {
						if(currentSeq.length() == 0) {
							System.err.println("Error. No sequence found for >" + currentTag);
							System.exit(1);
						}
						if(seqs.containsKey(currentTag)) {
							System.err.println("Error. Duplicate tag >" + currentTag);
							System.exit(1);
						}
						seqs.put(currentTag, currentSeq);
						currentSeq = "";
					}
					currentTag = line.trim().substring(1);	// tag is everything except the ">"
				} else {
					currentSeq += line.trim();
				}
			}
			if(currentSeq.length() == 0) {
				System.err.println("Error. No sequence found for >" + currentTag);
				System.exit(1);
			}
			if(seqs.containsKey(currentTag)) {
				System.err.println("Error. Duplicate tag >" + currentTag);
				System.exit(1);
			}
			seqs.put(currentTag, currentSeq);
			
		} catch (IOException e) {
			System.err.println("Error reading from Fasta file " + fastaFileName + ": " + e.getMessage());
			System.exit(1);
		}
		
		// remove console logger from root logger
		Logger rootLogger = Logger.getLogger("");
		Handler[] handlers = rootLogger.getHandlers();
		rootLogger.removeHandler(handlers[0]);
		
		// align all
		float sumSeqId = 0;
		float maxSeqId = 0;
		float minSeqId = 100;
		int numSeqId = 0;
		int x = 0;
		int y = 0;
		int n = seqs.size();
		float[][] ids = new float[n][n];
		for(String t1:seqs.keySet()) {
			for(String t2:seqs.keySet()) {
				if(y > x) {
					String s1 = seqs.get(t1);
					String s2 = seqs.get(t2);
					try {
						PairwiseSequenceAlignment al = new PairwiseSequenceAlignment(s1,s2,t1,t2);
						float id = al.getPercentIdentity();
						sumSeqId += id;
						maxSeqId = Math.max(maxSeqId, id);
						minSeqId = Math.min(minSeqId, id);
						ids[x][y] = id;
						ids[y][x] = id;
						numSeqId++;
					} catch (PairwiseSequenceAlignmentException e) {
						System.err.println("Error aligning sequences " + t1 + " and " + t2 + ": " + e.getMessage());
					}
				}
				y++;
			}
			x++;
			y = 0;
		}
		
		System.out.println("Pairwise percent identities:");
		if(seqs.size() <= 55) {
		System.out.print("  seq\t num\t");
		for(x = 0; x < n; x++) {
			System.out.printf("%3d", x);
		}
		System.out.println("\tavg\tmax\n");
		x = 0;
		for(String t1:seqs.keySet()) {
			System.out.printf("%5s\t%4d\t", t1, x);
			float rowAvg = 0;
			float rowMax = 0;
			float rowMin = 100;
			for(y = 0; y < n; y++) {
				if(y!=x) {
					System.out.printf("%3.0f",ids[x][y]);
					rowAvg+=ids[x][y];
					rowMax = Math.max(rowMax, ids[x][y]);
					rowMin = Math.min(rowMin, ids[x][y]);
				} else {
					System.out.printf("   ");
				}
			}
			System.out.printf("\t%3.0f\t%3.0f\n", rowAvg / (y-1), rowMax);
			x++;
		}
		} else {
			System.out.println("Output supressed for more than 55 sequences");
		}
		System.out.println();
		System.out.println("Number of sequences: " + seqs.size());
		System.out.println("Number of comparisons: " + numSeqId);
		System.out.println("Maximum percent identity: " + maxSeqId);
		System.out.println("Minimum percent identity: " + minSeqId);
		System.out.println("Average percent identity: " + (sumSeqId / numSeqId));

	}

}

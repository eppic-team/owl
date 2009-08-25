package proteinstructure;

import java.io.File;

/**
 * Class to store a decoy's score together with its file name and rmsd to native
 * Two objects of this class are comparable based on their scores.
 * @author duarte
 *
 */
public class DecoyScore implements Comparable<DecoyScore> {
	
	public File file;
	public double score;
	public double rmsd;
	
	public DecoyScore(File file, double score, double rmsd) {
		this.file = file;
		this.score = score;
		this.rmsd = rmsd;
	}
	
	public int compareTo(DecoyScore o) {
		return Double.compare(this.score, o.score);
	}
}

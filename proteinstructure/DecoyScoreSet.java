package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import tools.Statistics;

/**
 * A class to represent the scoring of a single decoy set (a native structure plus 
 * a number of decoys of that structure)
 * @author duarte
 *
 */
public class DecoyScoreSet implements Iterable<DecoyScore> {

	private String decoyName;
	private HashMap<String,DecoyScore> set; // decoy file name (no path) to DecoyScore

	/**
	 * Constructs an empty DecoyScoreSet 
	 */
	public DecoyScoreSet(String decoyName) {
		this.decoyName = decoyName;
		set = new HashMap<String, DecoyScore>();
	}
	
	/**
	 * Returns the decoy name of this set. The decoy name corresponds usually to the pdb 
	 * code of the native. 
	 * @return
	 */
	public String getDecoyName() {
		return decoyName;
	}
	
	/**
	 * Adds the given decoyScore to this set
	 * @param decoyScore
	 */
	public void addDecoyScore(DecoyScore decoyScore) {
		set.put(decoyScore.file.getName(),decoyScore);
	}
	
	/**
	 * Gets the DecoyScore (i.e. file, score, rmsd) for the given file name
	 * @param fileName the file name without the path
	 * @return
	 */
	public DecoyScore getDecoyScore(String fileName) {
		return set.get(fileName);
	}
	
	/**
	 * Tells whether the given file name is contained within this DecoyScoreSet
	 * @param fileName the file name without the path
	 * @return
	 */
	public boolean containsDecoyScore(String fileName) {
		return set.containsKey(fileName);
	}
	
	/**
	 * Returns the number of decoys in this decoy set
	 * @return
	 */
	public int size() {
		return set.size();
	}

	/**
	 * Gets the sorted list of decoys (based on score) as a new ArrayList 
	 * @return
	 */
	private ArrayList<DecoyScore> getSortedList() {
		ArrayList<DecoyScore> allScores = new ArrayList<DecoyScore>(set.values());
		Collections.sort(allScores);
		return allScores;
	}
	
	/**
	 * Writes to file the decoy set scores in three columns: file name, score and rmsd
	 * @param file
	 * @throws FileNotFoundException
	 */
	public void writeToFile(File file) throws FileNotFoundException {
		PrintWriter pw = new PrintWriter(file);
		ArrayList<DecoyScore> allScores = this.getSortedList();
		pw.printf("#%s\t%s\t%s\n","file","score","rmsd");
		for (DecoyScore fs:allScores) {
			pw.printf("%s\t%7.2f\t%6.3f\n",fs.file.getName(), fs.score, fs.rmsd);
		}
		pw.close();

	}
	
	/**
	 * Returns the z-score of the given structure (as a file name) with respect of the 
	 * whole set, i.e. how far the score of the given structure is from the mean of this set
	 * measured in units of standard deviation.  
	 * @param decoyFileName the file name without the path
	 * @return
	 * @throws NullPointerException if given file name not present in this set. Check it with 
	 * {@link #containsDecoyScore(String)}
	 */
	public double getZscore(String decoyFileName) {
		double nativeScore = this.getDecoyScore(decoyFileName).score;
		ArrayList<Double> scores = new ArrayList<Double>();
		for (DecoyScore fs:set.values()) {
			scores.add(fs.score);
		}
		return (nativeScore-Statistics.mean(scores))/Statistics.ssd(scores);
	}
	
	/**
	 * Tells whether given structure (as a file name) is ranked 1 (maximum score) in this set.
	 * @param decoyFileName the file name without the path
	 * @return true if given decoy file name is ranked first, false if not ranked first or given 
	 * file name not present in set
	 */
	public boolean isRank1(String decoyFileName) {
		ArrayList<DecoyScore> sorted = getSortedList();
		if (sorted.get(sorted.size()-1).file.getName().equals(decoyFileName)) {
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * Returns the spearman correlation of the scores vs the rmsds of this set.
	 * @return
	 */
	public double getSpearman() {
		double[] scores = new double[this.size()];
		double[] rmsds = new double[this.size()];
		int i = 0;
		for (DecoyScore ds:this){
			scores[i]=ds.score;
			rmsds[i]=ds.rmsd;
			i++;
		}
		return Statistics.spearman(scores, rmsds);
	}

	public Iterator<DecoyScore> iterator() {
		return set.values().iterator();
	}
	
	/**
	 * Read rmsds of a decoy set from a rmsd file in 'decoys r us' format.
	 * @param file
	 * @return map of file names to rmsd values
	 * @throws IOException
	 */
	public static HashMap<String,Double> readRmsds(File file) throws IOException {
		HashMap<String,Double> rmsds = new HashMap<String, Double>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		while ((line=br.readLine())!=null) {
			Pattern p = Pattern.compile("^cRMSD.*and\\s([^\\s]+\\.pdb)\\sis\\s+([\\d\\.]+)");
			Matcher m = p.matcher(line);
			if (m.matches()) {
				String filename = m.group(1);
				double rmsd = Double.parseDouble(m.group(2));
				rmsds.put(filename, rmsd);
			}
		}
		br.close();
		return rmsds;
	}

}

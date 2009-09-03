package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.Scorer.ScoringMethod;

import tools.MySQLConnection;
import tools.Statistics;

/**
 * A class to represent the scoring of a single decoy set (a native structure plus 
 * a number of decoys of that structure)
 * @author duarte
 *
 */
public class DecoyScoreSet implements Iterable<DecoyScore> {

	private String decoyName;
	private String decoyGroup;
	private String nativeFileName;
	
	private ScoringMethod scoringMethod;
	private String ct;
	private double cutoff;
	private int minSeqSep;
	private File trainingSetFile;
	private int numStructsTrainingSet;
	
	private HashMap<String,DecoyScore> set; // decoy file name (no path) to DecoyScore

	/**
	 * Constructs an empty DecoyScoreSet 
	 * @param decoyName
	 * @param scorer
	 */
	public DecoyScoreSet(String decoyName, String decoyGroup, Scorer scorer) {
		this.decoyName = decoyName;
		this.decoyGroup = decoyGroup;
		this.nativeFileName = decoyName+".pdb";
		set = new HashMap<String, DecoyScore>();
		
		this.scoringMethod = scorer.getScoringMethod();
		this.ct = scorer.getContactType();
		this.cutoff = scorer.getCutoff();
		this.minSeqSep = scorer.getMinSeqSep();
		this.trainingSetFile = scorer.getListFile();
		this.numStructsTrainingSet = scorer.sizeOfTrainingSet();
	}
	
	/**
	 * Constructs a DecoyScoreSet by reading the decoy scores from given file
	 * @param decoySetScoreFile
	 * @throws IOException
	 * @throws FileFormatError if native structure can't be found in file
	 */
	public DecoyScoreSet(File decoySetScoreFile) throws IOException, FileFormatError {
		set = new HashMap<String, DecoyScore>();
		readFromFile(decoySetScoreFile);
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
	 * Returns the decoy group to which this set belongs.
	 * @return
	 */
	public String getDecoyGroup() {
		return decoyGroup;
	}

	/**
	 * Returns the scoring method used to score this set.
	 * @return
	 */
	public ScoringMethod getScoringMethod() {
		return scoringMethod;
	}
	
	/**
	 * Returns the contact type upon which the scoring of
	 * this set was based (only vaid in residue-based scoring)
	 * @return
	 */
	public String getContactType() {
		return ct;
	}

	/**
	 * Returns the cutoff upon which the scoring of this set was based.
	 * @return
	 */
	public double getCutoff() {
		return cutoff;
	}

	/**
	 * Returns the minimum sequence separation used for filtering when compiling
	 * the scoring matrix and scoring.
	 * @return
	 */
	public int getMinSeqSep() {
		return minSeqSep;
	}

	/**
	 * Returns the name of the native PDB file for this set. The DecoyScore for it is not
	 * necessarily a member of this set (e.g. because it was missing in the decoy data).
	 * Check whether the native file is member of this set with {@link #containsNative()} 
	 * @return
	 */
	public String getNativeFileName() {
		return nativeFileName;
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
	 * Tells whether the native structure is contained within this DecoyScoreSet 
	 * @return true if native is in this set, false otherwise
	 */
	public boolean containsNative() {
		return containsDecoyScore(nativeFileName);
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
		
		pw.println("# SCORE METHOD: "+this.scoringMethod.getDescription());
		pw.println("# contact type: "+this.ct);
		pw.println("# cutoff: "+this.cutoff);
		pw.println("# min sequence separation: "+this.minSeqSep);
		pw.println("# structures: "+this.numStructsTrainingSet);
		pw.println("# list: "+this.trainingSetFile.toString());

		pw.printf("#%s\t%s\t%s\n","file","score","rmsd");
		for (DecoyScore fs:allScores) {
			pw.printf("%s\t%7.2f\t%6.3f\n",fs.file.getName(), fs.score, fs.rmsd);
		}
		pw.close();

	}
	
	/**
	 * Reads a decoy set scores file in our format. See {@link #writeToFile(File)}
	 * @param file
	 * @throws IOException
	 * @throws FileFormatError if native structure can't be found in file
	 */
	public void readFromFile(File file) throws IOException, FileFormatError {
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line;
		Pattern p = Pattern.compile("^(\\d\\w\\w\\w(-\\w)?)\\.pdb$");
		boolean nativeFound = false;
		while ((line=br.readLine())!=null) {
			if (line.startsWith("#")) {
				Pattern p1 = Pattern.compile("^# SCORE METHOD: (.*)$");
				Matcher m1 = p1.matcher(line);
				if (m1.matches()) {
					this.scoringMethod = ScoringMethod.getByDescription(m1.group(1));
				}				
				p1 = Pattern.compile("^# contact type: (.*)$");
				m1 = p1.matcher(line);
				if (m1.matches()) {
					this.ct = m1.group(1);
				}
				p1 = Pattern.compile("^# cutoff: (.*)$");
				m1 = p1.matcher(line);
				if (m1.matches()) {
					this.cutoff = Double.parseDouble(m1.group(1));
				}
				p1 = Pattern.compile("^# min sequence separation: (.*)$");
				m1 = p1.matcher(line);
				if (m1.matches()) {
					this.minSeqSep = Integer.parseInt(m1.group(1));
				}			
				p1 = Pattern.compile("^# structures: (.*)$");
				m1 = p1.matcher(line);
				if (m1.matches()) {
					this.numStructsTrainingSet = Integer.parseInt(m1.group(1));
				}
				p1 = Pattern.compile("^# list: (.*)$");
				m1 = p1.matcher(line);
				if (m1.matches()) {
					this.trainingSetFile = new File(m1.group(1));
				}			

			} else {
				String[] cols = line.split("\\s+");
				addDecoyScore(new DecoyScore(new File(cols[0]),Double.parseDouble(cols[1]),Double.parseDouble(cols[2])));
				Matcher m = p.matcher(cols[0]);
				if (m.matches()) {
					this.decoyName = m.group(1);
					this.nativeFileName = m.group(1)+".pdb";
					nativeFound = true;
				} else if (cols[2].equals("0.000")) {
					this.nativeFileName = cols[0];
					this.decoyName = cols[0].substring(0,cols[0].lastIndexOf(".pdb"));
					nativeFound = true;
				}
			}
		}
		br.close();
		if (!nativeFound)
			throw new FileFormatError("Couldn't find native decoy in decoy set file "+file);

	}
	
	/**
	 * Writes to db all data for this decoy score set
	 * @param conn
	 * @param db
	 * @param table
	 * @throws SQLException
	 */
	public void writeToDb(MySQLConnection conn, String db, String table) throws SQLException {
		Statement st = conn.createStatement();
		
		for (DecoyScore decoyScore:this) {
			String sql = "INSERT INTO "+db+"."+table+
			" (method, ct, cutoff, min_seq_sep, train_set_size, train_set_file, decoy_group, decoy_set_name, file, score, rmsd) " +
			" VALUES (" +
			"'"+scoringMethod.getId()+"', " +
			"'"+ct+"', "+cutoff+", "+minSeqSep+", "+numStructsTrainingSet+", '"+trainingSetFile.getName()+"', " +
			"'"+decoyGroup+"', '"+decoyName+"', '"+decoyScore.file.getName()+"', "
			+decoyScore.score+", "+decoyScore.rmsd+")";
			st.executeUpdate(sql);	
		}
		st.close();
		
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
	 * Returns the z-score of the native with respect to the whole set, i.e. how far the score 
	 * of the native is from the mean of this set measured in units of standard deviation.
	 * @return
	 * @throws NullPointerException if native not present in this set. Check with {@link #containsNative()}
	 */
	public double getNativeZscore() {
		return getZscore(nativeFileName);
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
	 * Tells whether native structure is ranked 1 (maximum score) in this set.
	 * @return true if native is ranked first, false if not ranked first or native not
	 * present in set
	 */
	public boolean isNativeRank1() {
		return isRank1(nativeFileName);
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

	/**
	 * Returns the minimum and maximum values of rmsds in this set as an array of size 2
	 * (excluding the native)
	 * @return array of size 2 with {min rmsd, max rmsd}
	 */
	public double[] getRmsdRange() {
		ArrayList<Double> rmsds = new ArrayList<Double>();
		for (DecoyScore ds:this){
			if (!ds.file.getName().equals(nativeFileName)) {
				rmsds.add(ds.rmsd);
			}
		}
		double[] minmax = {Collections.min(rmsds),Collections.max(rmsds)}; 
		return minmax;
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
	public static void writeStats(File file, ArrayList<DecoyScoreSet> resStats, ArrayList<DecoyScoreSet> atomStats) {
		try {
			PrintWriter pw = new PrintWriter(file);
			pw.printf("#%9s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\n",
					"decoy","nod","res_r1","resz","rescor","atom_r1","atomz","atomcor");
			double sumresz = 0,sumrescor = 0, sumatomz = 0, sumatomcor = 0;
			int countresr1 = 0, countatomr1 = 0;
			for (int i=0;i<resStats.size();i++) {
				String decoy = resStats.get(i).getDecoyName();
				int resNumScDec = resStats.get(i).size();
				double resz = resStats.get(i).getNativeZscore();
				double rescor = resStats.get(i).getSpearman();
				boolean resr1 = resStats.get(i).isNativeRank1();
				int atomNumScDec = atomStats.get(i).size();
				double atomz = atomStats.get(i).getNativeZscore();
				double atomcor = atomStats.get(i).getSpearman();
				boolean atomr1 = atomStats.get(i).isNativeRank1();
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

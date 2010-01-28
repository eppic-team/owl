package proteinstructure.DecoyScoring;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.FileFormatError;


/**
 * Class to score a PDB structure based on type frequencies (residue or atom types).
 * It contains code for compiling a scoring matrix from a set of non-redundant structures,
 * write/read the scoring matrix from/to file and scoring a given structure based on the 
 * scoring matrix.
 * 
 * See implementing subclasses {@link ResTypeScorer} and {@link AtomTypeScorer}
 * 
 * @author duarte
 *
 */
public abstract class TypeScorer extends Scorer {

	private static final int TOO_FEW_COUNTS_THRESHOLD = 20;
	private static final double TOO_FEW_COUNTS_SCORE = 0.0;
	
	protected int[] entityCounts; 		// the counts of the entities, i.e. of the residue or atom types
	protected int[][] pairCounts;		// the counts of the pairs, using always defined indices, see types2indices and indices2types maps
	protected double[][] scoringMat;	// the scoring matrix
	
	private int totalEntityCount;		// the total number of entities, i.e. of residues or atoms
	private long totalPairsCount; 		// the total number of pairs, to be on the safe side we make it 
										// a long (to avoid the dreaded overflow!) but it's unlikely to happen (>2E09 pairs)
	
	protected HashMap<String,Integer> types2indices;	// map of types to indices of the above arrays
	protected HashMap<Integer,String> indices2types;	// map of indices of the above arrays to types
	
	protected int numEntities;							// the number of entities, i.e. atom types (167) or residue types (20)
		
	protected TypeScorer() {

	}
	
	/**
	 * Performs the counts of the type pairs and stores them in the internal arrays.
	 * Use subsequently {@link #calcScoringMat()} to compute the scoring matrix from counts arrays.
	 * @throws SQLException if database server can't be accessed to get PDB data
	 * @throws IOException if list file can't be read  
	 */
	public abstract void countPairs() throws SQLException, IOException;
	
	private void countTotals() {
		totalEntityCount = 0;
		totalPairsCount = 0;
		for (int i=0;i<numEntities;i++) {
			totalEntityCount+=entityCounts[i];
			for (int j=0;j<numEntities;j++) {
				totalPairsCount+=pairCounts[i][j];
			}
		}		
	}
	
	protected void countPair(int i, int j) {
		if (j>i){
			pairCounts[i][j]++;
		} else {
			pairCounts[j][i]++;
		}
	}
	
	protected void countEntity(int i) {
		entityCounts[i]++;
	}
	
	/**
	 * Computes the scoring matrix from the counts arrays. The score is log2(p_obs/p_exp)
	 * Whenever the pair count of a certain pair is below the threshold ({@value #TOO_FEW_COUNTS_THRESHOLD}
	 * the score assigned for that pair is {@value #TOO_FEW_COUNTS_SCORE}
	 */
	public void calcScoringMat() {
		countTotals();
		scoringMat = new double[numEntities][numEntities];
		double logof2 = Math.log(2);
		for(int i=0;i<scoringMat.length;i++) {
			for(int j=i;j<scoringMat[i].length;j++) {
				if (pairCounts[i][j]<TOO_FEW_COUNTS_THRESHOLD) { // When counts are too small, we can't get significant statistics for them
					scoringMat[i][j] = TOO_FEW_COUNTS_SCORE;     // We use a score of 0 in this case, i.e. no information
				} else {
					scoringMat[i][j] = 
						Math.log(
								((double)pairCounts[i][j]/(double)totalPairsCount)/
								(((double)entityCounts[i]*(double)entityCounts[j])/((double)totalEntityCount*(double)totalEntityCount))
						)/logof2; // we explicitely cast every int/long to double to avoid overflows, java doesn't check for them!
				}
				
			}
		}
	}
	
	public int getEntityCount(int i) {
		return entityCounts[i];
	}
	
	public double getEntityFrequency(int i) {
		return (double)entityCounts[i]/(double)totalEntityCount;
	}
	
	public String getTypeFromIndex(int i) {
		return indices2types.get(i);
	}
	
	public int getIndexFromType(String type) {
		return types2indices.get(type);
	}
	
	/**
	 * Writes the scoring matrix to a file
	 * @param file the file to write to
	 * @param writeCounts whether to also write the matrix with the pair counts
	 * @throws FileNotFoundException if file not found
	 */
	protected void writeScMatToFile(File file, boolean writeCounts) throws FileNotFoundException {
		PrintWriter pw = new PrintWriter(file);

		pw.println("# SCORE METHOD: "+getScoringMethod().getDescription());
		pw.println("# contact type: "+ct);
		pw.println("# cutoff: "+cutoff);
		pw.println("# min sequence separation: "+minSeqSep);
		pw.println("# structures: "+sizeOfTrainingSet());
		pw.println("# list: "+getListFile().toString());
		pw.println("# nodes: "+totalEntityCount);
		pw.println("# pairs: "+totalPairsCount);
		pw.print("# type counts: ");
		for (int i=0;i<entityCounts.length;i++) {
			pw.print(entityCounts[i]+" ");
		}
		pw.println();
		
		if (writeCounts) {
			pw.printf("%7s","");
			for(int j=0;j<pairCounts.length;j++) {
				pw.printf("%7s",indices2types.get(j));
			}
			pw.println();

			for(int i=0;i<pairCounts.length;i++) {
				pw.printf("%7s",indices2types.get(i));
				for(int j=0;j<pairCounts[i].length;j++) {
					if (j<i) {
						pw.printf("%7s","");
					} else {
						pw.printf("%7d",pairCounts[i][j]);
					}
				}
				pw.println();
			}

			pw.println();
		}
		
		pw.printf("%7s","");
		for(int j=0;j<scoringMat.length;j++) {
			pw.printf("%7s",indices2types.get(j));
		}
		pw.println();
	
		for(int i=0;i<scoringMat.length;i++) {
			pw.printf("%7s",indices2types.get(i));
			for(int j=0;j<scoringMat[i].length;j++) {
				if (j<i) {
					pw.printf("%7s","");
				} else {
					pw.printf("%7.2f",scoringMat[i][j]);
				}
			}
			pw.println();
		}
		pw.close();
	}
	
	/**
	 * Reads the scoring matrix from a file.
	 * @param scMatFile
	 * @throws IOException
	 * @throws FileFormatError
	 */
	protected void readScMatFromFile(File scMatFile) throws IOException, FileFormatError {
		ArrayList<Integer> entityCountsAL = new ArrayList<Integer>();
		ArrayList<ArrayList<Double>> scoringMatAL = new ArrayList<ArrayList<Double>>();
		this.indices2types = new HashMap<Integer, String>();
		
		int rowCount = 0;
		BufferedReader br = new BufferedReader(new FileReader(scMatFile));
		String line;
		while ((line=br.readLine())!=null) {

			Pattern p = Pattern.compile("^# contact type: (.*)$");
			Matcher m = p.matcher(line);
			if (m.matches()) {
				this.ct = m.group(1);
			}
			p = Pattern.compile("^# cutoff: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.cutoff = Double.parseDouble(m.group(1));
			}
			p = Pattern.compile("^# min sequence separation: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.minSeqSep = Integer.parseInt(m.group(1));
			}			
			p = Pattern.compile("^# structures: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.totalStructures = Integer.parseInt(m.group(1));
			}
			p = Pattern.compile("^# list: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.listFile = new File(m.group(1));
			}			
			p = Pattern.compile("^# nodes: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.totalEntityCount = Integer.parseInt(m.group(1));
			}
			p = Pattern.compile("^# pairs: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.totalPairsCount = Long.parseLong(m.group(1));
			}
			p = Pattern.compile("^# type counts: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				String[] tokens = m.group(1).split("\\s+");
				for (int i=0;i<tokens.length;i++){
					entityCountsAL.add(Integer.parseInt(tokens[i]));
				}
			}
			if (!line.startsWith("#")) {
				// column names line
				p = Pattern.compile("^\\s*([\\w\\s]+)$");
				m = p.matcher(line);
				if (m.matches()) {
					String[] tokens = m.group(1).split("\\s+");
					for (int i=0;i<tokens.length;i++){
						indices2types.put(i, tokens[i]);
					}					
				}
				// matrix lines
				p = Pattern.compile("^\\s*\\w+\\s+([0-9.\\-\\s]+)$");
				m = p.matcher(line);
				if (m.matches()) {
					
					ArrayList<Double> thisRow = new ArrayList<Double>();
					scoringMatAL.add(thisRow);
					
					String[] tokens = m.group(1).split("\\s+");
					for (int i=0;i<tokens.length;i++){
						thisRow.add(Double.parseDouble(tokens[i]));
					}
					rowCount++;
				} 
			}
		}
		br.close();
		
		this.numEntities = entityCountsAL.size();
		this.entityCounts = new int[entityCountsAL.size()];
		for (int i=0;i<entityCountsAL.size();i++) {
			this.entityCounts[i] = entityCountsAL.get(i);
		}
		if (scoringMatAL.size()!=scoringMatAL.get(0).size()) {
			throw new FileFormatError("Size of matrix rows and columns don't match: "+scoringMatAL.size()+" "+scoringMatAL.get(0).size());
		}
		this.scoringMat = new double[scoringMatAL.size()][scoringMatAL.get(0).size()];
		for (int i=0;i<scoringMatAL.size();i++) {
			for (int j=0;j<scoringMatAL.get(i).size();j++) {
				this.scoringMat[i][j+i] = scoringMatAL.get(i).get(j);
			}
		}
		this.types2indices = new HashMap<String, Integer>();
		for (int i:indices2types.keySet()) {
			types2indices.put(indices2types.get(i),i);
		}
	}
	
}

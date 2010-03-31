package proteinstructure.decoyScoring;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import proteinstructure.FileFormatError;


/**
 * Class to score a PDB structure based on triplet type frequencies (residue or atom types).
 * It contains code for compiling a scoring matrix from a set of non-redundant structures,
 * write/read the scoring matrix from/to file and scoring a given structure based on the 
 * scoring matrix.
 * 
 * 
 * @author duarte
 *
 */
public abstract class TripletScorer extends Scorer {

	private static final int TOO_FEW_COUNTS_THRESHOLD = 10;
	private static final double TOO_FEW_COUNTS_SCORE = 0.0;
	
	protected int[] entityCounts; 		// the counts of the entities, i.e. of the residue or atom types
	protected int[][][] tripletCounts;	// the counts of the triplets, using always defined indices, see types2indices and indices2types maps
	protected double[][][] scoringMat;	// the scoring matrix
	
	private int totalEntityCount;		// the total number of entities, i.e. of residues or atoms
	private long totalTripletCount; 	// the total number of triplets, to be on the safe side we make it 
										// a long (to avoid the dreaded overflow!) but it's unlikely to happen (>2E09 triplets)
	
	protected HashMap<String,Integer> types2indices;	// map of types to indices of the above arrays
	protected HashMap<Integer,String> indices2types;	// map of indices of the above arrays to types
	
	protected int numEntities;							// the number of entities, i.e. atom types (167) or residue types (20)
		
	protected TripletScorer() {

	}
	
	/**
	 * Performs the counts of the type triplets and stores them in the internal arrays.
	 * Use subsequently {@link #calcScoringMat()} to compute the scoring matrix from counts arrays.
	 * @throws SQLException if database server can't be accessed to get PDB data
	 * @throws IOException if list file can't be read  
	 */
	public abstract void countTriplets() throws SQLException, IOException;
	
	private void countTotals() {
		totalEntityCount = 0;
		totalTripletCount = 0;
		for (int i=0;i<numEntities;i++) {
			totalEntityCount+=entityCounts[i];
			for (int j=i;j<numEntities;j++) {
				for (int k=j;k<numEntities;k++) {
					totalTripletCount+=tripletCounts[i][j][k];
				}
			}
		}		
	}
	
	protected void countTriplet(int i, int j, int k) {
		int[] ind = {i,j,k};
		Arrays.sort(ind);
		tripletCounts[ind[0]][ind[1]][ind[2]]++;
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
		scoringMat = new double[numEntities][numEntities][numEntities];
		double logof2 = Math.log(2);
		for(int i=0;i<numEntities;i++) {
			for(int j=i;j<numEntities;j++) {
				for(int k=j;k<numEntities;k++) {
					if (tripletCounts[i][j][k]<TOO_FEW_COUNTS_THRESHOLD) { // When counts are too small, we can't get significant statistics for them
						scoringMat[i][j][k] = TOO_FEW_COUNTS_SCORE;     // We use a score of 0 in this case, i.e. no information
					} else {
						scoringMat[i][j][k] = 
							Math.log(
									((double)tripletCounts[i][j][k]/
											(double)totalTripletCount)/
									(((double)entityCounts[i]*(double)entityCounts[j]*(double)entityCounts[k])/
											((double)totalEntityCount*(double)totalEntityCount*(double)totalEntityCount))
							)/logof2; // we explicitely cast every int/long to double to avoid overflows, java doesn't check for them!
					}
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
		pw.println("# triplets: "+totalTripletCount);
		pw.print("# type counts: ");
		for (int i=0;i<entityCounts.length;i++) {
			pw.print(entityCounts[i]+" ");
		}
		pw.println();
		
		if (writeCounts) {
			for(int i=0;i<numEntities;i++) {
				pw.println(i+":"+indices2types.get(i));
				printTypeHeaders(pw);
				for(int j=0;j<numEntities;j++) {
					pw.printf("%7s",indices2types.get(j));
					for(int k=0;k<numEntities;k++) {
						if (j<i || k<i || k<j) {
							pw.printf("%7s","");
						} else {
							pw.printf("%7d",tripletCounts[i][j][k]);
						}
					}
					pw.println();
				}
				pw.println();
				pw.println();
			}
			pw.println();
			pw.println();
			pw.println();
		}
		
		for(int i=0;i<numEntities;i++) {
			pw.println(i+":"+indices2types.get(i));
			printTypeHeaders(pw);
			for(int j=0;j<numEntities;j++) {
				pw.printf("%7s",indices2types.get(j));
				for(int k=0;k<numEntities;k++) {
					if (j<i || k<i || k<j) {
						pw.printf("%7s","");
					} else {
						pw.printf("%7.2f",scoringMat[i][j][k]);
					}
				}
				pw.println();
			}
			pw.println();
			pw.println();
		}
		pw.println();
		pw.println();
		pw.println();
		pw.close();
	}
	
	private void printTypeHeaders(PrintWriter pw) {
		pw.printf("%7s","");
		for(int j=0;j<numEntities;j++) {
			pw.printf("%7s",indices2types.get(j));
		}
		pw.println();
	}
	
	/**
	 * Reads the scoring matrix from a file.
	 * @param scMatFile
	 * @throws IOException
	 * @throws FileFormatError
	 */
	protected void readScMatFromFile(File scMatFile) throws IOException, FileFormatError {
		ArrayList<Integer> entityCountsAL = new ArrayList<Integer>();
		ArrayList<ArrayList<ArrayList<Double>>> scoringMatAL = new ArrayList<ArrayList<ArrayList<Double>>>();
		this.indices2types = new HashMap<Integer, String>();
		

		int iInd = 0;
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
			p = Pattern.compile("^# triplets: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				this.totalTripletCount = Long.parseLong(m.group(1));
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
				// headers
				p = Pattern.compile("^\\s*([\\w\\s]+)$");
				m = p.matcher(line);
				// the headers appear several times, we just want to fill the map the first time
				if (indices2types.isEmpty() && m.matches()) {
					String[] tokens = m.group(1).split("\\s+");
					for (int i=0;i<tokens.length;i++){
						indices2types.put(i, tokens[i]);
					}
				}

				// dimension i
				p = Pattern.compile("^(\\d+):\\w+$");
				m = p.matcher(line);
				if (m.matches()){
					iInd = Integer.parseInt(m.group(1));
					ArrayList<ArrayList<Double>> iMat = new ArrayList<ArrayList<Double>>();
					scoringMatAL.add(iMat);
				}
				
				// dimension j (matrix slices from the 3d matrix)
				p = Pattern.compile("^\\s*\\w+\\s+((?:[0-9.\\-]+\\s*)+)$");
				m = p.matcher(line);
				if (m.matches()) {
					
					ArrayList<Double> thisRow = new ArrayList<Double>();
					scoringMatAL.get(iInd).add(thisRow);
					
					String[] tokens = m.group(1).split("\\s+");
					for (int i=0;i<tokens.length;i++){
						thisRow.add(Double.parseDouble(tokens[i]));
					}
				} 
			}
		}
		br.close();
		
		this.numEntities = entityCountsAL.size();

		this.entityCounts = new int[entityCountsAL.size()];
		for (int i=0;i<entityCountsAL.size();i++) {
			this.entityCounts[i] = entityCountsAL.get(i);
		}
		if (numEntities!=scoringMatAL.get(0).size() || numEntities!=scoringMatAL.get(0).get(0).size()) {
			throw new FileFormatError("Dimensions of 3d matrix don't match: I "+numEntities+", J "+scoringMatAL.get(0).size()+", K "+scoringMatAL.get(0).get(0).size());
		}
		this.scoringMat = new double[numEntities][numEntities][numEntities];
		for (int i=0;i<numEntities;i++) {
			for (int j=0;j<scoringMatAL.get(i).size();j++) {
				for (int k=0;k<scoringMatAL.get(i).get(j).size();k++) {
					this.scoringMat[i][j+i][k+j+i] = scoringMatAL.get(i).get(j).get(k);
				}
			}
		}
		this.types2indices = new HashMap<String, Integer>();
		for (int i:indices2types.keySet()) {
			types2indices.put(indices2types.get(i),i);
		}
	}
	
	
}

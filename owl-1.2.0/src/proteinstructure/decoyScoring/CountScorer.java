package proteinstructure.decoyScoring;

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
 * Class to score a PDB structure based on contact count frequencies (residue or atom contacts).
 * It contains code for compiling a scoring matrix from a set of non-redundant structures,
 * write/read the scoring matrix from/to file and scoring a given structure based on the 
 * scoring matrix.
 * 
 * See implementing subclasses {@link ResCountScorer} and {@link AtomCountScorer}
 * 
 * @author duarte
 *
 */
public abstract class CountScorer extends Scorer {

	private static final int TOO_FEW_COUNTS_THRESHOLD = 10;
	protected static final double TOO_FEW_COUNTS_SCORE = 0.0;
	
	protected static final int NUM_COUNT_BINS= 120; // for high cut-offs it really goes very high
	
	protected int[] binCounts; 				// the counts of the members of each neighbour-count bin. Size: numCountBins
	protected int[] typeCounts;				// the counts of types
	protected int[][] binCountsPerType;		// the counts of the types per neighbour-count bin (thus array is not square). Sizes: numCountbins, numEntities
	protected double[][] scoringMat;		// the scoring matrix (not square in the case of counts). Sizes: numCountBins, numEntities	
	protected int totalEntityCount;			// the total number of entities, i.e. of residues or atoms
	protected int numCountBins;				// the number of neighbor-count bins (initialises to a constant)
	protected int numTypes;					// the number of types, i.e. atom types (167) or residue types (20)

	protected HashMap<String,Integer> types2indices;	// map of types to indices of the above arrays
	protected HashMap<Integer,String> indices2types;	// map of indices of the above arrays to types

	protected CountScorer() {
		
	}
		
	/**
	 * Performs the counts of the number of neighbours and stores them in the internal arrays.
	 * Use subsequently {@link #calcScoringMat()} to compute the scoring matrix from counts arrays.
	 * @throws SQLException if database server can't be access to get PDB data
	 * @throws IOException if list file can't be read
	 */
	public abstract void countNodes() throws SQLException, IOException;
	
	private void countTotals() {
		totalEntityCount = 0;
		binCounts = new int[numCountBins];
		typeCounts = new int[numTypes];
		for (int binIdx=0;binIdx<numCountBins;binIdx++) {
			for (int typeIdx=0;typeIdx<numTypes;typeIdx++) {
				typeCounts[typeIdx]+=binCountsPerType[binIdx][typeIdx];
				binCounts[binIdx]+=binCountsPerType[binIdx][typeIdx];
				totalEntityCount+=binCountsPerType[binIdx][typeIdx];
			}
		}		
	}
	
	protected void count(int binIdx, int typeIdx) {
		binCountsPerType[binIdx][typeIdx]++;
	}
		
	/**
	 * Computes the scoring matrix from the counts arrays. The score is log2(p_obs/p_exp)
	 * Whenever the count of a certain type in a neighbour-count bin is below the threshold ({@value #TOO_FEW_COUNTS_THRESHOLD}
	 * the score assigned for that matrix cell is {@value #TOO_FEW_COUNTS_SCORE}
	 */
	public void calcScoringMat() {
		countTotals();
		scoringMat = new double[numCountBins][numTypes];
		double logof2 = Math.log(2);
		for(int binIdx=0;binIdx<numCountBins;binIdx++) {
			for(int typeIdx=0;typeIdx<numTypes;typeIdx++) {
				if (binCountsPerType[binIdx][typeIdx]<TOO_FEW_COUNTS_THRESHOLD) {
					scoringMat[binIdx][typeIdx] = TOO_FEW_COUNTS_SCORE;
				} else {
					scoringMat[binIdx][typeIdx] = 
						Math.log(
								((double)binCountsPerType[binIdx][typeIdx]/(double)binCounts[binIdx])/
								(((double)typeCounts[typeIdx])/((double)totalEntityCount))
						)/logof2; // we explicitely cast every int/long to double to avoid overflows, java doesn't check for them!
				}
				
			}
		}
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
		pw.print("# neighbor-count bins: ");
		for (int i=0;i<numCountBins;i++) {
			if (binCounts[i]!=0) {
				pw.print(i+":"+binCounts[i]+" "); 
			}
		}
		pw.println();
		
		if (writeCounts) {
			pw.printf("%4s","");
			for(int j=0;j<numTypes;j++) {
				pw.printf("%7s",indices2types.get(j));
			}
			pw.println();

			for(int i=0;i<numCountBins;i++) {
				if (binCounts[i]!=0) {
					pw.printf("%4d",i);
					for(int j=0;j<numTypes;j++) {
						pw.printf("%7d",binCountsPerType[i][j]);
					}
					pw.println();
				}
			}

			pw.println();
		}
		
		pw.printf("%4s","");
		for(int j=0;j<numTypes;j++) {
			pw.printf("%7s",indices2types.get(j));
		}
		pw.println();
	
		for(int i=0;i<numCountBins;i++) {
			if (binCounts[i]!=0) {
				pw.printf("%4d",i);
				for(int j=0;j<numTypes;j++) {
					pw.printf("%7.2f",scoringMat[i][j]);
				}
				pw.println();
			}
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
		ArrayList<Integer> binCountsAL = new ArrayList<Integer>();
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
			p = Pattern.compile("^# neighbor-count bins: (.*)$");
			m = p.matcher(line);
			if (m.matches()) {
				String[] tokens = m.group(1).split("\\s+");
				for (int i=0;i<tokens.length;i++){
					binCountsAL.add(Integer.parseInt(tokens[i].substring(tokens[i].lastIndexOf(':')+1, tokens[i].length())));
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
		
		this.numTypes = indices2types.size();
		this.numCountBins = binCountsAL.size();
		this.binCounts = new int[binCountsAL.size()];
		for (int i=0;i<binCountsAL.size();i++) {
			this.binCounts[i] = binCountsAL.get(i);
		}
		this.scoringMat = new double[numCountBins][numTypes];
		for (int i=0;i<numCountBins;i++) {
			for (int j=0;j<numTypes;j++) {
				this.scoringMat[i][j] = scoringMatAL.get(i).get(j);
			}
		}
		this.types2indices = new HashMap<String, Integer>();
		for (int i:indices2types.keySet()) {
			types2indices.put(indices2types.get(i),i);
		}
	}

}

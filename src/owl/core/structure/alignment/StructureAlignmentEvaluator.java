package owl.core.structure.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import owl.core.runners.PolyposeRunner;
import owl.core.sequence.alignment.AlignmentConstructionException;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.structure.Pdb;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbfilePdb;
import owl.core.util.FileFormatException;

/**
 * Evaluates a multiple structure alignment in terms of conserved core size and RMSD.
 * @author stehr
 *
 */
public class StructureAlignmentEvaluator {

	/*------------------------------ constants ------------------------------*/
	//public static final String DEFAULT_CHAIN_CODE = "A";
	public static final String CCP4_PATH = "/project/StruPPi/Software/CCP4/ccp4-6.0.1";
	public static final String SHELL_PATH = "/bin/sh";
	
	
	/*--------------------------- member variables --------------------------*/
	MultipleSequenceAlignment al;
	TreeMap<String, Pdb> pdbs;
	PolyposeRunner pr;
	TreeSet<Integer> conservedCols;
	TreeSet<Integer> conservedCore;
	
	/*----------------------------- constructors ----------------------------*/
	
	public StructureAlignmentEvaluator(MultipleSequenceAlignment al, TreeMap<String, Pdb> pdbs) throws IOException {
		this.al = al;
		this.pdbs = pdbs;
		this.pr = new PolyposeRunner(CCP4_PATH, SHELL_PATH);		
	}
	
	/*---------------------------- private methods --------------------------*/
	
	/**
	 * Evaluates the RMSD for the given column set.
	 */
	private double getRMSD(TreeSet<Integer> cols) {
		double rmsd = -1.0;
		int[][] positions = new int[pdbs.size()][cols.size()];
		Pdb[] pdbArr = new Pdb[pdbs.size()];
		int idx = 0;
		for(String tag:pdbs.keySet()) {
			// add pdb
			Pdb pdb = pdbs.get(tag);
			pdbArr[idx] = pdb;
			// add positions
			TreeSet<Integer> newCols = alSet2SeqSet(cols, tag);
			int idx2 = 0;
			for(int col:newCols) {
				positions[idx][idx2] = col;
				idx2++;
			}
			idx++;
		}
		try {
			rmsd = pr.superimpose(pdbArr, positions);
		} catch (IOException e) {
			System.err.println(e.getMessage());
			System.err.println("Exiting.");
			System.exit(1);
			return -1.0;
		} 		
		return rmsd;
	}
	
	/**
	 * Evaluates the RMSD for the given column set.
	 * @returns the rmsd or null on error
	 */
	private Matrix3d[] getRotationMatrices(TreeSet<Integer> cols) {
		Matrix3d[] matrices = null;
		//double rmsd = -1.0;
		int[][] positions = new int[pdbs.size()][cols.size()];
		Pdb[] pdbArr = new Pdb[pdbs.size()];
		int idx = 0;
		for(String tag:pdbs.keySet()) {
			// add pdb
			Pdb pdb = pdbs.get(tag);
			pdbArr[idx] = pdb;
			// add positions
			TreeSet<Integer> newCols = alSet2SeqSet(cols, tag);
			int idx2 = 0;
			for(int col:newCols) {
				positions[idx][idx2] = col;
				idx2++;
			}
			idx++;
		}
		try {
			pr.superimpose(pdbArr, positions);
			matrices = pr.getRotationMatrices();
		} catch (IOException e) {
			System.err.println(e.getMessage());
			System.err.println("Exiting.");
			System.exit(1);
			return null;
		} 		
		return matrices;
	}
	
	/**
	 * Converts a set of positions in the alignment to a set of positions in the original sequence.
	 * TODO: Move this to alignment!
	 * @param cols
	 * @return
	 */
	private TreeSet<Integer> alSet2SeqSet(TreeSet<Integer> cols, String tag) {
		TreeSet<Integer> newCols = new TreeSet<Integer>();
		for(int col:cols) {
			newCols.add(this.al.al2seq(tag, col));
		}
		return newCols;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Transforms the given pdbs according to the minimum rmsd superposition on all conserved columns between start and end.
	 * @return the rmsd of the superposition
	 */
	public double superimposeAndTransform(int start, int end) {
		double rmsd = 0;
		Matrix3d[] matrices;
		conservedCore = al.getGaplessColumns(start, end);
		if(conservedCore.size() == 0) {
			System.err.println("Error calling superimposeAndTransform on empty column set.");
		}
		rmsd = getRMSD(conservedCore);	// executes the whole evaluation twice, could be skipped
		matrices = getRotationMatrices(conservedCore);
		if(matrices == null) {
			System.err.println("Error reading rotation matrix from Polypose output.");
		}
		
		// transform pdbs
		int c = 0;
		Vector3d zeros = new Vector3d();
		for (String tag:this.pdbs.keySet()) {	// this has to be the same oder as used in getRotationMatrices
			Pdb pdb = this.pdbs.get(tag);
			// get set of positions in sequence
			TreeSet<Integer> seqPositions = alSet2SeqSet(conservedCore, tag);			
			pdb.moveToOrigin(seqPositions);
			Matrix4d tm = new Matrix4d(matrices[c++], zeros, 1);
			pdb.transform(tm);
		}
		return rmsd;
	}
	
	/**
	 * Evaluates the optimal rmsd superposition for all core sizes down from the maximum number of conserved columns.
	 * Writes a log with the rmsd value for each core size.
	 */
	public void evaluate(PrintStream out, String alName) {
		// evaluate
		//al.printSimple();
		//Sequence.printSeqRuler(al.getAlignmentLength());
//		System.out.println("All non-gapped columns:");
//		conservedCols = al.getGaplessColumns();
//		printColSet(conservedCols);
//		System.out.println("Size: " + conservedCols.size());
//		double baseRMSD = getRMSD(conservedCols);
//		System.out.println("RMSD: " + baseRMSD);
		
//		System.out.println("Longest non-gapped region:");
//		conservedCore = al.getLongestNonGappedRegion();
//		printColSet(conservedCore);
//		System.out.println("Size: " + conservedCore.size());
//		System.out.println("RMSD: " + getRMSD(conservedCore));
		
		conservedCore = al.getGaplessColumns();
		//conservedCore = al.getLongestNonGappedRegion();
		double baseRMSD = getRMSD(conservedCore);
		double bestRMSD = baseRMSD;
		double newRMSD = baseRMSD;
		int bestCol = 999999;
		int coreSize = conservedCore.size();
		double mSAS = baseRMSD * 100 / coreSize;
		double lowestSAS = mSAS;
		out.printf("%s\t%d\t%d\t%.4f\t%.4f\n", alName, coreSize, -1, baseRMSD, lowestSAS);
		
		// greedy optimization
		while(bestCol > 0 && coreSize > 1) {
			baseRMSD = bestRMSD;	// baseRMSD=after last step
			//printColSet(conservedCore);
			//out.println("Core size: " + conservedCore.size());
			//out.println("Core RMSD: " + baseRMSD);
			//out.printf( "mSAS:      %.2f\n", mSAS);
			if(mSAS < lowestSAS) {
				lowestSAS = mSAS;
			}
			bestCol = -1;

			bestRMSD = 1000; // some very high value
			for(int col:conservedCore) {
				// try without this column
				TreeSet<Integer> colsCopy = new TreeSet<Integer>();
				colsCopy.addAll(conservedCore);
				colsCopy.remove(col);
				newRMSD = getRMSD(colsCopy);
				//System.out.println("Removing column " + col + ". New RMSD: " + newRMSD);
				//if(newRMSD > baseRMSD) System.out.println("Unexpected increase in RMSD when removing col " + col);
				if(newRMSD < bestRMSD) {
					bestRMSD = newRMSD;
					bestCol = col;	// always remove a column, even if RMSD increased
				}
			}
			// now remove the column which improved RMSD the most (new: or harmed it the least)
			if(bestCol > 0) {
				conservedCore.remove(bestCol);
				//out.println("Removal of column " + bestCol + " improved RMSD the most.");
			} else {
				//out.println("No improvement in RMSD by removing any columns.");
			}
			coreSize = conservedCore.size();
			mSAS = baseRMSD * 100 / coreSize;
			
			// output result line
			out.printf("%s\t%d\t%d\t%.4f\t%.4f\n", alName, coreSize, bestCol, bestRMSD, mSAS);
			
		}
		//out.println();
		//out.println("Final RMSD:  " + bestRMSD);
		//out.println("Lowest mSAS: " + lowestSAS);		
	}
	
	public int getConservedCoreSize() {
		return conservedCols.size();
	}
	
	public double getConservedCoreRMSD() {
		return 0.0;
	}
	
	public void printConservedCore() {
		printColSet(conservedCols);
	}
	
	public double getTotalRMSD() {
		return getRMSD(conservedCols);
	}
	
	/**
	 * Print a graphical overview of the given set of columns.
	 */
	public void printColSet(TreeSet<Integer> cols) {
		for (int i = 1; i <= cols.last(); i++) {
			if(cols.contains(i)) {
				System.out.print("*");
			} else {
				System.out.print(" ");
			}
		}
		System.out.println();
	}
	
	/** 
	 * Returns a string represenation of the given set 
	 */
	public String getColSetString(TreeSet<Integer> cols) {
		String s = "";
		for(int col:cols) {
			s += String.format("%d,", col);
		}
		return s.substring(0,  s.length()-1);
	}
	
	/*--------------------------------- main --------------------------------*/
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		if(args.length < 2) {
			System.out.println("Usage: AlignmentEvaluator <alignment_file> <pdb_list_file> [<log_file>]");
			System.exit(1);
		}
		
		String alignmentFileName; // = "/project/StruPPi/Software/multalign/results/globins/mustang.aln.fa";
		String listFileName; // = "/project/StruPPi/Software/multalign/data/malecon/globins/pdbs.list";
		String logFileName; // = "/project/StruPPi/Software/multalign/results/globins/mustang.log";
		
		alignmentFileName = args[0];
		listFileName = args[1];
		
		TreeMap<String,Pdb> tag2pdb = new TreeMap<String,Pdb>();
		
		// open log file
		PrintStream log = System.out;

		if(args.length > 2) {
			logFileName = args[2];
			File logFile = new File(logFileName);			
			try {
				log = new PrintStream(new FileOutputStream(logFile));
			} catch (FileNotFoundException e1) {
				System.err.println("Error writing to log file " + logFile.getAbsolutePath());
				System.exit(1);
			}
		}
		// open list file
		try {
		BufferedReader in = new BufferedReader(new FileReader(listFileName));
		String line;
		while((line = in.readLine()) != null) {
			String tag = line.substring(line.length()-9, line.length()-4);
			String chain = tag.substring(4,5);
			//System.out.println(tag);
			try {
			Pdb pdb = new PdbfilePdb(line);
			pdb.load(chain);
			tag2pdb.put(tag,pdb);
			} catch (PdbLoadException e) {
				System.err.println("Error reading from file " + line + ": " + e.getMessage());
			}
		}
		in.close();
		} catch (FileNotFoundException e) {
			System.err.println("File not found:" + listFileName);
			System.exit(1);
		} catch (IOException e) {
			System.err.println("Error reading from:" + listFileName);
			System.exit(1);
		} 

		// load alignment
		MultipleSequenceAlignment al = null;
		try {
			al = new MultipleSequenceAlignment(alignmentFileName, MultipleSequenceAlignment.FASTAFORMAT);
		} catch (FileFormatException e) {
			System.err.println("Error loading alignment " + alignmentFileName + ": " + e.getMessage());
			System.exit(1);
		} catch (AlignmentConstructionException e) {
			System.err.println("Error loading alignment " + alignmentFileName + ": " + e.getMessage());
			System.exit(1);
		} catch (IOException e) {
			System.err.println("Error loading alignment " + alignmentFileName + ": " + e.getMessage());
			System.exit(1);
		}
		
		// create Evaluator
		StructureAlignmentEvaluator alEv = null;
		try {
			alEv = new StructureAlignmentEvaluator(al, tag2pdb);
		} catch (IOException e) {
			System.err.println("Error creating PolyposeRunner: " + e.getMessage());
			System.exit(1);
		}
		// print results
		String alName = alignmentFileName.substring(alignmentFileName.lastIndexOf('/')+1, alignmentFileName.length()-3);
		log.println("#aln_file=" + alignmentFileName);
		log.println("#pdb_file=" + listFileName);
		log.println("#gapless_cols=" + alEv.getColSetString(al.getGaplessColumns()));
		log.println("#name\tncore\tremoved\trmsd\tmSAS");
		alEv.evaluate(log, alName);
		log.close();

	}

}

package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Evaluates a multiple structural in terms of conserved core size and RMSD.
 * @author stehr
 *
 */
public class AlignmentEvaluator {

	/*------------------------------ constants ------------------------------*/
	public static final String DEFAULT_CHAIN_CODE = "A";
	public static final String CCP4_PATH = "/project/StruPPi/Software/CCP4/ccp4-6.0.1";
	public static final String SHELL_PATH = "/bin/sh";
	
	
	/*--------------------------- member variables --------------------------*/
	Alignment al;
	TreeMap<String, Pdb> pdbs;
	PolyposeRunner pr;
	TreeSet<Integer> conservedCols;
	TreeSet<Integer> conservedCore;
	
	/*----------------------------- constructors ----------------------------*/
	
	public AlignmentEvaluator(Alignment al, TreeMap<String, Pdb> pdbs) throws IOException {
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
			return -1.0;
		} 		
		return rmsd;
	}
	
	/**
	 * Converts a set of positions in the alignment to a set of positions in the original sequence.
	 * TODO: Move this to alignment?
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
		int bestCol = 999999;
		int coreSize = conservedCore.size();
		double mSAS = baseRMSD * 100 / coreSize;
		double lowestSAS = mSAS;
		out.printf("%s\t%d\t%d\t%.4f\t%.4f\n", alName, coreSize, -1, baseRMSD, lowestSAS);
		
		// greedy optimization
		while(bestCol > 0 && coreSize > 1) {
			baseRMSD = bestRMSD;
			//printColSet(conservedCore);
			//out.println("Core size: " + conservedCore.size());
			//out.println("Core RMSD: " + baseRMSD);
			//out.printf( "mSAS:      %.2f\n", mSAS);
			if(mSAS < lowestSAS) {
				lowestSAS = mSAS;
			}
			bestCol = -1;

			for(int col:conservedCore) {
				// try without this column
				TreeSet<Integer> colsCopy = new TreeSet<Integer>();
				colsCopy.addAll(conservedCore);
				colsCopy.remove(col);
				double newRMSD = getRMSD(colsCopy);
				//System.out.println("Removing column " + col + ". New RMSD: " + newRMSD);
				//if(newRMSD > baseRMSD) System.out.println("Unexpected increase in RMSD when removing col " + col);
				if(newRMSD < bestRMSD) {
					bestRMSD = newRMSD;
					bestCol = col;
				}
			}
			// now remove the column which improved RMSD the most
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
	
	/*--------------------------------- main --------------------------------*/
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		if(args.length < 2) {
			System.out.println("Usage: AlignmentEvaluater <alignment_file> <pdb_list_file> [<log_file>]");
			System.exit(1);
		}
		
		String alignmentFileName = "/project/StruPPi/Software/multalign/results/globins/mustang.aln.fa";
		String listFileName = "/project/StruPPi/Software/multalign/data/malecon/globins/pdbs.list";
		String logFileName = "/project/StruPPi/Software/multalign/results/globins/mustang.log";
		
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
			//System.out.println(tag);
			try {
			Pdb pdb = new PdbfilePdb(line);
			pdb.load(DEFAULT_CHAIN_CODE);
			tag2pdb.put(tag,pdb);
			} catch (PdbLoadError e) {
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
		Alignment al = null;
		try {
			al = new Alignment(alignmentFileName, Alignment.FASTAFORMAT);
		} catch (FileFormatError e) {
			System.err.println("Error loading alignment " + alignmentFileName + ": " + e.getMessage());
			System.exit(1);
		} catch (AlignmentConstructionError e) {
			System.err.println("Error loading alignment " + alignmentFileName + ": " + e.getMessage());
			System.exit(1);
		} catch (IOException e) {
			System.err.println("Error loading alignment " + alignmentFileName + ": " + e.getMessage());
			System.exit(1);
		}
		
		// create Evaluator
		AlignmentEvaluator alEv = null;
		try {
			alEv = new AlignmentEvaluator(al, tag2pdb);
		} catch (IOException e) {
			System.err.println("Error creating PolyposeRunner: " + e.getMessage());
			System.exit(1);
		}
		// print results
		String alName = alignmentFileName.substring(alignmentFileName.lastIndexOf('/')+1, alignmentFileName.length()-7);
		log.println("#aln_file=" + alignmentFileName);
		log.println("#pdb_file=" + listFileName);
		log.println("#alignment_name\tncore\tremoved\trmsd\tmSAS");
		alEv.evaluate(log, alName);
		log.close();

	}

}

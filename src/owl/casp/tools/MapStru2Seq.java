package owl.casp.tools;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;

import owl.core.sequence.Sequence;
import owl.core.sequence.alignment.PairwiseSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.util.FileFormatException;


/**
 * This class is used to map Casp answers from the PDB to the original target sequences,
 * i.e. extract from the published structure the portion mathing the target sequence and
 * renumber accordingly.
 * @author stehr
 */
public class MapStru2Seq {
	
	public static final String PDBASEDB = "pdbase";
	public static final String DEFAULT_CHAIN_CODE = "A";
	public static final boolean WRITE_REMARKS = true;
	public static final String CIF_FILE_DIR = "/project/StruPPi/BiO/DBd/PDB-REMEDIATED/data/structures/unzipped/all/mmCIF/";
	
	public static void main(String[] args) throws FileNotFoundException {
		
		if(args.length < 3) {
			System.out.println("Usage: MapStru2Seq <seqfile> <pdbcode> <outfile>");
			System.exit(1);
		}
		
		String seqFileName = args[0];
		String rawPdbCode = args[1];
		String outFileName = args[2];
		//String debugPdbFile = "debug.pdb";
		String pdbCode = rawPdbCode;
		String chainCode = DEFAULT_CHAIN_CODE;
		
		if(rawPdbCode.length() == 5) {
			pdbCode = rawPdbCode.substring(0, 4);
			chainCode = rawPdbCode.substring(4,5);
		} else
		if (rawPdbCode.length() != 4) {
			System.err.println(rawPdbCode + " is not a pdb code.");
			System.exit(1);
		}

		// load pdb
		PdbChain pdb = null;
		try {
						
			
			// load from ciffile directory
			File cifFile = new File(CIF_FILE_DIR, pdbCode + ".cif");
			PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile);
			pdb = fullpdb.getChain(chainCode);				
		
		} catch (PdbLoadException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		} catch (IOException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		} catch (FileFormatException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}

		// load sequence
		Sequence seq = new Sequence();		
		try {
			seq.readFromFastaFile(new File(seqFileName));
		} catch (FileFormatException e) {
			System.err.println(e.getMessage());
			System.exit(1);

		} catch (IOException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
		
		// do alignment
		String s1 = seq.getSeq();
		String s2 = pdb.getSequence().getSeq();
		PairwiseSequenceAlignment al = null;
		try {
			al = new PairwiseSequenceAlignment(s1,s2,"target","structure");
		} catch (PairwiseSequenceAlignmentException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
		
		// print alignment
		System.out.println("Alignment of target sequence (top) against PDB structure (bottom):");
		System.out.println("");
		Sequence.printSeqRuler(al.getLength());
		for(String s:al.getAlignedSequences()) {
			System.out.println(s);
		}
		

		int[] mapping = al.getMapping2To1();

		// DEBUG: print mapping and write old structure
//		System.out.println("Structure:");
//		for (int i = 0; i < s2.length(); i++) {
//			System.out.printf("%3d",i+1);
//		}
//		System.out.println();
//		System.out.println("Target:");
//		for (int i = 0; i < s2.length(); i++) {
//			System.out.printf("%3d",mapping[i]+1);	
//		}
//		System.out.println();
//		pdb.writeAtomLines(new PrintStream(new File(debugPdbFile)), true);
		
		// write atom lines
		ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		PrintStream outBytes = new PrintStream(bytes);
		pdb.writeAtomLines(outBytes);
		outBytes.flush();
		String atomStr = bytes.toString();
		
		// update residue numbering
		int observedResidues = 0;
		int notInTarget = 0;
		String[] atomLines = atomStr.split("\n");
		int lastRes = -99999; // some ridiculously low number
		for (int i = 0; i < atomLines.length; i++) {
			String line = atomLines[i];
			int oldRes = Integer.parseInt(line.substring(22, 26).trim());
			int newRes = mapping[oldRes-1] + 1;
			if(lastRes != oldRes) {
				lastRes = oldRes;
				if(newRes < 1) {
					//System.err.println("Warning: Residue " + oldRes + " in structure does not have a corresponding residue in the target sequence. Please check manually.");
					notInTarget++;
				} else {
					observedResidues++;
				}
			}
			StringBuilder sb = new StringBuilder(line);
			sb.replace(22, 26, String.format("%4d", newRes));
			String newLine = sb.toString();
			atomLines[i] = newLine;
		}
		
		// write report
		System.out.println("");
		System.out.println("Residues in target but not in structure: " + (s1.length() - observedResidues));
		System.out.println("Residues in structure but not in target: " + notInTarget);		
		
		// write output file
		try {
			PrintWriter out = new PrintWriter(new File(outFileName));
			if(WRITE_REMARKS) {
				// header
				out.println("REMARK Original PDB: " + pdbCode + DEFAULT_CHAIN_CODE);			
				out.println("REMARK Alignment of target sequence (top) against PDB structure (bottom):");
				//Sequence.printSeqRuler(al.getLength());
				for(String s:al.getAlignedSequences()) {
					out.print("REMARK ");
					out.println(s);
				}
				out.println("REMARK Residues in target but not in structure: " + (s1.length() - observedResidues));
				out.println("REMARK Residues in structure but not in target: " + notInTarget);	
			}
			// atom lines
			for(String line:atomLines) {
				int res = Integer.parseInt(line.substring(22, 26).trim());
				if(res != 0) out.println(line);
			}
			out.close();
		} catch (FileNotFoundException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
	}
	
}

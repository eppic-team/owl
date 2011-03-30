import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.TreeMap;

import owl.core.sequence.alignment.AlignmentConstructionException;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.alignment.StructureAlignmentEvaluator;
import owl.core.util.FileFormatException;


/**
 * Executable for superimposing multiple PDB structures given a multiple alignment
 * @author stehr
 *
 */
public class superimposeMult {

	 public static void main(String[] args) {
		
		 // load pdb objects from file
		 // load alignment
		 // transform objects
		 // write objects to file

		 if(args.length < 3) { // TODO: Change back!
			 System.out.println("Usage: superimposeMult <alignment_file> <pdb_list_file> <out_dir> [<start> <end>]");
			 System.out.println("Superimposes the structure in list on all ungapped columns between start and end");
			 System.out.println("Note: pdb files need to be named pdbcode+chain, e.g. 1tdrB.pdb");
			 System.exit(1);
		 }

		 // Debug: 
//		 String alignmentFileName = "/project/StruPPi/Software/multalign/results/globins15/paul.fa";
//		 String listFileName = "/project/StruPPi/Software/multalign/data/globins15/globins15.pdbs.list";
//		 String outDirName = "/project/StruPPi/Software/multalign/results/globins15/super/paul/";

		 String alignmentFileName = args[0];
		 String listFileName = args[1];
		 String outDirName = args[2];
		 
		 File outDir = new File(outDirName);
		 if(!outDir.isDirectory() || !outDir.canWrite()) {
			 System.err.println("Error. Can not write to " + outDirName);
			 System.err.println("Exiting.");
			 System.exit(1);
		 }

		 TreeMap<String,PdbChain> tag2pdb = new TreeMap<String,PdbChain>();
		 
		 // open list file
		 try {
			 BufferedReader in = new BufferedReader(new FileReader(listFileName));
			 String line;
			 while((line = in.readLine()) != null) {
				 String tag = line.substring(line.length()-9, line.length()-4);
				 String chain = tag.substring(4,5);
				 //System.out.println(tag);
				 try {
					 PdbAsymUnit fullpdb = new PdbAsymUnit(new File(line));
					 PdbChain pdb = fullpdb.getChain(chain);
					 tag2pdb.put(tag,pdb);
				 } catch (PdbLoadException e) {
					 System.err.println("Error reading from file " + line + ": " + e.getMessage());
				 } catch (FileFormatException e) {
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

		 // parse start/end
		 int start = 1;
		 int end = al.getAlignmentLength();
		 if(args.length > 4) {
			 int newStart = Integer.parseInt(args[3]);
			 int newEnd = Integer.parseInt(args[4]);
			 if(newStart >= 1 && newEnd > newStart) {
				 start = newStart;
				 end = newEnd;
			 } else {
				 System.err.println("Alignment on start=" + newStart + " end=" + newEnd + " not possible. Exiting.");
				 System.exit(1);
			 }
		 }
		 System.out.println("start = " + start);
		 System.out.println("end = " + end);
		 
		 // create Evaluator
		 StructureAlignmentEvaluator alEv = null;
		 try {
			 alEv = new StructureAlignmentEvaluator(al, tag2pdb);
		 } catch (IOException e) {
			 System.err.println("Error creating PolyposeRunner: " + e.getMessage());
			 System.exit(1);
		 }
		 
		 System.out.println("Superimposing...");
		 // run transformation
		 double rmsd = alEv.superimposeAndTransform(start, end);
		 
		 // write pdbs to files
		 for(String tag:tag2pdb.keySet()) {
			 PdbChain pdb = tag2pdb.get(tag);
			 File outFile = new File(outDir,tag + ".pdb");
			 try {
				 System.out.println("Writing file " + outFile.getAbsolutePath());
				 PrintStream out = new PrintStream(new FileOutputStream(outFile));
				 pdb.writeAtomLines(out);
				 out.close();
			 } catch(IOException e) {
				 System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
			 }
		 }
		 System.out.println("RMSD=" + rmsd);
	 }
	
}

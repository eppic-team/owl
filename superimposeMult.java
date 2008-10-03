import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.TreeMap;

import proteinstructure.Alignment;
import proteinstructure.AlignmentConstructionError;
import proteinstructure.AlignmentEvaluator;
import proteinstructure.FileFormatError;
import proteinstructure.Pdb;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbfilePdb;

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

		 if(args.length < 0) { // TODO: Change back!
			 System.out.println("Usage: superimposeMult <alignment_file> <pdb_list_file> <out_dir> <start> <end>");
			 System.out.println("Superimposes the structure in list on all ungapped columns between start and end");
			 System.out.println("Note: pdb files need to be named pdbcode+chain, e.g. 1tdrB.pdb");
			 System.exit(1);
		 }

		 String alignmentFileName = "/project/StruPPi/Software/multalign/results/globins15/save/mcma_4.fa";
		 String listFileName = "/project/StruPPi/Software/multalign/data/globins15/globins15.pdbs.list";
		 String outDirName = "/project/StruPPi/Software/multalign/results/globins15/save/super/";
//		 String alignmentFileName = "/project/StruPPi/Software/multalign/results/sk_cl5_2/mcma_1.fa";
//		 String listFileName = "/project/StruPPi/Software/multalign/data/sk_cl5_2/sk_cl5_2.pdbs.list";
//		 String outDirName = "/project/StruPPi/Software/multalign/results/sk_cl5_2/super/";
		 int start = 1;
		 int end = 1000;

		 //alignmentFileName = args[0];
		 //listFileName = args[1];
		 //outDirName = args[2];
		 
		 File outDir = new File(outDirName);
		 if(!outDir.isDirectory() || !outDir.canWrite()) {
			 System.err.println("Error. Can not write to " + outDirName);
			 System.exit(1);
		 }

		 TreeMap<String,Pdb> tag2pdb = new TreeMap<String,Pdb>();
		 
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
		 // run transformation
		 double rmsd = alEv.superimposeAndTransform(start, end);
		 
		 // write pdbs to files
		 for(String tag:tag2pdb.keySet()) {
			 Pdb pdb = tag2pdb.get(tag);
			 File outFile = new File(outDir,tag + ".pdb");
			 try {
				 System.out.println("Writing file " + outFile.getAbsolutePath());
				 PrintStream out = new PrintStream(new FileOutputStream(outFile));
				 pdb.writeAtomLines(out, true);
				 out.close();
			 } catch(IOException e) {
				 System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
			 }
		 }
		 System.out.println("RMSD=" + rmsd);
	 }
	
}

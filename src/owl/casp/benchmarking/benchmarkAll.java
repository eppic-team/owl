package owl.casp.benchmarking;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.LinkedList;

/**
 * This class provides the main method for benchmarking all our Casp8/Casp9 comparative modelling predictions.
 * 
 * The target ids to be processed are passed as a list file (= text file with one target id per line).
 * 
 * The following files are created:
 * - text files with a ranking per target in text format
 * - img files with a comparsion of model vs. native for each of our models
 * - HTML report pages for each of the category and some variants (comparison with other groups, new targets)
 * 
 * The following external resources are being used:
 * - our submitted models from <base_dir>/submitted
 * - the renumbered native structures from <base_dir>/answers
 * - server models from <base_dir>/server_models
 * 
 * @author stehr
 * @date 2008-09-02
 * @date 2010-08-02	update for Casp9, added parameters basedir and outdir, group suffix still hardcoded
 */
public class benchmarkAll {

	public static String[] getTargetsFromDir(String dirName) {
		
		File targetDir = new File(dirName);
		
		String[] rawTargets = targetDir.list();
		LinkedList<String> targetList = new LinkedList<String>();
		
		for(String rawTarget:rawTargets) {
			if(rawTarget.endsWith("_1")) {
				targetList.add(rawTarget.substring(0, 5));
			}
		}
		String[] targets = new String[targetList.size()];
		targets = targetList.toArray(targets);
		return targets;
	}
	
	public static String[] getTargetsFromListFile(String listFileName) {
		LinkedList<String> targetList = new LinkedList<String>();
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(listFileName));
			String target;
			while((target = in.readLine()) != null) {
				targetList.add(target.trim());
			}
		} catch (IOException e) {
			System.err.println("Error reading list file: " + listFileName);
		}
		
		String[] targets = new String[targetList.size()];
		targets = targetList.toArray(targets);
		return targets;
	}
	
	public static void main(String[] args) {
	
		boolean writeTempBased = false;
		boolean writeCccp3D = false;
		boolean writeCccpCM = true;
		
		boolean writeImages = false;
		boolean writeDetails = false;
		boolean writeHTML = true;

		if(args.length < 2) {
			System.out.println("Usage: benchmarkAll <base_dir> <out_dir>");
			System.out.println("e.g. benchmarkAll /project/StruPPi/CASP8/ /home/web/lappe/casp8/");
			System.exit(1);
		}
		
		String baseDir = args[0];
		String outBaseDir = args[1];
		
		String allTargetsList = baseDir + "/targets/all_targets.list";
		String humanTargetsList = baseDir + "/targets/human_targets.list";
		String newTargetsList = baseDir + "/answers/last_update.list";
		String newHumanTargetsList = baseDir + "/answers/last_update_human.list";
		
		String[] allTargets = getTargetsFromListFile(allTargetsList);
		String[] humanTargets = getTargetsFromListFile(humanTargetsList);
		String[] newTargets = getTargetsFromListFile(newTargetsList);
		String[] newHumanTargets = getTargetsFromListFile(newHumanTargetsList);
		
		Benchmarking bm = new Benchmarking(null);
		String highlightGroup = null;
		File outDir = null;
		File predDir = null;
		String groupSuffix = null;
		boolean eval3D;
		PrintStream out;
		File outFile;
		int targetsPerRow;
				
		// 1. Template based
		if(writeTempBased) {
			System.out.println("Template based...");
			outDir = new File(outBaseDir);
			predDir = new File(baseDir + "/submitted");
			groupSuffix = "TS183";
			eval3D = true;
			targetsPerRow = 19;
			// Update images for template based
			if(writeImages) {
				System.out.println("Generating images...");
				try {
					bm.writeAllImages(predDir, groupSuffix, newTargets, 200, outDir);
					bm.writeAllImages(predDir, groupSuffix, newTargets, 600, outDir);
				} catch (IOException e) {
					System.err.println("Error creating images: " + e.getMessage());
				}
			}
			// Update txt files for template based
			if(writeDetails) {
				System.out.println("Writing txt files...");
				bm.writeAllTargetDetailFiles(predDir, groupSuffix, newTargets, outDir, eval3D);
			}
			// Write main HTML report
			if(writeHTML) {
				outFile = new File(outDir, "results.html");
				try {
					out = new PrintStream(new FileOutputStream(outFile));
					System.out.println("Writing " + outFile.getAbsolutePath());
					bm.printAllResultsHTML(predDir, groupSuffix, humanTargets, highlightGroup, outDir, out, eval3D,targetsPerRow);
				} catch (IOException e) {
					System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
				}
				// Write new targets HTML report
				outFile = new File(outDir, "new.html");
				try {
					out = new PrintStream(new FileOutputStream(outFile));
					System.out.println("Writing " + outFile.getAbsolutePath());
					bm.printAllResultsHTML(predDir, groupSuffix, newHumanTargets, highlightGroup, outDir, out, eval3D,targetsPerRow);
				} catch (IOException e) {
					System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
				}
				// Write comparison with Baker
				highlightGroup = "BAKER-ROBETTA";
				outFile = new File(outDir, "baker.html");
				try {
					out = new PrintStream(new FileOutputStream(outFile));
					System.out.println("Writing " + outFile.getAbsolutePath());
					bm.printAllResultsHTML(predDir, groupSuffix, humanTargets, highlightGroup, outDir, out, eval3D,targetsPerRow);
				} catch (IOException e) {
					System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
				}
				// Write comparison with Zhang
				highlightGroup = "Zhang";
				outFile = new File(outDir, "zhang.html");
				try {
					out = new PrintStream(new FileOutputStream(outFile));
					System.out.println("Writing " + outFile.getAbsolutePath());
					bm.printAllResultsHTML(predDir, groupSuffix, humanTargets, highlightGroup, outDir, out, eval3D,targetsPerRow);
				} catch (IOException e) {
					System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
				}
				// Write comparison with HHPred
				highlightGroup = "HHpred";
				outFile = new File(outDir, "HHpred.html");
				try {
					out = new PrintStream(new FileOutputStream(outFile));
					System.out.println("Writing " + outFile.getAbsolutePath());
					bm.printAllResultsHTML(predDir, groupSuffix, humanTargets, highlightGroup, outDir, out, eval3D,targetsPerRow);
				} catch (IOException e) {
					System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
				}
				// Write comparison with Torda
				highlightGroup = "torda";
				outFile = new File(outDir, "torda.html");
				try {
					out = new PrintStream(new FileOutputStream(outFile));
					System.out.println("Writing " + outFile.getAbsolutePath());
					bm.printAllResultsHTML(predDir, groupSuffix, humanTargets, highlightGroup, outDir, out, eval3D,targetsPerRow);
				} catch (IOException e) {
					System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
				}
			}
			highlightGroup = null;
		}

		// 2. Meta 3D
		if(writeCccp3D) {
			System.out.println("CCCP...");
			outDir = new File(outBaseDir + "/cccp/");
			predDir = new File(baseDir + "/submitted/cccp");
			groupSuffix = "TS014";
			targetsPerRow = 25;
			eval3D = true;
			// Update images for cccp-3D
			if(writeImages) {
				System.err.println("Generating images...");
				try {
					bm.writeAllImages(predDir, groupSuffix, newTargets, 200, outDir);
					bm.writeAllImages(predDir, groupSuffix, newTargets, 600, outDir);
				} catch (IOException e) {
					System.err.println("Error creating images: " + e.getMessage());
				}
			}
			// Update txt files for cccp-3D
			if(writeDetails) {
				System.out.println("Writing txt files...");
				bm.writeAllTargetDetailFiles(predDir, groupSuffix, newTargets, outDir, eval3D);
			}
			// Write main report
			if(writeHTML) {
				outFile = new File(outDir, "index.html");
				try {
					out = new PrintStream(new FileOutputStream(outFile));
					System.out.println("Writing " + outFile.getAbsolutePath());
					bm.printAllResultsHTML(predDir, groupSuffix, allTargets, highlightGroup, outDir, out, eval3D,targetsPerRow);
				} catch (IOException e) {
					System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
				}
				// write report for new targets
				outFile = new File(outDir, "new.html");
				try {
					out = new PrintStream(new FileOutputStream(outFile));
					System.out.println("Writing " + outFile.getAbsolutePath());
					bm.printAllResultsHTML(predDir, groupSuffix, newTargets, highlightGroup, outDir, out, eval3D,targetsPerRow);
				} catch (IOException e) {
					System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
				}
			}
		}

		// 3. Meta CM
		if(writeCccpCM) {
			outDir = new File(outBaseDir + "/cccp/cm/");
			predDir = new File(baseDir + "/submitted/cccp");
			groupSuffix = "TS014";
			eval3D = false;
			targetsPerRow = 25;
			// Update txt files for cccp-3D
			if(writeDetails) {
				System.out.println("Writing txt files...");
				bm.writeAllTargetDetailFiles(predDir, groupSuffix, newTargets, outDir, eval3D);
			}
			// Write main report
			if(writeHTML) {
				outFile = new File(outDir, "index.html");	
				try {
					out = new PrintStream(new FileOutputStream(outFile));
					System.out.println("Writing " + outFile.getAbsolutePath());
					bm.printAllResultsHTML(predDir, groupSuffix, allTargets, highlightGroup, outDir, out, eval3D,targetsPerRow);
				} catch (IOException e) {
					System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
				}
				// write report for new targets
				outFile = new File(outDir, "new.html");
				try {
					out = new PrintStream(new FileOutputStream(outFile));
					System.out.println("Writing " + outFile.getAbsolutePath());
					bm.printAllResultsHTML(predDir, groupSuffix, newTargets, highlightGroup, outDir, out, eval3D,targetsPerRow);
				} catch (IOException e) {
					System.err.println("Error writing to file " + outFile.getAbsolutePath() + ":" + e.getMessage());
				}
			}
		}

		// 4.
		// Write index page
		// Write summary page
		System.out.println("done.");
	}
}

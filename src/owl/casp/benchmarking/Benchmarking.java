package owl.casp.benchmarking;

import java.util.*;
import java.io.*;

import owl.core.runners.MaxClusterRunner;
import owl.core.runners.MaxClusterRunner.MaxClusterRow;
import owl.core.structure.Pdb;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbfilePdb;
import owl.core.structure.graphs.GraphComparisonResult;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.StreamGobbler;


public class Benchmarking {
	
	/*------------------------------ constants ------------------------------*/
	
	public static final boolean EVAL_3D = false;	// otherwise evaluate contacts (acc+cov)
	
	public static final String TMP_DIR = System.getProperty("java.io.tmpdir");
	public static final String maxClusterExecutable = "/project/StruPPi/bin/maxcluster";	
	public static final String tempFileName = "/tmp/casp.benchmark.tmp";

	// Template based:
//	public static final String GROUP_SUFFIX = "TS183"; // e.g. "T0xxx" + GROUP_SUFFIX + "_1"
//	public static final String predictionDir = "/project/StruPPi/CASP8/submitted";
//	public static final String targetList    = "/project/StruPPi/CASP8/targets/human_targets.list";
//	public static final int TARGETS_PER_ROW = 19; // for 57 human targets this gives exactly three rows
//	public static final String IMG_DIR = "/home/web/lappe/casp8/";
	
	// CCCP:
//	public static final String predictionDir = "/project/StruPPi/CASP8/submitted/cccp/";
//	public static final String GROUP_SUFFIX = "TS014"; // e.g. "T0xxx" + GROUP_SUFFIX + "_1"
//	public static final String targetList   = "/project/StruPPi/CASP8/targets/all_targets.list";
//	public static final int TARGETS_PER_ROW = 32; // for 128 server targets this gives exactly 8 rows
//	public static final String IMG_DIR = "/home/web/lappe/casp8/cccp/cm/";
		
	public static final int NO_ERROR = 0;
	public static final int PRED_NOT_FOUND = -1;
	public static final int ANSWER_NOT_FOUND = -2;
	public static final int OTHER_ERROR = -9;
	
	public static final String[] CmDescription = {"","Ensemble","Prediction","Native","Comparison",""};
	
	/*--------------------------- type definitions --------------------------*/	
	public class TargetResult {
		public String target;
		public int bestRank = 0;
		public int mod1Rank = 0;
		public int numRanks = 0;
		public double bestGdt = 0;
		public double mod1Gdt = 0;
		public double avgGdt = 0;
		public double medGdt = 0;
		
		public TargetResult(String target, int bestRank, int mod1Rank, int numRanks, double bestGdt, double mod1Gdt, double avgGdt, double medGdt) {
			this.target = target;
			this.bestRank = bestRank;
			this.mod1Rank = mod1Rank;
			this.numRanks = numRanks;
			this.bestGdt = bestGdt;
			this.mod1Gdt = mod1Gdt;
			this.avgGdt = avgGdt;
			this.medGdt = medGdt;
		}
	}
	
	/*--------------------------- member variables --------------------------*/
	File serverModelDir;
	File answerDir;
	ArrayList<MaxClusterRunner.MaxClusterRow> lastResultTable = null;
	int lastNumModels = -1;
	String lastTarget = null;
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Class to produce, hold and output benchmarking results for our Casp8 predictions
	 */
	public Benchmarking(String[] targets, File baseDir) {
		serverModelDir = new File(baseDir, "server_models");
		answerDir = new File(baseDir, "answers");
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Evaluates the given target and stores the results in the member variables
	 * lastResultTable and lastNumModels. Different outputs of the results can
	 * then be printed by the methods below (printLastResult...)
	 * TODO: This design is not optimal but works for the moment...
	 * @return An error code or 0 if everything was fine
	 */
	public int getTargetResult(File predDir, String groupSuffix, String target, boolean eval3D) {
		
		String predictionFileName = target + groupSuffix + "_1";
		String answerFileName = target + ".pdb";

		File predictionFile = new File(predDir, predictionFileName);
		File answerFile = new File(answerDir, answerFileName);
		File tempFile = new File(tempFileName);
		File ourModelDir = predDir;
		// For CCCP: File ourModelDir = new File(predictionDir, target);
		File otherModelDir = new File(serverModelDir, target);
		
		// check files
		if(!answerFile.canRead()) {
			//System.err.println("Could not read " + answerFile.getAbsolutePath());
			return ANSWER_NOT_FOUND;
		}
		if(!predictionFile.canRead()) {
			//System.err.println("Could not read " + predictionFile.getAbsolutePath());
			return PRED_NOT_FOUND;
		}
		
		// make list file
		PrintWriter out = null;
		try {
			out = new PrintWriter(tempFile);
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}
		String[] ourModels = ourModelDir.list();
		int numModels = 0;
		for(String modName: ourModels) {
			if(modName.startsWith(target)) {
			// For CCCP: if(modName.endsWith(".ts")) {
				File mod = new File(ourModelDir, modName);
				out.println(mod.getPath());
				numModels++;
			}
		}
		//out.println(predictionFile);
		String[] serverModels = otherModelDir.list();
		for(String modName: serverModels) {
			if(modName.endsWith("TS1")) {
				File mod = new File(otherModelDir, modName);
				out.println(mod.getPath());
			}
		}
		out.close();
		
		// run maxcluster
		MaxClusterRunner mcr;
		ArrayList<MaxClusterRunner.MaxClusterRow> resultTable = null;
		if(eval3D) {
			// evaluate by GDT
			try {
				mcr = new MaxClusterRunner(maxClusterExecutable);
				resultTable = mcr.calculateRanking(tempFile.getPath(), answerFile.getPath(), MaxClusterRunner.ScoreType.GDT);
			} catch (FileNotFoundException e) {
				System.err.println("Error running Maxcluster: " + e.getMessage());
				return OTHER_ERROR;
			} catch (IOException e) {
				System.err.println("Error running Maxcluster: " + e.getMessage());
				return OTHER_ERROR;
			}
		} else {
			// evaluate by contact map overlap
			resultTable = getContactMapOverlapResultTable(tempFile.getPath(), answerFile.getPath());
			if(resultTable == null) {
				return OTHER_ERROR;
			}
		}
		
		if(numModels >= resultTable.size() ) {
		    System.err.println("Error, resultTable.size() <= numModels:" + resultTable.size());
		    return OTHER_ERROR;
		}
		
		// order results
		Collections.sort(resultTable, new Comparator<MaxClusterRunner.MaxClusterRow>() {
			public int compare(MaxClusterRunner.MaxClusterRow arg0, MaxClusterRunner.MaxClusterRow arg1) {
				return new Integer(arg0.getRank()).compareTo(arg1.getRank());
			}
		});
		
		// store results in member variables
		this.lastNumModels = numModels;
		this.lastResultTable = resultTable;
		this.lastTarget = target;
		
		return NO_ERROR;
	}
	
	private ArrayList<MaxClusterRow> getContactMapOverlapResultTable(String listFileName, String answerFileName) {
		ArrayList<MaxClusterRow> resultTable = new ArrayList<MaxClusterRow>();

		// load answer structure
		Pdb answerPdb = new PdbfilePdb(answerFileName);
		try {
			answerPdb.load(answerPdb.getChains()[0]);
			// generate graph for answer file
			RIGraph answerGraph = answerPdb.getRIGraph("Cb", 8.0);
			// loop over list file
			BufferedReader in = null;
			try {
				in = new BufferedReader(new FileReader(listFileName));
				String line;
				int idx = 1;
				while((line=in.readLine()) != null) {
					String fileName = line.trim();
					// load structure
					Pdb pdb = new PdbfilePdb(fileName);
					try {
						pdb.load(pdb.getChains()[0]);
						// generate graph
						RIGraph rig = pdb.getRIGraph("Cb", 8.0);
						// compare graph to answer
						GraphComparisonResult eval = rig.evaluatePrediction(answerGraph);
						double score = 100.0 * (eval.accuracy + eval.coverage) / 2;
						// generate maxClusterRow
						MaxClusterRunner.MaxClusterRow resultRow = (new MaxClusterRunner(maxClusterExecutable)).new MaxClusterRow(fileName, idx, score);
						// add row to resultTable
						resultTable.add(resultRow);
						idx++;
					} catch (PdbLoadError e) {
						//e.printStackTrace();
					}
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		} catch (PdbLoadError e) {
			e.printStackTrace();
		}
		if(resultTable != null) {
			// order table by score
			Collections.sort(resultTable, new Comparator<MaxClusterRunner.MaxClusterRow>() {
				public int compare(MaxClusterRunner.MaxClusterRow arg0, MaxClusterRunner.MaxClusterRow arg1) {
					return new Double(arg1.getScore()).compareTo(arg0.getScore());
				}
			});
			// write ranks
			int rank=1;
			for(MaxClusterRow row:resultTable) {
				row.setRank(rank++);
			}
			// ordering by index not necessary here			
		} else {
			System.err.println("Error in getContactMapOverlapResultTable(): return null");
		}
		return resultTable;
	}

	public void printLastResultTable(PrintStream out) {
		// print summary
		printLastResultSummary(out);
		out.println();
		
		// print results
		for (MaxClusterRunner.MaxClusterRow row: lastResultTable) {
			if(row.getIndex() <= lastNumModels) out.println();
			out.println(row);
			if(row.getIndex() <= lastNumModels) out.println();
		}		
	}
	
	public void printLastResultSummary(PrintStream out) {
		
		TargetResult tr = getLastResultRow();
		
		// print results
	    out.println("Target: " + lastTarget);
//	    out.println("Our best rank: " + tr.bestRank + "/" + tr.numRanks);
//	    out.println("Our mod1 rank: " + tr.mod1Rank + "/" + tr.numRanks);
//	    out.println("Our best GDT: " + tr.bestGdt);
//	    out.println("Our mod1 GDT: " + tr.mod1Gdt);
//	    out.printf("Avg GDT: %.3f\n", tr.avgGdt);
//	    out.println("Med GDT: " + tr.medGdt);
	    
	    // for contact map overlap score
	    out.println("Our rank: " + tr.bestRank + "/" + tr.numRanks);
	    out.println("Our score: " + tr.bestGdt);
	    out.printf("Avg score: %.3f\n", tr.avgGdt);
	    out.println("Med score: " + tr.medGdt);
	}
	
	public void printLastResultSummaryHTML() {
		TargetResult tr = getLastResultRow();
		
		System.out.println("<table align=left>");
	    System.out.printf("<tr><td>Target:</td><td>%s</td></tr>\n", lastTarget);
	    System.out.printf("<tr><td>Our best rank:</td><td>%s/%s</td></tr>\n",tr.bestRank, tr.numRanks);
	    System.out.printf("<tr><td>Our mod1 rank:</td><td>%s/%s</td></tr>\n",tr.mod1Rank, tr.numRanks);
	    System.out.printf("<tr><td>Our best GDT:</td><td>%s</td></tr>\n",tr.bestGdt);
	    System.out.printf("<tr><td>Our mod1 GDT:</td><td>%s</td></tr>\n",tr.mod1Gdt);
	    System.out.printf("<tr><td>Avg GDT:</td><td>%.3f</td></tr>\n", tr.avgGdt);
	    System.out.printf("<tr><td>Med GDT:</td><td>%s</td></tr>\n", tr.medGdt);	
		System.out.println("</table>");
	}
	
	public TargetResult getLastResultRow() {

		// init results
		int ourBestRank = 0;
		int ourMod1Rank = 0;
		int numRanks = 0;
		double ourBestGdt = 0;
		double ourMod1Gdt = 0;
		double sumGdts = 0;
		int numGdts = 0;
		double avgGdt = 0;
		double medGdt = 0;
		double[] scores = new double[lastResultTable.size()-lastNumModels];
		
		// collect results
		int c = 0;
		for (MaxClusterRunner.MaxClusterRow row: lastResultTable) {
			if(row.getIndex() > lastNumModels) { // not our model
				numGdts++;
				sumGdts += row.getScore();
				scores[c++] = row.getScore();
			} else {	// our model
				if(row.getFileName().contains("_1")) { // model 1
					ourMod1Gdt = row.getScore();
					ourMod1Rank = row.getRank();
				}
				if(ourBestRank == 0) {	// no best rank yet
					ourBestRank = row.getRank();
					ourBestGdt = row.getScore();
				}
			}
		}
		avgGdt = sumGdts / numGdts;
		Arrays.sort(scores);
		medGdt = scores[scores.length/2];
		numRanks = lastResultTable.size();
		
		return new TargetResult(lastTarget,ourBestRank,ourMod1Rank,numRanks,ourBestGdt,ourMod1Gdt,avgGdt,medGdt);
	}
	
	public void printLastResultHTML(String highlightGroup, File imgFileDir, boolean eval3D, PrintStream out) {
		out.println("<a name=\"" + lastTarget + "\">");
		out.println("<table><tr><td>");
	    out.println("<a href=\"" + lastTarget + ".txt\"><img src=\"" 
	    		+ getGoogleChartURL(lastResultTable, lastNumModels, lastTarget, highlightGroup) 
	    		+ "\" border=0 /></a>");
	    out.println("</td><td>");
	    for (int mod = 1; mod <=5 ; mod++) {
	    //for (int mod = 1; mod <=1 ; mod++) {

	    	if(new File(imgFileDir, getImageFileName(lastTarget, mod, 200)).canRead()) {
		    	out.print("<table><caption align=\"bottom\"><font size=1 color=gray face=\"Verdana,Arial,Helvetica\">");
		    	if(eval3D) {
		    		out.print("Model " + mod);
		    	} else {
		    		out.print(CmDescription[mod]);
		    	}
		    	out.print("</font></caption><tr><td>");
			    String imgFileName = getImageFileName(lastTarget,mod,200);
			    String bigImage = getImageFileName(lastTarget,mod,600);
			    out.print("<a href=\"" + bigImage + "\"><img src=\"" + imgFileName + "\"/ border=0></a>");
			    out.println("</td></tr></table>");
			    out.println("</td><td>");
	    	}
		}
		out.println("</td></tr></table>");
	    out.println("<br>");
	    //out.println("<a href=\"" + lastTarget + ".txt\">[Details]</a><br>");
	}
	
	public void printLastResultRowHTML() {
		TargetResult tr = getLastResultRow();
		System.out.println("<tr>");
		System.out.printf("<td><a href=\"images.html#%s\">%s</a></td><td>%d/%d</td><td>%d/%d</td><td>%.3f</td><td>%.3f</td><td>%.3f</td><td>%.3f</td>\n", 
				tr.target, tr.target, tr.bestRank, tr.numRanks, tr.mod1Rank, tr.numRanks,
				tr.bestGdt, tr.mod1Gdt, tr.avgGdt, tr.medGdt);
		System.out.println("</tr>");
	}		
	
	
	public void printLastResultRow() {
		TargetResult tr = getLastResultRow();
	    System.out.printf("%s\t%d/%d\t%d/%d\t%.3f\t%.3f\t%.3f\t%.3f\n", 
			      tr.target, tr.bestRank, tr.numRanks, tr.mod1Rank, tr.numRanks,
			      tr.bestGdt, tr.mod1Gdt, tr.avgGdt, tr.medGdt);		
	}
	
	public void printAllResultsTable(File predDir, String groupSuffix, String[] targets, boolean eval3D) {
		System.out.printf("Target\tbestRnk\tmod1Rnk\tbestGdt\tmod1Gdt\tavgGdt\tmedGdt\n");
		
//		Benchmarking.TargetResult tr;
//		int sumBestRanks = 0;
//		int sumMod1Ranks = 0;
//		int numRanks = 0;
		
		for(String target:targets) {
			int ret = this.getTargetResult(predDir, groupSuffix, target, eval3D);
			if(ret==NO_ERROR) {
				this.printLastResultRow();
			}
		}
	}
	
	public void printAllResultsHTMLTable(File predDir, String groupSuffix, String[] targets, boolean eval3D) {
		System.out.println("<html>");
		System.out.println("<table>");
		System.out.printf("<tr><th>Target</th><th>bestRnk</th><th>mod1Rnk</th><th>bestGdt</th><th>mod1Gdt</th><th>avgGdt</th><th>medGdt</th></tr>\n");
		for(String target:targets) {
			int ret = this.getTargetResult(predDir, groupSuffix, target, eval3D);
			if(ret==NO_ERROR) {
				this.printLastResultRowHTML();
			} else
				if(ret==PRED_NOT_FOUND) {
					System.out.printf("<tr><td>%s</td><td colspan=6>No predictions found</td></tr>\n", target);
				} else
				if(ret==ANSWER_NOT_FOUND) {
					System.out.printf("<tr><td>%s</td><td colspan=6>Structure not released yet</td></tr>\n", target);				
				} else {
					System.out.printf("<tr><td>%s</td><td colspan=6>Processing error</td></tr>\n", target);
			}
		}			
		System.out.println("</table>");		
		System.out.println("</html>");
		
	}
	
	/**
	 * 
	 * @param targets the list of targets to be processed in the form T0xxx
	 * @param highlightGroup the name of a group to hightlight in the results or null if none
	 */
	public void printAllResultsHTML(File predDir, String groupSuffix, String[] targets, String highlightGroup, File outDir, PrintStream out, boolean eval3D, int targetsPerRow) {
		String msg;
		if(eval3D) {
			out.println("<html>");
			out.printf("<font face=\"Verdana,Arial,Helvetica\" size=\"1\" color=\"gray\">");
			out.println("Bar chart: GDT score of all server models (blue) compared to our models (orange)<br>");
			out.println("Images: Our model (purple) compared to native (green)<br>");
			out.println("Click on bar chart to view the actual numbers<br>");
			out.println("Click on the image to get a larger version<br>");
			out.println("Targets without a link below have not been released yet<br>");
			out.println("<p>");
		} else {
			out.println("<html>");
			out.printf("<font face=\"Verdana,Arial,Helvetica\" size=\"1\" color=\"gray\">");
			out.println("Bar chart: 100*(acc+cov)/2 score of server models converted to contact maps (blue) and our predicted contact map (orange)<br>");
			out.println("acc = number of correctly predicted contacts (TP) / number of predicted contacts (TP+FP)<br>");
			out.println("cov = number of correctly predicted contacts (TP) / number of native contacts (TP+FN)<br>");
			out.println("Ensemble contact map: Overlay of all 3D server models converted to contact maps<br>");
			out.println("Comparison contact map: black=true positives (TP), red=false positives (FP), green=false negatives (FN)<br>");
			out.println("Click on bar chart to view the actual numbers<br>");
			out.println("Click on contact maps to get a larger version<br>");
			out.println("Targets without a link below have not been released yet<br>");
			out.println("<p>");
		}
		int c = 1;
		for(String target:targets) {
			if(new File(answerDir, target + ".pdb").canRead()) {
				out.print("<a href=\"#" + target + "\">" + target + "</a> ");
			} else {
				out.println(target + " ");
			}
			if(c++ % targetsPerRow == 0) out.print("<br>");
		}
		out.printf("</font><br><p>");
		for(String target:targets) {
			int ret = this.getTargetResult(predDir, groupSuffix, target, eval3D);
			if(ret==NO_ERROR) {
				msg = "Ok";
				this.printLastResultHTML(highlightGroup, outDir, eval3D, out);
			} else {
				if(ret==PRED_NOT_FOUND) {
					msg = "No predictions found";
				} else
					if(ret==ANSWER_NOT_FOUND) {
						msg = "Structure not released yet";
					} else {
						msg = "Processing error";
					}
				
				//out.printf("<font face=\"Verdana,Arial,Helvetica\" size=\"1\" color=\"gray\">Target %s: %s</font><br>", target, msg);
				
			}
			System.out.printf("%s: %s\n", target, msg);
		}		
		
		out.println("</html>");
	}
	
	/**
	 * Generate image files for the given list of targets.
	 * @param targets
	 * @param size
	 * @throws IOException
	 */
	public void writeAllImages(File predDir, String groupSuffix, String[] targets, int size, File outDir) throws IOException {
		for(String target:targets) {
			for (int mod = 1; mod <= 5; mod++) {
				String imgFileName = getImageFileName(target, mod, size);
				File pred = new File(predDir, target + groupSuffix + "_" + mod);
				File answer = new File(answerDir, target + ".pdb");
				writePymolImage(pred, answer, new File(outDir, imgFileName), size);				
			}

		}
	}
	
	/**
	 * Return the default name for an image of the given target with the given size.
	 * @param target
	 * @param the model number (1-5)
	 * @param size
	 * @return
	 */
	public String getImageFileName(String target, int model, int size) {
		return target + "_" + model + "_" + size + ".png";
	}
	
	/**
	 * Writes a pymol script to compare prediction and native structure for the given target and write to an image file
	 * @param target
	 * @param size the size of the output image in pixels (height=width)
	 */
	public boolean writePymolImage(File pred, File answer, File imgFile, int size) throws IOException {
		
		String PYMOL_EXECUTABLE = "/project/StruPPi/bin/pymol-1.0 -c";
		String TMP_SCRIPT_FILE = Math.abs(new Random().nextInt()) + ".pml";
		
		File script = new File(TMP_DIR, TMP_SCRIPT_FILE);
		script.deleteOnExit();
		
		if(!pred.canRead()) {
			System.err.println("File not found: " + pred.getAbsolutePath());
			return false;
		}
		
		if(!answer.canRead()) {
			System.err.println("File not found: " + answer.getAbsolutePath());
			return false;
		}
		

		// write pymol script
		PrintWriter out = new PrintWriter(script);
		int width = size;
		int height = size;
		int dpi = -1;
		int ray = 1;
		String display = "cartoon";
		String background = "white";

		// load both structures
		out.println("bg_color " + background);
		out.println("load " + pred.getAbsolutePath() + ", pred");
		out.println("load " + answer.getAbsolutePath() + ", answer");
		out.println("super pred, answer");
		out.println("hide all");
		out.println("show " + display);
		out.println("orient");
		out.println("color magenta, pred");
		out.println("color green, answer");
		out.println("png " + imgFile.getAbsolutePath() + "," + width + "," + height + "," + dpi + "," + ray);
		out.println("quit");
		out.close();

		// execute script (writes image file)
		//System.out.println("Script written to " + script.getAbsolutePath());
		String cmd = PYMOL_EXECUTABLE + " " + script.getAbsolutePath();
		Process p = Runtime.getRuntime().exec(cmd);
		StreamGobbler s1 = new StreamGobbler ("stdout", p.getInputStream ());
		StreamGobbler s2 = new StreamGobbler ("stderr", p.getErrorStream ());
		s1.start();
		s2.start();
		try {
			p.waitFor();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
		if(imgFile.canRead()) {
			System.err.println("Image file written to " + imgFile.getAbsolutePath());
		} else {
			System.err.println("Problems with PyMol, image file not created.");
		}
		return true;
	}
	
	/**
	 * @return The URL to produce a google chart for the given result table.
	 * @param hightlightGroup the name of a group to highlight or null if none
	 */
	public String getGoogleChartURL(ArrayList<MaxClusterRunner.MaxClusterRow> resultTable, int numModels, String title, String highlightGroup) {
		int width = 580;
		int height = 200;
		int barWidth = 5;
		int barSpace = 2;
		String ourModelColor = "FFAA00";
		String othersColor = "76A4FB";
		String bakerColor = "3664BB";
		
		// chart type and title
		String url = "http://chart.apis.google.com/chart?chtt=" + title + "&cht=bvs";
		
		// sizes
		url += "&chbh="+barWidth+","+barSpace+"&chs=" + width + "x" + height;
		
		// numeric values
		 url += "&chd=t:";
		for (MaxClusterRunner.MaxClusterRow row: resultTable) {
			url += String.format("%.1f,", row.getScore());
		}
		url = url.substring(0, url.length()-1); // remove trailing comma
		
		// labels
		url += "&chxt=y,x";
		url += "&chxl=0:|0|50|100|";
		url += "1:|";
		for (MaxClusterRunner.MaxClusterRow row: resultTable) {
			if(row.getIndex() <= numModels) {
				String file = row.getFileName().trim();
				url += file.substring(file.length()-1);
			}
			url += "|";
		}		
		
		// colors
		String color;
		url += "&chco=";		
		for (MaxClusterRunner.MaxClusterRow row: resultTable) {
			if(row.getIndex() <= numModels) {
				color = ourModelColor;
			} else {
				if(highlightGroup != null && row.getFileName().contains(highlightGroup)) {
					color = bakerColor;
				} else
					color = othersColor;
			}
			url += String.format("%s|", color);
		}
		url = url.substring(0, url.length()-1); // remove trailing comma
		
		return url;
	}
	
	public void writeAllTargetDetailFiles(File predDir, String groupSuffix, String[] targets, File outDir, boolean eval3D) {
		for(String target:targets) {
			int ret = this.getTargetResult(predDir, groupSuffix, target, eval3D);
			if(ret==NO_ERROR) {
				String fileName = target + ".txt";
				File outFile = new File(outDir, fileName);
				try {
					PrintStream out = new PrintStream(new FileOutputStream(outFile));
					this.printLastResultTable(out);
					out.close();
				} catch (IOException e) {
					System.err.println("Error creating file " + fileName + ":" + e.getMessage());
				}
			}
		}	
	}
		
}

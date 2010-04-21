package owl.core.connections;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Connection class to query the JPred server from the Barton group in Dundee for secondary structure-, burial-
 * and coiled coil prediction. For secondary structure prediction you might also want to consider {@link PsipredRunner}.
 * The main difference (apart from the different method) is that this connection obviously depends on the availability
 * of the online resource. PsipredConnection works offline but requires different installed executables and an up-to-date
 * sequence database file.
 * @author stehr
 */
public class JPredConnection {

	// constants
	static final String SUBMIT_URL	  = "http://www.compbio.dundee.ac.uk/www-jpred/cgi-bin/jpred_form?seq=%s&input=seq&pdb=on"; //&queryName=%s";
	static final String CHECK_URL 	  = "http://www.compbio.dundee.ac.uk/www-jpred/cgi-bin/chklog?%s";
	static final String CHECK_URL2 	  = "http://www.compbio.dundee.ac.uk/www-jpred/cgi-bin/chklog?keywords=%s";
	static final String RESULT_URL 	  = "http://www.compbio.dundee.ac.uk/www-jpred/results/jp_UX8aVue/jp_UX8aVue.concise";
	static final String REGEX_ERROR   = "(Error - \\w+)";
	static final String REGEX_STATUS  = "Job Status";
	static final String REGEX_STAT_1  = "Your job is next to be submitted";
	static final String REGEX_STAT_2  = "Your job [(]<b>\\w+</b>[)] started";
	static final String REGEX_JOB_ID  = "http://www.compbio.dundee.ac.uk/www-jpred/cgi-bin/chklog[?](\\w+)";
	static final String REGEX_RESULT  = "http://www.compbio.dundee.ac.uk/www-jpred/results/%s/%s.results.html";
	static final File TEMP_DIR 	  	  = new File(System.getProperty("java.io.tmpdir"));
	static final File DEBUG_FILE_1    = new File(TEMP_DIR, "jpred_submit.log");	// file created during online query in debug mode
	static final File DEBUG_FILE_2    = new File(TEMP_DIR, "jpred_status.log"); // file created during online query in debug mode
	static final File DEBUG_FILE_3    = new File(TEMP_DIR, "jpred_result.log"); // file created during online query in debug mode
	
	static final int TIMEOUT = 60;		// waiting for results timeout in seconds
	static final int INTERVAL = 10; 	// result checking interval in seconds
	
	// members
	HashMap<String,String[]> resultMap;		// stores the results of the last query
	boolean debug;						// if true, output debug information to stdout
	
	/**
	 * Create a new JPredConnection.
	 */
	public JPredConnection() {
		resultMap = null;
		debug = false;
		
	}
	
	/**
	 * Sets the debug mode flag to the given value.
	 * In debug mode, some status information is printed to stdout and
	 * log files will not be deleted on exit.
	 * @param b
	 */
	public void setDebugMode(boolean b) {
		this.debug = b;
	}
	
	/**
	 * Loads precomputed results from a text file in the JPred specific 'concise' format overwriting any
	 * previously loaded results. The file is not checked for the correct format. If the format is incorrect,
	 * the results of the method call are undefined. This is mainly used for testing with result files
	 * downloaded from JPred, so that the server does not need to be queried every time during debugging.
	 * It can also be used to retrieve results saved using local method saveToConciseFile.
	 * @param inFile the concise file to be read from
	 * @throws IOException if an i/o error occured during reading
	 */
	public void loadFromConciseFile(File inFile) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(inFile));
		String line;
		resultMap = new LinkedHashMap<String,String[]>();
		while((line = in.readLine()) != null) {
			String[] fields = line.split(":");
			String resulTag = fields[0];
			String[] resultSeq = fields[1].split(",");
			resultMap.put(resulTag, resultSeq);
		}
		in.close();
	}
	
	/**
	 * Stores the currently loaded results (e.g. obtained by submitQuery()) to a
	 * text file in the JPred specific concise format.
	 * @param outFile the file to be written to
	 * @throws IOException if an i/o error occured during writing or no results have been loaded
	 */
	public void saveToConciseFile(File outFile) throws IOException {
		if(resultMap == null) {
			throw new IOException("JPred connection error: Trying to save when no results have been loaded.");
		}
		FileWriter out = new FileWriter(outFile);
		for(String tag:resultMap.keySet()) {
			out.write(tag + ":");
			for(String s: resultMap.get(tag)) {
				out.write(s + ",");
			}
			out.write("\n");
		}
		out.close();
	}
	
	/**
	 * Helper method for printResults for printing a single result line
	 * @param key
	 */
	private void printResultLine(String key) {
		System.out.println(">" + key);
		for(String s: resultMap.get(key)) {
			System.out.print(s);
		}
		System.out.println();
	}
	
	/**
	 * Prints the currently loaded results in Fasta format.
	 */
	public void printResults() {
		if(resultMap == null) {
			System.err.println("Error: no results loaded.");
			return;
		}

		printResultLine("align1;QUERY");
		printResultLine("jnetpred");
		printResultLine("JNETCONF");
		printResultLine("JNETSOL0");
		printResultLine("JNETSOL5");
		printResultLine("JNETSOL25");
		printResultLine("Lupas_14");
		printResultLine("Lupas_21");
		printResultLine("Lupas_28");
	}
	
	/**
	 * Submit a JPred query for the given sequence. JPred returns secondary structure-,
	 * burial and coiled coil prediction results.
	 * After successful completion the results are stored in this object and can be retrieved
	 * with the different get... methods.
	 * @param seq the sequence to submit to the server
	 * @throws IOException
	 */
	public void submitQuery(String seq) throws IOException {
		
		if(!debug) {
			DEBUG_FILE_1.deleteOnExit();
			DEBUG_FILE_2.deleteOnExit();
			DEBUG_FILE_3.deleteOnExit();
		}		
		
		// submit query
		URL submitUrl;
		submitUrl = new URL(String.format(SUBMIT_URL, seq));
		if(debug) System.out.println("Sending query: " + submitUrl);
		BufferedReader in = new BufferedReader(new InputStreamReader(submitUrl.openStream()));
		FileWriter out = null;
		if(debug) out = new FileWriter(DEBUG_FILE_1);
		
		// check output for Error or JobId		
		String line;
		Pattern pErr = Pattern.compile(REGEX_ERROR);
		Pattern pJobId = Pattern.compile(REGEX_JOB_ID);
		Matcher m;
		String error = null;
		String jobId = null;
		while((line = in.readLine()) != null) {
			if(debug) out.write(line + "\n");
			m = pErr.matcher(line);
			if(m.find()) {
				error = m.group(1);
			}
			m = pJobId.matcher(line);
			if(m.find()) {
				jobId = m.group(1);
			}
		}
		in.close();
		if(debug) out.close();
		if(error != null) {
			if(debug) System.out.println(error);
			throw new IOException("JPred " + error);
		} else
		if(jobId != null) {
			if(debug) System.out.println("JobId="+jobId);
		} else {
			if(debug) System.out.println("Unexpected output.");
			throw new IOException("JPred Error - Unexpected output");
		}
				
		// check status
		long startTime = System.currentTimeMillis();
		boolean timeout = false;
		boolean result = false;
		if(debug) out =	new FileWriter(DEBUG_FILE_2);
		while(!timeout && !result) {
			
			// wait a moment
			try {
				Thread.sleep(1000 * INTERVAL);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			
			// check job status
			URL checkUrl = new URL(String.format(CHECK_URL, jobId));
			if(debug) System.out.println("Checking job status...");
			BufferedReader in2 = new BufferedReader(new InputStreamReader(checkUrl.openStream()));
			String line2;
			Pattern pStat1 = Pattern.compile(REGEX_STAT_1);
			Pattern pStat2 = Pattern.compile(REGEX_STAT_2);
			Pattern pResult = Pattern.compile(String.format(REGEX_RESULT, jobId, jobId));
			Matcher m2;
			while((line2 = in2.readLine()) != null) {
				if(debug) out.write(line2 + "\n");
				m2 = pStat1.matcher(line2);
				if(m2.find()) {
					if(debug) System.out.println("Job waiting");
				}
				m2 = pStat2.matcher(line2);
				if(m2.find()) {
					if(debug) System.out.println("Job running");
				}
				m2 = pResult.matcher(line2);
				if(m2.find()) {
					result = true;
				}
			}			
			in2.close();
			
			timeout = System.currentTimeMillis() - startTime > 1000 * TIMEOUT; 
		}
		if(debug) out.close();
		
		// parse results
		if(!result) {
			if(debug) System.out.println("Connection time out. Could not connect to JPred server.");
			throw new IOException("JPred connection timeout. Server is currently not available.");
		} else {
			if(debug) System.out.println("Job finished. Reading results...");
			URL resultUrl = new URL(String.format(RESULT_URL, jobId, jobId));
			BufferedReader in3 = new BufferedReader(new InputStreamReader(resultUrl.openStream()));
			resultMap = new LinkedHashMap<String,String[]>();
			if(debug) out = new FileWriter(DEBUG_FILE_3);
			String line3;
			while((line3 = in3.readLine()) != null) {
				if(debug) out.write(line3 + "\n");
				if(line3.length() > 0) {
					String[] fields = line3.split(":");
					String resulTag = fields[0];
					String[] resultSeq = fields[1].split(",");
					resultMap.put(resulTag, resultSeq);
				}
			}
			in3.close();		
			if(debug) out.close();
			if(debug) System.out.println("Result file written to " + DEBUG_FILE_3);
		}
	}

	/**
	 * Returns the query sequence stored in the results. It may be useful to compare this
	 * to the original query sequence passed to submitQuery to check for discrepancies.
	 * @return the query sequence for which results are stored or null if no results are loaded.
	 */
	public String getQuerySequence() {
		return resultMap==null?null:compressString(resultMap.get("align1;QUERY"));
	}
	
	/**
	 * Returns the JNet secondary structure prediction for the query sequence using
	 * symbols 'H' for helix and 'E' for extended conformation.
	 * @return the JNet secondary structure prediction as a string or null if no results are loaded.
	 */
	public String getSecondaryStructurePrediction() {
		return resultMap==null?null:compressString(resultMap.get("jnetpred"));
	}
	
	/**
	 * Returns confidence values from 0 to 9 for each position of the JNet secondary structure prediction.
	 * @return the JNet secondary structure prediction confidence values as a string or null if no results are loaded.
	 */
	public String getSecondaryStructureConfidence() {
		return resultMap==null?null:compressString(resultMap.get("JNETCONF"));
	}

	/**
	 * Returns the JNet burial prediction for the query sequence using symbols 'B' for buried
	 * at 0% relative surface accessibility and '-' for exposed.
	 * @return the JNet burial prediction as a string or null if no results are loaded.
	 */
	public String getBurialPrediction0() {
		return resultMap==null?null:compressString(resultMap.get("JNETSOL0"));
	}
	
	/**
	 * Returns the JNet burial prediction for the query sequence using symbols 'B' for buried
	 * at 5% relative surface accessibility and '-' for exposed.
	 * @return the JNet burial prediction as a string or null if no results are loaded.
	 */
	public String getBurialPrediction5() {
		return resultMap==null?null:compressString(resultMap.get("JNETSOL5"));
	}
	
	/**
	 * Returns the JNet burial prediction for the query sequence using symbols 'B' for buried
	 * at 25% relative surface accessibility and '-' for exposed.
	 * @return the JNet burial prediction as a string or null if no results are loaded.
	 */
	public String getBurialPrediction25() {
		return resultMap==null?null:compressString(resultMap.get("JNETSOL25"));
	}
	
	/**
	 * Returns the Lupas coiled coil prediction for the query sequence using a window length of 14 and symbols
	 * 'C' for high confidence coiled coil, 'c' for low confidence coiled coil and '-' for no coiled coil.
	 * @return the Lupas coiled coil prediction as a string or null if no results are loaded.
	 */
	public String getCoiledCoilPrediction14() {
		return resultMap==null?null:compressString(resultMap.get("Lupas_14"));
	}
	
	/**
	 * Returns the Lupas coiled coil prediction for the query sequence using a window length of 21 and symbols
	 * 'C' for high confidence coiled coil, 'c' for low confidence coiled coil and '-' for no coiled coil.
	 * @return the Lupas coiled coil prediction as a string or null if no results are loaded.
	 */
	public String getCoiledCoilPrediction21() {
		return resultMap==null?null:compressString(resultMap.get("Lupas_21"));
	}
	
	/**
	 * Returns the Lupas coiled coil prediction for the query sequence using a window length of 28 and symbols
	 * 'C' for high confidence coiled coil, 'c' for low confidence coiled coil and '-' for no coiled coil.
	 * @return the Lupas coiled coil prediction as a string or null if no results are loaded.
	 */
	public String getCoiledCoilPrediction28() {
		return resultMap==null?null:compressString(resultMap.get("Lupas_28"));
	}
	
	/**
	 * Helper function to convert a string array which contains only char entries to a string.
	 * For every entry the first character is taken such that if the input array contains only
	 * single letter strings, the desired result is achieved. Otherwise the result may not be
	 * useful.
	 */
	private String compressString(String[] arr) {
		StringBuffer b = new StringBuffer(arr.length);
		for(String s:arr) {
			b.append(s.charAt(0));
		}
		return b.toString();
	}
	
	public static void main(String[] args) throws IOException {
		
		// submit request
		if(args.length < 1) {
			System.out.println("Usage: jpred <sequence>");
			System.exit(1);
		}
		
		// submit query
		String seq = args[0];
		JPredConnection conn = new JPredConnection();
		conn.setDebugMode(true);
		conn.submitQuery(seq);
		
		// print results
		conn.printResults();
		
		// test
		File testFile = new File("jpred.out");
		File testFile2 = new File("jpred.test");
		conn.saveToConciseFile(testFile);
		conn.loadFromConciseFile(testFile);
		conn.saveToConciseFile(testFile2);	// if everything is ok, then these two files and DEBUG_FILE_3 should be the same
		
	}
	

	
}

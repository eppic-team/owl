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
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.connections.JPredProgressRetriever.Status;
import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;

/**
 * Package:		owl.core.connections
 * Class: 		JPredConnection
 * Author:		Henning Stehr, stehr@molgen.mpg.de
 * Date:		21/Apr/2010
 * 
 * Connection class to query the JPred server of the Barton group in Dundee for secondary structure-, burial-
 * and coiled coil prediction. Provides methods to set query parameters, submit the query and retrieve the
 * results. Results of the last query are stored locally and can be obtained with the get... methods or saved
 * to a file. 
 * 
 * Provides a main method to submit a JPred query from the command line.
 *  
 * For secondary structure prediction see also {@link PsipredRunner}.
 * The main difference (apart from the different method) is that this connection obviously depends on the availability
 * of the online resource. PsipredConnection works offline but requires different installed executables and an up-to-date
 * sequence database file.
 * 
 * Changelog:
 * 2010/04/21 first created by HS
 * 2010/04/26 new method getSecondaryStructurePredictionObject()
 * 2010/07/05 new method setTimeout()
 */
public class JPredConnection {

	/*------------------------------ constants ------------------------------*/
	
	// server URLs for sending queries and retrieving results
	static final String SUBMIT_URL	  = "http://www.compbio.dundee.ac.uk/www-jpred/cgi-bin/jpred_form?seq=%s&input=seq&pdb=on";
	static final String CHECK_URL_1   = "http://www.compbio.dundee.ac.uk/www-jpred/cgi-bin/chklog?%s";
	static final String CHECK_URL_2   = "http://www.compbio.dundee.ac.uk/www-jpred/cgi-bin/chklog?keywords=%s";
	static final String RESULT_URL 	  = "http://www.compbio.dundee.ac.uk/www-jpred/results/%s/%s.concise";
	
	// regular expressions to parse result pages
	static final String REGEX_ERROR   = "(Error - \\w+)";
	static final String REGEX_STATUS  = "Job Status";
	static final String REGEX_STAT_1  = "Your job is next to be submitted";
	static final String REGEX_STAT_2  = "Your job [(]<b>\\w+</b>[)] started";
	static final String REGEX_JOB_ID  = "http://www.compbio.dundee.ac.uk/www-jpred/cgi-bin/chklog[?](\\w+)";
	static final String REGEX_RESULT  = "http://www.compbio.dundee.ac.uk/www-jpred/results/%s/%s.results.html";
	
	// files created in debug mode
	static final File TEMP_DIR 	  	  = new File(System.getProperty("java.io.tmpdir"));
	static final File DEBUG_FILE_1    = new File(TEMP_DIR, "jpred_submit.log"); // server answer for submission
	static final File DEBUG_FILE_2    = new File(TEMP_DIR, "jpred_status.log"); // server answer for status query
	static final File DEBUG_FILE_3    = new File(TEMP_DIR, "jpred_result.log"); // server answer for result retrieval
	
	static final int TIMEOUT = 180;		// default waiting-for-results timeout in seconds, can be overriden with setTimeout()
	static final int INTERVAL = 10; 	// result checking interval in seconds
	
	public static final String SAMPLE_QUERY  = "APAFSVSPASGASDGQSVSVSVAAAGETYYIAQCAPVGGQDACNPATATSFTTDASGA";
	
	/*--------------------------- member variables --------------------------*/
	
	HashMap<String,String[]> resultMap;			// stores the results of the last query
	boolean debug;								// if true, output debug information to stdout
	JPredProgressRetriever progressRetriever;	// retriever for progress updates, may also be null
	JPredStopNotifier stopNotifier;				// holds a flag which indicates that executation should stop
												// this flag can be set by another thread (e.g. in response
												// to a cancel button being pressed and is queried by
												// submitQuery in regular intervals to stop execuation if requested
												// if this is null, it will be ignored
	private int timeout;						// the timeout after which connection to the server will be closed
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Creates a new JPredConnection. Use {@link submitQuery()} to send the query.
	 */
	public JPredConnection() {
		resultMap = null;
		debug = false;
		progressRetriever = null;
		stopNotifier = null;
		this.timeout = TIMEOUT;
	}
		
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Sets the debug mode flag to the given value.
	 * In debug mode, some status information is printed to stdout and
	 * log files will not be deleted on exit.
	 * @param b the value to be set
	 */
	public void setDebugMode(boolean b) {
		this.debug = b;
	}
	
	/**
	 * Sets the timeout after which connection to the server will be closed
	 * if we are still waiting for results. In this case an IOException will
	 * be thrown by submitQuery().
	 * @param seconds the timeout in seonds
	 */
	public void setTimeout(int seconds) {
		this.timeout = seconds;
	}
	
	/**
	 * Registers the progress retriever which will be notified of the status while the submitQuery method is running.
	 * @param r the progress retriever to register
	 */
	public void setProgressRetriever(JPredProgressRetriever r) {
		this.progressRetriever = r;
	}
	
	/**
	 * Registers the stop notifier which contains a flag indicating whether execution of submitQuery
	 * is supposed to stop.
	 * @param s
	 */
	public void setStopNotifier(JPredStopNotifier s) {
		this.stopNotifier = s;
	}
	
	/**
	 * Loads precomputed results from a text file in the JPred specific 'concise' format overwriting any
	 * previously loaded results. The file is not checked for the correct format. If the format is incorrect,
	 * the results of the method calls are undefined. This is mainly used for testing with result files
	 * downloaded from JPred, so that the server does not need to be queried every time during debugging.
	 * It can also be used to retrieve results saved using local method {@link saveToConciseFile()}.
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
	 * text file in the JPred specific 'concise' format.
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
	 * After successful completion, the results are stored in this object and can be retrieved
	 * with the different get... methods.
	 * @param seq the sequence to submit to the server
	 * @throws IOException
	 */
	public void submitQuery(String seq) throws IOException {
		
		if(progressRetriever != null) progressRetriever.setStatus(Status.START);
		
		if(!debug) {
			DEBUG_FILE_1.deleteOnExit();
			DEBUG_FILE_2.deleteOnExit();
			DEBUG_FILE_3.deleteOnExit();
		}		
		
		// submit query
		URL submitUrl;
		submitUrl = new URL(String.format(SUBMIT_URL, seq));
		if(debug) System.out.println("Sending query: " + submitUrl);
		if(progressRetriever != null) progressRetriever.setStatus(Status.SENDING);
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
			if(progressRetriever != null) progressRetriever.setStatus(Status.JOBID);
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
			
			// check whether thread was notified to stop
			if(stopNotifier != null && stopNotifier.getStop()) {
				throw new IOException("Execution was interrupted");
			}
			
			// wait a moment
			try {
				Thread.sleep(1000 * INTERVAL);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			
			// check job status
			URL checkUrl = new URL(String.format(CHECK_URL_1, jobId));
			if(debug) System.out.println("Checking job status...");
			if(progressRetriever != null) progressRetriever.setStatus(Status.CHECKING);
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
					if(progressRetriever != null) progressRetriever.setStatus(Status.WAITING);
				}
				m2 = pStat2.matcher(line2);
				if(m2.find()) {
					if(debug) System.out.println("Job running");
					if(progressRetriever != null) progressRetriever.setStatus(Status.RUNNING);
				}
				m2 = pResult.matcher(line2);
				if(m2.find()) {
					result = true;
				}
			}			
			in2.close();
			
			timeout = System.currentTimeMillis() - startTime > 1000 * this.timeout; 
		}
		if(debug) out.close();
		
		// parse results
		if(!result) {
			if(debug) System.out.println("Connection time out. Could not connect to JPred server.");
			throw new IOException("JPred connection timeout. Server is currently not available.");
		} else {
			if(debug) System.out.println("Job finished. Reading results...");
			if(progressRetriever != null) progressRetriever.setStatus(Status.RESULT);
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
			if(progressRetriever != null) progressRetriever.setStatus(Status.DONE);
			if(progressRetriever != null) progressRetriever.setResult(this.getSecondaryStructurePredictionObject());
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
	 * Returns a three state secondary structure object with the JNet prediction.
	 * @return a secondary structure object or null if no results were loaded.
	 */
	public SecondaryStructure getSecondaryStructurePredictionObject() {
		if(resultMap == null) return null;
		TreeMap<Integer, Character> ssTypes = new TreeMap<Integer,Character>();
		TreeMap<Integer, Double> ssConfs = new TreeMap<Integer,Double>();
		String ssPred = this.getSecondaryStructurePrediction();
		String ssConf = this.getSecondaryStructureConfidence();
		for(int i = 0; i < ssPred.length(); i++) {
			char ssType;
			switch(ssPred.charAt(i)) {
			case 'H' : ssType = SecStrucElement.HELIX; break;
			case 'E' : ssType = SecStrucElement.EXTENDED; break;
			default  : ssType = SecStrucElement.OTHER;
			}
			ssTypes.put(i+1, ssType);
			double confidence = Double.parseDouble("0."+String.valueOf(ssConf.charAt(i)));
			ssConfs.put(i+1, confidence);
		}
		SecondaryStructure secondaryStructure = new SecondaryStructure(this.getQuerySequence());
		char lastType = ssTypes.get(ssTypes.firstKey());
		int lastResSer = ssTypes.firstKey();
		char thisType;
		int start = 1;
		int elementCount = 0;
		SecStrucElement ssElem;
		String ssId;
		for(int resSer:ssTypes.keySet()) {
			thisType = ssTypes.get(resSer);
			if(thisType != lastType || resSer > lastResSer+1) {
				// finish previous element, start new one
				elementCount++;
				ssId = new Character(lastType).toString() + new Integer(elementCount).toString();
				ssElem = new SecStrucElement(lastType,start,lastResSer,ssId);
				secondaryStructure.add(ssElem);
				start = resSer;
				lastType = thisType;
			}
			lastResSer = resSer;
		}
		// finish last element
		elementCount++;
		ssId = new Character(lastType).toString() + new Integer(elementCount).toString();
		ssElem = new SecStrucElement(lastType, start,ssTypes.lastKey(),ssId);
		secondaryStructure.add(ssElem);

		secondaryStructure.setComment("JPred");
		
		return secondaryStructure;
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

	/*---------------------------- private methods --------------------------*/
	
	/**
	 * Helper function for printResults, prints a single result line for the given key
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
	 * Helper function to convert a string array which contains only char entries to a string.
	 * For every entry the first character is taken such that if the input array contains only
	 * single letter strings, the desired result is achieved. Otherwise the result may not be
	 * useful.
	 * @param arr the array of strings to be compressed
	 * @return the compressed string
	 */
	private String compressString(String[] arr) {
		StringBuffer b = new StringBuffer(arr.length);
		for(String s:arr) {
			b.append(s.charAt(0));
		}
		return b.toString();
	}
	
	/*--------------------------------- main --------------------------------*/
	
	/**
	 * Queries the JPred server for secondary structure-, burial- and coiled-coil
	 * predictions for the given protein sequence.
	 */
	public static void main(String[] args) {
		
		if(args.length < 1) {
			System.out.println("Usage: JPredConnection <sequence>");
			System.exit(1);
		}
	
		String seq = args[0];
		JPredConnection conn = new JPredConnection();
		conn.setDebugMode(true);
		try {
			
			conn.submitQuery(seq);
			conn.printResults();
			
		} catch (IOException e) {
			System.err.println(e.getMessage());
			System.exit(1);
		}
	}
}

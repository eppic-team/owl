package owl.core.connections;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.SQLException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.regex.*;

import owl.core.util.MySQLConnection;


/**
 * An interface for querying Prosite using their REST-like webservice.
 * See http://www.expasy.ch/tools/scanprosite/ScanPrositeREST.html
 * @author stehr
 */
public class PrositeScanner {

	/*------------------------------ constants ------------------------------*/
	
	public static final String 	PROSITE_URL = "http://www.expasy.org/cgi-bin/prosite/PSScan.cgi";
	public static final String 	QUERY_STR = PROSITE_URL + "?seq=%s&output=xml&skip=true";
	public static final String 	MATCH_PATTERN = "signature_ac";	// parsing hack: lines containing this pattern are considered hit lines
	
	public static final String 	REGEXP_START = "<start>(\\d+)</start>";
	public static final String 	REGEXP_STOP = "<stop>(\\d+)</stop>";
	public static final String 	REGEXP_SIGNATURE = "<signature_ac>(PS\\d+)</signature_ac>";
	public static final String 	REGEXP_SCORE = "<score>(\\d*\\.\\d+)</score>";
	public static final String 	REGEXP_LEVEL = "<level>(\\d+)</level>";
	
	/*--------------------------- member variables --------------------------*/
	
	LinkedList<String> latestResult;	// xml/text output of prosite query
	LinkedList<PrositeHit> latestHits;	// prosite hits parsed from output
	
	/*----------------------------- constructors ----------------------------*/
	
	public PrositeScanner() {
		latestResult = null;
		latestHits = null;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Searches for all prosite patterns/signatures in the given sequence.
	 * Results can be retrieved by retrieveLastResult() or printLatestHits().
	 * @param seq sequence in raw format or uniprotId to scan
	 */
	public void scanSeq(String seq) {
		try {
			URL url = new URL(String.format(QUERY_STR, seq));
			BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));
			this.latestResult = new LinkedList<String>();
			this.latestHits = new LinkedList<PrositeHit>();
			String line;
			while((line = in.readLine()) != null) {
				this.latestResult.add(line);
				// instead of parsing proper xml, we look for a pattern identifying a 'hit' line
				if(line.indexOf("signature_ac") > 0) {
					Pattern p;
					Matcher m;
					
					int start = PrositeHit.UNDEF_START;
					int stop = PrositeHit.UNDEF_STOP;
					String signatureAc = PrositeHit.UNDEF_SIGNATURE;
					double score = PrositeHit.UNDEF_SCORE;
					int level = PrositeHit.UNDEF_LEVEL;
					
					p = Pattern.compile(REGEXP_START);
					m = p.matcher(line);
					if(m.find()) {
						start = Integer.parseInt(m.group(1));
					}
					
					p = Pattern.compile(REGEXP_STOP);
					m = p.matcher(line);
					if(m.find()) {
						stop = Integer.parseInt(m.group(1));
					}	
					
					p = Pattern.compile(REGEXP_SIGNATURE);
					m = p.matcher(line);
					if(m.find()) {
						signatureAc = m.group(1);
					}
					
					p = Pattern.compile(REGEXP_SCORE);
					m = p.matcher(line);
					if(m.find()) {
						score = Double.parseDouble(m.group(1));
					}
					
					p = Pattern.compile(REGEXP_LEVEL);
					m = p.matcher(line);
					if(m.find()) {
						level = Integer.parseInt(m.group(1));
					}
					PrositeHit hit = new PrositeHit(start, stop, signatureAc, score, level);
					this.latestHits.add(hit);
				}
			}
			in.close();
		} catch (MalformedURLException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Returns the hits returned by the last Prosite search (performed by scanSeq)
	 * @return
	 */
	public Collection<PrositeHit> getLatestHits() {
		return this.latestHits;
	}
	
	/**
	 * Prints the output of the last Prosite search (which may contain an error message)
	 */
	public void printLatestOutput() {
		for(String s:this.latestResult) {
			System.out.println(s);
		}
	}
	
	/**
	 * Prints a table of the hits returned by the last Prosite search (performed by scanSeq)
	 */
	public void printLatesHits() {
		for(PrositeHit hit:this.latestHits) {
			System.out.print(hit);
			System.out.println(" " + getPatternDescription(hit.signatureAc));
		}
	}
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * TODO: This method is prone to SQL injection. Revise.
	 */
	public static String getPatternDescription(String prositeAc) {
		String sql = "SELECT descr FROM prosite.motif WHERE ac=\"%s\"";
		try {
			MySQLConnection conn = new MySQLConnection();
			String query = String.format(sql, prositeAc);
			String result = conn.getStringFromDb(query);
			conn.close();
			return result;
		} catch(SQLException e) {
			System.err.println("Error while reading from database: " + e.getMessage());
		}
		return null;
	}
	
	/*--------------------------------- main --------------------------------*/
	
	/**
	 * Example for using this class
	 */
	public static void main(String[] args) {
		String uniprotId = "ERBB2_HUMAN";		
		PrositeScanner ps = new PrositeScanner();
		System.out.println("Contacting Prosite...");
		ps.scanSeq(uniprotId);
		ps.printLatesHits();
		System.out.println("done.");	
	}
}

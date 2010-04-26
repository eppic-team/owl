package owl.core.connections;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.List;

import owl.core.sequence.Sequence;

/**
 * Connection class with static methods to download EMBL data with EMBL's DB fetch
 * web service.
 * See http://www.ebi.ac.uk/Tools/webservices/services/dbfetch_rest
 * 
 * Entries are retrieved with URLs of the form:
 * 
 * http://www.ebi.ac.uk/Tools/webservices/rest/dbfetch/{db}/{id}/{format}
 * 
 * 
 * @author duarte_j
 *
 */
public class EmblWSDBfetchConnection {
	
	private static final String BASE_URL = "http://www.ebi.ac.uk/Tools/webservices/rest/dbfetch";
	private static final int MAX_ENTRIES_PER_REQUEST = 200; // see the embl ws dbfetch docs (url above) 
	
	public enum Db {
		EMBLCDS("emblcds"),
		UNIPROTKB("uniprotkb");
		
		private String dbfetchStr;
		private Db(String dbfetchStr) {
			this.dbfetchStr = dbfetchStr;
		}
		public String getDBfetchStr() {
			return dbfetchStr;
		}
	}
	
	public enum Format {
		FASTA("fasta");
		
		private String dbfetchStr;
		private Format(String dbfetchStr) {
			this.dbfetchStr = dbfetchStr;
		}
		public String getDBfetchStr() {
			return dbfetchStr;
		}
	}
	
	/**
	 * Fetches a list of entries from EMBL DB fetch web service of the given db type  
	 * in FASTA format
	 * @param db
	 * @param ids
	 * @return
	 * @throws IOException
	 * @throws NoMatchFoundException
	 */
	private static List<Sequence> fetch(Db db, List<String> ids) throws IOException, NoMatchFoundException {
		List<Sequence> allSeqs = new ArrayList<Sequence>();
		// we do batches of MAX_ENTRIES_PER_REQUEST as the EMBL web service has a limit of entries per request
		for (int i=0;i<ids.size();i+=MAX_ENTRIES_PER_REQUEST) {
			String commaSepList = "";
			for (int c=i;c<i+MAX_ENTRIES_PER_REQUEST && c<ids.size();c++) {
				if (c!=i) commaSepList+=",";
				commaSepList+=ids.get(c);
			}
			allSeqs.addAll(fetch(db, commaSepList));
		}
		return allSeqs;
	}
	
	/**
	 * Retrieves the sequence data in FASTA format from EMBL DB fetch web service
	 * given a comma separated list of entry identifiers.
	 * @param db
	 * @param commaSepList
	 * @return
	 * @throws IOException
	 * @throws NoMatchFoundException
	 */
	private static List<Sequence> fetch(Db db, String commaSepList) throws IOException, NoMatchFoundException {
		URL url = new URL(BASE_URL+"/"+db.getDBfetchStr()+"/"+commaSepList+"/"+Format.FASTA.getDBfetchStr());
		URLConnection urlc = url.openConnection();
		InputStream is = urlc.getInputStream();
	    BufferedReader reader = new BufferedReader(new InputStreamReader(is));    
	    List<Sequence> out = new ArrayList<Sequence>();
		String 	nextLine = "",
		currentSeq = "",
		lastSeqTag = "";
	    while ((nextLine = reader.readLine())!=null) {
		    if (nextLine.startsWith("No entries found")) 
		    	throw new NoMatchFoundException("No "+db.getDBfetchStr()+" match found for ids "+commaSepList);
			nextLine = nextLine.trim();					    // remove whitespace
			if(nextLine.length() > 0) {						// ignore empty lines
				if (nextLine.startsWith(">")){
					if (!lastSeqTag.equals("")) {
						out.add(new Sequence(lastSeqTag,currentSeq));
						currentSeq = "";
					}
					lastSeqTag=nextLine.substring(1, nextLine.length()).trim();
				} else {
					currentSeq += nextLine;
				}
			}

	    }
	    is.close();
	    reader.close();
	    // adding last sequence (missed by loop above)
	    out.add(new Sequence(lastSeqTag,currentSeq));
	    
	    return out;		
	}
	
	/**
	 * Fetches EMBLCDS entries from EMBL DB fetch web service in FASTA format
	 * @param emblcdsIds
	 * @return
	 * @throws IOException
	 * @throws NoMatchFoundException
	 */
	public static List<Sequence> fetchEMBLCDS(List<String> emblcdsIds) throws IOException, NoMatchFoundException {
		return fetch(Db.EMBLCDS, emblcdsIds);
	}
	
	/**
	 * Fetches UNIPROT entries from EMBL DB fetch web service in FASTA format 
	 * @param uniprotIds
	 * @return
	 * @throws IOException
	 * @throws NoMatchFoundException
	 */
	public static List<Sequence> fetchUniprot(List<String> uniprotIds) throws IOException, NoMatchFoundException {
		return fetch(Db.UNIPROTKB, uniprotIds);
	}

	// testing
	public static void main(String[] args) throws Exception {
		System.out.println("Fetching EMBLCDS:");
		List<String> emblcdsids = new ArrayList<String>();
		emblcdsids.add("CAA84586");
		emblcdsids.add("ACB36185");
		List<Sequence> emblcdsSeqs = fetchEMBLCDS(emblcdsids);
		for (Sequence seq:emblcdsSeqs) {
			seq.writeToPrintStream(System.out);
		}
		System.out.println();
		System.out.println("Fetching a Uniprot:");
		List<String> uniIds = new ArrayList<String>();
		uniIds.add("P12830");
		List<Sequence> uniSeqs = fetchUniprot(uniIds);
		for (Sequence seq:uniSeqs) {
			seq.writeToPrintStream(System.out);
		}

	}
	
}

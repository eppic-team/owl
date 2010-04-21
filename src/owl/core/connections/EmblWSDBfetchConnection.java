package owl.core.connections;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;

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
	 * 
	 * @param db
	 * @param id
	 * @param format
	 * @return
	 * @throws IOException
	 * @throws NoMatchFoundException
	 */
	public static String fetch(Db db, String id, Format format) throws IOException, NoMatchFoundException {
		URL url = new URL(BASE_URL+"/"+db.getDBfetchStr()+"/"+id+"/"+format.getDBfetchStr());
		URLConnection urlc = url.openConnection();
		InputStream is = urlc.getInputStream();
	    BufferedReader reader = new BufferedReader(new InputStreamReader(is));
	    StringBuilder sb = new StringBuilder();
	    String line;
	    while ((line = reader.readLine())!=null) {
	      sb.append(line+"\n");
	    }
	    is.close();
	    reader.close();
	    String out = sb.toString();
	    if (out.startsWith("No entries found")) 
	    	throw new NoMatchFoundException("No "+db.getDBfetchStr()+" match found for id "+id);
	    return out;
	}
	
	/**
	 * Fetches a EMBLCDS entry from EMBL DB fetch web service in FASTA format
	 * @param emblcdsId
	 * @return
	 * @throws IOException
	 * @throws NoMatchFoundException
	 */
	public static String fetchEMBLCDS(String emblcdsId) throws IOException, NoMatchFoundException {
		return fetch(Db.EMBLCDS, emblcdsId, Format.FASTA);
	}
	
	/**
	 * Fetches a UNIPROT entry from EMBL DB fetch web service in FASTA format 
	 * @param uniprotId
	 * @return
	 * @throws IOException
	 * @throws NoMatchFoundException
	 */
	public static String fetchUniprot(String uniprotId) throws IOException, NoMatchFoundException {
		return fetch(Db.UNIPROTKB, uniprotId, Format.FASTA);
	}

	// testing
	public static void main(String[] args) throws Exception {
		System.out.println("Fetching an EMBLCDS:");
		System.out.println(fetchEMBLCDS("CAA84586"));
		System.out.println();
		System.out.println("Fetching a Uniprot:");
		System.out.println(fetchUniprot("P12830"));
	}
	
}

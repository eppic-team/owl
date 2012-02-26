package owl.core.connections;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import owl.core.sequence.UnirefEntry;
import owl.core.util.MySQLConnection;

/**
 * Class to get Uniprot Knowledge Base info from a local MySQL database. 
 * 
 * Two databases are needed: 
 * 1) a uniref database downloaded from 
 * ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/ 
 * or for achived ones from
 * ftp://ftp.uniprot.org/pub/databases/uniprot/previous_releases/
 * and then parsed with {@link UnirefXMLParser}, containing a minimal subset of the 
 * UniProt Knowledge Base.
 * The database must be named like uniprot_yyyy_mm or for the older style versions (e.g. 15.3)
 * the name convention must be uniprot_xx_y. The two tables of the database must be named "uniprot" 
 * and "uniprot_clusters"
 * 2) a taxonomy database as downloaded from http://www.uniprot.org/taxonomy
 * 
 * @see UniProtConnection
 * @see UnirefXMLParser
 * @author duarte_j
 *
 */
public class UniprotLocalConnection {

	private static final String DATA_TABLE = "uniprot"; 
	private static final String CLUSTERS_TABLE = "uniprot_clusters";
	private static final String TAX_TABLE = "taxonomy";
	
	private static final Log LOGGER = LogFactory.getLog(UniprotLocalConnection.class);
	
	private class TaxonomyRecord {
		
		@SuppressWarnings("unused")
		public String scientificName;
		public List<String> taxons;
		public TaxonomyRecord(String scientificName, List<String> taxons) {
			this.scientificName = scientificName;
			this.taxons = taxons;
		}
	}
	
	
	
	private MySQLConnection conn;
	private String dbName;
	private String taxonomyDbName;
	private String uniprotVer;
	
	private HashSet<String> nonReturnedIdsLastMultipleRequest;
	
	public UniprotLocalConnection(String dbName, String taxonomyDbName) throws SQLException {
		
		conn = new MySQLConnection();
		
		this.dbName = dbName;
		this.taxonomyDbName = taxonomyDbName;
		
		this.uniprotVer = dbName.substring(dbName.indexOf('_')+1, dbName.length());
		if (this.uniprotVer.length()<7) { // for old style (pre 2010) version numbers
			// version numbers are 7 characters (2010_01) but older style uniprot version numbers (15.3) are only 3 to 4 chars
			this.uniprotVer = uniprotVer.replace("_", "."); 
		}
		
	}
	
	public String getVersion() {
		return uniprotVer;
	}
	
	/**
	 * Given a uniprot or uniparc id returns the corresponding UnirefEntry record
	 * @param uniId
	 * @return
	 * @throws SQLException if something goes wrong while querying or if multiple matches are returned
	 * @throws NoMatchFoundException if no match found
	 */
	public UnirefEntry getUnirefEntry(String uniId) throws SQLException, NoMatchFoundException {
		String repId = null;
		String idColumn = null;
		if (uniId.startsWith("UPI")) {
			// uniparc id
			repId = uniId;
			idColumn = "uniparc_id";
		} else {
			repId = getRepresentative(uniId);
			idColumn = "uniprot_id";
		}
		
		Statement st = conn.createStatement();
		String sql = "SELECT id, uniprot_id, uniparc_id, tax_id, sequence FROM "+dbName+"."+DATA_TABLE+" WHERE "+idColumn+"='"+repId+"'";
		ResultSet rs = st.executeQuery(sql);
		UnirefEntry uniref = null;
		int count = 0;
		while (rs.next()) {
			uniref = new UnirefEntry();
			uniref.setId(rs.getString(1));
			uniref.setUniprotId(rs.getString(2));
			uniref.setUniparcId(rs.getString(3));
			uniref.setNcbiTaxId(rs.getInt(4));
			uniref.setSequence(rs.getString(5));
			count++;
		}
		rs.close();
		st.close();
		if (uniref==null) 
			throw new NoMatchFoundException("No match in table "+dbName+"."+DATA_TABLE+" for id "+uniId);
		if (count>1) 
			throw new SQLException("Multiple matches in table "+dbName+"."+DATA_TABLE+" for id "+repId);
		
		
		if (taxonomyDbName!=null) {
			TaxonomyRecord tax = getTaxonomy(uniref.getNcbiTaxId());
			if (tax!=null) {
				uniref.setTaxons(tax.taxons);
			} else {
				LOGGER.warn("No taxonomy information could be found for uniprot/uniparc id "+uniref.getUniId()+"(tax_id="+uniref.getNcbiTaxId()+")");
			}
		} else {
			LOGGER.warn("No taxonomy database specified, no taxonomy information will be available");
		}
		
		return uniref;
	}
	
	private String getRepresentative(String uniId) throws SQLException, NoMatchFoundException {
		Statement st = conn.createStatement();
		String sql = "SELECT representative FROM "+dbName+"."+CLUSTERS_TABLE+" WHERE member='"+uniId+"'";
		ResultSet rs = st.executeQuery(sql);
		String repUniId = null;
		int count = 0;
		while (rs.next()) {
			repUniId = rs.getString(1);
			count++;
		}
		rs.close();
		st.close();
		if (repUniId==null) 
			throw new NoMatchFoundException("No match in clusters table "+dbName+"."+CLUSTERS_TABLE+" for uniprot id "+uniId);
		if (count>1) 
			throw new SQLException("Multiple matches in clusters table "+dbName+"."+CLUSTERS_TABLE+" for uniprot id "+uniId);
			
		return repUniId;
	}
	
	/**
	 * Given a list of uniprot or uniparc ids returns the corresponding UnirefEntry records.
	 * If the query does not return all requested ids a warning is logged and the list of non-returned 
	 * ids can be retrieved through {@link #getNonReturnedIdsLastMultipleRequest()}  
	 * @param uniIds
	 * @return
	 * @throws SQLException
	 */
	public List<UnirefEntry> getMultipleUnirefEntries(List<String> uniIds) throws SQLException {
		nonReturnedIdsLastMultipleRequest = new HashSet<String>();
		
		List<UnirefEntry> entries = new ArrayList<UnirefEntry>();
		for (String uniId:uniIds) {			
			try {
				entries.add(getUnirefEntry(uniId));
			} catch (NoMatchFoundException e) {
				nonReturnedIdsLastMultipleRequest.add(uniId);
				LOGGER.warn("Information for uniprot/uniparc ID "+uniId+" could not be retrieved from local Uniprot.");
			}
			
		}
		return entries;
	}
	
	public HashSet<String> getNonReturnedIdsLastMultipleRequest() {
		return nonReturnedIdsLastMultipleRequest;
	}
	
	public TaxonomyRecord getTaxonomy(int taxId) throws SQLException {
		if (taxId==0) return null;
		
		Statement st = conn.createStatement();
		String sql = "SELECT scientific,lineage FROM "+taxonomyDbName+"."+TAX_TABLE+" WHERE tax_id="+taxId; 
		ResultSet rs = st.executeQuery(sql);

		String scientific = null;
		String lineage = null;
		while (rs.next()) {
			 scientific = rs.getString(1);
			 lineage = rs.getString(2);
		}
		rs.close();
		st.close();
		
		if (scientific==null) 
			return null;

		String[] taxons = lineage.split("; ");
		List<String> taxonsAL = Arrays.asList(taxons);
		
		return new TaxonomyRecord(scientific, taxonsAL);
	}
}

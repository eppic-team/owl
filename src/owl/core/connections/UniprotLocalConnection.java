package owl.core.connections;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import owl.core.sequence.UnirefEntry;
import owl.core.util.MySQLConnection;

/**
 * Class to get Uniprot Knowledge Base info from a local MySQL database
 * as parsed with {@link UnirefXMLParser} containing a minimal subset of the 
 * UniProt Knowledge Base.
 * 
 * @see UniProtConnection
 * @author duarte_j
 *
 */
public class UniprotLocalConnection {

	private static final String TABLES_PREFIX = "uniprot_"; 
	private static final Log LOGGER = LogFactory.getLog(UniprotLocalConnection.class);
	
	private MySQLConnection conn;
	private String dbName;
	private String uniprotTable;
	private String clustersTable;
	
	private HashSet<String> nonReturnedIdsLastMultipleRequest;
	
	public UniprotLocalConnection(String dbName, String uniprotVersion) throws SQLException {
		
		conn = new MySQLConnection();
		
		this.dbName = dbName;
		if (uniprotVersion.contains(".")) {
			// for older style uniprot version numbers
			uniprotVersion = uniprotVersion.replace(".", "_");
		}
		uniprotTable = TABLES_PREFIX+uniprotVersion;
		clustersTable = TABLES_PREFIX+uniprotVersion+"_clusters";
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
		String sql = "SELECT id, uniprot_id, uniparc_id, tax_id, sequence FROM "+dbName+"."+uniprotTable+" WHERE "+idColumn+"='"+repId+"'";
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
			throw new NoMatchFoundException("No match in table "+dbName+"."+uniprotTable+" for id "+uniId);
		if (count>1) 
			throw new SQLException("Multiple matches in table "+dbName+"."+uniprotTable+" for id "+repId);
		
		return uniref;
	}
	
	private String getRepresentative(String uniId) throws SQLException, NoMatchFoundException {
		Statement st = conn.createStatement();
		String sql = "SELECT representative FROM "+dbName+"."+clustersTable+" WHERE member='"+uniId+"'";
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
			throw new NoMatchFoundException("No match in clusters table "+dbName+"."+clustersTable+" for uniprot id "+uniId);
		if (count>1) 
			throw new SQLException("Multiple matches in clusters table "+dbName+"."+clustersTable+" for uniprot id "+uniId);
			
		return repUniId;
	}
	
	/**
	 * Given a list of uniprot or uniparc ids returns the corresponding UnirefEntry records.
	 * If the query does not return all requested ids a warning is logged and the list of non-returned 
	 * ids can be retrieved through {@link #getNonReturnedIdsLastMultipleRequest()}  
	 * @param uniprotIds
	 * @return
	 * @throws SQLException
	 */
	public List<UnirefEntry> getMultipleUnirefEntries(List<String> uniprotIds) throws SQLException {
		nonReturnedIdsLastMultipleRequest = new HashSet<String>();
		
		List<UnirefEntry> entries = new ArrayList<UnirefEntry>();
		for (String uniprotId:uniprotIds) {			
			try {
				entries.add(getUnirefEntry(uniprotId));
			} catch (NoMatchFoundException e) {
				nonReturnedIdsLastMultipleRequest.add(uniprotId);
				LOGGER.warn("Information for uniprot ID "+uniprotId+" could not be retrieved from local Uniprot.");
			}
			
		}
		return entries;
	}
	
	public HashSet<String> getNonReturnedIdsLastMultipleRequest() {
		return nonReturnedIdsLastMultipleRequest;
	}
}

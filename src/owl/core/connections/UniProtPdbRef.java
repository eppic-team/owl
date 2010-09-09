package owl.core.connections;

import java.sql.SQLException;

import owl.core.util.MySQLConnection;

/**
 * Simple storage class to hold a Uniprot x-reference to a PDB chain.
 * Also provides methods to write records of this type to database.
 * @author stehr
 */
public class UniProtPdbRef {
	
	/*--------------------------- member variables --------------------------*/
	
	// making all these public for easy access, sorry Niklaus Wirth ;)
	public String geneName;			// name of the gene, may be null
	public String uniprotId;		// uniprot ID, may be null
	public String pdbCode;			// pdb four letter code
	public String chainCode;		// pdb one letter chain code
	char type;						// X=xray, N=nmr, O=other 
	public double resolution;		// crystallographic resolution (in case of xray)
	public int begPos;				// starting position in uniprot sequence
	public int endPos;				// end position in uniprot sequence

	/*----------------------------- constructors ----------------------------*/
	
	public UniProtPdbRef(String geneName, String uniprotId, String pdbCode,
			String chainCode, char type, double resolution, int begPos, int endPos) {
		this.geneName = geneName;
		this.uniprotId = uniprotId;
		this.pdbCode = pdbCode;
		this.chainCode = chainCode;
		this.type = type;
		this.resolution = resolution;
		this.begPos = begPos;
		this.endPos = endPos;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Creates a database table to hold records of this class type.
	 * Hibernate can suck it!
	 * @throws SQLException 
	 */
	public void createDbTable(MySQLConnection conn, String table) throws SQLException {
		String sql = "CREATE TABLE %s (gene_name varchar(30), sp_id varchar(6), pos_beg int, pos_end int, pdb_code char(4), chain_code char(1), type char(1), resol float)";
		conn.executeSql(String.format(sql, table));
	}
	
	/**
	 * Writes this record to database.
	 * @param conn an active database connection
	 * @param table the target database and table
	 * @throws SQLException 
	 */
	public void writeToDb(MySQLConnection conn, String table) throws SQLException {
		String sql = "INSERT INTO %s (gene_name,sp_id,pos_beg,pos_end,pdb_code,chain_code,type,resol) VALUES (%s,%s,%d,%d,%s,%s,%f)";
		conn.executeSql(String.format(sql,table,this.geneName,this.uniprotId,this.begPos,this.endPos,this.pdbCode,this.chainCode,this.type,this.resolution));
	}
	
}

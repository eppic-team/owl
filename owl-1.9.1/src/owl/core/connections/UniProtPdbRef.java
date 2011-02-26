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
	public String[] chains;			// pdb one letter chain code(s)
	public String method;			// Xray, NMR, ... 
	public double resolution;		// crystallographic resolution (in case of xray) - otherwise Double.NaN
	public int begPos;				// starting position in uniprot sequence
	public int endPos;				// end position in uniprot sequence

	/*----------------------------- constructors ----------------------------*/
	
	public UniProtPdbRef(String geneName, String uniprotId, String pdbCode, 
			String[] chains, String method, double resolution, int begPos, int endPos) {
		this.geneName = geneName;
		this.uniprotId = uniprotId;
		this.pdbCode = pdbCode;
		this.chains = chains;
		this.method = method;
		this.resolution = resolution;
		this.begPos = begPos;
		this.endPos = endPos;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Writes this record to database. 
	 * Notes: 
	 * - In case of multiple chains, only the first one in written
	 * - Only the first character of the method is written (X=xray, N=nmr,...)
	 * @param conn an active database connection
	 * @param table the target database and table
	 * @throws SQLException 
	 */
	public void writeToDb(MySQLConnection conn, String table) throws SQLException {
		double resol = Double.isNaN(this.resolution)?-0:this.resolution;
		String sql = "INSERT INTO %s (gene_name,sp_id,pos_beg,pos_end,pdb_code,chain_code,type,resol) VALUES ('%s','%s',%d,%d,'%s','%s','%s',%.2f)";
		conn.executeSql(String.format(sql,table,this.geneName,this.uniprotId,this.begPos,this.endPos,this.pdbCode,this.chains[0].charAt(0),this.method.charAt(0),resol));
	}
	
	/**
	 * Returns a string represenation of this PDB cross reference.
	 */
	public String toString() {
		String chainsStr = "";
		for(String chain:chains) {
			chainsStr += chain;
			chainsStr += "/";
		}
		chainsStr = chainsStr.substring(0,chainsStr.length()-1);
		return String.format("%s\t%s\t%d-%d\t%s\t%s\t%s\t%.2f", this.geneName,this.uniprotId,this.begPos,this.endPos,this.pdbCode,chainsStr,this.method,this.resolution);
	}
	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Creates a database table to hold records of this class type (if not exists).
	 * Hibernate can suck it ;)
	 * @throws SQLException 
	 */
	public static void createDbTable(MySQLConnection conn, String table) throws SQLException {
		String sql = "CREATE TABLE IF NOT EXISTS %s (gene_name varchar(30), sp_id varchar(6), pos_beg int, pos_end int, pdb_code char(4), chain_code char(1), type char(1), resol float)";
		conn.executeSql(String.format(sql, table));
	}
	
}

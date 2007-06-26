package proteinstructure;

import tools.MySQLConnection;

import java.sql.Statement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;

public class MsdsdInfo {
	private final static String MYSQLSERVER="white";
	private final static String MYSQLUSER=getUserName();
	private final static String MYSQLPWD="nieve";
	private final static String DEFAULT_MYMSDSD_DB="my_msdsd_00_07_a";
	private final static String DEFAULT_MSDSD_DB="msdsd_00_07_a";
	
	private static final int DEFAULT_MODEL=1;

	private MySQLConnection conn;
	private String pdbCode="";
	private String pdbChainCode="";
	private int model=DEFAULT_MODEL;
	private int chainid;
	private int modelid;
	private String chainCode=""; // the internal msdsd/my_msdsd chain identifier (pchain_code)
	
	

	public MsdsdInfo (String pdbCode, String pdbChainCode, int model_serial, String db) throws MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException {
		this.pdbCode=pdbCode;
		this.pdbChainCode=pdbChainCode;
		this.model=model_serial;
		this.conn = new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD,db);
		getchainid();// initialises chainid and modelid
        if (check_inconsistent_res_numbering()){
            throw new MsdsdInconsistentResidueNumbersError("Inconsistent residue numbering in msdsd for accession_code "+this.pdbCode+", chain_pdb_code "+this.pdbChainCode);
        }

	}

	public MsdsdInfo (String pdbCode, String pdbChainCode, String db) throws MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException  {
		this(pdbCode,pdbChainCode,DEFAULT_MODEL,db);
	}
	public MsdsdInfo (String pdbCode, String pdbChainCode, int model_serial) throws MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException  {
		this(pdbCode,pdbChainCode,model_serial,DEFAULT_MSDSD_DB);
	}
	
	public MsdsdInfo (String pdbCode, String pdbChainCode) throws MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException  {
		this(pdbCode,pdbChainCode,DEFAULT_MODEL,DEFAULT_MSDSD_DB);
	}

	/** get user name from operating system (for use as database username) */
	private static String getUserName() {
		String user = null;
		user = System.getProperty("user.name");
		if(user == null) {
			System.err.println("Could not get user name from operating system. Exiting");
			System.exit(1);
		}
		return user;
	}
	
	public void close() {
		conn.close();
	}

	private void getchainid() throws MsdsdAcCodeNotFoundError {
		chainid=0;
		String chaincodestr="='"+pdbChainCode+"'";
		if (pdbChainCode.equals("NULL")){
			chaincodestr="IS NULL";
		}
		String sql = "SELECT chain_id, model_id, pchain_code " +
				" FROM "+DEFAULT_MYMSDSD_DB+".mmol_chain_info " +
				" WHERE accession_code='"+pdbCode+"' " +
				" AND chain_pdb_code "+chaincodestr +
				" AND chain_type='C' " +
				" AND asu_chain=1 " +
				" AND model_serial="+model;
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			if (rsst.next()) {
				chainid = rsst.getInt(1);
				modelid = rsst.getInt(2);
				chainCode=rsst.getString(3);
				if (! rsst.isLast()) {
					System.err.println("More than 1 chain_id match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode);
					throw new MsdsdAcCodeNotFoundError("More than 1 chain_id match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode);					
				}
			} else {
				System.err.println("No chain_id match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode);
				throw new MsdsdAcCodeNotFoundError("No chain_id could be matched for accession_code "+pdbCode+", chain_pdb_code "+pdbChainCode);
			} 
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}

	}
	
	private boolean check_inconsistent_res_numbering(){
		int count=0;
		int numserial=0;
		try {
			String sql="SELECT count(*) " +
			" FROM "+DEFAULT_MYMSDSD_DB+".problem_serial_chain " +
			" WHERE chain_id="+chainid +
			" AND (min_serial!=1 OR num_serial!=num_dist_serial OR num_serial!=max_serial-min_serial+1)";
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			while (rsst.next()) {
				count = rsst.getInt(1);
				if (count>0){
					return true;
				}
			}
			sql="SELECT num_serial FROM "+DEFAULT_MYMSDSD_DB+".problem_serial_chain WHERE chain_id="+chainid;
			rsst = stmt.executeQuery(sql);
			int check = 0;
			while (rsst.next()){
				check++;
				numserial=rsst.getInt(1);
			}
			if (check!=1){
				System.err.println("No num_serial match or more than 1 match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode);
			}
			String allresseq = read_seq();
			if (allresseq.length()!=numserial){
				System.err.println("num_serial and length of all_res_seq don't match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode);
				return true;
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return false;
	}
	
	public ArrayList<ArrayList> read_atomData(){
		ArrayList<ArrayList> resultset = new ArrayList<ArrayList>();
		String sql = "SELECT serial,chem_atom_name,code_3_letter,residue_serial,x,y,z " +
				" FROM atom_data " +
				" WHERE	(model_id = "+modelid+") " +
				" AND (chain_id = "+chainid+") " +
				" AND (graph_alt_code_used = 1) " +
				" AND (graph_standard_aa=1) " +
				" AND (pdb_group = 'A')" +
				" ORDER BY chain_code, residue_serial, serial";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			int count=0;
			while (rsst.next()){
				count++;
				ArrayList<Comparable> thisrecord = new ArrayList<Comparable>();
				thisrecord.add(rsst.getInt(1)); //atomserial
				thisrecord.add(rsst.getString(2).trim()); // atom
				thisrecord.add(rsst.getString(3).trim()); // res_type
				thisrecord.add(rsst.getInt(4)); //res_serial
				thisrecord.add(rsst.getDouble(5)); // x
				thisrecord.add(rsst.getDouble(6)); // y
				thisrecord.add(rsst.getDouble(7)); // z
				
				resultset.add(thisrecord);
			}
			if (count==0){
				System.err.println("atom data query returned no data at all for model_id="+modelid+", model_id="+modelid);
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return resultset;
	}
	
	public String read_seq(){
		String allresseq="";
		String sql="SELECT all_res_seq FROM "+DEFAULT_MYMSDSD_DB+".chain_seq WHERE chain_id="+chainid;
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			int check = 0;
			if (rsst.next()) {
				check++;
				allresseq=rsst.getString(1);
			} 
			if (check!=1) {
				System.err.println("No all_res_seq match or more than 1 match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode+", chain_id="+chainid);
			} 
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}

		return allresseq;
	}
	
	public HashMap<String,Integer> get_ressers_mapping() {
		HashMap<String,Integer> map = new HashMap<String, Integer>();
		String sql="SELECT serial, concat(pdb_seq,IF(pdb_insert_code IS NULL,'',pdb_insert_code)) " +
				" FROM residue " +
				" WHERE chain_id="+chainid+
				" AND pdb_seq IS NOT NULL";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			int count=0;
			while (rsst.next()) {
				count++;
				int resser = rsst.getInt(1);
				String pdbresser = rsst.getString(2);
				map.put(pdbresser, resser);
			} 
			if (count==0) {
				System.err.println("No residue serials mapping data match for chain_id="+chainid);
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}

		return map;
	}
	
	public String getChainCode(){
		return this.chainCode;
	}


}

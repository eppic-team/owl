package proteinstructure;

import tools.MySQLConnection;

import java.sql.Statement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;

public class MsdsdInfo {
	public final static String MYSQLSERVER="white";
	public final static String MYSQLUSER=getUserName();
	public final static String MYSQLPWD="nieve";
	public final static String mymsdsdDB="my_msdsd_00_07_a";
	public final static String msdsdDB="msdsd_00_07_a";
	public final static String pdbaseDB="pdbase";

	MySQLConnection conn;
	String accode="";
	String chaincode="";
	int model=DEFAULT_MODEL;
	int chainid;
	int modelid;
	String chain=""; // the internal msdsd/my_msdsd chain identifier (pchain_code)
	
	static int DEFAULT_MODEL=1;

	MsdsdInfo (String accode, String chaincode, int model_serial, String db) throws MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError {
		this.accode=accode;
		this.chaincode=chaincode;
		this.model=model_serial;
		this.conn = new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD,db);
		getchainid();// initialises chainid and modelid
        if (check_inconsistent_res_numbering()){
            throw new MsdsdInconsistentResidueNumbersError("Inconsistent residue numbering in msdsd for accession_code "+this.accode+", chain_pdb_code "+this.chaincode);
        }

	}

	MsdsdInfo (String accode, String chaincode, String db) throws MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError  {
		this(accode,chaincode,DEFAULT_MODEL,db);
	}
	MsdsdInfo (String accode, String chaincode, int model_serial) throws MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError  {
		this(accode,chaincode,model_serial,msdsdDB);
	}
	
	MsdsdInfo (String accode, String chaincode) throws MsdsdAcCodeNotFoundError, MsdsdInconsistentResidueNumbersError  {
		this(accode,chaincode,DEFAULT_MODEL,msdsdDB);
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

	public void getchainid() throws MsdsdAcCodeNotFoundError {
		chainid=0;
		String chaincodestr="='"+chaincode+"'";
		if (chaincode.equals("NULL")){
			chaincodestr="IS NULL";
		}
		String sql = "SELECT chain_id, model_id, pchain_code " +
				" FROM "+mymsdsdDB+".mmol_chain_info " +
				" WHERE accession_code='"+accode+"' " +
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
				chain=rsst.getString(3);
				if (! rsst.isLast()) {
					System.err.println("More than 1 chain_id match for accession_code="+accode+", chain_pdb_code="+chaincode);
					throw new MsdsdAcCodeNotFoundError("More than 1 chain_id match for accession_code="+accode+", chain_pdb_code="+chaincode);					
				}
			} else {
				System.err.println("No chain_id match for accession_code="+accode+", chain_pdb_code="+chaincode);
				throw new MsdsdAcCodeNotFoundError("No chain_id could be matched for accession_code "+accode+", chain_pdb_code "+chaincode);
			} 
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}

	}
	
	public boolean check_inconsistent_res_numbering(){
		int count=0;
		int numserial=0;
		try {
			String sql="SELECT count(*) " +
			" FROM "+mymsdsdDB+".problem_serial_chain " +
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
			sql="SELECT num_serial FROM "+mymsdsdDB+".problem_serial_chain WHERE chain_id="+chainid;
			rsst = stmt.executeQuery(sql);
			int check = 0;
			while (rsst.next()){
				check++;
				numserial=rsst.getInt(1);
			}
			if (check!=1){
				System.err.println("No num_serial match or more than 1 match for accession_code="+accode+", chain_pdb_code="+chaincode);
			}
			String allresseq = read_seq();
			if (allresseq.length()!=numserial){
				System.err.println("num_serial and length of all_res_seq don't match for accession_code="+accode+", chain_pdb_code="+chaincode);
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
				ArrayList thisrecord = new ArrayList();
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
		String sql="SELECT all_res_seq FROM "+mymsdsdDB+".chain_seq WHERE chain_id="+chainid;
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			int check = 0;
			if (rsst.next()) {
				check++;
				allresseq=rsst.getString(1);
			} 
			if (check!=1) {
				System.err.println("No all_res_seq match or more than 1 match for accession_code="+accode+", chain_pdb_code="+chaincode+", chain_id="+chainid);
			} 
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}

		return allresseq;
	}
}

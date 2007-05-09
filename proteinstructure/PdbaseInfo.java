package proteinstructure;

import tools.MySQLConnection;
import java.sql.Statement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Collections;

public class PdbaseInfo {
	public final static String MYSQLSERVER="white";
	public final static String MYSQLUSER="duarte"; //TODO implement getuser method
	public final static String MYSQLPWD="nieve";
	public final static String mymsdsdDB="my_msdsd_00_07_a";
	public final static String msdsdDB="msdsd_00_07_a";
	public final static String pdbaseDB="pdbase";
	
	MySQLConnection conn;
	String accode="";
	String chaincode="";
	int model=1;
	int entrykey;
	String asymid;
	int entitykey;
	String alt_locs_sql_str;
	
	PdbaseInfo (String accode, String chaincode, int model_serial, String db) {
		this.accode=accode;
		this.chaincode=chaincode;
		this.model=model_serial;
		this.conn = new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD,db);
        this.entrykey=get_entry_key();
        this.asymid=get_asym_id();
        this.entitykey=get_entity_key();
        this.alt_locs_sql_str=get_atom_alt_locs();

	}
	
	PdbaseInfo (String accode, String chaincode, String db) {
		this(accode,chaincode,1,db);
	}
	PdbaseInfo (String accode, String chaincode, int model_serial) {
		this(accode,chaincode,model_serial,pdbaseDB);
	}
	
	PdbaseInfo (String accode, String chaincode) {
		this(accode,chaincode,1,pdbaseDB);
	}
	
	public void close() {
		conn.close();
	}
	
	public int get_entry_key() {
		String sql="SELECT entry_key FROM struct WHERE entry_id='"+accode.toUpperCase()+"'";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			if (rsst.next()) {
				entrykey = rsst.getInt(1);
				if (! rsst.isLast()) {
					System.err.println("More than 1 entry_key match for accession_code="+accode+", chain_pdb_code="+chaincode);
					//TODO implement exception class
					//throw PdbaseAcCodeNotFound;					
				}
			} else {
				System.err.println("No entry_key match for accession_code="+accode+", chain_pdb_code="+chaincode);
				//TODO implement exception class
				//throw PdbaseAcCodeNotFoundError;
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return entrykey;
	}
	
	public String get_asym_id(){
		String pdbstrandid=chaincode;
		if (chaincode.equals("NULL")){
			pdbstrandid="A";
		}
		String sql="SELECT asym_id " +
				" FROM pdbx_poly_seq_scheme " +
				" WHERE entry_key=" + entrykey +
				" AND pdb_strand_id='"+pdbstrandid+"' " +
				" LIMIT 1";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			if (rsst.next()) {
				asymid = rsst.getString(1);
			} else {
				System.err.println("No asym_id match for entry_key="+entrykey+", pdb_strand_id="+chaincode);
				// TODO implement exceoption class
				//throw PdbaseInconsistencyError;
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return asymid;	
	}
	
	public int get_entity_key() {
		String sql="SELECT entity_key " +
				" FROM struct_asym " +
				" WHERE entry_key="+ entrykey +
				" AND id='"+asymid+"'";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			if (rsst.next()) {
				entitykey = rsst.getInt(1);
				if (! rsst.isLast()) {
					System.err.println("More than 1 entity_key match for entry_key="+entrykey+", asym_id="+asymid);
					//TODO implement exception class
					//throw PdbaseInconsistencyError;					
				}
			} else {
				System.err.println("No entity_key match for entry_key="+entrykey+", asym_id="+asymid);
				//TODO implement exception class
				//throw PdbaseInconsistencyError;
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return entitykey;
	}
	
	public String get_atom_alt_locs(){
		ArrayList<String> alt_ids = new ArrayList<String>();
		HashMap<String,Integer> alt_ids2keys = new HashMap<String,Integer>();
		String alt_loc_field="label_alt_key";
		String sql="SELECT id, atom_sites_alt_key FROM atom_sites_alt WHERE entry_key="+entrykey;
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			int count=0;
			while (rsst.next()) {
				count++;
				alt_ids.add(rsst.getString(1));
				alt_ids2keys.put(rsst.getString(1), rsst.getInt(2));
			}
			if (count!=0){
				if ((! alt_ids.contains(".")) || alt_ids.indexOf(".")!=alt_ids.lastIndexOf(".")){ // second term is a way of finding out if there is more than 1 ocurrence of "." in the ArrayList 
					System.err.println("alt_codes exist for entry_key "+entrykey+" but there is either no default value '.' or more than 1 '.'. Something wrong with this entry_key or with "+pdbaseDB+" db!");
					//TODO implement Exception class
					//throw PdbaseInconsistencyError;
				}
				alt_ids.remove(".");
				Collections.sort(alt_ids);
				String lowest_alt_id = alt_ids.get(0);
				alt_locs_sql_str = "("+alt_loc_field+"="+alt_ids2keys.get(".")+" OR "+alt_loc_field+"="+alt_ids2keys.get(lowest_alt_id)+")";
			} else {
				alt_locs_sql_str=alt_loc_field+" IS NULL";
			} 

			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
		return alt_locs_sql_str;
	}
	
	public ArrayList<ArrayList> read_atomData() {
		ArrayList<ArrayList> resultset = new ArrayList<ArrayList>();
		String sql = "SELECT id, label_atom_id, label_comp_id, label_asym_id, label_seq_id, Cartn_x, Cartn_y, Cartn_z " +
				" FROM atom_site " +
				" WHERE entry_key="+entrykey +
				" AND label_asym_id='"+asymid+"' " +
				" AND label_entity_key="+ entitykey +
				" AND model_num="+ model +
				" AND "+alt_locs_sql_str;
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			//TODO throw exception if rsst is empty
			while (rsst.next()){
				ArrayList thisrecord = new ArrayList();
				thisrecord.add(rsst.getInt(1)); //atomserial
				thisrecord.add(rsst.getString(2).trim()); // atom
				thisrecord.add(rsst.getString(3).trim()); // res_type
				thisrecord.add(rsst.getString(4).trim()); // chain
				thisrecord.add(rsst.getInt(5)); //res_serial
				thisrecord.add(rsst.getDouble(6)); // x
				thisrecord.add(rsst.getDouble(7)); // y
				thisrecord.add(rsst.getDouble(8)); // z
				
				resultset.add(thisrecord);
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return resultset;
	}
	
	public String read_seq(){
		String sequence="";
		String pdbstrandid=chaincode;
		if (chaincode.equals("NULL")){
			pdbstrandid="A";
		}
        // we use seq_id+0 (implicitly converts to int) in ORDER BY because seq_id is varchar!!
        String sql="SELECT mon_id" +
        		" FROM pdbx_poly_seq_scheme " +
        		" WHERE entry_key=" + entrykey +
        		" AND asym_id='"+asymid+"' " +
        		" AND pdb_strand_id='"+pdbstrandid+"' " +
        		" ORDER BY seq_id+0";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			ArrayList<String> aalist=AA.aas();
			int count=0;
			while (rsst.next()) {
				count++;
				String res_type = rsst.getString(1);
				if (aalist.contains(res_type)){
					sequence+=AA.threeletter2oneletter(res_type);
				} else {
					sequence+="X";
				}
			} 
			if (count==0) {
				System.err.println("No sequence data match for entry_key="+entrykey+", asym_id="+asymid+", pdb_strand_id="+pdbstrandid);
				//TODO implement exception class
				//throw PdbaseInconsistencyError;
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}

		return sequence;
	}
}

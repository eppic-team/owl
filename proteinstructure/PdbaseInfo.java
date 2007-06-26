package proteinstructure;

import tools.MySQLConnection;

import java.sql.Statement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Collections;

public class PdbaseInfo {
	private final static String MYSQLSERVER="white";
	private final static String MYSQLUSER=getUserName();
	private final static String MYSQLPWD="nieve";
	public final static String DEFAULT_PDBASE_DB="pdbase";
	
	private static final int DEFAULT_MODEL=1;
	
	private MySQLConnection conn;
	private String pdbCode="";
	private String pdbChainCode="";
	private int model=DEFAULT_MODEL;
	private String chainCode=""; // the internal pdbase chain identifier (asym_id)
	private int entrykey;
	private String asymid;
	private int entitykey;
	private String alt_locs_sql_str;
	
	
	
	public PdbaseInfo (String pdbCode, String pdbChainCode, int model_serial, String db) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, SQLException{
		this.pdbCode=pdbCode;
		this.pdbChainCode=pdbChainCode;
		this.model=model_serial;
		this.conn = new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD,db);
        this.entrykey=get_entry_key();
        this.asymid=get_asym_id();
        this.entitykey=get_entity_key();
        this.alt_locs_sql_str=get_atom_alt_locs();

	}
	
	public PdbaseInfo (String pdbCode, String pdbChainCode, String db) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, SQLException {
		this(pdbCode,pdbChainCode,DEFAULT_MODEL,db);
	}
	public PdbaseInfo (String pdbCode, String pdbChainCode, int model_serial) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, SQLException {
		this(pdbCode,pdbChainCode,model_serial,DEFAULT_PDBASE_DB);
	}
	
	public PdbaseInfo (String pdbCode, String pdbChainCode) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, SQLException {
		this(pdbCode,pdbChainCode,DEFAULT_MODEL,DEFAULT_PDBASE_DB);
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
	
	private int get_entry_key() throws PdbaseAcCodeNotFoundError {
		String sql="SELECT entry_key FROM struct WHERE entry_id='"+pdbCode.toUpperCase()+"'";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			if (rsst.next()) {
				entrykey = rsst.getInt(1);
				if (! rsst.isLast()) {
					System.err.println("More than 1 entry_key match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode);
					throw new PdbaseAcCodeNotFoundError();					
				}
			} else {
				System.err.println("No entry_key match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode);
				throw new PdbaseAcCodeNotFoundError();
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return entrykey;
	}
	
	private String get_asym_id() throws PdbaseInconsistencyError {
		String pdbstrandid=pdbChainCode;
		if (pdbChainCode.equals("NULL")){
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
				System.err.println("No asym_id match for entry_key="+entrykey+", pdb_strand_id="+pdbChainCode);
				throw new PdbaseInconsistencyError("No asym_id match for entry_key="+entrykey+", pdb_strand_id="+pdbChainCode);
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		// we set the internal chain identifier chainCode from asymid
		chainCode = asymid;
		return asymid;	
	}
	
	private int get_entity_key() throws PdbaseInconsistencyError {
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
					throw new PdbaseInconsistencyError("More than 1 entity_key match for entry_key="+entrykey+", asym_id="+asymid);					
				}
			} else {
				System.err.println("No entity_key match for entry_key="+entrykey+", asym_id="+asymid);
				throw new PdbaseInconsistencyError("No entity_key match for entry_key="+entrykey+", asym_id="+asymid);
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return entitykey;
	}
	
	private String get_atom_alt_locs() throws PdbaseInconsistencyError{
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
					System.err.println("alt_codes exist for entry_key "+entrykey+" but there is either no default value '.' or more than 1 '.'. Something wrong with this entry_key or with "+DEFAULT_PDBASE_DB+" db!");
					throw new PdbaseInconsistencyError("alt_codes exist for entry_key "+entrykey+" but there is either no default value '.' or more than 1 '.'. Something wrong with this entry_key or with "+DEFAULT_PDBASE_DB+" db!");
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
	
	public ArrayList<ArrayList> read_atomData() throws PdbaseInconsistencyError{
		ArrayList<ArrayList> resultset = new ArrayList<ArrayList>();
		String sql = "SELECT id, label_atom_id, label_comp_id, label_seq_id, Cartn_x, Cartn_y, Cartn_z " +
				" FROM atom_site " +
				" WHERE entry_key="+entrykey +
				" AND label_asym_id='"+asymid+"' " +
				" AND label_entity_key="+ entitykey +
				" AND model_num="+ model +
				" AND "+alt_locs_sql_str;
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
				throw new PdbaseInconsistencyError("atom data query returned no data at all for entry_key="+entrykey+", asym_id="+asymid+", entity_key="+entitykey+", model_num="+model+", alt_locs_sql_str='"+alt_locs_sql_str+"'");
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return resultset;
	}
	
	public String read_seq() throws PdbaseInconsistencyError{
		String sequence="";
		String pdbstrandid=pdbChainCode;
		if (pdbChainCode.equals("NULL")){
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
				throw new PdbaseInconsistencyError("No sequence data match for entry_key="+entrykey+", asym_id="+asymid+", pdb_strand_id="+pdbstrandid);
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}

		return sequence;
	}
	
	public HashMap<String,Integer> get_ressers_mapping() throws PdbaseInconsistencyError{
		String pdbstrandid=pdbChainCode;
		if (pdbChainCode.equals("NULL")){
			pdbstrandid="A";
		}

		HashMap<String,Integer> map = new HashMap<String, Integer>();
		String sql="SELECT seq_id, concat(auth_seq_num,IF(pdb_ins_code='.','',pdb_ins_code))" +
					" FROM pdbx_poly_seq_scheme " +
					" WHERE entry_key=" + entrykey +
					" AND asym_id='"+asymid+"' " +
					" AND pdb_strand_id='"+pdbstrandid+"' " +
					" AND auth_seq_num!='?'" +
					" ORDER BY seq_id+0";
		try {
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			int count=0;
			while (rsst.next()) {
				count++;
				int resser = Integer.parseInt(rsst.getString(1));
				String pdbresser = rsst.getString(2);
				map.put(pdbresser, resser);
			} 
			if (count==0) {
				System.err.println("No residue serials mapping data match for entry_key="+entrykey+", asym_id="+asymid+", pdb_strand_id="+pdbstrandid);
				throw new PdbaseInconsistencyError("No residue serials mapping data match for entry_key="+entrykey+", asym_id="+asymid+", pdb_strand_id="+pdbstrandid);
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

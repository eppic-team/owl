package proteinstructure;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import tools.MySQLConnection;

/**
 * A single chain pdb protein structure loaded from a PDBASE database
 * See http://openmms.sdsc.edu/OpenMMS-1.5.1_Std/openmms/docs/guides/PDBase.html to know what PDBASE is
 * 
 * @author		Jose Duarte
 * Class:		PdbasePdb
 * Package:		proteinstructure
 */
public class PdbasePdb extends Pdb {

	private final static String MYSQLSERVER="white";
	private final static String MYSQLUSER=MySQLConnection.getUserName();
	private final static String MYSQLPWD="nieve";
	private final static String DEFAULT_PDBASE_DB="pdbase";

	private MySQLConnection conn;
	
	private int entrykey;
	private String asymid;
	private int entitykey;
	private String alt_locs_sql_str;

	/**
	 * Constructs Pdb object given pdb code and pdb chain code. 
	 * Model will be DEFAULT_MODEL
	 * MySQLConnection is taken from defaults in PdbasePdb class: MYSQLSERVER, MYSQLUSER, MYSQLPWD
	 * Database is taken from default pdbase database in PdbasePdb class: DEFAULT_PDBASE_DB
	 * @param pdbCode
	 * @param pdbChainCode
	 * @throws PdbaseInconsistencyError
	 * @throws PdbaseAcCodeNotFoundError
	 * @throws SQLException 
	 */
	public PdbasePdb (String pdbCode, String pdbChainCode) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, SQLException {
		this(pdbCode, pdbChainCode, DEFAULT_MODEL, DEFAULT_PDBASE_DB, new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD));
	}

	/**
	 * Constructs Pdb object given pdb code, pdb chain code, source db and a MySQLConnection 
	 * Model will be DEFAULT_MODEL
	 * The db must be a pdbase database
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param db
	 * @param conn
	 * @throws PdbaseInconsistencyError
	 * @throws PdbaseAcCodeNotFoundError 
	 * @throws SQLException 
	 */
	public PdbasePdb (String pdbCode, String pdbChainCode, String db, MySQLConnection conn) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, SQLException {
		this(pdbCode,pdbChainCode,DEFAULT_MODEL,db, conn);		
	}

	/**
	 * Constructs Pdb object given pdb code, pdb chain code and model serial.
	 * MySQLConnection is taken from defaults in PdbasePdb class: MYSQLSERVER, MYSQLUSER, MYSQLPWD
	 * Database is taken from default pdbase database in PdbasePdb class: DEFAULT_PDBASE_DB
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param model_serial
	 * @throws PdbaseInconsistencyError
	 * @throws PdbaseAcCodeNotFoundError
	 * @throws SQLException
	 */
	public PdbasePdb (String pdbCode, String pdbChainCode, int model_serial) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, SQLException {
		this(pdbCode, pdbChainCode, model_serial, DEFAULT_PDBASE_DB, new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD));
	}
	
	/**
	 * Constructs Pdb object given pdb code, pdb chain code, model serial, source db and a MySQLConnection.
	 * The db must be a pdbase database 
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param model_serial
	 * @param db
	 * @param conn
	 * @throws PdbaseInconsistencyError
	 * @throws PdbaseAcCodeNotFoundError 
	 * @throws SQLException 
	 */
	public PdbasePdb (String pdbCode, String pdbChainCode, int model_serial, String db, MySQLConnection conn) throws PdbaseInconsistencyError, PdbaseAcCodeNotFoundError, SQLException {
		this.pdbCode=pdbCode;
		this.pdbChainCode=pdbChainCode;
		this.model=model_serial;
		this.db=db;
		
		this.conn = conn;
        this.entrykey=get_entry_key();
        this.asymid=get_asym_id();		// sets asymid and chainCode
        this.entitykey=get_entity_key();
        this.alt_locs_sql_str=get_atom_alt_locs();
		
		this.chainCode = getChainCode();
		this.sequence = read_seq();
		this.pdbresser2resser = get_ressers_mapping();
    
		this.read_atomData(); // populates resser_atom2atomserial, resser2restype, atomser2coord, atomser2resser

		// we initialise resser2pdbresser from the pdbresser2resser HashMap
		this.resser2pdbresser = new HashMap<Integer, String>();
		for (String pdbresser:pdbresser2resser.keySet()){
			resser2pdbresser.put(pdbresser2resser.get(pdbresser), pdbresser);
		}
	}

	private int get_entry_key() throws PdbaseAcCodeNotFoundError, SQLException {
		String sql="SELECT entry_key FROM "+db+".struct WHERE entry_id='"+pdbCode.toUpperCase()+"'";
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
		return entrykey;
	}
	
	private String get_asym_id() throws PdbaseInconsistencyError, SQLException {
		String pdbstrandid=pdbChainCode;
		if (pdbChainCode.equals("NULL")){
			pdbstrandid="A";
		}
		String sql="SELECT asym_id " +
				" FROM "+db+".pdbx_poly_seq_scheme " +
				" WHERE entry_key=" + entrykey +
				" AND pdb_strand_id='"+pdbstrandid+"' " +
				" LIMIT 1";

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
		// we set the internal chain identifier chainCode from asymid
		chainCode = asymid;
		return asymid;	
	}
	
	private int get_entity_key() throws PdbaseInconsistencyError, SQLException {
		String sql="SELECT entity_key " +
				" FROM "+db+".struct_asym " +
				" WHERE entry_key="+ entrykey +
				" AND id='"+asymid+"'";

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
		return entitykey;
	}
	
	private String get_atom_alt_locs() throws PdbaseInconsistencyError, SQLException{
		ArrayList<String> alt_ids = new ArrayList<String>();
		HashMap<String,Integer> alt_ids2keys = new HashMap<String,Integer>();
		String alt_loc_field="label_alt_key";
		String sql="SELECT id, atom_sites_alt_key FROM "+db+".atom_sites_alt WHERE entry_key="+entrykey;

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

		return alt_locs_sql_str;
	}
	
	private void read_atomData() throws PdbaseInconsistencyError, SQLException{
		resser_atom2atomserial = new HashMap<String,Integer>();
		resser2restype = new HashMap<Integer,String>();
		atomser2coord = new HashMap<Integer,Double[]>();
		atomser2resser = new HashMap<Integer,Integer>();

		
		String sql = "SELECT id, label_atom_id, label_comp_id, label_seq_id, Cartn_x, Cartn_y, Cartn_z " +
				" FROM "+db+".atom_site " +
				" WHERE entry_key="+entrykey +
				" AND label_asym_id='"+asymid+"' " +
				" AND label_entity_key="+ entitykey +
				" AND model_num="+ model +
				" AND "+alt_locs_sql_str;

		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		int count=0;
		while (rsst.next()){
			count++;

			int atomserial = rsst.getInt(1); 			// atomserial
			String atom = rsst.getString(2).trim();  	// atom
			String res_type = rsst.getString(3).trim();	// res_type
			int res_serial = rsst.getInt(4);			// res_serial
			double x = rsst.getDouble(5);				// x
			double y = rsst.getDouble(6);				// y
			double z = rsst.getDouble(7);				// z
			Double[] coords = {x, y, z};
			ArrayList<String> aalist=AA.aas();
			if (aalist.contains(res_type)) {
				atomser2coord.put(atomserial, coords);
				atomser2resser.put(atomserial, res_serial);
				resser2restype.put(res_serial, res_type);
				ArrayList<String> atomlist = aas2atoms.get(res_type);
				if (atomlist.contains(atom)){
					resser_atom2atomserial.put(res_serial+"_"+atom, atomserial);
				}
			}

		}
		if (count==0){
			throw new PdbaseInconsistencyError("atom data query returned no data at all for entry_key="+entrykey+", asym_id="+asymid+", entity_key="+entitykey+", model_num="+model+", alt_locs_sql_str='"+alt_locs_sql_str+"'");
		}
		rsst.close();
		stmt.close();
	}
	
	private String read_seq() throws PdbaseInconsistencyError, SQLException{
		String sequence="";
		String pdbstrandid=pdbChainCode;
		if (pdbChainCode.equals("NULL")){
			pdbstrandid="A";
		}
        // we use seq_id+0 (implicitly converts to int) in ORDER BY because seq_id is varchar!!
        String sql="SELECT mon_id" +
        		" FROM "+db+".pdbx_poly_seq_scheme " +
        		" WHERE entry_key=" + entrykey +
        		" AND asym_id='"+asymid+"' " +
        		" AND pdb_strand_id='"+pdbstrandid+"' " +
        		" ORDER BY seq_id+0";

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
        		sequence+=NONSTANDARD_AA_LETTER;
        	}
        } 
        if (count==0) {
        	System.err.println("No sequence data match for entry_key="+entrykey+", asym_id="+asymid+", pdb_strand_id="+pdbstrandid);
        	throw new PdbaseInconsistencyError("No sequence data match for entry_key="+entrykey+", asym_id="+asymid+", pdb_strand_id="+pdbstrandid);
        }
        rsst.close();
        stmt.close();

		return sequence;
	}
	
	private HashMap<String,Integer> get_ressers_mapping() throws PdbaseInconsistencyError, SQLException{
		String pdbstrandid=pdbChainCode;
		if (pdbChainCode.equals("NULL")){
			pdbstrandid="A";
		}

		HashMap<String,Integer> map = new HashMap<String, Integer>();
		String sql="SELECT seq_id, concat(auth_seq_num,IF(pdb_ins_code='.','',pdb_ins_code))" +
					" FROM "+db+".pdbx_poly_seq_scheme " +
					" WHERE entry_key=" + entrykey +
					" AND asym_id='"+asymid+"' " +
					" AND pdb_strand_id='"+pdbstrandid+"' " +
					" AND auth_seq_num!='?'" +
					" ORDER BY seq_id+0";

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

		return map;
	}

}

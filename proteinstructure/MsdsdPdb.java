package proteinstructure;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

import tools.MySQLConnection;

/**
 * A single chain pdb protein structure loaded from a MSDSD database
 * See http://www.ebi.ac.uk/msd-srv/docs/dbdoc/refaindex.html to know what MSDSD is
 * 
 * @author		Jose Duarte
 * Class:		MsdsdPdb
 * Package:		proteinstructure
 */
public class MsdsdPdb extends Pdb {
	
	private final static String MYSQLSERVER="white";
	private final static String MYSQLUSER=MySQLConnection.getUserName();
	private final static String MYSQLPWD="nieve";
	//private final static String DEFAULT_MYMSDSD_DB="my_msdsd_00_07_a";
	private final static String DEFAULT_MSDSD_DB="msdsd_00_07_a";

	private MySQLConnection conn;
	
	private int chainid;
	private int modelid;

	// TODO for this to be able to be used by other people we need to do things without a myMsdsdDb (or also distribute our fixes database)
	private String myMsdsdDb; // our database with add-ons and fixes to msdsd

	/**
	 * Constructs Pdb object given pdb code and pdb chain code. 
	 * Model will be DEFAULT_MODEL
	 * MySQLConnection is taken from defaults in MsdsdPdb class: MYSQLSERVER, MYSQLUSER, MYSQLPWD
	 * Database is taken from default msdsd database in MsdsdPdb class: DEFAULT_MSDSD_DB
	 * @param pdbCode
	 * @param pdbChainCode
	 * @throws PdbCodeNotFoundError
	 * @throws MsdsdInconsistentResidueNumbersError
	 * @throws SQLException 
	 */
	public MsdsdPdb (String pdbCode, String pdbChainCode) throws PdbCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException {
		this(pdbCode,pdbChainCode,DEFAULT_MODEL,DEFAULT_MSDSD_DB,new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD));		
	}

	/**
	 * Constructs Pdb object given pdb code, pdb chain code, db and MySQLConnection
	 * Model will be DEFAULT_MODEL
	 * db must be a msdsd database
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param db
	 * @param conn
	 * @throws PdbCodeNotFoundError
	 * @throws MsdsdInconsistentResidueNumbersError 
	 * @throws SQLException 
	 */
	public MsdsdPdb (String pdbCode, String pdbChainCode, String db, MySQLConnection conn) throws PdbCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException {		
		this(pdbCode,pdbChainCode,DEFAULT_MODEL,db,conn);		
	}
	
	/**
	 * Constructs Pdb object given pdb code, pdb chain code and a model serial
	 * MySQLConnection is taken from defaults in MsdsdPdb class: MYSQLSERVER, MYSQLUSER, MYSQLPWD
	 * Database is taken from default msdsd database in MsdsdPdb class: DEFAULT_MSDSD_DB
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param model_serial
	 * @throws PdbCodeNotFoundError
	 * @throws MsdsdInconsistentResidueNumbersError
	 * @throws SQLException 
	 */
	public MsdsdPdb (String pdbCode, String pdbChainCode, int model_serial) throws PdbCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException {
		this(pdbCode,pdbChainCode,model_serial,DEFAULT_MSDSD_DB,new MySQLConnection(MYSQLSERVER,MYSQLUSER,MYSQLPWD));		
	}

	/**
	 * Constructs Pdb object given pdb code, pdb chain code, model serial, a source db and a MySQLConnection.
	 * db must be a msdsd database
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param model_serial
	 * @param db
	 * @param conn
	 * @throws PdbCodeNotFoundError
	 * @throws MsdsdInconsistentResidueNumbersError 
	 * @throws SQLException 
	 */
	public MsdsdPdb (String pdbCode, String pdbChainCode, int model_serial, String db, MySQLConnection conn) throws PdbCodeNotFoundError, MsdsdInconsistentResidueNumbersError, SQLException {
		this.pdbCode=pdbCode.toLowerCase();				// our convention: pdb codes are lower case
		this.pdbChainCode=pdbChainCode.toUpperCase();	// our convention: chain codes are upper case
		this.model=model_serial;
		this.db=db;
		this.myMsdsdDb="my_"+db;  // i.e. for db=msdsd_00_07_a then myMsdsdDb=my_msdsd_00_07_a

		this.conn = conn;
		
		this.getchainid();// initialises chainid, modelid and chainCode
        
		if (check_inconsistent_res_numbering()){
            throw new MsdsdInconsistentResidueNumbersError("Inconsistent residue numbering in msdsd for accession_code "+this.pdbCode+", chain_pdb_code "+this.pdbChainCode);
        }
		
		this.sequence = read_seq();
		this.pdbresser2resser = get_ressers_mapping();

		this.read_atomData();
		
		// we initialise resser2pdbresser from the pdbresser2resser HashMap
		this.resser2pdbresser = new HashMap<Integer, String>();
		for (String pdbresser:pdbresser2resser.keySet()){
			resser2pdbresser.put(pdbresser2resser.get(pdbresser), pdbresser);
		}
		
		this.readSecStructure();
	}

	private void getchainid() throws PdbCodeNotFoundError, SQLException {
		chainid=0;
		String chaincodestr="='"+pdbChainCode+"'";
		if (pdbChainCode.equals("NULL")){
			chaincodestr="IS NULL";
		}
		String sql = "SELECT chain_id, model_id, pchain_code " +
				" FROM "+myMsdsdDb+".mmol_chain_info " +
				" WHERE accession_code='"+pdbCode+"' " +
				" AND chain_pdb_code "+chaincodestr +
				" AND chain_type='C' " +
				" AND asu_chain=1 " +
				" AND model_serial="+model;

		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		if (rsst.next()) {
			chainid = rsst.getInt(1);
			modelid = rsst.getInt(2);
			chainCode=rsst.getString(3);
			if (! rsst.isLast()) {
				//System.err.println("More than 1 chain_id match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode);
				throw new PdbCodeNotFoundError("More than 1 chain_id match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode);					
			}
		} else {
			//System.err.println("No chain_id match for accession_code="+pdbCode+", chain_pdb_code="+pdbChainCode);
			throw new PdbCodeNotFoundError("No chain_id could be matched for accession_code "+pdbCode+", chain_pdb_code "+pdbChainCode);
		} 
		rsst.close();
		stmt.close();
	}
	
	private boolean check_inconsistent_res_numbering() throws SQLException{
		int count=0;
		int numserial=0;

		String sql="SELECT count(*) " +
		" FROM "+myMsdsdDb+".problem_serial_chain " +
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
		sql="SELECT num_serial FROM "+myMsdsdDb+".problem_serial_chain WHERE chain_id="+chainid;
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
		return false;
	}
	
	private void read_atomData() throws SQLException{
		resser_atom2atomserial = new HashMap<String,Integer>();
		resser2restype = new HashMap<Integer,String>();
		atomser2coord = new HashMap<Integer,Double[]>();
		atomser2resser = new HashMap<Integer,Integer>();

		String sql = "SELECT serial,chem_atom_name,code_3_letter,residue_serial,x,y,z " +
				" FROM "+db+".atom_data " +
				" WHERE	(model_id = "+modelid+") " +
				" AND (chain_id = "+chainid+") " +
				" AND (graph_alt_code_used = 1) " +
				" AND (graph_standard_aa=1) " +
				" AND (pdb_group = 'A')" +
				" ORDER BY chain_code, residue_serial, serial";

		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		int count=0;
		while (rsst.next()){
			count++;

			int atomserial = rsst.getInt(1); 			// atomserial
			String atom = rsst.getString(2).trim();		// atom
			String res_type = rsst.getString(3).trim(); // res_type
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
			System.err.println("atom data query returned no data at all for model_id="+modelid+", model_id="+modelid);
		}
		rsst.close();
		stmt.close();
	}
	
	private String read_seq() throws SQLException{
		String allresseq="";
		String sql="SELECT all_res_seq FROM "+myMsdsdDb+".chain_seq WHERE chain_id="+chainid;

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

		return allresseq;
	}
	
	private HashMap<String,Integer> get_ressers_mapping() throws SQLException {
		HashMap<String,Integer> map = new HashMap<String, Integer>();
		String sql="SELECT serial, concat(pdb_seq,IF(pdb_insert_code IS NULL,'',pdb_insert_code)) " +
				" FROM "+db+".residue " +
				" WHERE chain_id="+chainid+
				" AND pdb_seq IS NOT NULL";

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

		return map;
	}
	
	private void readSecStructure() throws SQLException{
		resser2secstruct = new HashMap<Integer, String>();
		
		//HELIX -- helix table
		String sql = "SELECT helix_serial, beg_residue_serial, end_residue_serial " +
				" FROM "+db+".helix " +
				" WHERE	(model_id = "+modelid+") " +
				" AND (chain_id = "+chainid+") ";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		int count=0;
		while (rsst.next()) {
			count++;
			int serial = rsst.getInt(1);
			int beg = rsst.getInt(2);
			int end =rsst.getInt(3);
			String ssId = "H"+serial;
			for (int i=beg;i<=end;i++){
				if (resser2secstruct.containsKey(i)){ // if already assigned we print a warning and then assign it
					System.err.println("Inconsistency in secondary structure assignment. " +
							"Residue "+i+" is getting reassigned from "+resser2secstruct.get(i)+" to "+ssId);
				}
				resser2secstruct.put(i,ssId);	
			}
		} 
		rsst.close();
		stmt.close();
		//SHEET -- strand table
		sql = "SELECT sheet_serial, strand_serial, strand_beg_residue_serial, strand_end_residue_serial " +
				" FROM "+db+".strand " +
				" WHERE	(model_id = "+modelid+") " +
				" AND (chain_id = "+chainid+") ";
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		// we store everything in these 2 maps to assign later to resser2secstruct based on our own ids (ids are not very consistent in msdsd)
		HashMap<Integer,int[]> strands2begEnd = new HashMap<Integer, int[]>(); 
		TreeMap<Integer,ArrayList<Integer>> sheets2strands = new TreeMap<Integer, ArrayList<Integer>>();
		count=0;
		while (rsst.next()) {
			count++;
			int sheetSerial = rsst.getInt(1);
			int strandSerial = rsst.getInt(2);
			int beg = rsst.getInt(3);
			int end =rsst.getInt(4);
			int[] begEnd = {beg, end};
			strands2begEnd.put(strandSerial, begEnd);
			if (sheets2strands.containsKey(sheetSerial)){
				sheets2strands.get(sheetSerial).add(strandSerial);
			} else {
				ArrayList<Integer> strands = new ArrayList<Integer>();
				strands.add(strandSerial);
				sheets2strands.put(sheetSerial, strands);
			}
		} 
		rsst.close();
		stmt.close();
		char sheet='A';
		for (int sheetSerial:sheets2strands.keySet()){
			int strand=1;
			for (int strandSerial:sheets2strands.get(sheetSerial)){
				int[] begEnd = strands2begEnd.get(strandSerial);
				for (int i=begEnd[0];i<=begEnd[1];i++){
					String ssId = "S"+sheet+strand;
					if (resser2secstruct.containsKey(i)){ // if already assigned we print a warning and then assign it
						System.err.println("Inconsistency in secondary structure assignment. " +
								"Residue "+i+" is getting reassigned from "+resser2secstruct.get(i)+" to "+ssId);
					} 
					resser2secstruct.put(i,ssId);
				}
				strand++;
			}
			sheet++;
		}

		//TURN -- turn table	
		// they forgot to fill up the turn_serial field so we have to use turn_id and get a serial from it that is unique within the chain only
		sql = "SELECT turn_id, res_1_residue_serial, res_2_residue_serial, res_3_residue_serial, res_4_residue_serial " +
				" FROM "+db+".turn " +
				" WHERE	(model_id = "+modelid+") " +
				" AND (chain_id = "+chainid+") ";
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		TreeMap<Integer,int[]> turns = new TreeMap<Integer, int[]>();
		count=0;
		while (rsst.next()) {
			count++;
			int dbId = rsst.getInt(1);
			int res1 = rsst.getInt(2);
			int res2 = rsst.getInt(3);
			int res3 = rsst.getInt(4);
			int res4 = rsst.getInt(5);
			int[] residues = {res1,res2,res3,res4};
			turns.put(dbId, residues);				
		} 
		rsst.close();
		stmt.close();
		int serial=1;
		for (int dbId:turns.keySet()){
			for (int i:turns.get(dbId)){
				if (i!=0){ // if res1, res2, res3 or res4 where null in database then they are 0 here, we skip them
					String ssId="T"+serial;
					if (resser2secstruct.containsKey(i)){ // if already assigned we print a warning and then assign it
						System.err.println("Inconsistency in secondary structure assignment. " +
								"Residue "+i+" is getting reassigned from "+resser2secstruct.get(i)+" to "+ssId);
					}
					resser2secstruct.put(i,ssId);
				}
			}
			serial++;
		}

	}
}


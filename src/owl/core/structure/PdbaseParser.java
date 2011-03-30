package owl.core.structure;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Point3d;

import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;
import owl.core.util.MySQLConnection;


/**
 * A PDB data reader from PDBASE database format. 
 * See http://openmms.sdsc.edu/OpenMMS-1.5.1_Std/openmms/docs/guides/PDBase.html 
 * 
 */
public class PdbaseParser {

	public final static String DEFAULT_PDBASE_DB="pdbase";

	private String db;				// the pdbase db from which we have taken the data
	private MySQLConnection conn;
	
	//public enum ChainCodeType {PDB_CHAIN_CODE, CIF_CHAIN_CODE};
	
	private int model;
	private String pdbChainCode;
	
	private String pdbCode;
	private int entrykey;
	private String asymid;
	private int entitykey;
	private String alt_locs_sql_str;
	
	private String[] chainsArray;
	private Integer[] modelsArray;


	/**
	 * Constructs a pdbase reader given pdb code, source db and a MySQLConnection.
	 * Data will be loaded from database upon call of load(pdbChainCode, modelSerial)  
	 * @param pdbCode
	 * @param db a pdbase database
	 * @param conn
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundException 
	 */
	public PdbaseParser (String pdbCode, String db, MySQLConnection conn) throws PdbCodeNotFoundException, SQLException {
		this.db=db;
		
		this.conn = conn;
		
		this.pdbCode = pdbCode;
		this.entrykey = getEntryKey(); 
	}

	/**
	 * Reads PDB data from pdbase given a chain code of type ccType.
	 * @param pdbChainCode
	 * @param modelSerial
	 * @throws PdbLoadException
	 */
	public PdbChain readChain(String pdbChainCode, int modelSerial) throws PdbLoadException {
		PdbChain pdb = new PdbChain();
		try {
			this.model = modelSerial;
			this.pdbChainCode=pdbChainCode;	// NOTE! pdb chain code are case sensitive!
			pdb.setPdbChainCode(pdbChainCode);
			
			this.asymid=getAsymId();		
			pdb.setChainCode(asymid);
 
			this.entitykey=getEntityKey();
			this.alt_locs_sql_str=getAtomAltLocs(pdb);
			
			pdb.setSequence(readSeq());
			
			TreeMap<String,Integer> pdbresser2resser = getRessersMapping(); 
			pdb.setPdbresser2resserMap(pdbresser2resser);
			// we initialise resser2pdbresser from the pdbresser2resser TreeMap
			TreeMap<Integer,String> resser2pdbresser = new TreeMap<Integer, String>();
			for (String pdbresser:pdbresser2resser.keySet()){
				resser2pdbresser.put(pdbresser2resser.get(pdbresser), pdbresser);
			}
			pdb.setResser2pdbresserMap(resser2pdbresser);

			// this needs the info in the pdbress2resser maps
			this.readAtomData(pdb);
			
			SecondaryStructure secondaryStructure = new SecondaryStructure(pdb.getSequence().getSeq());	// create empty secondary structure first to make sure object is not null
			readSecStructure(secondaryStructure);
			pdb.setSecondaryStructure(secondaryStructure);
			if(!secondaryStructure.isEmpty()) {
				secondaryStructure.setComment("Pdbase");
				pdb.initialiseResiduesSecStruct();
			}
			
			pdb.initialiseMaps();
			
		} catch (SQLException e) {
			throw new PdbLoadException(e);
		} 

		return pdb;
	}
	
	/**
	 * Returns all alphabetically sorted PDB chain codes for given pdbCode entry in pdbase
	 * It caches the result so that next time called no loading has to be done.  
	 * @param pdbCode
	 * @param db
	 * @param conn
	 * @return array with all pdb chain codes
	 * @throws PdbLoadException
	 */
	public String[] getChains() throws PdbLoadException {
		if (chainsArray!=null) {
			return chainsArray;
		}
		TreeSet<String> chains = new TreeSet<String>();
		try {
			String sql = "SELECT DISTINCT pdb_strand_id FROM "+db+".pdbx_poly_seq_scheme WHERE entry_key="+entrykey;
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			while (rsst.next()) {
				chains.add(rsst.getString(1));
			}
			rsst.close();
			stmt.close();
		} catch (SQLException e) {
			throw new PdbLoadException(e);
		}
		
		if (chains.isEmpty()) return null;
		
		chainsArray = new String[chains.size()];
		chains.toArray(chainsArray);
		return chainsArray;
	}
	
	/**
	 * Returns all model serials for this entry in pdbase
	 * It caches the result so that next time called no loading has to be done. 
	 * @param pdbCode
	 * @param db
	 * @param conn
	 * @return array with all model serials
	 * @throws PdbLoadException 
	 */
	public Integer[] getModels() throws PdbLoadException {
		if (modelsArray!=null) {
			return modelsArray;
		}
		TreeSet<Integer> models = new TreeSet<Integer>();
		try {
			String sql = "SELECT DISTINCT model_num FROM "+db+".atom_site WHERE entry_key="+entrykey;
			Statement stmt = conn.createStatement();
			ResultSet rsst = stmt.executeQuery(sql);
			while (rsst.next()) {
				models.add(rsst.getInt(1));
			}
			rsst.close();
			stmt.close();
	
		} catch (SQLException e) {
			throw new PdbLoadException(e);
		}
		
		if (models.isEmpty()) return null;		
		modelsArray = new Integer[models.size()];
		models.toArray(modelsArray);
		return modelsArray;
	}
	
	/**
	 * Finds the entry key for this structure in the database. 
	 * If not found, throws an exception.
	 * @param pdbCode
	 * @param db
	 * @param conn
	 * @return the entry key
	 * @throws PdbCodeNotFoundException
	 * @throws SQLException
	 */
	private int getEntryKey() throws PdbCodeNotFoundException, SQLException {
		int entrykey = -1;
		String sql="SELECT entry_key FROM "+db+".struct WHERE entry_id='"+pdbCode.toUpperCase()+"'";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		if (rsst.next()) {
			entrykey = rsst.getInt(1);
			if (! rsst.isLast()) {
				throw new PdbCodeNotFoundException("More than 1 entry_key match for accession_code="+pdbCode);					
			}
		} else {
			throw new PdbCodeNotFoundException("No entry_key match for accession_code="+pdbCode);
		}
		rsst.close();
		stmt.close();
		return entrykey;
	}
	
	/**
	 * Gets the asym_id given a pdb_strand_id (from field pdbChainCode) and entry_key (from field entrykey)
	 * @return
	 * @throws PdbLoadException
	 * @throws SQLException
	 */
	private String getAsymId() throws PdbLoadException, SQLException {
		String asymid = null;
		String pdbstrandid=pdbChainCode;
		if (pdbChainCode.equals(PdbAsymUnit.NULL_CHAIN_CODE)){
			pdbstrandid="A";
		}
		// NOTE: as pdbx_poly_seq_scheme contains a record per residue, this query returns many (identical) records
		// We use the 'LIMIT 1' simply to be able to easily catch the error of no matches with an if (rsst.next()), see below
		// NOTE2: pdb_strand_id case sensitive!
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
			//System.err.println("No asym_id match for entry_key="+entrykey+", pdb_strand_id="+pdbChainCode);
			throw new PdbLoadException("No asym_id match for entry_key="+entrykey+", pdb_strand_id="+pdbChainCode);
		}
		rsst.close();
		stmt.close();

		return asymid;	
	}

	/**
	 * Gets the pdb_strand_id given an asym_id (from field asymid) and entry_key (from field entrykey)
	 * This is useful to get the pdbChainCode given the cif chain code. We don't use it anymore, thus is now commented out.
	 * @return
	 * @throws PdbLoadException
	 * @throws SQLException
	 */
//	private String getPdbStrandId() throws PdbLoadException, SQLException {
//		String pdbstrandid = null;
//		// NOTE: as pdbx_poly_seq_scheme contains a record per residue, this query returns many (identical) records
//		// We use the 'LIMIT 1' simply to be able to easily catch the error of no matches with an if (rsst.next()), see below
//		String sql="SELECT pdb_strand_id " +
//				" FROM "+db+".pdbx_poly_seq_scheme " +
//				" WHERE entry_key=" + entrykey +
//				" AND asym_id='"+asymid+"' " +
//				" LIMIT 1";
//
//		Statement stmt = conn.createStatement();
//		ResultSet rsst = stmt.executeQuery(sql);
//		if (rsst.next()) {
//			pdbstrandid = rsst.getString(1);
//		} else {
//			throw new PdbLoadException("No pdb_strand_id match for entry_key="+entrykey+", asym_id="+asymid);
//		}
//		rsst.close();
//		stmt.close();
//
//		return pdbstrandid;	
//	}
	
	// NOTE: Entity key not really needed since there can be only one entity_key
	// per entry_key,asym_id combination 
	private int getEntityKey() throws PdbLoadException, SQLException {
		String sql="SELECT entity_key " +
				" FROM "+db+".struct_asym " +
				" WHERE entry_key="+ entrykey +
				" AND id='"+asymid+"'";

		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		if (rsst.next()) {
			entitykey = rsst.getInt(1);
			if (! rsst.isLast()) {
				throw new PdbLoadException("More than 1 entity_key match for entry_key="+entrykey+", asym_id="+asymid);					
			}
		} else {
			throw new PdbLoadException("No entity_key match for entry_key="+entrykey+", asym_id="+asymid);
		}
		rsst.close();
		stmt.close();
		return entitykey;
	}
	
	private String getAtomAltLocs(PdbChain pdb) throws PdbLoadException, SQLException{
		ArrayList<String> alt_ids = new ArrayList<String>();
		String alt_loc_field="label_alt_id";
		String sql = "SELECT DISTINCT " + alt_loc_field + 
					" FROM "+db+".atom_site " +
					" WHERE entry_key="+entrykey +
					" AND label_asym_id='"+asymid+"' " +
					" AND label_entity_key="+ entitykey +
					" AND model_num="+ model;
		
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		int count=0;
		while (rsst.next()) {
			count++;
			alt_ids.add(rsst.getString(1));
		}
		if (count>0){
			if (! alt_ids.contains(".")){ 
				throw new PdbLoadException("alt_codes exist for entry_key "+entrykey+" but there is no default value '.'. Something wrong with this entry_key or with "+DEFAULT_PDBASE_DB+" db!");
			}
			if (count==1) {
				alt_locs_sql_str = alt_loc_field+"='.'";
			} else {
				pdb.setHasAltCodes(true);
				alt_ids.remove(".");
				Collections.sort(alt_ids);
				String lowest_alt_id = alt_ids.get(0);
				alt_locs_sql_str = "("+alt_loc_field+"='.' OR "+alt_loc_field+"='"+lowest_alt_id+"')";
			}
		} else {
			throw new PdbLoadException("No records returned from atom_site table for entry_key="+entrykey+", entity_key="+entitykey+", asym_id="+asymid+", model_num="+model);
		} 

		rsst.close();
		stmt.close();

		return alt_locs_sql_str;
	}
	
	protected SpaceGroup readSpaceGroup() throws SQLException, PdbLoadException {
		String sql = "SELECT space_group_name_h_m " +
				" FROM "+db+".symmetry " +
				" WHERE entry_key="+entrykey;
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		
		SpaceGroup spaceGroup = null;
		String sg = null;
		if (rsst.next()) {
			sg = rsst.getString(1);
		}
		if (sg!=null) {
			// for some pdb entries (e.g. NMRs) there's no crystal information at all
			spaceGroup = SymoplibParser.getSpaceGroup(sg);
			if (spaceGroup==null) {
				throw new PdbLoadException("The space group found '"+sg+"' is not recognised as a standard space group");
			}
		}

		rsst.close();
		stmt.close();
		return spaceGroup;
	}
	
	protected CrystalCell readCrystalCell () throws SQLException {
		
		String sql = "SELECT cell_length_a, cell_length_b, cell_length_c, cell_angle_alpha, cell_angle_beta, cell_angle_gamma " +
			  " FROM "+db+".cell " +
			  " WHERE entry_key="+entrykey;
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);

		CrystalCell crystalCell = null;
		if (rsst.next()) {
			double a = rsst.getFloat(1);
			double b = rsst.getFloat(2);
			double c = rsst.getFloat(3);
			double alpha = rsst.getFloat(4);
			double beta = rsst.getFloat(5);
			double gamma = rsst.getFloat(6);
			crystalCell = new CrystalCell(a, b, c, alpha, beta, gamma);
		}
		rsst.close();
		stmt.close();
		return crystalCell;
	}
	
	protected String readExpMethod() throws SQLException {
		String expMethod = null;
		String sql = "SELECT method FROM "+db+".exptl WHERE entry_key="+entrykey;
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		if (rsst.next()) {
			expMethod = rsst.getString(1).trim();
		}
		rsst.close();
		stmt.close();
		return expMethod;
	}
	
	/**
	 * Returns an array of size 3 with the quality parameters for a crystal structure: resolution, rFree and rSym
	 * @return
	 * @throws SQLException
	 */
	protected double[] readQparams() throws SQLException {
		double[] qParams = {-1,-1,-1};
		// NOTE that in case of non-xray structures no records will be present in refine or reflns and nothing will be set
		
		String sql = "SELECT ls_d_res_high, ls_r_factor_r_free FROM "+db+".refine WHERE entry_key="+entrykey;
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		if (rsst.next()) {
			double resol = rsst.getFloat(1);
			if (resol<1000) { // for some reason there are a few 3.4e+38 in the db, which is basically a null
				qParams[0] = resol;
			}
			double val =rsst.getFloat(2);
			if (val<1000) { // for some reason there are lots of 3.4e+38 in the db, which is basically a null
				qParams[1] = val;
			}
		}
		
		// if both Rsym/Rmerge are present, we don't compare them but take the Rsym value to be 
		// the right one (there's not much consensus in the field as to what's the 
		// right thing to do anyway!)
		sql = "SELECT pdbx_Rsym_value, pdbx_Rmerge_I_obs FROM "+db+".reflns WHERE entry_key="+entrykey;
		rsst = stmt.executeQuery(sql);
		if (rsst.next()) {
			double rsymval = rsst.getFloat(1);
			double rmergeval = rsst.getFloat(2);
			if (rsymval<1000) {
				qParams[2]=rsymval;
			} else if (rmergeval<1000) {
				qParams[2]=rmergeval;
			}
		}
		rsst.close();
		stmt.close();
		return qParams;
	}
	
	protected String readTitle() throws SQLException {
		String sql="SELECT title FROM "+db+".struct WHERE entry_key="+entrykey+"";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		String title = "";
		if (rsst.next()) {
			title = rsst.getString(1);
		}
		rsst.close();
		stmt.close();
		return title;
	}
	
	private void readAtomData(PdbChain pdb) throws PdbLoadException, SQLException{
		// NOTE: label_entity_key not really needed since there can be only one entity_key
		// per entry_key,asym_id combination
		String sql = "SELECT id, label_atom_id, type_symbol_id, label_comp_id, label_seq_id, Cartn_x, Cartn_y, Cartn_z, occupancy, b_iso_or_equiv " +
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
			String element = rsst.getString(3).trim();  //element
			String res_type = rsst.getString(4).trim();	// res_type
			int res_serial = rsst.getInt(5);			// res_serial
			double x = rsst.getDouble(6);				// x
			double y = rsst.getDouble(7);				// y
			double z = rsst.getDouble(8);				// z
			Point3d coords = new Point3d(x, y, z);
			double occupancy = rsst.getDouble(9);       // occupancy
			double bfactor = rsst.getDouble(10);        // bfactor
			if (AminoAcid.isStandardAA(res_type)) {
				if (!pdb.containsResidue(res_serial)) { 
					Residue residue = new Residue(AminoAcid.getByThreeLetterCode(res_type), res_serial, pdb);
					residue.setPdbSerial(pdb.getPdbResSerFromResSer(res_serial));
					pdb.addResidue(residue);
				}
				if (AminoAcid.isValidAtomWithOXT(res_type,atom)){
					Residue residue = pdb.getResidue(res_serial);
					residue.addAtom(new Atom(atomserial, atom, element, coords, residue, occupancy, bfactor));
				}
			}

		}
		if (count==0){
			throw new PdbLoadException("Atom data query returned no data at all for entry_key="+entrykey+", asym_id="+asymid+", entity_key="+entitykey+", model_num="+model+", alt_locs_sql_str='"+alt_locs_sql_str+"'");
		}
		rsst.close();
		stmt.close();
	}
	
	private String readSeq() throws PdbLoadException, SQLException{
		String sequence="";

        // we use seq_id+0 (implicitly converts to int) in ORDER BY because seq_id is varchar!!
        String sql="SELECT mon_id" +
        		" FROM "+db+".pdbx_poly_seq_scheme " +
        		" WHERE entry_key=" + entrykey +
        		" AND asym_id='"+asymid+"' " +
        		" ORDER BY seq_id+0";

        Statement stmt = conn.createStatement();
        ResultSet rsst = stmt.executeQuery(sql);
        int count=0;
        while (rsst.next()) {
        	count++;
        	String res_type = rsst.getString(1);
        	if (AminoAcid.isStandardAA(res_type)){
        		sequence+=AminoAcid.three2one(res_type);
        	} else {
        		sequence+=AminoAcid.XXX.getOneLetterCode();
        	}
        } 
        if (count==0) {
        	throw new PdbLoadException("No sequence data match for entry_key="+entrykey+", asym_id="+asymid);
        }
        rsst.close();
        stmt.close();

		return sequence;
	}
	
	private TreeMap<String,Integer> getRessersMapping() throws PdbLoadException, SQLException{
		String pdbstrandid=pdbChainCode;
		if (pdbChainCode.equals(PdbAsymUnit.NULL_CHAIN_CODE)){
			pdbstrandid="A";
		}

		TreeMap<String,Integer> map = new TreeMap<String, Integer>();
		String sql="SELECT seq_id, concat(pdb_seq_num,IF(pdb_ins_code='.','',pdb_ins_code))" +
					" FROM "+db+".pdbx_poly_seq_scheme " +
					" WHERE entry_key=" + entrykey +
					" AND asym_id='"+asymid+"' " +
					" AND pdb_strand_id='"+pdbstrandid+"' " +
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
			throw new PdbLoadException("No residue serials mapping data match for entry_key="+entrykey+", asym_id="+asymid+", pdb_strand_id="+pdbstrandid);
		}
		rsst.close();
		stmt.close();

		return map;
	}

	private void readSecStructure(SecondaryStructure secondaryStructure) throws SQLException {
		
		// HELIX AND TURN -- struct_conf table
		String sql = "SELECT id,beg_label_seq_id,end_label_seq_id " +
				" FROM "+db+".struct_conf " +
				" WHERE entry_key="+entrykey+
				" AND beg_label_asym_id='"+asymid+"'";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		int count=0;
		while (rsst.next()) {
			count++;
			String id = rsst.getString(1).trim(); // id is either HELIX_Pnn or TURN_Pnn
			Pattern p = Pattern.compile("^(\\w).+_P(\\d)+$");
			Matcher m = p.matcher(id);
			String ssId="Unknown";
			if (m.find()){
				ssId = m.group(1)+m.group(2); // e.g.: Hnn (helices) or Tnn (turns) 				
			}
			int beg = rsst.getInt(2);
			int end =rsst.getInt(3);
			char ssType = SecStrucElement.OTHER;
			if(id.startsWith("H")) {
				ssType = SecStrucElement.HELIX;
			} else if(id.startsWith("T")) {
				ssType = SecStrucElement.TURN;
			} else {
				System.err.println("Unknown secondary structure type " + id + " encountered when reading from Pdbase. Skipping.");
			}
			if(ssType != SecStrucElement.OTHER) {
				SecStrucElement ssElem = new SecStrucElement(ssType, beg, end, ssId);
				secondaryStructure.add(ssElem);
			}
		} 
		rsst.close();
		stmt.close();
		
		// SHEET -- struct_sheet_range table
		sql = "SELECT sheet_id, id, beg_label_seq_id, end_label_seq_id " +
				" FROM "+db+".struct_sheet_range " +
				" WHERE entry_key="+entrykey+
				" AND beg_label_asym_id='"+asymid+"'";
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		count=0;
		while (rsst.next()) {
			count++;
			String sheetid = rsst.getString(1).trim();
			int id = rsst.getInt(2);
			int beg = rsst.getInt(3);
			int end =rsst.getInt(4);
			String ssId=SecStrucElement.STRAND+sheetid+id; // e.g.: SA1, SA2..., SB1, SB2,...
			SecStrucElement ssElem = new SecStrucElement(SecStrucElement.STRAND, beg, end, ssId);
			secondaryStructure.add(ssElem);
		} 
		rsst.close();
		stmt.close();

	}
	
	/**
	 * Gets a mapping of observed residue serials (numbered from 1 to last 
	 * observed residue, with no gaps) to internal residue serials (cif)
	 * This method is designed specifically to be used with GTGParser 
	 * NOTE: the mapping might not be 100% perfect: we are not sure whether unobserved residues 
	 * 		 always have auth_seq_num!='?'. So far it looks fine (tested with a few hundred PDB entries) 
	 * @return
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param conn
	 * @param pdbaseDb
	 * @throws SQLException
	 * @throws PdbCodeNotFoundException
	 * @throws PdbLoadException
	 */
	public static TreeMap<Integer,Integer> getObservedResMapping(String pdbCode, String pdbChainCode, MySQLConnection conn, String pdbaseDb) 
	throws SQLException, PdbCodeNotFoundException, PdbLoadException {
		
		// finding entry key
		int entrykey = -1;
		String sql="SELECT entry_key FROM "+pdbaseDb+".struct WHERE entry_id='"+pdbCode.toUpperCase()+"'";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		if (rsst.next()) {
			entrykey = rsst.getInt(1);
			if (! rsst.isLast()) {
				throw new PdbCodeNotFoundException("More than 1 entry_key match for accession_code="+pdbCode);					
			}
		} else {
			throw new PdbCodeNotFoundException("No entry_key match for accession_code="+pdbCode);
		}
		rsst.close();
		stmt.close();
		
		// finding asymid
		String asymid = null;
		String pdbstrandid=pdbChainCode;
		if (pdbChainCode.equals(PdbAsymUnit.NULL_CHAIN_CODE)){
			pdbstrandid="A";
		}
		// NOTE: as pdbx_poly_seq_scheme contains a record per residue, this query returns many (identical) records
		// We use the 'LIMIT 1' simply to be able to easily catch the error of no matches with an if (rsst.next()), see below
		// NOTE2: pdb_strand_id case sensitive!
		sql="SELECT asym_id " +
				" FROM "+pdbaseDb+".pdbx_poly_seq_scheme " +
				" WHERE entry_key=" + entrykey +
				" AND pdb_strand_id='"+pdbstrandid+"' " +
				" LIMIT 1";

		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		if (rsst.next()) {
			asymid = rsst.getString(1);
		} else {
			//System.err.println("No asym_id match for entry_key="+entrykey+", pdb_strand_id="+pdbChainCode);
			throw new PdbLoadException("No asym_id match for entry_key="+entrykey+", pdb_strand_id="+pdbChainCode);
		}
		rsst.close();
		stmt.close();

		// and finally getting mapping
		TreeMap<Integer,Integer> map = new TreeMap<Integer, Integer>();
		sql = "SELECT seq_id " +
					" FROM "+pdbaseDb+".pdbx_poly_seq_scheme"+
					" WHERE entry_key=" + entrykey +
					" AND asym_id='"+asymid+"' " +
					" AND auth_seq_num!='?' " + // this gives only observed residues
					" ORDER BY seq_id+0";
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		int serial = 0;
		while (rsst.next()) {
			serial++;
			map.put(serial, Integer.parseInt(rsst.getString(1)));
		}
		rsst.close();
		stmt.close();
		return map;
	}
	
	/**
	 * Returns the database name from which the PDB data has
	 * been read.
	 * @return
	 */
	public String getDb() {
		return this.db;
	}
}

package owl.core.structure;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
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
	
	private String pdbCode;
	private int entrykey;
	
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
	 * Reads PDB data from pdbase given a PDB chain code
	 * @param pdbAsymUnit
	 * @param modelSerial
	 * @throws PdbLoadException
	 */
	public void readChains(PdbAsymUnit pdbAsymUnit, int modelSerial) throws PdbLoadException {

		try {
 			
			readPdbxPolySeq(pdbAsymUnit);
			
			// this needs the info in the pdbress2resser maps
			this.readAtomSite(pdbAsymUnit,modelSerial);
			
			for (PdbChain pdb:pdbAsymUnit.getAllChains()) {
				pdb.initialiseMaps();
			}
			
			readSecStructure(pdbAsymUnit);
			
		} catch (SQLException e) {
			throw new PdbLoadException(e);
		} 

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
	
	private void readAtomSite(PdbAsymUnit pdbAsymUnit, int model) throws PdbLoadException, SQLException{
		
		AtomLineList atomLines = new AtomLineList();
		String sql = 
			"SELECT label_asym_id, label_alt_id, id, label_atom_id, type_symbol_id, label_comp_id, label_seq_id, auth_seq_id, Cartn_x, Cartn_y, Cartn_z, occupancy, b_iso_or_equiv, auth_asym_id " +
			" FROM "+db+".atom_site " +
			" WHERE entry_key="+entrykey +
			" AND model_num="+ model+
			" ORDER BY label_asym_id, id+0";

		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		while (rsst.next()){
			String labelAsymId = rsst.getString(1).trim();	// asym_id
			String labelAltId = rsst.getString(2).trim();	// altCode
			int atomserial = rsst.getInt(3); 			// atomserial
			String atom = rsst.getString(4).trim();  	// atom
			String element = rsst.getString(5).trim();  //element
			String res_type = rsst.getString(6).trim();	// res_type
			String field7 = rsst.getString(7);			// res_serial
			int resSerial = -1;
			if (!field7.equals(".")) {
				resSerial = Integer.parseInt(field7);
			}
			int pdbResSerial = rsst.getInt(8);			// pdb res serial (needed for het chains where res_serial is not set)
			double x = rsst.getDouble(9);				// x
			double y = rsst.getDouble(10);				// y
			double z = rsst.getDouble(11);				// z
			double occupancy = rsst.getDouble(12);      // occupancy
			double bfactor = rsst.getDouble(13);        // bfactor
			String authAsymId = rsst.getString(14);		// auth_asym_id
			
			if (!res_type.equals("HOH")) {
				// note we don't really use the insCode and nonPoly fields of AtomLine (we use them only in pdb file parser), we fill them with null and false
				atomLines.addAtomLine(new AtomLine(labelAsymId,labelAltId,atomserial,atom,element,res_type,resSerial,pdbResSerial,null,new Point3d(x,y,z),occupancy,bfactor,authAsymId,false,false));
			}
		}
		rsst.close();
		stmt.close();
		if (atomLines.isEmpty()){
			throw new PdbLoadException("Atom data query returned no data at all for entry_key="+entrykey+", model_num="+model);
		}
		
		String altLoc = atomLines.getAtomAltLoc();
		String lastChainCode = null;
		
		for (AtomLine atomLine:atomLines) {
			// we read only the alt locs we want
			if (altLoc!=null && !atomLine.labelAltId.equals(altLoc) && !atomLine.labelAltId.equals(".")) continue;
			
			if (lastChainCode!=null && !lastChainCode.equals(atomLine.labelAsymId)) {
				if (!pdbAsymUnit.containsChainCode(atomLine.labelAsymId)) {
					// in readPdbxPolySeq we already added the polymer chains, now we only need to add the missing ones: non-polymer chains
					PdbChain pdb = new PdbChain();
					pdb.setChainCode(atomLine.labelAsymId);
					pdb.setPdbChainCode(atomLine.authAsymId);
					pdbAsymUnit.setNonPolyChain(atomLine.labelAsymId,pdb);
					pdb.setIsNonPolyChain(true);
				}
			} 
			PdbChain pdb = pdbAsymUnit.getChainForChainCode(atomLine.labelAsymId);
			if ((atomLine.resSerial!=-1 && !pdb.containsResidue(atomLine.resSerial)) ||
					(atomLine.resSerial==-1 && !pdb.containsResidue(atomLine.pdbResSerial))) {
				Residue residue = null;
				if (AminoAcid.isStandardAA(atomLine.res_type)) {
					residue = new AaResidue(AminoAcid.getByThreeLetterCode(atomLine.res_type), atomLine.resSerial, pdb);
				} else if (Nucleotide.isStandardNuc(atomLine.res_type)) {
					Nucleotide nuc = Nucleotide.getByCode(atomLine.res_type);
					residue = new NucResidue(nuc,atomLine.resSerial,pdb);
				} else {
					// this check is valid for both protein or nucleotide sequences
					if (!pdb.isNonPolyChain() && pdb.getSequence().getSeq().charAt(atomLine.resSerial-1)!=AminoAcid.XXX.getOneLetterCode()) {
						throw new PdbLoadException("HET residue with residue serial "+atomLine.resSerial+" and type "+atomLine.res_type+" does not match an X in the SEQRES sequence");
					}
					if (pdb.isNonPolyChain()) {
						// in het chains the res serial is not set in cif, we've got to use the pdb res serial 
						residue = new HetResidue(atomLine.res_type,atomLine.pdbResSerial,pdb);

					} else {
						residue = new HetResidue(atomLine.res_type,atomLine.resSerial,pdb);
					}
				}
				if (pdb.isNonPolyChain()) {
					residue.setPdbSerial(String.valueOf(atomLine.pdbResSerial));
					residue.setSerial(atomLine.pdbResSerial);
				} else {
					residue.setPdbSerial(pdb.getPdbResSerFromResSer(atomLine.resSerial));
				}
				pdb.addResidue(residue);
			}

			Residue residue = null;
			if (pdb.isNonPolyChain()) {
				residue = pdb.getResidue(atomLine.pdbResSerial);
			} else {
				residue = pdb.getResidue(atomLine.resSerial);
			}

			residue.addAtom(new Atom(atomLine.atomserial, atomLine.atom, atomLine.element, atomLine.coords, residue, atomLine.occupancy, atomLine.bfactor));

			lastChainCode = atomLine.labelAsymId;
		}
		
		if (altLoc!=null) {
			for (PdbChain pdb:pdbAsymUnit.getAllChains()) {
				pdb.setHasAltCodes(true);
			}
		}

		for (PdbChain pdb:pdbAsymUnit.getPolyChains()) {
			if (pdb.getObsLength()>pdb.getFullLength()) {
				throw new PdbLoadException("Length of observed (atom_site) sequence longer than pdbx_poly_seq_scheme sequence for CIF chain "+pdb.getChainCode()+". Inconsistent PDB entry. Report to the PDB.");
			}
			for (Residue residue:pdb) {
				if (residue.getShortCode()!=pdb.getSequence().getSeq().charAt(residue.getSerial()-1)){
					throw new PdbLoadException("atom_site sequence does not match sequence in pdbx_poly_seq_scheme for CIF chain "+pdb.getChainCode()+". Inconsistent PDB entry. Report to the PDB.");
				}	
			}
		}

		
	}
	
	private void readPdbxPolySeq(PdbAsymUnit pdbAsymUnit) throws PdbLoadException, SQLException{
		
		PdbxPolySeqLineList list  = new PdbxPolySeqLineList();
        
		// we use seq_id+0 (implicitly converts to int) in ORDER BY because seq_id is varchar!!
		// we order by asym_id, thus our list will have that order (we then read atoms in the same way)
        String sql="SELECT asym_id, pdb_strand_id, mon_id, seq_id, pdb_seq_num, pdb_ins_code " +
        		" FROM "+db+".pdbx_poly_seq_scheme " +
        		" WHERE entry_key=" + entrykey +
        		" ORDER BY asym_id, seq_id+0"; 
		
        Statement stmt = conn.createStatement();
        ResultSet rsst = stmt.executeQuery(sql);

        while (rsst.next()) {
        	String asymId = rsst.getString(1).trim();
        	String pdbChainCode = rsst.getString(2).trim();
        	String res_type = rsst.getString(3).trim();
			int resser = Integer.parseInt(rsst.getString(4).trim());
			int pdb_seq_num = Integer.parseInt(rsst.getString(5).trim());
			String pdb_ins_code = rsst.getString(6).trim();
			
			list.add(new PdbxPolySeqLine(asymId, resser, res_type, pdb_seq_num, pdbChainCode, pdb_ins_code));
        }
        rsst.close();
        stmt.close();
        if (list.isEmpty()) {
        	throw new PdbLoadException("No sequence data match for entry_key="+entrykey);
        }
        
        ArrayList<PdbxPolySeqGroup> groups = new ArrayList<PdbxPolySeqGroup>();

        String lastAsymId = null;
        for (PdbxPolySeqLine line:list) {
        	if (lastAsymId==null || !lastAsymId.equals(line.asym_id)) {
        		PdbxPolySeqGroup group = new PdbxPolySeqGroup();
        		group.add(line);
        		groups.add(group);
        	} else {
        		PdbxPolySeqGroup group = groups.get(groups.size()-1);
        		group.add(line);
        	}
        	lastAsymId = line.asym_id;
        } 
        
        TreeMap<String,String> pdbchaincode2chaincode = new TreeMap<String,String>();
        pdbAsymUnit.setPdbchaincode2chaincode(pdbchaincode2chaincode);

        for (PdbxPolySeqGroup group:groups) {

        	pdbchaincode2chaincode.put(group.getPdbChainCode(),group.getChainCode());

        	PdbChain pdb = new PdbChain();
            pdb.setSequence(group.getSequence(), group.isProtein());
 
            pdb.setPdbresser2resserMap(group.getPdbresser2resserMap());
			pdb.setResser2pdbresserMap(group.getResser2pdbresserMap());

			pdb.setPdbChainCode(group.getPdbChainCode());			
			pdb.setChainCode(group.getChainCode());

			pdb.setIsNonPolyChain(false);
			
			pdbAsymUnit.setPolyChain(group.getChainCode(),pdb);
			
        }


	}

	private void readSecStructure(PdbAsymUnit pdbAsymUnit) throws SQLException {
		for (PdbChain pdb:pdbAsymUnit.getPolyChains()) {
			SecondaryStructure secondaryStructure = new SecondaryStructure(pdb.getSequence().getSeq());	// create empty secondary structure first to make sure object is not null
			pdb.setSecondaryStructure(secondaryStructure);
		}
		// HELIX AND TURN -- struct_conf table
		String sql = "SELECT id,beg_label_seq_id,end_label_seq_id,beg_label_asym_id " +
				" FROM "+db+".struct_conf " +
				" WHERE entry_key="+entrykey+
				" ORDER BY beg_label_asym_id";
		Statement stmt = conn.createStatement();
		ResultSet rsst = stmt.executeQuery(sql);
		int count=0;
		while (rsst.next()) {
			count++;
			String id = rsst.getString(1).trim(); // id is either HELIX_Pnn or TURN_Pnn
			if (id.startsWith("TURN")) continue; // we don't parse turns anymore as they were dropped from PDB files and anyway the annotation is VERY inconsistent
			Pattern p = Pattern.compile("^(\\w).+_P(\\d+)$");
			Matcher m = p.matcher(id);
			String ssId="Unknown";
			if (m.find()){
				ssId = m.group(1)+m.group(2); // e.g.: Hnn (helices) or Tnn (turns) 				
			}
			int beg = rsst.getInt(2);
			int end =rsst.getInt(3);
			String asymId = rsst.getString(4);
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
				pdbAsymUnit.getChainForChainCode(asymId).getSecondaryStructure().add(ssElem);
			}
		} 
		rsst.close();
		stmt.close();
		
		// SHEET -- struct_sheet_range table
		sql = "SELECT sheet_id, id, beg_label_seq_id, end_label_seq_id, beg_label_asym_id " +
				" FROM "+db+".struct_sheet_range " +
				" WHERE entry_key="+entrykey+
				" ORDER BY beg_label_asym_id";
		stmt = conn.createStatement();
		rsst = stmt.executeQuery(sql);
		count=0;
		while (rsst.next()) {
			count++;
			String sheetid = rsst.getString(1).trim();
			int id = rsst.getInt(2);
			int beg = rsst.getInt(3);
			int end =rsst.getInt(4);
			String asymId = rsst.getString(5);
			String ssId=SecStrucElement.STRAND+sheetid+id; // e.g.: SA1, SA2..., SB1, SB2,...
			SecStrucElement ssElem = new SecStrucElement(SecStrucElement.STRAND, beg, end, ssId);
			pdbAsymUnit.getChainForChainCode(asymId).getSecondaryStructure().add(ssElem);
		} 
		rsst.close();
		stmt.close();

		for (PdbChain pdb:pdbAsymUnit.getPolyChains()) {
			if(!pdb.getSecondaryStructure().isEmpty()) {
				pdb.getSecondaryStructure().setComment("Pdbase");
				pdb.initialiseResiduesSecStruct();
			}
		}
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
	
}

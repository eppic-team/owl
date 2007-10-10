package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;

import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import javax.vecmath.Point3d;


/**
 * A single chain pdb protein structure loaded from an mmCIF file or downloaded from the PDB FTP site 
 * 
 * @author		Jose Duarte
 * Class:		CiffilePdb
 * Package:		proteinstructure
 */
public class CiffilePdb extends Pdb {

	/*------------------------------ constants ------------------------------*/
	public static final String PDB_FTP_URL = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/mmCIF/";
	public static final String CIF_FILE_EXTENSION = ".cif.gz";
	
	/*--------------------------- member variables --------------------------*/
	
	// input file
	private File cifFile;

	// fields we will read
	private static final String entryId = "_entry";
	private static final String atomSiteId = "_atom_site";
	private static final String atomSitesAltId = "_atom_sites_alt";
	private static final String pdbxPolySeqId = "_pdbx_poly_seq_scheme";
	private static final String structConfId = "_struct_conf";
	private static final String structSheetId = "_struct_sheet_range"; 
	private static final String[] ids = {entryId,atomSitesAltId,atomSiteId,pdbxPolySeqId,structConfId,structSheetId};
	
	private TreeMap<String,Integer> ids2elements;					// map of ids to element serials
	private TreeMap<String,String> fields2values;					// map of field names (id.field) to values (for non-loop elements)
	private TreeMap<String,Integer> fields2indices;					// map of field names (id.field) to index (for loop elements)
	private TreeMap<String,Integer> ids2fieldsIdx;					// map of element ids to field index counter (after parseCifFile method done it contains the total number of fields per element id)
	private TreeSet<Integer> loopElements; 							// contains list of elements that are of loop type
	private TreeMap<Integer,Interval> loopelements2contentIndex;    // begin and end line index of each loop element
	
	private String altLoc;
 
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Constructs Pdb object from online PDB given pdb code and pdb chain code. 
	 * The DEFAULT_MODEL (see superclass) and default PDB_FTP_URL are used.
	 * @param pdbCode
	 * @param pdbChainCode
	 * @throws PdbChainCodeNotFoundError
	 * @throws IOException 
	 * @throws CiffileFormatError 
	 */
	public CiffilePdb(String pdbCode, String pdbChainCode) throws PdbChainCodeNotFoundError, IOException, CiffileFormatError {
		this(pdbCode, pdbChainCode, DEFAULT_MODEL, PDB_FTP_URL);
	}	

	/**
	 * Constructs Pdb object from online PDB given pdb code, pdb chain code and model serial.
	 * The default PDB_FTP_URL is used.
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param model_serial
	 * @throws PdbChainCodeNotFoundError
	 * @throws IOException 
	 * @throws CiffileFormatError 
	 */
	public CiffilePdb(String pdbCode, String pdbChainCode, int model_serial) throws PdbChainCodeNotFoundError, IOException, CiffileFormatError {
		this(pdbCode, pdbChainCode, model_serial, PDB_FTP_URL);
	}
	
	/**
	 * Constructs Pdb object from online PDB given pdb code, pdb chain code and pdbFtpUrl. 
	 * Model will be DEFAULT_MODEL (see superclass).
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param pdbFtpUrl
	 * @throws PdbChainCodeNotFoundError
	 * @throws IOException 
	 * @throws CiffileFormatError 
	 */
	public CiffilePdb (String pdbCode, String pdbChainCode, String pdbFtpUrl) throws PdbChainCodeNotFoundError, IOException, CiffileFormatError {
		this(pdbCode, pdbChainCode, DEFAULT_MODEL, pdbFtpUrl);
	}
	
	/**
	 * Constructs Pdb object from online PDB given pdb code, pdb chain code, model serial and pdbFtpUrl 
	 * @param pdbCode
	 * @param pdbChainCode
	 * @param model_serial
	 * @param pdbFtpUrl
	 * @throws PdbChainCodeNotFoundError
	 * @throws IOException 
	 * @throws CiffileFormatError 
	 */
	public CiffilePdb (String pdbCode, String pdbChainCode, int model_serial, String pdbFtpUrl) throws PdbChainCodeNotFoundError, IOException, CiffileFormatError {
		String tempDir = System.getProperty("java.io.tmpdir");	// TODO: read from FTP directly
		String gzCifFileName = pdbCode+CIF_FILE_EXTENSION;
		File gzCifFile = new File(tempDir,gzCifFileName);
		gzCifFile.deleteOnExit();
		this.cifFile = new File(tempDir,pdbCode + ".cif");
		this.cifFile.deleteOnExit();
		
		// getting gzipped cif file from ftp
		URL url = new URL(pdbFtpUrl+gzCifFileName);
		URLConnection urlc = url.openConnection();
		InputStream is = urlc.getInputStream();
		FileOutputStream os = new FileOutputStream(gzCifFile);
		int b;
		while ( (b=is.read())!=-1) {
			os.write(b);
		}
		is.close();
		os.close();
		
		// unzipping downloaded file
		GZIPInputStream zis = new GZIPInputStream(new FileInputStream(gzCifFile));
		os = new FileOutputStream(cifFile);
		while ( (b=zis.read())!=-1) {
			os.write(b);
		}
		zis.close();
		os.close();
		
		// here we would like to call the constructor this(ciffile, pdbChainCode, model_serial); which does not work, so we use copy/paste:
		
		// load from temp file
		this.pdbChainCode=pdbChainCode.toUpperCase();	// our convention: chain codes are upper case
		this.model=model_serial;
		
		parseCifFile();			
		this.pdbCode = readPdbCode();
		readAtomAltLocs(); // sets altLoc String (needed in readAtomSite to get the right alt atom locations)		
		readPdbxPolySeq(); // sets chainCode, sequence, pdbresser2resser		
		readAtomSite(); // populates resser_atom2atomserial, resser2restype, atomser2coord, atomser2resser 		
		secondaryStructure = new SecondaryStructure();	// create empty secondary structure first to make sure object is not null		
		readSecStructure(); // populates secondaryStructure	
		this.fullLength = sequence.length();		
		this.obsLength = resser2restype.size();		
		if(!secondaryStructure.isEmpty()) {
			secondaryStructure.setComment("CIFfile");
		}
		
		// we initialise resser2pdbresser from the pdbresser2resser HashMap
		this.resser2pdbresser = new HashMap<Integer, String>();
		for (String pdbresser:pdbresser2resser.keySet()){
			resser2pdbresser.put(pdbresser2resser.get(pdbresser), pdbresser);
		}
		
		// initialising atomser2atom from resser_atom2atomserial
		atomser2atom = new HashMap<Integer, String>();
		for (String resser_atom:resser_atom2atomserial.keySet()){
			int atomserial = resser_atom2atomserial.get(resser_atom);
			String atom = resser_atom.split("_")[1];
			atomser2atom.put(atomserial,atom);
		}
	}
	
	/**
	 * Constructs Pdb object given cif file and pdb chain code. 
	 * Model will be DEFAULT_MODEL
	 * @param ciffile
	 * @param pdbChainCode
	 * @throws PdbChainCodeNotFoundError
	 * @throws IOException 
	 * @throws CiffileFormatError 
	 */
	public CiffilePdb (File ciffile, String pdbChainCode) throws PdbChainCodeNotFoundError, IOException, CiffileFormatError {
		this(ciffile, pdbChainCode, DEFAULT_MODEL);
	}

	/**
	 * Constructs Pdb object given cif file, pdb chain code and model serial 
	 * @param ciffile
	 * @param pdbChainCode
	 * @param model_serial
	 * @throws PdbChainCodeNotFoundError
	 * @throws IOException 
	 * @throws CiffileFormatError 
	 */
	public CiffilePdb (File ciffile, String pdbChainCode, int model_serial) throws PdbChainCodeNotFoundError, IOException, CiffileFormatError {
		this.cifFile = ciffile;
		this.pdbChainCode=pdbChainCode.toUpperCase();	// our convention: chain codes are upper case
		this.model=model_serial;
		
		parseCifFile();			
		this.pdbCode = readPdbCode();
		readAtomAltLocs(); // sets altLoc String (needed in readAtomSite to get the right alt atom locations)		
		readPdbxPolySeq(); // sets chainCode, sequence, pdbresser2resser		
		readAtomSite(); // populates resser_atom2atomserial, resser2restype, atomser2coord, atomser2resser 		
		secondaryStructure = new SecondaryStructure();	// create empty secondary structure first to make sure object is not null		
		readSecStructure(); // populates secondaryStructure	
		this.fullLength = sequence.length();		
		this.obsLength = resser2restype.size();		
		if(!secondaryStructure.isEmpty()) {
			secondaryStructure.setComment("CIFfile");
		}
		
		// we initialise resser2pdbresser from the pdbresser2resser HashMap
		this.resser2pdbresser = new HashMap<Integer, String>();
		for (String pdbresser:pdbresser2resser.keySet()){
			resser2pdbresser.put(pdbresser2resser.get(pdbresser), pdbresser);
		}
		
		// initialising atomser2atom from resser_atom2atomserial
		atomser2atom = new HashMap<Integer, String>();
		for (String resser_atom:resser_atom2atomserial.keySet()){
			int atomserial = resser_atom2atomserial.get(resser_atom);
			String atom = resser_atom.split("_")[1];
			atomser2atom.put(atomserial,atom);
		}
	}
	
	/*---------------------------- private methods --------------------------*/
	
	private void parseCifFile() throws IOException, CiffileFormatError{
		// data structures to store the parsed fields
		ids2elements = new TreeMap<String, Integer>();
		fields2indices = new TreeMap<String,Integer>();
		fields2values = new TreeMap<String, String>();
		loopElements = new TreeSet<Integer>(); // contains list of elements that are of loop type
		loopelements2contentIndex = new TreeMap<Integer,Interval>();
		ids2fieldsIdx = new TreeMap<String,Integer>(); // this map holds the field index counters for each element id
		
		BufferedReader fcif = new BufferedReader(new FileReader(cifFile));
		int element = 0;
		String line;
		line = fcif.readLine(); // read first line
		Pattern p = Pattern.compile("^data_\\d\\w\\w\\w");
		if (!p.matcher(line).find()){
			throw new CiffileFormatError("The file doesn't seem to be a cif file");
		}
		int linecount = 1; // we have read one line already, we initialise count to 1
		while((line = fcif.readLine()) != null ) {
			linecount++;
			if (line.startsWith("#")) {
				element++;
				continue;
			}
			if (line.startsWith("loop_")) {
				loopElements.add(element);
				continue;
			}
			
			for (String id:ids){
				if (!ids2fieldsIdx.containsKey(id)) ids2fieldsIdx.put(id,0);
				p = Pattern.compile("^"+id+"\\.(\\w+)(?:\\s+(.*))?$");
				Matcher m = p.matcher(line);
				if (m.find()){
					ids2elements.put(id,element);
					String field = id + "." + m.group(1);
					if (!loopElements.contains(element)) { // if not a loop element
						fields2values.put(field, m.group(2)); // 2nd capture group only matches for non-loops where the value of the field is in same line as field name
					} else { // for loop elements we fill the fields2indices TreeMap
						fields2indices.put(field,ids2fieldsIdx.get(id));
					}
					ids2fieldsIdx.put(id,ids2fieldsIdx.get(id)+1); 
					continue;
				}
			}
			if (!line.startsWith("_") && !line.startsWith("#")){ // not in field definition, we are in values of a loop element 
				if (ids2elements.containsValue(element)) { // if this is one of the fields we want to parse (members of String[] ids)
					if (!loopelements2contentIndex.containsKey(element)) {
						//loopelements2content.put(element,line+"\n");
						Interval interval = new Interval(linecount,linecount);
						loopelements2contentIndex.put(element,interval);
					} else {
						//loopelements2content.put(element,loopelements2content.get(element)+line+"\n");
						loopelements2contentIndex.get(element).end=linecount;
					}
				}
			}			
		} // end scanning lines
		
		fcif.close();
	}
	
	private String readPdbCode(){
		return fields2values.get(entryId+".id").trim();
	}
	
	private void readAtomAltLocs() throws IOException, CiffileFormatError {
		// The read of the atom_sites_alt element must be done in a separate scan of the file, previous to scanning the atom_site element
		// This is because the order of the different elements in the cif files is not guaranteed, so atom_sites_alt can come before or after atom_site
		// (and altLoc needs to be set before starting reading the atom_site element)

		ArrayList<String> altLocs = new ArrayList<String>();
		// we initialise to ".", this is the default value in the cif files for the alt loc field. If no atom_sites_alt is present it's ok to stay with this value
		altLoc = ".";  
		
		// atom_sites_alt element is optional
		Interval intAtomSitesAlt = null;
		if (ids2elements.containsKey(atomSitesAltId)){
			intAtomSitesAlt = loopelements2contentIndex.get(ids2elements.get(atomSitesAltId));
		}

		BufferedReader fcif = new BufferedReader(new FileReader(cifFile));
		String line;
		int linecount=0;
		while((line = fcif.readLine()) != null ) {
			linecount++; 
			// atom_sites_alt (optional element)
			if (intAtomSitesAlt!=null && linecount>=intAtomSitesAlt.beg && linecount<=intAtomSitesAlt.end){
				int idIdx = fields2indices.get(atomSitesAltId+".id");
				// id=0
				// A ?
				String[] tokens = tokeniseFields(line);
				if (tokens.length!=ids2fieldsIdx.get(atomSitesAltId)) {
					throw new CiffileFormatError("Line "+linecount+" doesn't have the right number of fields for loop element "+atomSitesAltId);
				}
				if (!tokens[idIdx].equals(".")) {
					altLocs.add(tokens[idIdx]);
				}
			} 
		}
		fcif.close();
		if (!altLocs.isEmpty()){
			altLoc = Collections.min(altLocs);
		}
	}
	
	private void readAtomSite() throws IOException, PdbChainCodeNotFoundError, CiffileFormatError {
		resser_atom2atomserial = new HashMap<String,Integer>();
		resser2restype = new HashMap<Integer,String>();
		atomser2coord = new HashMap<Integer,Point3d>();
		atomser2resser = new HashMap<Integer,Integer>();
		
		Interval intAtomSite = loopelements2contentIndex.get(ids2elements.get(atomSiteId));
		
		boolean empty = true;
		BufferedReader fcif = new BufferedReader(new FileReader(cifFile));
		String line;
		int linecount=0;
		while((line = fcif.readLine()) != null ) {
			linecount++; 
			// atom_site
			if (linecount>=intAtomSite.beg && linecount<=intAtomSite.end){ 
				int groupPdbIdx = fields2indices.get(atomSiteId+".group_PDB");
				int idIdx = fields2indices.get(atomSiteId+".id");
				int labelAtomIdIdx = fields2indices.get(atomSiteId+".label_atom_id");
				int labelAltIdIdx = fields2indices.get(atomSiteId+".label_alt_id");
				int labelCompIdIdx = fields2indices.get(atomSiteId+".label_comp_id");
				int labelAsymIdIdx = fields2indices.get(atomSiteId+".label_asym_id");
				int labelSeqIdIdx = fields2indices.get(atomSiteId+".label_seq_id");
				int cartnXIdx = fields2indices.get(atomSiteId+".Cartn_x");
				int cartnYIdx = fields2indices.get(atomSiteId+".Cartn_y");
				int cartnZIdx = fields2indices.get(atomSiteId+".Cartn_z");
				int pdbxPDBModelNumIdx = fields2indices.get(atomSiteId+".pdbx_PDB_model_num");
				// group_PDB=0, auth_asym_id=22, pdbx_PDB_model_num=24, label_alt_id=4, id=1, label_atom_id=3, label_comp_id=5, label_asym_id=6, label_seq_id=8, Cartn_x=10, Cartn_y=11, Cartn_z=12
				//   0   1    2  3  4   5 6 7 8  9     10    11       12    13    14 151617181920   2122 23 24
				//ATOM   2    C CA  . MET A 1 1  ? 38.591 8.543   15.660  1.00 77.79  ? ? ? ? ? 1  MET A CA  1
				String[] tokens = tokeniseFields(line);
				if (tokens.length!=ids2fieldsIdx.get(atomSiteId)) {
					throw new CiffileFormatError("Line "+linecount+" doesn't have the right number of fields for loop element "+atomSiteId);
				}
				if (tokens[groupPdbIdx].equals("ATOM") && tokens[labelAsymIdIdx].equals(chainCode) && Integer.parseInt(tokens[pdbxPDBModelNumIdx])==model) { // match our given chain and model 
					empty = false;
					if (tokens[labelAltIdIdx].equals(".") || tokens[labelAltIdIdx].equals(altLoc)) { // don't read lines with something else as "." or altLoc
						int atomserial=Integer.parseInt(tokens[idIdx]); // id
						String atom = tokens[labelAtomIdIdx]; // label_atom_id
						String res_type = tokens[labelCompIdIdx]; // label_comp_id
						int res_serial = Integer.parseInt(tokens[labelSeqIdIdx]); // label_seq_id
						double x = Double.parseDouble(tokens[cartnXIdx]); // Cartn_x
						double y = Double.parseDouble(tokens[cartnYIdx]); // Cartn_y
						double z = Double.parseDouble(tokens[cartnZIdx]); // Cartn_z
						Point3d coords = new Point3d(x,y,z);
						if (AAinfo.isValidAA(res_type)) {
							atomser2coord.put(atomserial, coords);
							atomser2resser.put(atomserial, res_serial);
							resser2restype.put(res_serial, res_type);
							if (AAinfo.isValidAtomWithOXT(res_type,atom)){
								resser_atom2atomserial.put(res_serial+"_"+atom, atomserial);
							}
						}
					}
				}
				continue;
			}
		}
		fcif.close();
		if (empty) { // no atom data was found for given pdb chain code and model
			throw new PdbChainCodeNotFoundError("Couldn't find _atom_site data for given pdbChainCode: "+pdbChainCode+", model: "+model);
		}
	}
	
	private void readPdbxPolySeq() throws IOException, CiffileFormatError {
		pdbresser2resser = new HashMap<String, Integer>();
		sequence = "";
		
		String chainCodeStr=pdbChainCode;
		if (pdbChainCode.equals(Pdb.NULL_CHAIN_CODE)) chainCodeStr="A";
		
		Interval intPdbxPoly = loopelements2contentIndex.get(ids2elements.get(pdbxPolySeqId));
		
		BufferedReader fcif = new BufferedReader(new FileReader(cifFile));
		String line;
		int linecount=0;
		while((line = fcif.readLine()) != null ) {
			linecount++; 
			// pdbx_poly_seq_scheme
			if (linecount>=intPdbxPoly.beg && linecount<=intPdbxPoly.end){
				int asymIdIdx = fields2indices.get(pdbxPolySeqId+".asym_id");
				int seqIdIdx = fields2indices.get(pdbxPolySeqId+".seq_id");
				int authSeqNumIdx = fields2indices.get(pdbxPolySeqId+".auth_seq_num");
				int pdbInsCodeIdx = fields2indices.get(pdbxPolySeqId+".pdb_ins_code");
				int monIdIdx = fields2indices.get(pdbxPolySeqId+".mon_id");
				int pdbStrandIdIdx = fields2indices.get(pdbxPolySeqId+".pdb_strand_id");
				// asym_id=0, seq_id=2, auth_seq_num=6, pdb_ins_code=10, mon_id=3 
				// 0 1 2     3 4   5   6     7   8 910
				// A 1 1   ASP 1   1   1   ASP ASP A .
				String[] tokens = tokeniseFields(line);
				if (tokens.length!=ids2fieldsIdx.get(pdbxPolySeqId)) {
					throw new CiffileFormatError("Line "+linecount+" doesn't have the right number of fields for loop element "+pdbxPolySeqId);
				}
				if (tokens[pdbStrandIdIdx].equals(chainCodeStr)) { // we can't rely on using chainCode, because the order of elements is not guranteed (pdbx_poly_seq_scheme doesn't always come after atom_site)
					int res_serial = Integer.parseInt(tokens[seqIdIdx]); // seq_id
					chainCode = tokens[asymIdIdx];
					//TODO revise: do we want auth_seq_num or pdb_seq_num here??
					String pdb_res_serial = tokens[authSeqNumIdx]; // auth_seq_num
					String pdb_ins_code = tokens[pdbInsCodeIdx]; // pdb_ins_code
					String pdb_res_serial_with_icode = pdb_res_serial;
					if (!pdb_ins_code.equals(".")) {
						pdb_res_serial_with_icode=pdb_res_serial+pdb_ins_code;
					}
					String res_type = tokens[monIdIdx]; // mon_id
					// sequence
					if (AAinfo.isValidAA(res_type)){
		        		sequence+=AAinfo.threeletter2oneletter(res_type);
		        	} else {
		        		sequence+=NONSTANDARD_AA_LETTER;
		        	}
					// pdbresser2resser
					if (!pdb_res_serial_with_icode.startsWith("?")) { // question marks are author missing serials, we don't want them in the map
						pdbresser2resser.put(pdb_res_serial_with_icode,res_serial);
					}
				}
				continue;
			}

		}
		fcif.close();
	}
	
	private void readSecStructure() throws IOException, CiffileFormatError {
		secondaryStructure = new SecondaryStructure();
		
		// struct_conf element is optional
		Interval intStructConf = null;
		if (ids2elements.containsKey(structConfId)) {
			// if not a loop element then intStructConf stays null (because loopelements2contentIndex will return null)
			intStructConf = loopelements2contentIndex.get(ids2elements.get(structConfId));
		} 
		// taking care of cases where struct_conf is not a loop element but a one value field
		if (ids2elements.containsKey(structConfId) && !loopElements.contains(ids2elements.get(structConfId))){  
			String begChainCode = fields2values.get(structConfId+".beg_label_asym_id").trim();
			if (begChainCode.equals(chainCode)) { // chainCode has been set already in reading pdbx_poly_seq_scheme			
				String id = fields2values.get(structConfId+".id").trim();
				int beg = Integer.parseInt(fields2values.get(structConfId+".beg_label_seq_id").trim());
				int end = Integer.parseInt(fields2values.get(structConfId+".end_label_seq_id").trim());
				Pattern p = Pattern.compile("^(\\w).+_P(\\d)+$");
				Matcher m = p.matcher(id);
				String ssId="Unknown";
				if (m.find()){
					ssId = m.group(1)+m.group(2); // e.g.: Hnn (helices) or Tnn (turns) 				
				}
				char ssType = SecStrucElement.OTHER;
				if(id.startsWith("H")) {
					ssType = SecStrucElement.HELIX;
				} else if(id.startsWith("T")) {
					ssType = SecStrucElement.TURN;
				} else {
					System.err.println("Unknown secondary structure type " + id + " encountered when reading from ciffile. Skipping.");
				}
				if(ssType != SecStrucElement.OTHER) {
					SecStrucElement ssElem = new SecStrucElement(ssType, beg, end, ssId);
					secondaryStructure.add(ssElem);
				}
			}
		}
		// struct_sheet_range element is optional
		Interval intStructSheet = null; 
		if (ids2elements.containsKey(structSheetId)) {
			// if not a loop element intStructSheet stays null (because loopelements2contentIndex will return null)
			intStructSheet = loopelements2contentIndex.get(ids2elements.get(structSheetId));
		}
		// taking care of cases where struct_sheet_range is not a loop element but a one value field
		if (ids2elements.containsKey(structSheetId) && !loopElements.contains(ids2elements.get(structSheetId))){
			String begChainCode = fields2values.get(structSheetId+".beg_label_asym_id").trim();
			if (begChainCode.equals(chainCode)){ // chainCode has been set already in reading pdbx_poly_seq_scheme
				String sheetid = fields2values.get(structSheetId+".sheet_id").trim(); //tokens[sheetIdIdx];
				int id = Integer.parseInt(fields2values.get(structSheetId+".id").trim()); //Integer.parseInt(tokens[idIdx]);
				int beg = Integer.parseInt(fields2values.get(structSheetId+".beg_label_seq_id").trim()); //tokens[begLabelSeqIdIdx]);
				int end = Integer.parseInt(fields2values.get(structSheetId+".end_label_seq_id").trim()); //tokens[endLabelSeqIdIdx]);
				String ssId=SecStrucElement.STRAND+sheetid+id; // e.g.: SA1, SA2..., SB1, SB2,...
				SecStrucElement ssElem = new SecStrucElement(SecStrucElement.STRAND, beg, end, ssId);
				secondaryStructure.add(ssElem);
			}

		}
		
		BufferedReader fcif = new BufferedReader(new FileReader(cifFile));
		String line;
		int linecount=0;
		while((line = fcif.readLine()) != null ) {
			linecount++;
			// struct_conf (optional element), HELIX and TURN secondary structure
			if (intStructConf!=null && linecount>=intStructConf.beg && linecount<=intStructConf.end){
				int idIdx = fields2indices.get(structConfId+".id");
				int begLabelAsymIdIdx = fields2indices.get(structConfId+".beg_label_asym_id");
				int begLabelSeqIdIdx = fields2indices.get(structConfId+".beg_label_seq_id");
				int endLabelSeqIdIdx = fields2indices.get(structConfId+".end_label_seq_id");
				//id=1, beg_label_seq_id=5, end_label_seq_id=9, beg_label_asym_id=4
				//     0       1  2    3 4 5   6   7 8  9 10  111213    1415 16 1718 19
				//HELX_P HELX_P1  1  ASN A 2   ? GLY A 12  ? ASN A 2   GLY A 12  1 ? 11
				String[] tokens = tokeniseFields(line);
				if (tokens.length!=ids2fieldsIdx.get(structConfId)) {
					throw new CiffileFormatError("Line "+linecount+" doesn't have the right number of fields for loop element "+structConfId);
				}
				if (tokens[begLabelAsymIdIdx].equals(chainCode)) { // chainCode has been set already in reading pdbx_poly_seq_scheme
					String id = tokens[idIdx];
					Pattern p = Pattern.compile("^(\\w).+_P(\\d)+$");
					Matcher m = p.matcher(id);
					String ssId="Unknown";
					if (m.find()){
						ssId = m.group(1)+m.group(2); // e.g.: Hnn (helices) or Tnn (turns) 				
					}
					int beg = Integer.parseInt(tokens[begLabelSeqIdIdx]);
					int end = Integer.parseInt(tokens[endLabelSeqIdIdx]);
					char ssType = SecStrucElement.OTHER;
					if(id.startsWith("H")) {
						ssType = SecStrucElement.HELIX;
					} else if(id.startsWith("T")) {
						ssType = SecStrucElement.TURN;
					} else {
						System.err.println("Unknown secondary structure type " + id + " encountered when reading from ciffile. Skipping.");
					}
					if(ssType != SecStrucElement.OTHER) {
						SecStrucElement ssElem = new SecStrucElement(ssType, beg, end, ssId);
						secondaryStructure.add(ssElem);
					}
				}
				continue;
			}
			// struct_sheet_range (optional element), SHEETs
			if (intStructSheet!=null && linecount>=intStructSheet.beg && linecount<=intStructSheet.end){
				int sheetIdIdx = fields2indices.get(structSheetId+".sheet_id");
				int idIdx = fields2indices.get(structSheetId+".id");
				int begLabelAsymIdIdx = fields2indices.get(structSheetId+".beg_label_asym_id");
				int begLabelSeqIdIdx = fields2indices.get(structSheetId+".beg_label_seq_id");
				int endLabelSeqIdIdx = fields2indices.get(structSheetId+".end_label_seq_id");
				//sheet_id=0, id=1, beg_label_seq_id=4, end_label_seq_id=8, beg_label_asym_id=3
				//0 1   2 3  4 5   6 7  8 910  1112 13  1415 16
				//A 1 ARG A 14 ? LYS A 19 ? ? ARG A 14 LYS A 19
				String[] tokens = tokeniseFields(line);
				if (tokens.length!=ids2fieldsIdx.get(structSheetId)) {
					throw new CiffileFormatError("Line "+linecount+" doesn't have the right number of fields for loop element "+structSheetId);
				}
				if (tokens[begLabelAsymIdIdx].equals(chainCode)){ // chainCode has been set already in reading pdbx_poly_seq_scheme
					String sheetid = tokens[sheetIdIdx];
					int id = Integer.parseInt(tokens[idIdx]);
					int beg = Integer.parseInt(tokens[begLabelSeqIdIdx]);
					int end = Integer.parseInt(tokens[endLabelSeqIdIdx]);
					String ssId=SecStrucElement.STRAND+sheetid+id; // e.g.: SA1, SA2..., SB1, SB2,...
					SecStrucElement ssElem = new SecStrucElement(SecStrucElement.STRAND, beg, end, ssId);
					secondaryStructure.add(ssElem);
				}
				continue;
			}

		}
		fcif.close();
	}
	
	/**
	 * Splits a space separated line into its individual tokens returning an array with all tokens
	 * Takes care of quoted fields that contain spaces
	 * e.g. HELX_P HELX_P2 H4 GLY A 111 ? GLU A 127 ? GLY A 112 GLU A 128 1 'SEE REMARK 650' 17
	 * 
	 * TODO maybe StreamTokenizer can do all this with less pain
	 * @param line
	 * @return
	 */
	private String[] tokeniseFields(String line) {
		String[] tokens;
		if (line.contains("'")) { // if there are single quotes in the line
			ArrayList<String> tokensAL = new ArrayList<String>();
			Pattern p = Pattern.compile("'[^']*'|[^ \\t]+"); // note: regex doesn't work inverting the order of expressions in the 'OR' (in python it does!)
			Matcher m = p.matcher(line);
			while (m.find()){
				tokensAL.add(m.group());
			}
			tokens = new String[tokensAL.size()];
			tokensAL.toArray(tokens);
		} else { // if no quotes we simply split by columns using spaces as delimiters
			tokens = line.split("\\s+");
		}
		if (line.contains("\"")){ // in some rare cases some fields are quoted with double quotes, this seems to be to escape single quotes within them (used normally as a "prime" symbol) 
			for (int i=0;i<tokens.length;i++){ // we get rid of the double quoting
				tokens[i] = tokens[i].replaceAll("\"", "");
			}
		}
		return tokens;
	}
	
}
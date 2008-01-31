package proteinstructure;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.RandomAccessFile;

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
	private TreeMap<Integer,Long[]> loopelements2contentOffset;    // begin and end line index of each loop element
	
	private String altLoc;
	
	private RandomAccessFile fcif;
 
	private boolean fieldsTitlesRead;
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Constructs an empty Pdb object from online PDB given pdb code
	 * Data will be downloaded an stored in local file
	 * but will only be loaded from local file upon call of load(pdbChainCode, modelSerial)  
	 * The default PDB_FTP_URL is used.
	 * @param pdbCode
	 * @throws IOException  
	 */
	public CiffilePdb(String pdbCode) throws IOException {
		this(pdbCode, PDB_FTP_URL);
	}
	
	/**
	 * Constructs an empty Pdb object from online PDB given pdb code and pdbFtpUrl
	 * Data will be downloaded an stored in local file
	 * but will only be loaded from local file upon call of load(pdbChainCode, modelSerial)  
	 * @param pdbCode
	 * @param pdbFtpUrl
	 * @throws IOException  
	 */
	public CiffilePdb (String pdbCode, String pdbFtpUrl) throws IOException {
		this.dataLoaded = false;
		this.fieldsTitlesRead = false;
		
		// we store the file locally instead of reading directly from the ftp stream, so that the file can be cached locally in applications like CMView
		String tempDir = System.getProperty("java.io.tmpdir");	
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
		
	}
	
	/**
	 * Constructs an empty Pdb object given cif file
	 * Data will be loaded from file upon call of load(pdbChainCode, modelSerial) 
	 * @param ciffile
	 */
	public CiffilePdb (File ciffile) {
		this.cifFile = ciffile;
		this.dataLoaded = false;
		this.fieldsTitlesRead = false;
		
	}
	
	public File getCifFile() {
		return cifFile;
	}
	
	public void load(String pdbChainCode, int modelSerial) throws PdbLoadError{
		try {
			this.model = modelSerial;
			this.pdbChainCode=pdbChainCode;				// NOTE! pdb chain codes are case sensitive
			fcif = new RandomAccessFile(cifFile,"r");
			parseCifFile();
			fcif.close();

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
			
			dataLoaded = true;
			
		} catch (CiffileFormatError e) {
			throw new PdbLoadError(e);
		} catch (IOException e) {
			throw new PdbLoadError(e);
		} catch (PdbChainCodeNotFoundError e) {
			throw new PdbLoadError(e);
		}

	}
	
	public String[] getChains() throws PdbLoadError {
		TreeSet<String> chains = new TreeSet<String>();
	
		try {
			fcif = new RandomAccessFile(cifFile,"r");
			if (!fieldsTitlesRead) {
				readFieldsTitles();
			}
			Long[] intPdbxPoly = loopelements2contentOffset.get(ids2elements.get(pdbxPolySeqId));
			
			int recordCount=0;
			
			fcif.seek(intPdbxPoly[0]);
			while(fcif.getFilePointer()<intPdbxPoly[1]) {
				recordCount++; 

				int pdbStrandIdIdx = fields2indices.get(pdbxPolySeqId+".pdb_strand_id");
				int numberFields = ids2fieldsIdx.get(pdbxPolySeqId);
				String[] tokens = tokeniseFields(numberFields);
				if (tokens.length!=numberFields) {
					throw new CiffileFormatError("Incorrect number of fields for record "+recordCount+" in loop element "+pdbxPolySeqId);
				}
				chains.add(tokens[pdbStrandIdIdx]);
			}
			fcif.close();
			
		} catch (IOException e) {
			throw new PdbLoadError(e);
		} catch (CiffileFormatError e) {
			throw new PdbLoadError(e);
		}
		
		if (chains.isEmpty()) return null;
		
		String[] chainsArray = new String[chains.size()];
		chains.toArray(chainsArray);
		return chainsArray;
	}
	
	public Integer[] getModels() throws PdbLoadError {
		TreeSet<Integer> models = new TreeSet<Integer>();
		try {
			fcif = new RandomAccessFile(cifFile,"r");
			if (!fieldsTitlesRead) {
				readFieldsTitles();
			}
			Long[] intAtomSite = loopelements2contentOffset.get(ids2elements.get(atomSiteId));
			
			int recordCount=0;
			
			fcif.seek(intAtomSite[0]);
			while(fcif.getFilePointer()<intAtomSite[1]) {
				recordCount++; 

				int pdbxPDBModelNumIdx = fields2indices.get(atomSiteId+".pdbx_PDB_model_num");
				int numberFields = ids2fieldsIdx.get(atomSiteId);
				String[] tokens = tokeniseFields(numberFields);
				if (tokens.length!=numberFields) {
					throw new CiffileFormatError("Incorrect number of fields for record "+recordCount+" in loop element "+atomSiteId);
				}
				models.add(Integer.parseInt(tokens[pdbxPDBModelNumIdx]));
			}
			fcif.close();
			
		} catch (IOException e) {
			throw new PdbLoadError(e);
		} catch (CiffileFormatError e) {
			throw new PdbLoadError(e);
		}
		
		if (models.isEmpty()) return null;
		
		Integer[] modelsArray = new Integer[models.size()];
		models.toArray(modelsArray);
		return modelsArray;
	}
	
	/*---------------------------- private methods --------------------------*/
	
	private void parseCifFile() throws IOException, CiffileFormatError, PdbChainCodeNotFoundError{
		
		if (!fieldsTitlesRead) {
			readFieldsTitles();
		}		
		// now reading separate elements separately using private methods
		// the order in the elements in the file is not guaranteed, that's why (among other reasons) we have to use RandomAccessFile
		this.pdbCode = readPdbCode();
		readAtomAltLocs(); // sets altLoc String (needed in readAtomSite to get the right alt atom locations)		
		readPdbxPolySeq(); // sets chainCode, sequence, pdbresser2resser		
		readAtomSite(); // populates resser_atom2atomserial, resser2restype, atomser2coord, atomser2resser 		
		secondaryStructure = new SecondaryStructure();	// create empty secondary structure first to make sure object is not null		
		readSecStructure(); // populates secondaryStructure	

	}
	
	private void readFieldsTitles() throws IOException, CiffileFormatError {
		// data structures to store the parsed fields
		ids2elements = new TreeMap<String, Integer>();
		fields2indices = new TreeMap<String,Integer>();
		fields2values = new TreeMap<String, String>();
		loopElements = new TreeSet<Integer>(); // contains list of elements that are of loop type
		loopelements2contentOffset = new TreeMap<Integer,Long[]>();
		ids2fieldsIdx = new TreeMap<String,Integer>(); // this map holds the field index counters for each element id
		
		int element = 0;
		String line;
		line = fcif.readLine(); // read first line
		Pattern p = Pattern.compile("^data_\\d\\w\\w\\w");
		if (!p.matcher(line).find()){
			throw new CiffileFormatError("The file doesn't seem to be a cif file");
		}
		int linecount = 1; // we have read one line already, we initialise count to 1
		// we need to store the last line's byte offset (which indicates the beginning of this line) 
		long lastLineOffset=fcif.getFilePointer();
		while((line = fcif.readLine()) != null ) {
			long currentOffset = fcif.getFilePointer(); //this gets byte offset at end of line
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
					if (!loopelements2contentOffset.containsKey(element)) {
						//loopelements2content.put(element,line+"\n");
						Long[] interval = {lastLineOffset, currentOffset};
						loopelements2contentOffset.put(element,interval);
					} else {
						//loopelements2content.put(element,loopelements2content.get(element)+line+"\n");
						loopelements2contentOffset.get(element)[1]=currentOffset;
					}
				}
			}
			lastLineOffset = currentOffset; //we store this line's offset to have it for next iteration
		} // end scanning lines
		
		fieldsTitlesRead = true;
	}
	
	private String readPdbCode(){
		return fields2values.get(entryId+".id").trim().toLowerCase();
	}
	
	private void readAtomAltLocs() throws IOException, CiffileFormatError {
		// The read of the atom_sites_alt element must be done previously to scanning the atom_site element
		// This is because the order of the different elements in the cif files is not guaranteed, so atom_sites_alt can come before or after atom_site
		// (and altLoc needs to be set before starting reading the atom_site element)

		ArrayList<String> altLocs = new ArrayList<String>();
		// we initialise to ".", this is the default value in the cif files for the alt loc field. If no atom_sites_alt is present it's ok to stay with this value
		altLoc = ".";  
		
		// atom_sites_alt element is optional
		Long[] intAtomSitesAlt = null;
		if (ids2elements.containsKey(atomSitesAltId)){
			intAtomSitesAlt = loopelements2contentOffset.get(ids2elements.get(atomSitesAltId));
		}

		int recordCount = 0;
		// atom_sites_alt (optional element)
		if (intAtomSitesAlt!=null) {
			
			fcif.seek(intAtomSitesAlt[0]);			
			while(fcif.getFilePointer()<intAtomSitesAlt[1]) {
				recordCount++;
				 
				int idIdx = fields2indices.get(atomSitesAltId+".id");
				// id=0
				// A ?
				int numberFields = ids2fieldsIdx.get(atomSitesAltId);
				String[] tokens = tokeniseFields(numberFields);
				if (tokens.length!=numberFields) {
					throw new CiffileFormatError("Incorrect number of fields for record "+recordCount+" in loop element "+atomSitesAltId);
				}
				if (!tokens[idIdx].equals(".")) {
					altLocs.add(tokens[idIdx]);
				}
			}
		}
		if (!altLocs.isEmpty()){
			altLoc = Collections.min(altLocs);
		}
	}
		
	private void readAtomSite() throws IOException, PdbChainCodeNotFoundError, CiffileFormatError {
		resser_atom2atomserial = new HashMap<String,Integer>();
		resser2restype = new HashMap<Integer,String>();
		atomser2coord = new HashMap<Integer,Point3d>();
		atomser2resser = new HashMap<Integer,Integer>();
		
		Long[] intAtomSite = loopelements2contentOffset.get(ids2elements.get(atomSiteId));
		
		boolean empty = true;
		int recordCount = 0;
		
		fcif.seek(intAtomSite[0]);
		while(fcif.getFilePointer()<intAtomSite[1]) {
			recordCount++;

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
			int numberFields = ids2fieldsIdx.get(atomSiteId);
			String[] tokens = tokeniseFields(numberFields);
			if (tokens.length!=numberFields) {
				throw new CiffileFormatError("Incorrect number of fields for record "+recordCount+" in loop element "+atomSiteId);
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
		}
		if (empty) { // no atom data was found for given pdb chain code and model
			throw new PdbChainCodeNotFoundError("Couldn't find _atom_site data for given pdbChainCode: "+pdbChainCode+", model: "+model);
		}
	}
	
	private void readPdbxPolySeq() throws IOException, CiffileFormatError {
		pdbresser2resser = new HashMap<String, Integer>();
		sequence = "";
		
		String chainCodeStr=pdbChainCode;
		if (pdbChainCode.equals(Pdb.NULL_CHAIN_CODE)) chainCodeStr="A";
		
		Long[] intPdbxPoly = loopelements2contentOffset.get(ids2elements.get(pdbxPolySeqId));
		
		int recordCount=0;
		
		fcif.seek(intPdbxPoly[0]);
		while(fcif.getFilePointer()<intPdbxPoly[1]) {
			recordCount++; 

			int asymIdIdx = fields2indices.get(pdbxPolySeqId+".asym_id");
			int seqIdIdx = fields2indices.get(pdbxPolySeqId+".seq_id");
			//int authSeqNumIdx = fields2indices.get(pdbxPolySeqId+".auth_seq_num");
			int pdbSeqNumIdx = fields2indices.get(pdbxPolySeqId+".pdb_seq_num");
			int pdbInsCodeIdx = fields2indices.get(pdbxPolySeqId+".pdb_ins_code");
			int monIdIdx = fields2indices.get(pdbxPolySeqId+".mon_id");
			int pdbStrandIdIdx = fields2indices.get(pdbxPolySeqId+".pdb_strand_id");
			// asym_id=0, seq_id=2, auth_seq_num=6, pdb_ins_code=10, mon_id=3 
			// 0 1 2     3 4   5   6     7   8 910
			// A 1 1   ASP 1   1   1   ASP ASP A .
			int numberFields = ids2fieldsIdx.get(pdbxPolySeqId);
			String[] tokens = tokeniseFields(numberFields);
			if (tokens.length!=numberFields) {
				throw new CiffileFormatError("Incorrect number of fields for record "+recordCount+" in loop element "+pdbxPolySeqId);
			}
			if (tokens[pdbStrandIdIdx].equals(chainCodeStr)) { // we can't rely on using chainCode, because the order of elements is not guranteed (pdbx_poly_seq_scheme doesn't always come after atom_site)
				int res_serial = Integer.parseInt(tokens[seqIdIdx]); // seq_id
				chainCode = tokens[asymIdIdx];
				String pdb_res_serial = tokens[pdbSeqNumIdx]; // pdb_seq_num
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
					sequence+=AAinfo.NONSTANDARD_AA_ONE_LETTER;
				}
				// pdbresser2resser
				pdbresser2resser.put(pdb_res_serial_with_icode,res_serial);

			}
		}
	}
	
	private void readSecStructure() throws IOException, CiffileFormatError {
		secondaryStructure = new SecondaryStructure();
		
		// struct_conf element is optional
		Long[] intStructConf = null;
		if (ids2elements.containsKey(structConfId)) {
			// if not a loop element then intStructConf stays null (because loopelements2contentIndex will return null)
			intStructConf = loopelements2contentOffset.get(ids2elements.get(structConfId));
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
		Long[] intStructSheet = null; 
		if (ids2elements.containsKey(structSheetId)) {
			// if not a loop element intStructSheet stays null (because loopelements2contentIndex will return null)
			intStructSheet = loopelements2contentOffset.get(ids2elements.get(structSheetId));
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
				
		if (intStructConf!=null) {
			int recordCount=0;
			
			fcif.seek(intStructConf[0]);
			while(fcif.getFilePointer()<intStructConf[1]) {
				recordCount++;
				// struct_conf (optional element), HELIX and TURN secondary structure

				int idIdx = fields2indices.get(structConfId+".id");
				int begLabelAsymIdIdx = fields2indices.get(structConfId+".beg_label_asym_id");
				int begLabelSeqIdIdx = fields2indices.get(structConfId+".beg_label_seq_id");
				int endLabelSeqIdIdx = fields2indices.get(structConfId+".end_label_seq_id");
				//id=1, beg_label_seq_id=5, end_label_seq_id=9, beg_label_asym_id=4
				//     0       1  2    3 4 5   6   7 8  9 10  111213    1415 16 1718 19
				//HELX_P HELX_P1  1  ASN A 2   ? GLY A 12  ? ASN A 2   GLY A 12  1 ? 11
				int numFields = ids2fieldsIdx.get(structConfId);
				String[] tokens = tokeniseFields(numFields);
				if (tokens.length!=numFields) {
					throw new CiffileFormatError("Incorrect number of fields for record "+recordCount+" in loop element "+structConfId);
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
			}
		}
		if (intStructSheet!=null) {
			int recordCount=0;
			
			fcif.seek(intStructSheet[0]);
			while(fcif.getFilePointer()<intStructSheet[1]) {
				recordCount++;
				// struct_sheet_range (optional element), SHEETs
				int sheetIdIdx = fields2indices.get(structSheetId+".sheet_id");
				int idIdx = fields2indices.get(structSheetId+".id");
				int begLabelAsymIdIdx = fields2indices.get(structSheetId+".beg_label_asym_id");
				int begLabelSeqIdIdx = fields2indices.get(structSheetId+".beg_label_seq_id");
				int endLabelSeqIdIdx = fields2indices.get(structSheetId+".end_label_seq_id");
				//sheet_id=0, id=1, beg_label_seq_id=4, end_label_seq_id=8, beg_label_asym_id=3
				//0 1   2 3  4 5   6 7  8 910  1112 13  1415 16
				//A 1 ARG A 14 ? LYS A 19 ? ? ARG A 14 LYS A 19
				int numFields = ids2fieldsIdx.get(structSheetId);
				String[] tokens = tokeniseFields(numFields);
				if (tokens.length!=numFields) {
					throw new CiffileFormatError("Incorrect number of fields for record "+recordCount+" in loop element "+structSheetId);
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
			}
		}
	}

	/**
	 * Splits a space separated line into its individual tokens returning an array with all tokens
	 * Takes care of all particularities of the format of a record in the ciffiles:
	 *  - fields within records are separated by spaces
	 *  - spaces can be used within quoted strings (at the moment this only supports single quotes, not double)
	 *  - free style with all characters allowed if something is quoted with \n; ;\n 
	 * The java class StreamTokenizer could have done all this, but it was limited to do all that we needed to do
	 *  
	 *  
	 * This method is black magic. I don't fully understand it myself as I write it.
	 * If you need to come back to this and read it, good luck!!
	 * 
	 * @param numberTokens
	 * @return
	 */
	private String[] tokeniseFields(int numberTokens) throws IOException {
		String[] tokens = new String[numberTokens];
		// initialise tokens to empty strings
		for (int i=0; i<numberTokens;i++){
			tokens[i]="";
		}
		
		int i = 0;
		char lastChar=' ';
		char quoteChar = 0;
		while (true) {
			char currentChar = (char)fcif.readByte();
			
			// '' quoting
			if (quoteChar!=';' && currentChar=='\'' && (lastChar==' ' || lastChar=='\n')){
				quoteChar = '\'';
			}
			else if (quoteChar!=';' && currentChar==' ' && lastChar=='\''){
				quoteChar = 0;
			}
			// "" quoting
			if (quoteChar!=';' && currentChar=='"' && (lastChar==' ' || lastChar=='\n')){
				quoteChar = '"';
			}
			else if (quoteChar!=';' && currentChar==' ' && lastChar=='"'){
				quoteChar = 0;
			}			
			// ;; quoting (multi-line quoting)
			if (quoteChar!=';' && currentChar==';' && lastChar=='\n'){
				quoteChar = ';';
			}
			else if (quoteChar==';' && currentChar==';' && lastChar=='\n'){
				quoteChar = 0;
			}
			
			// reading field
			if (quoteChar==0) { // not within quotes
				if (currentChar==' ' || currentChar=='\n') { 
					if (currentChar!=lastChar && !(currentChar=='\n' && lastChar==' ')) i++; // we only increment when we move from a non-space to a space or from non-space to \n
				} else {
					tokens[i]+=currentChar;
					// if we are adding the last ; of a ;;-quoted string then strip the starting ';' and ending "\n;" out 
					if (currentChar==';' && lastChar=='\n' && tokens[i].startsWith(";") && tokens[i].endsWith("\n;")) {
						tokens[i]=tokens[i].replaceFirst("^;", "");
						tokens[i]=tokens[i].replaceFirst("\n;","");
					}									
				} 
			} else {			// within quotes (of type '', ""  or  ;;)
				tokens[i]+=currentChar;
				// if string is surrounded by '' or "" then strip them out (except when string is length 1 and thus beginning and end are quoteChar)
				if (tokens[i].length()!=1 && tokens[i].startsWith(Character.toString(quoteChar)) && tokens[i].endsWith(Character.toString(quoteChar))) tokens[i]=tokens[i].replaceAll(Character.toString(quoteChar), "");

			}
			
			lastChar = currentChar;
			 
			if (i==numberTokens) {
				// for the last record of an element it is important to have read up to the end of the line (including the '\n'), 
				// otherwise the condition : "while (current_pointer<max_pointer_of_this_element)" won't work
				// we read one more character at a time: test whether it is a ' ' or a '\n', if not then we have overread so we need to rewind back
				while (true) {
					long currentPos = fcif.getFilePointer(); // get current position to rewind back to it if needed
					currentChar = (char) fcif.readByte();
					if (currentChar!='\n' && currentChar!=' '){ 
						fcif.seek(currentPos);
						break;
					}
				}
				return tokens;
			}
		}
	}

}

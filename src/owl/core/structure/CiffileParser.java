package owl.core.structure;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.RandomAccessFile;

import java.io.IOException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;
import owl.core.util.FileFormatException;


/**
 * A mmCIF file format parser. 
 * 
 * @author		Jose Duarte
 */
public class CiffileParser {

	/*------------------------------ constants ------------------------------*/
	public static final String PDB_FTP_URL = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/mmCIF/";
	public static final String CIF_FILE_EXTENSION = ".cif.gz";
	
	private static final String TMP_DIR = System.getProperty("java.io.tmpdir");
	
	/*--------------------------- member variables --------------------------*/
	
	// input file
	private File cifFile;

	// fields we will read
	private static final String entryId = "_entry";
	private static final String structId = "_struct";
	private static final String atomSiteId = "_atom_site";
	private static final String pdbxPolySeqId = "_pdbx_poly_seq_scheme";
	private static final String structConfId = "_struct_conf";
	private static final String structSheetId = "_struct_sheet_range";
	private static final String cell = "_cell";
	private static final String symmetry = "_symmetry";
	private static final String exptl = "_exptl";
	private static final String reflns = "_reflns";
	private static final String refine = "_refine";
	private static final String atomSitesId = "_atom_sites";
	private static final String[] ids = 
		{entryId,structId,atomSiteId,pdbxPolySeqId,structConfId,structSheetId,cell,symmetry,exptl,reflns,refine,atomSitesId};
	
	private TreeMap<String,Integer> ids2elements;					// map of ids to element serials
	private TreeMap<String,String> fields2values;					// map of field names (id.field) to values (for non-loop elements)
	private TreeMap<String,Integer> fields2indices;					// map of field names (id.field) to index (for loop elements)
	private TreeMap<String,Integer> ids2fieldsIdx;					// map of element ids to field index counter (after parseCifFile method done it contains the total number of fields per element id)
	private TreeSet<Integer> loopElements; 							// contains list of elements that are of loop type
	private TreeMap<Integer,Long[]> loopelements2contentOffset;    // begin and end line index of each loop element, will also contain begin and end of multi-line quoted (';' quoting) values in non-loop elements
	
	private transient RandomAccessFile fcif;
 
	private boolean fieldsTitlesRead;
	
	private String[] chainsArray;
	private Integer[] modelsArray;


	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Constructs a cif file parser from online PDB given pdb code
	 * Data will be downloaded an stored in local file
	 * but will only be loaded from local file upon call of readChain(pdbChainCode, modelSerial)  
	 * The default PDB_FTP_URL is used.
	 * @param pdbCode
	 * @throws IOException  
	 * @throws FileFormatException 
	 */
	public CiffileParser(String pdbCode) throws IOException, FileFormatException {
		this(pdbCode, PDB_FTP_URL);
	}
	
	/**
	 * Constructs a cif file parser from online PDB given pdb code and pdbFtpUrl
	 * Data will be downloaded an stored in local file
	 * but will only be loaded from local file upon call of readChain(pdbChainCode, modelSerial)  
	 * @param pdbCode
	 * @param pdbFtpUrl
	 * @throws IOException  
	 * @throws FileFormatException 
	 */
	public CiffileParser (String pdbCode, String pdbFtpUrl) throws IOException, FileFormatException {
		this.fieldsTitlesRead = false;
		
		// we store the file locally instead of reading directly from the ftp stream, so that the file can be cached locally in applications like CMView	
		String gzCifFileName = pdbCode+CIF_FILE_EXTENSION;
		File gzCifFile = new File(TMP_DIR,gzCifFileName);
		gzCifFile.deleteOnExit();
		this.cifFile = new File(TMP_DIR,pdbCode + ".cif");
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
		
		fcif = new RandomAccessFile(cifFile,"r");
		readFieldsTitles();
	}
	
	/**
	 * Constructs a cif file parser object given cif file
	 * Data will be loaded from file upon call of {@link #readChain(String, int)} 
	 * @param ciffile
	 * @throws FileFormatException 
	 * @throws IOException 
	 */
	public CiffileParser (File ciffile) throws IOException, FileFormatException {
		this.cifFile = ciffile;
		this.fieldsTitlesRead = false;
		fcif = new RandomAccessFile(cifFile,"r");
		readFieldsTitles();
	}
	
	public File getCifFile() {
		return cifFile;
	}
	
	/**
	 * Reads PDB data (coordinates, sequence, etc.) from the cif file
	 * for given pdbChainCode and modelSerial
	 * @param pdbAsymUnit
	 * @param modelSerial
	 * @throws PdbLoadException
	 */
	public void readChains(PdbAsymUnit pdbAsymUnit, int modelSerial) throws PdbLoadException{

		try {

			if (!fieldsTitlesRead) {
				readFieldsTitles();
			}

			readPdbxPolySeq(pdbAsymUnit);

			// this needs the info in the pdbress2resser maps
			this.readAtomSite(pdbAsymUnit,modelSerial);

			for (PdbChain pdb:pdbAsymUnit.getAllChains()) {
				pdb.initialiseMaps();				
			}

			readSecStructure(pdbAsymUnit);
		
		} catch (FileFormatException e) {
			throw new PdbLoadException(e);
		} catch (IOException e) {
			throw new PdbLoadException(e);
		} 

	}
	
	public void closeFile() throws IOException {
		fcif.close();
	}
	
	/**
	 * Returns all alphabetically sorted PDB chain codes found in the cif file.
	 * It caches the result so that next time called no parsing has to be done. 
	 * @return array with all pdb chain codes
	 */
	public String[] getChains() throws PdbLoadException {
		if (chainsArray!=null) {
			return chainsArray;
		}
		TreeSet<String> chains = new TreeSet<String>();
	
		try {
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
					throw new FileFormatException("Incorrect number of fields for record "+recordCount+" in loop element "+pdbxPolySeqId+" of CIF file "+cifFile);
				}
				chains.add(tokens[pdbStrandIdIdx]);
			}
			
		} catch (IOException e) {
			throw new PdbLoadException(e);
		} catch (FileFormatException e) {
			throw new PdbLoadException(e);
		}
		
		if (chains.isEmpty()) return null;
		
		chainsArray = new String[chains.size()];
		chains.toArray(chainsArray);
		return chainsArray;
	}
	
	/**
	 * Returns all model serials found in the cif file.
	 * It caches the result so that next time called no parsing has to be done. 
	 * @return array with all model serials
	 */
	public Integer[] getModels() throws PdbLoadException {
		if (modelsArray!=null) {
			return modelsArray;
		}
		TreeSet<Integer> models = new TreeSet<Integer>();
		try {
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
					throw new FileFormatException("Incorrect number of fields for record "+recordCount+" in loop element "+atomSiteId+" of CIF file "+cifFile);
				}
				models.add(Integer.parseInt(tokens[pdbxPDBModelNumIdx]));
			}
			
		} catch (IOException e) {
			throw new PdbLoadException(e);
		} catch (FileFormatException e) {
			throw new PdbLoadException(e);
		}
		
		if (models.isEmpty()) return null;
		
		modelsArray = new Integer[models.size()];
		models.toArray(modelsArray);
		return modelsArray;
	}
	
	/*---------------------------- private methods --------------------------*/
	
	private void readFieldsTitles() throws IOException, FileFormatException {
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
			throw new FileFormatException("The file "+cifFile+" doesn't seem to be a cif file");
		}
		//int linecount = 1; // we have read one line already, we initialise count to 1
		// we need to store the last line's byte offset (which indicates the beginning of this line) 
		long lastLineOffset=fcif.getFilePointer();
		while((line = fcif.readLine()) != null ) {
			long currentOffset = fcif.getFilePointer(); //this gets byte offset at end of line
			//linecount++;
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
				p = Pattern.compile("^"+id+"\\.([\\w\\-\\[\\]]+)(?:\\s+(.*))?$");
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
			if (!line.startsWith("_") && !line.startsWith("#")){ // not in field definition, we are in values of a loop element or in a multi-line quoted value (';' quoting)
				if (ids2elements.containsValue(element)) { // if this is one of the fields we want to parse (members of String[] ids)
					if (!loopelements2contentOffset.containsKey(element)) {
						// this happens for both values in loop-elements and for multi-line quoted values in non-loop elements
						Long[] interval = {lastLineOffset, currentOffset};
						loopelements2contentOffset.put(element,interval);
					} else {
						// condition in next line is a bad hack: needed to be able to also parse 
						// multi-line quoted (';' quoting) values in non-loop elements with our 
						// tokeniser (for instance it happens often for _struct.title) 
						if (loopElements.contains(element)) { // i.e. we only do this if we really are in a loop-element (we could be in a non-loop with a multi-line quoted value) 
							loopelements2contentOffset.get(element)[1]=currentOffset; // i.e. if we are in a loop-element we keep extending the offset to the end of the last value
						}
					}
				}
			}
			lastLineOffset = currentOffset; //we store this line's offset to have it for next iteration
		} // end scanning lines
		
		fieldsTitlesRead = true;
	}
	
	protected String readPdbCode(){
		return fields2values.get(entryId+".id").trim().toLowerCase();
	}
	
	protected String readTitle() throws IOException, PdbLoadException {
		String title = fields2values.get(structId+".title").trim();
		if (title.isEmpty()) { 
			Long[] intStruct = loopelements2contentOffset.get(ids2elements.get(structId));
			fcif.seek(intStruct[0]);
			while(fcif.getFilePointer()<intStruct[1]) { 
				String[] tokens = tokeniseFields(1);
				if (tokens.length!=1) {
					throw new PdbLoadException("More than 1 field in element "+structId+" of CIF file "+cifFile);
				}
				title = tokens[0];
			}
		}

		if (title.charAt(0)=='\'' && title.charAt(title.length()-1)=='\'') {
			title = title.substring(1,title.length()-1);
		}

		return title;
	}
	
	protected CrystalCell readCrystalCell() throws PdbLoadException {
		// cell and symmetry fields are optional, e.g. NMR entries don't have them at all
		// we first see if the fields are there at all 
		if (!fields2values.containsKey(cell+".length_a")) {
			return null;
		}
		if (!isDouble(fields2values.get(cell+".length_a").trim()) ||
				!isDouble(fields2values.get(cell+".length_b").trim()) ||
				!isDouble(fields2values.get(cell+".length_c").trim()) ||
				!isDouble(fields2values.get(cell+".angle_alpha").trim()) ||
				!isDouble(fields2values.get(cell+".angle_beta").trim()) ||
				!isDouble(fields2values.get(cell+".angle_gamma").trim()) ) {
			// some NMR entries (e.g. 1fmm, 3iyx) do have cell params with a question mark value! we love the PDB!!
			return null;
		}
		
		CrystalCell crystalCell = null;
		double a = Double.parseDouble(fields2values.get(cell+".length_a").trim());
		double b = Double.parseDouble(fields2values.get(cell+".length_b").trim());
		double c = Double.parseDouble(fields2values.get(cell+".length_c").trim());
		double alpha = Double.parseDouble(fields2values.get(cell+".angle_alpha").trim());
		double beta = Double.parseDouble(fields2values.get(cell+".angle_beta").trim());
		double gamma = Double.parseDouble(fields2values.get(cell+".angle_gamma").trim());
		crystalCell = new CrystalCell(a, b, c, alpha, beta, gamma);
		return crystalCell;
	}
	
	protected SpaceGroup readSpaceGroup() throws PdbLoadException {
		// cell and symmetry fields are optional, e.g. NMR entries don't have them at all
		// we first see if the fields are there at all 
		if (!fields2values.containsKey(symmetry+".space_group_name_H-M")) {
			return null;
		}		
		SpaceGroup spaceGroup = null;
		// fields can be either single or double quoted: we have to remove 
		// both types of quotes (most pdbs use single, e.g. of double is 1hhu)
		String sg = fields2values.get(symmetry+".space_group_name_H-M").replaceAll("['\"]", "").trim();
		spaceGroup = SymoplibParser.getSpaceGroup(sg);
		if (spaceGroup==null) {
			throw new PdbLoadException("The space group found '"+sg+"' is not recognised as a standard space group");
		}
		return spaceGroup;
	}
	
	protected Matrix4d readScaleMatrix() throws PdbLoadException {

		if (!fields2values.containsKey(atomSitesId+".fract_transf_matrix[1][1]")) {
			return null;
		}
		
		Matrix4d scaleMatrix = new Matrix4d();
		scaleMatrix.m00 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_matrix[1][1]").trim());
		scaleMatrix.m01 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_matrix[1][2]").trim());
		scaleMatrix.m02 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_matrix[1][3]").trim());
		scaleMatrix.m10 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_matrix[2][1]").trim());
		scaleMatrix.m11 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_matrix[2][2]").trim());
		scaleMatrix.m12 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_matrix[2][3]").trim());
		scaleMatrix.m20 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_matrix[3][1]").trim());
		scaleMatrix.m21 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_matrix[3][2]").trim());
		scaleMatrix.m22 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_matrix[3][3]").trim());
		scaleMatrix.m03 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_vector[1]").trim());
		scaleMatrix.m13 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_vector[2]").trim());
		scaleMatrix.m23 = Double.parseDouble(fields2values.get(atomSitesId+".fract_transf_vector[3]").trim());

		return scaleMatrix;
	}
	
	protected String readExpMethod() throws IOException, PdbLoadException {
		String expMethod = null;
		if (!loopElements.contains(ids2elements.get(exptl))){  
			String expMethWithQuotes = fields2values.get(exptl+".method").trim();
			expMethod = expMethWithQuotes.substring(1,expMethWithQuotes.length()-1);
		} else {
			// normally _exptl.method is a single value element, but in some cases, e.g. 2krl, the _exptl.method is a loop element
			// in those cases we simply take the first value appearing (condition recordCount==1)
			Long[] intExptl = loopelements2contentOffset.get(ids2elements.get(exptl));
			int methodIdx = fields2indices.get(exptl+".method");
			int numberFields = ids2fieldsIdx.get(exptl);
			int recordCount = 0;
			fcif.seek(intExptl[0]);
			while(fcif.getFilePointer()<intExptl[1]) {
				recordCount++; 
				String[] tokens = tokeniseFields(numberFields);
				if (tokens.length!=numberFields) {
					throw new PdbLoadException("Incorrect number of fields for record "+recordCount+" in loop element "+exptl+" of CIF file "+cifFile);
				}
				if (recordCount==1) {
					expMethod = tokens[methodIdx];
				}
			}
		}
		return expMethod;
	}
	
	/**
	 * Returns an array of size 3 with the quality parameters for a crystal structure: resolution, rFree and rSym
	 * @return
	 */
	protected double[] readQparams() {
		double[] qParams = {-1,-1,-1};
		// refine only present in xray structures
		if (fields2values.containsKey(refine+".ls_d_res_high")) {
			String resolStr = fields2values.get(refine+".ls_d_res_high").trim();
			if ( isDouble(resolStr) ) {
				qParams[0] = Double.parseDouble(resolStr);
			}
		}
		
		if (fields2values.containsKey(refine+".ls_R_factor_R_free")) {
			String rfreeStr = fields2values.get(refine+".ls_R_factor_R_free").trim();
			if ( isDouble(rfreeStr) ) {
				qParams[1] = Double.parseDouble(rfreeStr);
			}
		}
		
		if (fields2values.containsKey(reflns+".pdbx_Rsym_value")){
			String rsymvalStr = fields2values.get(reflns+".pdbx_Rsym_value").trim();
			String rmergevalStr = null;
			if (fields2values.containsKey(reflns+".pdbx_Rmerge_I_obs")) {
				rmergevalStr = fields2values.get(reflns+".pdbx_Rmerge_I_obs").trim();}
			// if both are present, we don't compare them but take the Rsym value to be 
			// the right one (there's not much consensus in the field as to what's the 
			// right thing to do anyway!)
			
			if ( isDouble(rsymvalStr) && isDouble(rmergevalStr) ) {
				qParams[2] = Double.parseDouble(rsymvalStr);
			}
			else if(isDouble(rsymvalStr)){
				qParams[2] = Double.parseDouble(rsymvalStr);
			}
			else if(isDouble(rmergevalStr)){
				qParams[2] = Double.parseDouble(rmergevalStr);
			}
			
		}
		return qParams;
	}
	
	private void readAtomSite(PdbAsymUnit pdbAsymUnit, int model) throws IOException, FileFormatException, PdbLoadException {
		
		AtomLineList atomLines = new AtomLineList();
		
		Long[] intAtomSite = loopelements2contentOffset.get(ids2elements.get(atomSiteId));
		int groupPdbIdx = fields2indices.get(atomSiteId+".group_PDB");
		int idIdx = fields2indices.get(atomSiteId+".id");
		int typeSymbolIdx = fields2indices.get(atomSiteId+".type_symbol");
		int labelAtomIdIdx = fields2indices.get(atomSiteId+".label_atom_id");
		int labelAltIdIdx = fields2indices.get(atomSiteId+".label_alt_id");
		int labelCompIdIdx = fields2indices.get(atomSiteId+".label_comp_id");
		int labelAsymIdIdx = fields2indices.get(atomSiteId+".label_asym_id");
		int labelSeqIdIdx = fields2indices.get(atomSiteId+".label_seq_id");
		int authSeqIdIdx = fields2indices.get(atomSiteId+".auth_seq_id");
		int cartnXIdx = fields2indices.get(atomSiteId+".Cartn_x");
		int cartnYIdx = fields2indices.get(atomSiteId+".Cartn_y");
		int cartnZIdx = fields2indices.get(atomSiteId+".Cartn_z");
		int occupancyIdx = fields2indices.get(atomSiteId+".occupancy");
		int bIsoOrEquivIdx = fields2indices.get(atomSiteId+".B_iso_or_equiv");
		int pdbxPDBModelNumIdx = fields2indices.get(atomSiteId+".pdbx_PDB_model_num");
		int authAsymIdIdx = fields2indices.get(atomSiteId+".auth_asym_id");
		int numberFields = ids2fieldsIdx.get(atomSiteId);
		
		int recordCount = 0;
		
		fcif.seek(intAtomSite[0]);
		while(fcif.getFilePointer()<intAtomSite[1]) {
			recordCount++;

			// group_PDB=0, auth_asym_id=22, pdbx_PDB_model_num=24, label_alt_id=4, id=1, label_atom_id=3, label_comp_id=5, label_asym_id=6, label_seq_id=8, Cartn_x=10, Cartn_y=11, Cartn_z=12, occupancy=13, B_iso_or_equiv=14
			//   0   1    2  3  4   5 6 7 8  9     10    11       12    13    14 151617181920   2122 23 24
			//ATOM   2    C CA  . MET A 1 1  ? 38.591 8.543   15.660  1.00 77.79  ? ? ? ? ? 1  MET A CA  1
			String[] tokens = tokeniseFields(numberFields);
			if (tokens.length!=numberFields) {
				throw new FileFormatException("Incorrect number of fields for record "+recordCount+" in loop element "+atomSiteId+ " of CIF file "+cifFile);
			}
			if ((tokens[groupPdbIdx].equals("ATOM") || tokens[groupPdbIdx].equals("HETATM")) &&  
				Integer.parseInt(tokens[pdbxPDBModelNumIdx])==model) { // match our given chain and model 
				
				String asymId = tokens[labelAsymIdIdx];
				String labelAltId = tokens[labelAltIdIdx];
				int atomserial=Integer.parseInt(tokens[idIdx]); // id
				String element = tokens[typeSymbolIdx]; // type_symbol
				String atom = tokens[labelAtomIdIdx]; // label_atom_id
				String res_type = tokens[labelCompIdIdx]; // label_comp_id
				int resSerial = -1;
				if (!tokens[labelSeqIdIdx].equals(".")) {
					resSerial = Integer.parseInt(tokens[labelSeqIdIdx]); // label_seq_id
				}
				int pdbResSerial = Integer.parseInt(tokens[authSeqIdIdx]); 
				double x = Double.parseDouble(tokens[cartnXIdx]); // Cartn_x
				double y = Double.parseDouble(tokens[cartnYIdx]); // Cartn_y
				double z = Double.parseDouble(tokens[cartnZIdx]); // Cartn_z
				Point3d coords = new Point3d(x,y,z);
				double occupancy = Double.parseDouble(tokens[occupancyIdx]); // occupancy
				double bfactor = Double.parseDouble(tokens[bIsoOrEquivIdx]); // bfactor
				String authAsymId = tokens[authAsymIdIdx];
				if (!res_type.equals(HetResidue.WATER) && !res_type.equals(HetResidue.DEUT_WATER)) {
					// note we don't really use the insCode and nonPoly fields of AtomLine (we use them only in pdb file parser), we fill them with null and false
					atomLines.addAtomLine(new AtomLine(asymId, labelAltId, atomserial, atom, element, res_type, resSerial, pdbResSerial, null, coords, occupancy, bfactor, authAsymId, false, false));
				}
			}
		}
		if (atomLines.isEmpty()) { // no atom data was found for given pdb chain code and model
			throw new PdbLoadException("Couldn't find _atom_site data for given model: "+model);
		}
		
		String altLoc = atomLines.getAtomAltLoc();
		String lastChainCode = null;

		for (AtomLine atomLine:atomLines) {
			// we read only the alt locs we want
			if (altLoc!=null && !atomLine.labelAltId.equals(altLoc) && !atomLine.labelAltId.equals(".")) continue;
			
			if (lastChainCode!=null && !lastChainCode.equals(atomLine.labelAsymId)) {
				// in readPdbxPolySeq we already added the polymer chains, now we only need to add the missing ones: non-polymer chains
				if (!pdbAsymUnit.containsChainCode(atomLine.labelAsymId)) {
					PdbChain pdb = new PdbChain();
					pdb.setChainCode(atomLine.labelAsymId);
					pdb.setPdbChainCode(atomLine.authAsymId);
					pdbAsymUnit.setNonPolyChain(atomLine.labelAsymId, pdb);
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
					throw new PdbLoadException("atom_site sequence does not match sequence in pdbx_poly_seq_scheme for CIF chain "+pdb.getChainCode()+" at residue "+residue.getSerial()+" "+residue.getLongCode()+" (pdb serial: "+residue.getPdbSerial()+"). Inconsistent PDB entry. Report to the PDB.");
				}	
			}
		}

	}

	
	
	private void readPdbxPolySeq(PdbAsymUnit pdbAsymUnit) throws IOException, FileFormatException, PdbLoadException {
		
		PdbxPolySeqLineList list  = new PdbxPolySeqLineList();
		
		if(!ids2elements.containsKey(pdbxPolySeqId)) {
			throw new PdbLoadException("Missing pdbx_poly_seq_scheme field in mmCIF file. Is there no protein or nucleotide chains in this structure?");
		}
		
        	
		Long[] intPdbxPoly = loopelements2contentOffset.get(ids2elements.get(pdbxPolySeqId));
		int asymIdIdx = fields2indices.get(pdbxPolySeqId+".asym_id");
		int seqIdIdx = fields2indices.get(pdbxPolySeqId+".seq_id");
		//int authSeqNumIdx = fields2indices.get(pdbxPolySeqId+".auth_seq_num");
		int pdbSeqNumIdx = fields2indices.get(pdbxPolySeqId+".pdb_seq_num");
		int pdbInsCodeIdx = fields2indices.get(pdbxPolySeqId+".pdb_ins_code");
		int monIdIdx = fields2indices.get(pdbxPolySeqId+".mon_id");
		int pdbStrandIdIdx = fields2indices.get(pdbxPolySeqId+".pdb_strand_id");
		int numberFields = ids2fieldsIdx.get(pdbxPolySeqId);
		
		int recordCount=0;
		
		fcif.seek(intPdbxPoly[0]);
		while(fcif.getFilePointer()<intPdbxPoly[1]) {
			recordCount++; 

			// asym_id=0, seq_id=2, auth_seq_num=6, pdb_ins_code=10, mon_id=3 
			// 0 1 2     3 4   5   6     7   8 910
			// A 1 1   ASP 1   1   1   ASP ASP A .
			String[] tokens = tokeniseFields(numberFields);
			if (tokens.length!=numberFields) {
				throw new FileFormatException("Incorrect number of fields for record "+recordCount+" in loop element "+pdbxPolySeqId+" of CIF file "+cifFile);
			}			
			
			String asymId = tokens[asymIdIdx];
			String pdbChainCode = tokens[pdbStrandIdIdx];
			int resser = Integer.parseInt(tokens[seqIdIdx]); // seq_id
			String res_type = tokens[monIdIdx]; // mon_id
			int pdbSeqNum = Integer.parseInt(tokens[pdbSeqNumIdx]); // pdb_seq_num
			String pdb_ins_code = tokens[pdbInsCodeIdx]; // pdb_ins_code

			list.add(new PdbxPolySeqLine(asymId, resser, res_type, pdbSeqNum, pdbChainCode, pdb_ins_code));
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
        	pdb.setParent(pdbAsymUnit);
			pdb.setPdbChainCode(group.getPdbChainCode());			
			pdb.setChainCode(group.getChainCode());

            pdb.setSequence(group.getSequence(), group.isProtein());
 
            pdb.setPdbresser2resserMap(group.getPdbresser2resserMap());
			pdb.setResser2pdbresserMap(group.getResser2pdbresserMap());

			pdb.setIsNonPolyChain(false);
			
			pdbAsymUnit.setPolyChain(group.getChainCode(), pdb);
			
        }
	}
	
	private void readSecStructure(PdbAsymUnit pdbAsymUnit) throws IOException, FileFormatException, PdbLoadException {
		for (PdbChain pdb:pdbAsymUnit.getPolyChains()) {
			SecondaryStructure secondaryStructure = new SecondaryStructure(pdb.getSequence().getSeq());	// create empty secondary structure first to make sure object is not null
			pdb.setSecondaryStructure(secondaryStructure);
		}
		
		// struct_conf element is optional
		Long[] intStructConf = null;
		if (ids2elements.containsKey(structConfId)) {
			// if not a loop element then intStructConf stays null 
			if (loopElements.contains(ids2elements.get(structConfId))) {
				intStructConf = loopelements2contentOffset.get(ids2elements.get(structConfId));
			}
		} 
		// taking care of cases where struct_conf is not a loop element but a one value field
		if (ids2elements.containsKey(structConfId) && !loopElements.contains(ids2elements.get(structConfId))){  
			String begChainCode = fields2values.get(structConfId+".beg_label_asym_id").trim();
			String id = fields2values.get(structConfId+".id").trim();
			int beg = Integer.parseInt(fields2values.get(structConfId+".beg_label_seq_id").trim());
			int end = Integer.parseInt(fields2values.get(structConfId+".end_label_seq_id").trim());
			// we don't parse turns anymore as they were dropped from PDB files and anyway the annotation is VERY inconsistent
			if (!id.startsWith("TURN")) {
				Pattern p = Pattern.compile("^(\\w).+_P(\\d+)$");
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
					pdbAsymUnit.getChainForChainCode(begChainCode).getSecondaryStructure().add(ssElem);
				}
			}
		}

		// struct_sheet_range element is optional
		Long[] intStructSheet = null; 
		if (ids2elements.containsKey(structSheetId)) {
			// if not a loop element intStructSheet stays null 
			if (loopElements.contains(ids2elements.get(structSheetId))) {
				intStructSheet = loopelements2contentOffset.get(ids2elements.get(structSheetId));
			}
		}
		// taking care of cases where struct_sheet_range is not a loop element but a one value field
		if (ids2elements.containsKey(structSheetId) && !loopElements.contains(ids2elements.get(structSheetId))){
			String begChainCode = fields2values.get(structSheetId+".beg_label_asym_id").trim();
			String sheetid = fields2values.get(structSheetId+".sheet_id").trim(); //tokens[sheetIdIdx];
			int id = Integer.parseInt(fields2values.get(structSheetId+".id").trim()); //Integer.parseInt(tokens[idIdx]);
			int beg = Integer.parseInt(fields2values.get(structSheetId+".beg_label_seq_id").trim()); //tokens[begLabelSeqIdIdx]);
			int end = Integer.parseInt(fields2values.get(structSheetId+".end_label_seq_id").trim()); //tokens[endLabelSeqIdIdx]);
			String ssId=SecStrucElement.STRAND+sheetid+id; // e.g.: SA1, SA2..., SB1, SB2,...
			SecStrucElement ssElem = new SecStrucElement(SecStrucElement.STRAND, beg, end, ssId);
			pdbAsymUnit.getChainForChainCode(begChainCode).getSecondaryStructure().add(ssElem);


		}
				
		if (intStructConf!=null) {
			int idIdx = fields2indices.get(structConfId+".id");
			int begLabelAsymIdIdx = fields2indices.get(structConfId+".beg_label_asym_id");
			int begLabelSeqIdIdx = fields2indices.get(structConfId+".beg_label_seq_id");
			int endLabelSeqIdIdx = fields2indices.get(structConfId+".end_label_seq_id");
			int numFields = ids2fieldsIdx.get(structConfId);
			
			int recordCount=0;
			
			fcif.seek(intStructConf[0]);
			while(fcif.getFilePointer()<intStructConf[1]) {
				recordCount++;
				// struct_conf (optional element), HELIX and TURN secondary structure

				//id=1, beg_label_seq_id=5, end_label_seq_id=9, beg_label_asym_id=4
				//     0       1  2    3 4 5   6   7 8  9 10  111213    1415 16 1718 19
				//HELX_P HELX_P1  1  ASN A 2   ? GLY A 12  ? ASN A 2   GLY A 12  1 ? 11
				String[] tokens = tokeniseFields(numFields);
				if (tokens.length!=numFields) {
					throw new FileFormatException("Incorrect number of fields for record "+recordCount+" in loop element "+structConfId+" of CIF file "+cifFile);
				}
				String begChainCode = tokens[begLabelAsymIdIdx];
				String id = tokens[idIdx];
				if (id.startsWith("TURN")) continue; // we don't parse turns anymore as they were dropped from PDB files and anyway the annotation is VERY inconsistent
				Pattern p = Pattern.compile("^(\\w).+_P(\\d+)$");
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
					pdbAsymUnit.getChainForChainCode(begChainCode).getSecondaryStructure().add(ssElem);
				}
			}
		}
		if (intStructSheet!=null) {
			int sheetIdIdx = fields2indices.get(structSheetId+".sheet_id");
			int idIdx = fields2indices.get(structSheetId+".id");
			int begLabelAsymIdIdx = fields2indices.get(structSheetId+".beg_label_asym_id");
			int begLabelSeqIdIdx = fields2indices.get(structSheetId+".beg_label_seq_id");
			int endLabelSeqIdIdx = fields2indices.get(structSheetId+".end_label_seq_id");
			int numFields = ids2fieldsIdx.get(structSheetId);
			
			int recordCount=0;
			
			fcif.seek(intStructSheet[0]);
			while(fcif.getFilePointer()<intStructSheet[1]) {
				recordCount++;
				// struct_sheet_range (optional element), SHEETs
				//sheet_id=0, id=1, beg_label_seq_id=4, end_label_seq_id=8, beg_label_asym_id=3
				//0 1   2 3  4 5   6 7  8 910  1112 13  1415 16
				//A 1 ARG A 14 ? LYS A 19 ? ? ARG A 14 LYS A 19
				String[] tokens = tokeniseFields(numFields);
				if (tokens.length!=numFields) {
					throw new FileFormatException("Incorrect number of fields for record "+recordCount+" in loop element "+structSheetId+" of CIF file "+cifFile);
				}
				String begChainCode = tokens[begLabelAsymIdIdx];
				String sheetid = tokens[sheetIdIdx];
				int id = Integer.parseInt(tokens[idIdx]);
				int beg = Integer.parseInt(tokens[begLabelSeqIdIdx]);
				int end = Integer.parseInt(tokens[endLabelSeqIdIdx]);
				String ssId=SecStrucElement.STRAND+sheetid+id; // e.g.: SA1, SA2..., SB1, SB2,...
				SecStrucElement ssElem = new SecStrucElement(SecStrucElement.STRAND, beg, end, ssId);
				pdbAsymUnit.getChainForChainCode(begChainCode).getSecondaryStructure().add(ssElem);
	
			}
		}
		
		for (PdbChain pdb:pdbAsymUnit.getPolyChains()) {
			if(!pdb.getSecondaryStructure().isEmpty()) {
				if (!pdb.getSequence().isProtein()) throw new PdbLoadException("Secondary structure records present for a non-protein chain");
				pdb.getSecondaryStructure().setComment("CIFfile");
				pdb.initialiseResiduesSecStruct();
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
		char lastChar = 0; // ' '
		char quoteChar = 0;
		while (true) {
			char currentChar = (char)fcif.readByte();
			
			// '' quoting
			if (quoteChar!=';' && currentChar=='\'' && (lastChar==' ' || lastChar=='\n' || lastChar==0)){
				quoteChar = '\'';
			}
			else if (quoteChar!=';' && currentChar==' ' && lastChar=='\''){
				quoteChar = 0;
			}
			// "" quoting
			if (quoteChar!=';' && currentChar=='"' && (lastChar==' ' || lastChar=='\n' || lastChar==0)){
				quoteChar = '"';
			}
			else if (quoteChar!=';' && currentChar==' ' && lastChar=='"'){
				quoteChar = 0;
			}			
			// ;; quoting (multi-line quoting)
			if (quoteChar!=';' && currentChar==';' && (lastChar=='\n' || lastChar==0)){ // was && lastChar=='\n'
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
	
	private static boolean isDouble(String value)
	{	if(value==null) return false;
	    try{
	    	Double.parseDouble(value);
	    	return true;
	    }catch (NumberFormatException ex)
	    {	return false; }
	}
}

package owl.core.structure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.IOException;
import java.net.URL;
import java.net.URLConnection;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Locale;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;
import owl.core.structure.io.AtomLine;
import owl.core.structure.io.AtomLineList;
import owl.core.structure.io.BioUnitAssembly;
import owl.core.structure.io.BioUnitAssemblyGen;
import owl.core.structure.io.BioUnitOperation;
import owl.core.structure.io.CifFieldInfo;
import owl.core.structure.io.PdbxPolySeqGroup;
import owl.core.structure.io.PdbxPolySeqLine;
import owl.core.structure.io.PdbxPolySeqLineList;
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
	private static final String entry = "_entry";
	private static final String struct = "_struct";
	private static final String atomSite = "_atom_site";
	private static final String pdbxPolySeq = "_pdbx_poly_seq_scheme";
	private static final String structConf = "_struct_conf";
	private static final String structSheet = "_struct_sheet_range";
	private static final String cell = "_cell";
	private static final String symmetry = "_symmetry";
	private static final String exptl = "_exptl";
	private static final String reflns = "_reflns";
	private static final String refine = "_refine";
	private static final String atomSites = "_atom_sites";
	private static final String pdbxStructAssembly = "_pdbx_struct_assembly";
	private static final String pdbxStructAssemblyGen = "_pdbx_struct_assembly_gen";
	private static final String pdbxStructOperList = "_pdbx_struct_oper_list";
	private static final String databasePdbRev = "_database_PDB_rev";
	
	private static final String[] ids = 
		{entry, struct, atomSite, pdbxPolySeq,structConf, 
		structSheet, cell, symmetry, exptl, reflns, refine, atomSites, 
		pdbxStructAssembly, pdbxStructAssemblyGen, pdbxStructOperList, databasePdbRev};
	
	
	private HashMap<String, CifFieldInfo> fields;
	
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
		
		scanFile();
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
		scanFile();
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
				scanFile();
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
				scanFile();
			}
			
			CifFieldInfo pdbxPolySeqField = fields.get(pdbxPolySeq);
			
			while (pdbxPolySeqField.hasMoreData()) {
				String[] tokens = pdbxPolySeqField.getNextTokens();
				chains.add(tokens[pdbxPolySeqField.getIndexForSubField("pdb_strand_id")]);				
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
				scanFile();
			}
			
			CifFieldInfo atomSiteField = fields.get(atomSite);
			int modelIdx = atomSiteField.getIndexForSubField("pdbx_PDB_model_num");
			
			while(atomSiteField.hasMoreData()) {
				String[] tokens = atomSiteField.getNextTokens();
				models.add(Integer.parseInt(tokens[modelIdx]));
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
	
	/**
	 * Scanning of data structure in file, everything is stored into the fields map.
	 * Only data in ids array is read.
	 * @throws IOException
	 * @throws FileFormatException
	 */
	private void scanFile() throws IOException, FileFormatException {
		// initialising data structure to store the parsed fields
		fields = new HashMap<String, CifFieldInfo>();
		for (String id:ids){
			fields.put(id, new CifFieldInfo(id));
		}
		
		BufferedReader bf = new BufferedReader(new FileReader(cifFile));
		
		String line;
		line = bf.readLine(); // read first line
		Pattern p = Pattern.compile("^data_\\w+");
		if (!p.matcher(line).find()){
			bf.close();
			throw new FileFormatException("The file "+cifFile+" doesn't seem to be a cif file");
		}
		
		boolean inLoop = false;
		boolean inWantedData = false;
		String lastFieldId = null;
		String currentFieldId = null;

		while((line = bf.readLine()) != null ) {
			
			if (line.startsWith("#")) continue; // hashes are comments
			
			if (line.startsWith("loop_")) {
				inLoop = true;
				continue;
			}
			
			boolean isDataLine = false;
			
			p = Pattern.compile("^\\s*(_\\w+)\\..*$");
			Matcher m = p.matcher(line);
			if (m.find()) {
				currentFieldId = m.group(1);  				
			} else {
				// otherwise the line has values of loop element or is in a multi-line quoted value
				isDataLine = true;
			}

			if (!isDataLine) { // id line: we scan all the ids that we are interested in
				
				if (fields.containsKey(currentFieldId)) { // we only scan those we are interested in
					CifFieldInfo currentField = fields.get(currentFieldId);
					
					// setting inWantedData to true to flag that we want to parse any data lines coming within this field 
					inWantedData = true;

					// setting field to loop if seen in last line
					if (inLoop) {
						currentField.setLoop(inLoop);
						inLoop = false;
					}

					p = Pattern.compile("^\\s*"+currentFieldId+"\\.([\\w\\-\\[\\]]+)(?:\\s+(.*))?$");
					m = p.matcher(line);
					if (m.find()){
						currentField.addSubField(m.group(1));
						
						if (m.group(2)!=null) { // id there is data in same line we add it
							// note that we add the end of line so that strings are full, needed for the tokeniser method
							currentField.addDataToLastSubFieldDataBuffer(m.group(2)+"\n"); 
						}
					}
				} else {
					inLoop = false; // we need to reset inLoop in case we are in an unwanted field
				}

			} else { // data line 
			
				if (inWantedData) { // if this is one of the fields we want to parse 
					
					CifFieldInfo currentField = fields.get(currentFieldId);
					// note that we add the end of line so that strings are full, needed for the tokeniser method
					currentField.addDataToLastSubFieldDataBuffer(line+"\n");

				}
			}

			if (lastFieldId!=null && !lastFieldId.equals(currentFieldId)) {
				// we are in a new field (id):
				// resetting inWantedData to false
				inWantedData = false;
			}
			lastFieldId = currentFieldId;
		} 
		
		bf.close();
		
		for (CifFieldInfo field:fields.values()) {
			if (!field.checkConsistency()) throw new FileFormatException("Inconsistent data structure while scanning mmCIF file");
		}
		
		
		fieldsTitlesRead = true;
	}
	
	protected String readPdbCode() {

		CifFieldInfo entryField = fields.get(entry); 
		if (!entryField.isSubFieldPresent("id")) return PdbAsymUnit.NO_PDB_CODE;  
		
		return entryField.getSubFieldData("id").toLowerCase();
	}
	
	protected Date readReleaseDate() throws FileFormatException {
		CifFieldInfo databasePdbRevField = fields.get(databasePdbRev);
		if (databasePdbRevField.isEmpty()) return null;

		Date releaseDate = null;
		//for non-loop cases
		if (!databasePdbRevField.isLoop()) {
			
			if (!databasePdbRevField.isSubFieldPresent("date")) return null;
			
			String date = databasePdbRevField.getSubFieldData("date");
			try {
				releaseDate = new SimpleDateFormat("yyyy-MM-dd", Locale.ENGLISH).parse(date);
			} catch(ParseException pe){
				throw new FileFormatException("Release date is in unexpected format.");
			}
		} else {

			int numIdx = databasePdbRevField.getIndexForSubField("num");
			int dateIdx = databasePdbRevField.getIndexForSubField("date");
			
			while(databasePdbRevField.hasMoreData()) {  

				String[] tokens = databasePdbRevField.getNextTokens();
				
				String num = tokens[numIdx].trim();
				
				if (isInteger(num) && (Integer.parseInt(num) == 1)) {
					String date = tokens[dateIdx].trim();
					try {
						releaseDate = new SimpleDateFormat("yyyy-MM-dd", Locale.ENGLISH).parse(date);
					} catch(ParseException pe){
						throw new FileFormatException("Release date is in unexpected format.");
					}
				}
			}
		}
		
		return releaseDate;
	}
	
	protected String readTitle() {
		CifFieldInfo structField = fields.get(struct);
		if (structField.isEmpty()) return null;
		
		if (!structField.isSubFieldPresent("title")) return null; // subfield is missing
		
		String title = structField.getSubFieldData("title");
		
		return title;
	}
	
	protected CrystalCell readCrystalCell() {
		CifFieldInfo cellField = fields.get(cell);
		// cell and symmetry fields are optional, e.g. NMR entries don't have them at all
		// we first see if the fields are there at all 
		if (cellField.isEmpty()) return null;


		String aStr = cellField.getSubFieldData("length_a");
		String bStr = cellField.getSubFieldData("length_b");
		String cStr = cellField.getSubFieldData("length_c");
		String alphaStr = cellField.getSubFieldData("angle_alpha");
		String betaStr = cellField.getSubFieldData("angle_beta");
		String gammaStr = cellField.getSubFieldData("angle_gamma");

		
		if (!isDouble(aStr) ||
			!isDouble(bStr) ||
			!isDouble(cStr) ||
			!isDouble(alphaStr) ||
			!isDouble(betaStr) ||
			!isDouble(gammaStr)) {
			// some NMR entries (e.g. 1fmm, 3iyx) do have cell params with a question mark value! we love the PDB!!
			return null;
		}
		
		CrystalCell crystalCell = null;
		double a = Double.parseDouble(aStr);
		double b = Double.parseDouble(bStr);
		double c = Double.parseDouble(cStr);
		double alpha = Double.parseDouble(alphaStr);
		double beta = Double.parseDouble(betaStr);
		double gamma = Double.parseDouble(gammaStr);
		crystalCell = new CrystalCell(a, b, c, alpha, beta, gamma);
		return crystalCell;
	}
	
	protected SpaceGroup readSpaceGroup() throws PdbLoadException {
		CifFieldInfo symmetryField = fields.get(symmetry);
		// cell and symmetry fields are optional, e.g. NMR entries don't have them at all
		// we first see if the fields are there at all 
		if (symmetryField.isEmpty()) return null;
		
		SpaceGroup spaceGroup = null;

		String sg = symmetryField.getSubFieldData("space_group_name_H-M");
		spaceGroup = SymoplibParser.getSpaceGroup(sg);
		if (spaceGroup==null) {
			throw new PdbLoadException("The space group found '"+sg+"' is not recognised as a standard space group");
		}
		return spaceGroup;
	}
	
	protected Matrix4d readScaleMatrix() {
		CifFieldInfo atomSitesField = fields.get(atomSites);
		if (atomSitesField.isEmpty()) return null;
		
		Matrix4d scaleMatrix = new Matrix4d();
		scaleMatrix.m00 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_matrix[1][1]"));
		scaleMatrix.m01 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_matrix[1][2]"));
		scaleMatrix.m02 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_matrix[1][3]"));
		scaleMatrix.m10 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_matrix[2][1]"));
		scaleMatrix.m11 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_matrix[2][2]"));
		scaleMatrix.m12 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_matrix[2][3]"));
		scaleMatrix.m20 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_matrix[3][1]"));
		scaleMatrix.m21 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_matrix[3][2]"));
		scaleMatrix.m22 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_matrix[3][3]"));
		scaleMatrix.m03 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_vector[1]"));
		scaleMatrix.m13 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_vector[2]"));
		scaleMatrix.m23 = Double.parseDouble(atomSitesField.getSubFieldData("fract_transf_vector[3]"));

		return scaleMatrix;
	}
	
	protected String readExpMethod() {
		CifFieldInfo exptlField = fields.get(exptl);
		
		String expMethod = null;
		
		if (!exptlField.isLoop()) {
			expMethod = exptlField.getSubFieldData("method");
			
		} else {
			// normally _exptl.method is a single value element, but in some cases, e.g. 2krl, the _exptl.method is a loop element
			// in those cases we simply take the first value appearing (condition recordCount==1)

			int methodIdx = exptlField.getIndexForSubField("method");
			
			while(exptlField.hasMoreData()) {

				String[] tokens = exptlField.getNextTokens();
				
				expMethod = tokens[methodIdx];
				
				// this means we read only first method and ignore the rest!
				break;
			}
		}
		
		return expMethod;
	}
	
	/**
	 * Returns an array of size 3 with the quality parameters for a crystal structure: resolution, rFree and rSym
	 * @return
	 */
	protected double[] readQparams() {
		CifFieldInfo refineField = fields.get(refine);
		CifFieldInfo reflnsField = fields.get(reflns);
		
		double[] qParams = {-1,-1,-1};
		
		// refine only present in xray structures
		if (!refineField.isEmpty()) {
			
			// note: refine can be a loop in some cases e.g. 3otj
			if (!refineField.isLoop()) {
				if (refineField.isSubFieldPresent("ls_d_res_high")) {
					String resolStr = refineField.getSubFieldData("ls_d_res_high");
					if ( isDouble(resolStr) ) {
						qParams[0] = Double.parseDouble(resolStr);
					}
				}

				if (refineField.isSubFieldPresent("ls_R_factor_R_free")) {
					String rfreeStr = refineField.getSubFieldData("ls_R_factor_R_free");
					if ( isDouble(rfreeStr) ) {
						qParams[1] = Double.parseDouble(rfreeStr);
					}
				}
			} else {
				int reshighIdx = refineField.getIndexForSubField("ls_d_res_high");
				int rfactorrfreeIdx = refineField.getIndexForSubField("ls_R_factor_R_free");
				
				while(refineField.hasMoreData()) {
					// NOTE we take first only and ignore the rest
					String[] tokens = refineField.getNextTokens();
					String resolStr = tokens[reshighIdx];
					String rfreeStr = tokens[rfactorrfreeIdx];
					if ( isDouble(resolStr) ) {
						qParams[0] = Double.parseDouble(resolStr);
					}
					if ( isDouble(rfreeStr) ) {
						qParams[1] = Double.parseDouble(rfreeStr);
					}
					break;
				}
			}
		}
		
		if (!reflnsField.isEmpty()) {
			String rsymvalStr = null;
			if (reflnsField.isSubFieldPresent("pdbx_Rsym_value")) 
				rsymvalStr = reflnsField.getSubFieldData("pdbx_Rsym_value");
			String rmergevalStr = null;
			if (reflnsField.isSubFieldPresent("pdbx_Rmerge_I_obs")) 
				rmergevalStr = reflnsField.getSubFieldData("pdbx_Rmerge_I_obs");

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
	
	/**
	 * Method to read in the values of temporary biounit classes after parsing
	 * @param pdbAsymUnit
	 * @throws IOException
	 * @throws FileFormatException
	 */
	protected void readBioUnit(PdbAsymUnit pdbAsymUnit) {
		ArrayList<BioUnitAssembly> assemblies = new ArrayList<BioUnitAssembly>();
		ArrayList<BioUnitAssemblyGen> generators = new ArrayList<BioUnitAssemblyGen>();
		ArrayList<BioUnitOperation> operations = new ArrayList<BioUnitOperation>();
		
		//**********************************************
		// _pdb_struct_assembly fields
		//**********************************************
		CifFieldInfo pdbxStructAssemblyField = fields.get(pdbxStructAssembly);
		
		if (pdbxStructAssemblyField.isEmpty()) return;
		
		//for non-loop cases
		if (!pdbxStructAssemblyField.isLoop()){
			BioUnitAssembly assembly = new BioUnitAssembly();
			//Set entry number
			String id = pdbxStructAssemblyField.getSubFieldData("id");
			if(isInteger(id)) {
				int i = Integer.parseInt(id);
				assembly.setId(i);
			}
			
			boolean isSoftwarePresent  = pdbxStructAssemblyField.isSubFieldPresent("method_details");
			boolean isCountPresent = pdbxStructAssemblyField.isSubFieldPresent("oligomeric_count");
			boolean isOligomerPresent = pdbxStructAssemblyField.isSubFieldPresent("oligomeric_details");
			String details=null; String method=null; String sizeStr=null;
			
			//Add assembly details
			details = pdbxStructAssemblyField.getSubFieldData("details").toLowerCase();
			if(isSoftwarePresent) method = pdbxStructAssemblyField.getSubFieldData("method_details").toLowerCase();
			if(isCountPresent) sizeStr = pdbxStructAssemblyField.getSubFieldData("oligomeric_count");
			
			//Check if the size contained in .oligomeric_count is a valid size; otherwise take it from .oligomeric_details field
			if(isCountPresent && isInteger(sizeStr)){
				int size = Integer.parseInt(sizeStr);
				assembly.setSize(size);
			}else if(isOligomerPresent){
				sizeStr = pdbxStructAssemblyField.getSubFieldData("oligomeric_details").toLowerCase();
				assembly.setSize(sizeStr);
			}
			//Add details on the type of assignments
			if(details.contains("author")) assembly.addType("authors");
			if(isSoftwarePresent && !method.contains("?")){
				String[] software = method.split(",");
				for(String soft:software) assembly.addType(soft.trim().toLowerCase());
			}
						
			
			//Add to the list of entries
			assemblies.add(assembly);
			
		} else { //for loop cases
			
			//boolean to check if fields are present
			boolean isSoftwarePresent  = pdbxStructAssemblyField.isSubFieldPresent("method_details");
			boolean isCountPresent = pdbxStructAssemblyField.isSubFieldPresent("oligomeric_count");
			boolean isOligomerPresent = pdbxStructAssemblyField.isSubFieldPresent("oligomeric_details");
			int softwareIdx=0;
			int sizeCountIdx=0;
			int sizeCountStrIdx=0;
			
			int idIdx = pdbxStructAssemblyField.getIndexForSubField("id");
			int detailsIdx = pdbxStructAssemblyField.getIndexForSubField("details");
			if(isSoftwarePresent) softwareIdx = pdbxStructAssemblyField.getIndexForSubField("method_details");
			if(isCountPresent) sizeCountIdx = pdbxStructAssemblyField.getIndexForSubField("oligomeric_count");
			if(isOligomerPresent) sizeCountStrIdx = pdbxStructAssemblyField.getIndexForSubField("oligomeric_details");

			//Read _pdbx_struct_assembly block values
			
			while(pdbxStructAssemblyField.hasMoreData()) {
				
				BioUnitAssembly assembly = new BioUnitAssembly();
				String[] tokens = pdbxStructAssemblyField.getNextTokens();
				
				//Add assembly details
				
				//1. Set Id
				String id = tokens[idIdx].trim();
				if(isInteger(id)) {
					int i = Integer.parseInt(id);
					assembly.setId(i);
				}
				
				//2. Set size
				if(isCountPresent){
					String sizeStr = tokens[sizeCountIdx].trim();
					if(isInteger(sizeStr)){
						int size = Integer.parseInt(sizeStr);
						assembly.setSize(size);
					}else if(isOligomerPresent){
						sizeStr = tokens[sizeCountStrIdx].trim().toLowerCase();
						assembly.setSize(sizeStr);
					}
				} else if(isOligomerPresent){
					String sizeStr = tokens[sizeCountStrIdx].trim().toLowerCase();
					assembly.setSize(sizeStr);
				}
				//3. set details on type of assignments
				String details = tokens[detailsIdx].trim().toLowerCase();
				if(details.contains("author")) assembly.addType("authors");
				
				if(isSoftwarePresent){
					String method = tokens[softwareIdx].trim().toLowerCase();
					if(!method.contains("?")){
						String[] software = method.split(",");
						for(String soft:software) assembly.addType(soft.replaceAll("\\s+","").toLowerCase());
					}
				}
				//Finally add the assembly
				assemblies.add(assembly);
			}
			
			
		}
		
		//**********************************************
		// _pdb_struct_assembly_gen fields
		//**********************************************
		CifFieldInfo pdbxStructAssemblyGenField = fields.get(pdbxStructAssemblyGen);
		
		if (pdbxStructAssemblyGenField.isEmpty()) return;
				
		//Check if the _pdb_struct_assembly_gen fields are in loop or not
		//for non-loop cases
		if(!pdbxStructAssemblyGenField.isLoop()){
			BioUnitAssemblyGen gen = new BioUnitAssemblyGen(1);
			
			String assemblyId = pdbxStructAssemblyGenField.getSubFieldData("assembly_id");
			String operatorIds = pdbxStructAssemblyGenField.getSubFieldData("oper_expression").toLowerCase();
			
			String chainIds = pdbxStructAssemblyGenField.getSubFieldData("asym_id_list"); 
			
			//Add details
			
			//1. Add the generator Id to the corresponding assembly
			for(BioUnitAssembly assembly:assemblies)
				if(isInteger(assemblyId) && Integer.parseInt(assemblyId) == assembly.getId()){
					assembly.addGeneratorId(gen.getId());
					break;
				}
				
			//2. Add the corresponding operator Id's
			String[] oids = operatorIds.split(",");
			for(String oid:oids)
				if(isInteger(oid))
					gen.addOperationId(Integer.parseInt(oid.trim()));
			
			//3. Add corresponding chain id's
			String[] chainCodes = chainIds.split(",");
			gen.addPdbChainCodes(chainCodes);
			
			//Finally add the generator
			generators.add(gen);
			
		} else { //for loop cases
			
			int genidIdx = pdbxStructAssemblyGenField.getIndexForSubField("assembly_id");
			int operidIdx = pdbxStructAssemblyGenField.getIndexForSubField("oper_expression");
			int chainIdx = pdbxStructAssemblyGenField.getIndexForSubField("asym_id_list");
			
			BioUnitAssemblyGen gen = new BioUnitAssemblyGen();

			//Read _pdbx_struct_assembly_gen values

			while(pdbxStructAssemblyGenField.hasMoreData()) {

				String[] tokens = pdbxStructAssemblyGenField.getNextTokens();
				
				gen = new BioUnitAssemblyGen(gen.getId()+1);
				
				
				String assemblyId = tokens[genidIdx].trim();
				String operatorIds = tokens[operidIdx].trim().toLowerCase();
				String chainIds = tokens[chainIdx].trim();
				
				//Add details
				
				//1. Add the generator Id to the corresponding assembly
				for(BioUnitAssembly assembly:assemblies)
					if(isInteger(assemblyId) && Integer.parseInt(assemblyId) == assembly.getId()){
						assembly.addGeneratorId(gen.getId());
						break;
					}
					
					
				//2. Add the corresponding operator Id's
				String[] oids = operatorIds.split(",");
				for(String oid:oids)
					if(isInteger(oid))
						gen.addOperationId(Integer.parseInt(oid.trim()));
				
				//3. Add corresponding chain id's
				String[] chainCodes = chainIds.split(",");
				gen.addPdbChainCodes(chainCodes);
				
				//Finally add the generator
				generators.add(gen);
				
			}
		}
		
		
		//**********************************************
		// _pdb_struct_oper_list fields
		//**********************************************
		CifFieldInfo pdbxStructOperListField = fields.get(pdbxStructOperList);
		
		if (pdbxStructOperListField.isEmpty()) return;
		
		//Check if the fields are in loop or not
		//for non-loop cases
		if(!pdbxStructOperListField.isLoop()){
			BioUnitOperation operation = new BioUnitOperation();

			String idStr = pdbxStructOperListField.getSubFieldData("id");
			String m00 = pdbxStructOperListField.getSubFieldData("matrix[1][1]");
			String m01 = pdbxStructOperListField.getSubFieldData("matrix[1][2]");
			String m02 = pdbxStructOperListField.getSubFieldData("matrix[1][3]");
			String m03 = pdbxStructOperListField.getSubFieldData("vector[1]");
			String m10 = pdbxStructOperListField.getSubFieldData("matrix[2][1]");
			String m11 = pdbxStructOperListField.getSubFieldData("matrix[2][2]");
			String m12 = pdbxStructOperListField.getSubFieldData("matrix[2][3]");
			String m13 = pdbxStructOperListField.getSubFieldData("vector[2]");
			String m20 = pdbxStructOperListField.getSubFieldData("matrix[3][1]");
			String m21 = pdbxStructOperListField.getSubFieldData("matrix[3][2]");
			String m22 = pdbxStructOperListField.getSubFieldData("matrix[3][3]");
			String m23 = pdbxStructOperListField.getSubFieldData("vector[3]");


			if(isInteger(idStr))
				//System.err.println("Warning: Encountered a non-integer id in "+pdbxStructOperList+" ; will not add it");
			//else
				if(isDouble(m00) && isDouble(m01) && isDouble(m02) && isDouble(m03) &&
					isDouble(m10) && isDouble(m11) && isDouble(m12) && isDouble(m13) &&
					isDouble(m20) && isDouble(m21) && isDouble(m22) && isDouble(m23) ){
				int id = Integer.parseInt(idStr);
				operation.setId(id);

				operation.setOperatorValue(0, Double.parseDouble(m00));
				operation.setOperatorValue(1, Double.parseDouble(m01));
				operation.setOperatorValue(2, Double.parseDouble(m02));
				operation.setOperatorValue(3, Double.parseDouble(m03));
				operation.setOperatorValue(4, Double.parseDouble(m10));
				operation.setOperatorValue(5, Double.parseDouble(m11));
				operation.setOperatorValue(6, Double.parseDouble(m12));
				operation.setOperatorValue(7, Double.parseDouble(m13));
				operation.setOperatorValue(8, Double.parseDouble(m20));
				operation.setOperatorValue(9, Double.parseDouble(m21));
				operation.setOperatorValue(10, Double.parseDouble(m22));
				operation.setOperatorValue(11, Double.parseDouble(m23));

				operations.add(operation);
			}//else
				//System.err.println("Warning: Encountered a non-numeric field in "+pdbxStructOperList+" matrix list; will not add it");
			
			
		} else { //for loop-elements

			int idIdx = pdbxStructOperListField.getIndexForSubField("id");
			int m00Idx = pdbxStructOperListField.getIndexForSubField("matrix[1][1]");
			int m01Idx = pdbxStructOperListField.getIndexForSubField("matrix[1][2]");
			int m02Idx = pdbxStructOperListField.getIndexForSubField("matrix[1][3]");
			int m03Idx = pdbxStructOperListField.getIndexForSubField("vector[1]");
			int m10Idx = pdbxStructOperListField.getIndexForSubField("matrix[2][1]");
			int m11Idx = pdbxStructOperListField.getIndexForSubField("matrix[2][2]");
			int m12Idx = pdbxStructOperListField.getIndexForSubField("matrix[2][3]");
			int m13Idx = pdbxStructOperListField.getIndexForSubField("vector[2]");
			int m20Idx = pdbxStructOperListField.getIndexForSubField("matrix[3][1]");
			int m21Idx = pdbxStructOperListField.getIndexForSubField("matrix[3][2]");
			int m22Idx = pdbxStructOperListField.getIndexForSubField("matrix[3][3]");
			int m23Idx = pdbxStructOperListField.getIndexForSubField("vector[3]");


			//Read _pdbx_struct_oper_list block values
			while(pdbxStructOperListField.hasMoreData()) {
				
				BioUnitOperation operation = new BioUnitOperation();
				String[] tokens = pdbxStructOperListField.getNextTokens();
				
				String idStr = tokens[idIdx].trim();
				String m00 = tokens[m00Idx].trim();
				String m01 = tokens[m01Idx].trim();
				String m02 = tokens[m02Idx].trim();
				String m03 = tokens[m03Idx].trim();
				String m10 = tokens[m10Idx].trim();
				String m11 = tokens[m11Idx].trim();
				String m12 = tokens[m12Idx].trim();
				String m13 = tokens[m13Idx].trim();
				String m20 = tokens[m20Idx].trim();
				String m21 = tokens[m21Idx].trim();
				String m22 = tokens[m22Idx].trim();
				String m23 = tokens[m23Idx].trim();
				
				//Add data
				if(isInteger(idStr))
					//System.err.println("Warning: Encountered a non-integer id in "+pdbxStructOperList+" ; will not add it");
				//else 
					if(isDouble(m00) && isDouble(m01) && isDouble(m02) && isDouble(m03) &&
						isDouble(m10) && isDouble(m11) && isDouble(m12) && isDouble(m13) &&
						isDouble(m20) && isDouble(m21) && isDouble(m22) && isDouble(m23) ){
					int id = Integer.parseInt(idStr);
					operation.setId(id);

					operation.setOperatorValue(0, Double.parseDouble(m00));
					operation.setOperatorValue(1, Double.parseDouble(m01));
					operation.setOperatorValue(2, Double.parseDouble(m02));
					operation.setOperatorValue(3, Double.parseDouble(m03));
					operation.setOperatorValue(4, Double.parseDouble(m10));
					operation.setOperatorValue(5, Double.parseDouble(m11));
					operation.setOperatorValue(6, Double.parseDouble(m12));
					operation.setOperatorValue(7, Double.parseDouble(m13));
					operation.setOperatorValue(8, Double.parseDouble(m20));
					operation.setOperatorValue(9, Double.parseDouble(m21));
					operation.setOperatorValue(10, Double.parseDouble(m22));
					operation.setOperatorValue(11, Double.parseDouble(m23));

					operations.add(operation);
				}//else
				//	System.err.println("Warning: Encountered a non-numeric field in "+pdbxStructOperList+" matrix list; will not add it");
				
			}
		}
		
		pdbAsymUnit.setPdbBioUnitList(new PdbBioUnitList(pdbAsymUnit, assemblies, generators, operations,"cif"));
		
	}
	
	private void readAtomSite(PdbAsymUnit pdbAsymUnit, int model) throws PdbLoadException {
		
		AtomLineList atomLines = new AtomLineList();
		
		CifFieldInfo atomSiteField = fields.get(atomSite);
		
		int groupPdbIdx = atomSiteField.getIndexForSubField("group_PDB");
		int idIdx = atomSiteField.getIndexForSubField("id");
		int typeSymbolIdx = atomSiteField.getIndexForSubField("type_symbol");
		int labelAtomIdIdx = atomSiteField.getIndexForSubField("label_atom_id");
		int labelAltIdIdx = atomSiteField.getIndexForSubField("label_alt_id");
		int labelCompIdIdx = atomSiteField.getIndexForSubField("label_comp_id");
		int labelAsymIdIdx = atomSiteField.getIndexForSubField("label_asym_id");
		int labelSeqIdIdx = atomSiteField.getIndexForSubField("label_seq_id");
		int authSeqIdIdx = atomSiteField.getIndexForSubField("auth_seq_id");
		int cartnXIdx = atomSiteField.getIndexForSubField("Cartn_x");
		int cartnYIdx = atomSiteField.getIndexForSubField("Cartn_y");
		int cartnZIdx = atomSiteField.getIndexForSubField("Cartn_z");
		int occupancyIdx = atomSiteField.getIndexForSubField("occupancy");
		int bIsoOrEquivIdx = atomSiteField.getIndexForSubField("B_iso_or_equiv");
		int pdbxPDBModelNumIdx = atomSiteField.getIndexForSubField("pdbx_PDB_model_num");
		int authAsymIdIdx = atomSiteField.getIndexForSubField("auth_asym_id");
		
		while(atomSiteField.hasMoreData()) {

			// group_PDB=0, auth_asym_id=22, pdbx_PDB_model_num=24, label_alt_id=4, id=1, label_atom_id=3, label_comp_id=5, label_asym_id=6, label_seq_id=8, Cartn_x=10, Cartn_y=11, Cartn_z=12, occupancy=13, B_iso_or_equiv=14
			//   0   1    2  3  4   5 6 7 8  9     10    11       12    13    14 151617181920   2122 23 24
			//ATOM   2    C CA  . MET A 1 1  ? 38.591 8.543   15.660  1.00 77.79  ? ? ? ? ? 1  MET A CA  1
			String[] tokens = atomSiteField.getNextTokens();

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

	
	
	private void readPdbxPolySeq(PdbAsymUnit pdbAsymUnit) throws PdbLoadException {
		
		PdbxPolySeqLineList list  = new PdbxPolySeqLineList();
		
		CifFieldInfo pdbxPolySeqField = fields.get(pdbxPolySeq);
		if (pdbxPolySeqField.isEmpty()) {
			throw new PdbLoadException("Missing pdbx_poly_seq_scheme field in mmCIF file. Is there no protein or nucleotide chains in this structure?");
		}
        	
		int asymIdIdx = pdbxPolySeqField.getIndexForSubField("asym_id");
		int seqIdIdx = pdbxPolySeqField.getIndexForSubField("seq_id");
		//int authSeqNumIdx = pdbxPolySeqField.getIndexForSubField("auth_seq_num");
		int pdbSeqNumIdx = pdbxPolySeqField.getIndexForSubField("pdb_seq_num");
		int pdbInsCodeIdx = pdbxPolySeqField.getIndexForSubField("pdb_ins_code");
		int monIdIdx = pdbxPolySeqField.getIndexForSubField("mon_id");
		int pdbStrandIdIdx = pdbxPolySeqField.getIndexForSubField("pdb_strand_id");

		while(pdbxPolySeqField.hasMoreData()) {

			// asym_id=0, seq_id=2, auth_seq_num=6, pdb_ins_code=10, mon_id=3 
			// 0 1 2     3 4   5   6     7   8 910
			// A 1 1   ASP 1   1   1   ASP ASP A .
			String[] tokens = pdbxPolySeqField.getNextTokens();
			
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
	
	private void readSecStructure(PdbAsymUnit pdbAsymUnit) throws PdbLoadException {
		for (PdbChain pdb:pdbAsymUnit.getPolyChains()) {
			SecondaryStructure secondaryStructure = new SecondaryStructure(pdb.getSequence().getSeq());	// create empty secondary structure first to make sure object is not null
			pdb.setSecondaryStructure(secondaryStructure);
		}
		
		CifFieldInfo structConfField = fields.get(structConf);
		
		// struct_conf element is optional
		if (!structConfField.isEmpty()) {
			
			// if present it can be loop or non-loop
			
			// non-loop
			if (!structConfField.isLoop()) {
				String begChainCode = structConfField.getSubFieldData("beg_label_asym_id");
				String id = structConfField.getSubFieldData("id");
				int beg = Integer.parseInt(structConfField.getSubFieldData("beg_label_seq_id"));
				int end = Integer.parseInt(structConfField.getSubFieldData("end_label_seq_id"));
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

			// loop
			} else {
				int idIdx = structConfField.getIndexForSubField("id");
				int begLabelAsymIdIdx = structConfField.getIndexForSubField("beg_label_asym_id");
				int begLabelSeqIdIdx = structConfField.getIndexForSubField("beg_label_seq_id");
				int endLabelSeqIdIdx = structConfField.getIndexForSubField("end_label_seq_id");
				
				while(structConfField.hasMoreData()) {

					// struct_conf (optional element), HELIX and TURN secondary structure

					//id=1, beg_label_seq_id=5, end_label_seq_id=9, beg_label_asym_id=4
					//     0       1  2    3 4 5   6   7 8  9 10  111213    1415 16 1718 19
					//HELX_P HELX_P1  1  ASN A 2   ? GLY A 12  ? ASN A 2   GLY A 12  1 ? 11
					String[] tokens = structConfField.getNextTokens();

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
		}
		
		
		CifFieldInfo structSheetField = fields.get(structSheet);

		// struct_sheet_range element is optional, when present it can be loop or non-loop		
		if (!structSheetField.isEmpty()) {
			
			// non-loop
			if (!structSheetField.isLoop()) {
				String begChainCode = structSheetField.getSubFieldData("beg_label_asym_id");
				String sheetid = structSheetField.getSubFieldData("sheet_id");
				int id = Integer.parseInt(structSheetField.getSubFieldData("id"));
				int beg = Integer.parseInt(structSheetField.getSubFieldData("beg_label_seq_id"));
				int end = Integer.parseInt(structSheetField.getSubFieldData("end_label_seq_id"));
				String ssId=SecStrucElement.STRAND+sheetid+id; // e.g.: SA1, SA2..., SB1, SB2,...
				SecStrucElement ssElem = new SecStrucElement(SecStrucElement.STRAND, beg, end, ssId);
				pdbAsymUnit.getChainForChainCode(begChainCode).getSecondaryStructure().add(ssElem);
				
			// loop
			} else {
				int sheetIdIdx = structSheetField.getIndexForSubField("sheet_id");
				int idIdx = structSheetField.getIndexForSubField("id");
				int begLabelAsymIdIdx = structSheetField.getIndexForSubField("beg_label_asym_id");
				int begLabelSeqIdIdx = structSheetField.getIndexForSubField("beg_label_seq_id");
				int endLabelSeqIdIdx = structSheetField.getIndexForSubField("end_label_seq_id");
				
				while(structSheetField.hasMoreData()) {
					
					// struct_sheet_range (optional element), SHEETs
					//sheet_id=0, id=1, beg_label_seq_id=4, end_label_seq_id=8, beg_label_asym_id=3
					//0 1   2 3  4 5   6 7  8 910  1112 13  1415 16
					//A 1 ARG A 14 ? LYS A 19 ? ? ARG A 14 LYS A 19
					String[] tokens = structSheetField.getNextTokens();

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
		}
		
		for (PdbChain pdb:pdbAsymUnit.getPolyChains()) {
			if(!pdb.getSecondaryStructure().isEmpty()) {
				if (!pdb.getSequence().isProtein()) 
					throw new PdbLoadException("Secondary structure records present for a non-protein chain");
				pdb.getSecondaryStructure().setComment("CIFfile");
				pdb.initialiseResiduesSecStruct();
			}
		}
	}
	
	private static boolean isDouble(String value) {	
		if(value==null) return false;
		try{
			Double.parseDouble(value);
			return true;
		} catch (NumberFormatException ex) {	
			return false; 
		}
	}
	
	private static boolean isInteger(String value) {
		if (value==null) return false;
		
	    try { 
	    	Integer.parseInt(value);
	    	return true;
	    }
	    catch (NumberFormatException ex) {	
	    	return false; 
	    }
	}
	
	public static void main(String[] args) throws Exception {
		
		// 4AX3, 1BXY for multi-line title
		CiffileParser parser = new CiffileParser(new File("/home/duarte_j/1bxy.cif"));
		
		for (CifFieldInfo field:parser.fields.values()) {
			System.out.println(field.getFieldName()+" "+field.getNumSubFields()+" "+(field.isLoop()?"loop":""));
		}
		
		String[] chains = parser.getChains();
		for (String chain:chains) {
			System.out.print(chain+" ");
		}		
		System.out.println();

		Integer[] models = parser.getModels();
		for (int model:models) {
			System.out.print(model+" ");
		}		
		System.out.println();

		
		System.out.println("-->"+parser.readPdbCode()+"<--");
		System.out.println("-->"+parser.readReleaseDate()+"<--");
		System.out.println("-->"+parser.readTitle()+"<--");
		CrystalCell cell = parser.readCrystalCell();
		if (cell==null) System.out.println("null cell");
		else System.out.println(cell.getA()+" "+cell.getB()+" "+cell.getC()+" "+cell.getAlpha()+" "+cell.getBeta()+" "+cell.getGamma());
		SpaceGroup sg = parser.readSpaceGroup();
		if (sg==null) System.out.println("null SG");
		else System.out.println(parser.readSpaceGroup().getShortSymbol());		
		System.out.println(parser.readScaleMatrix());
		System.out.println("-->"+parser.readExpMethod()+"<--");
		System.out.println("q params: "+parser.readQparams()[0]+" "+parser.readQparams()[1]+" "+parser.readQparams()[2]);
		
		System.out.println("**done");
	}
}

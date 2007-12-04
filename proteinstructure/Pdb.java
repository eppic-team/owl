package proteinstructure;

import java.io.*;
import java.net.*;
import java.util.Collections;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.ArrayList;
import java.util.TreeSet;
import java.util.regex.*;
import java.util.zip.GZIPInputStream;
import java.util.Arrays;

import javax.vecmath.Point3d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

import tools.MySQLConnection;
import java.sql.SQLException;
import java.sql.Statement;

/**
 * A single chain pdb protein structure
 * 
 */
public abstract class Pdb {
	
	protected final static int DEFAULT_MODEL=1;				// default model serial (NMR structures)
	public final static String NONSTANDARD_AA_LETTER="X";   // letter we assign to nonstandard aas to use in sequence
	public static final String NULL_CHAIN_CODE = "NULL";	// to specify no chain code
	
	protected HashMap<String,Integer> resser_atom2atomserial; // residue serial+atom name (separated by underscore) to atom serials
	protected HashMap<Integer,String> resser2restype;   	// residue serial to 3 letter residue type 
	protected HashMap<Integer,Point3d> atomser2coord;  		// atom serials to 3D coordinates
	protected HashMap<Integer,Integer> atomser2resser;  	// atom serials to residue serials
	protected HashMap<Integer,String> atomser2atom;     	// atom serials to atom names
	protected HashMap<String,Integer> pdbresser2resser; 	// pdb (author) residue serials (can include insetion codes so they are strings) to internal residue serials
	protected HashMap<Integer,String> resser2pdbresser; 	// internal residue serials to pdb (author) residue serials (can include insertion codes so they are strings)
	
	protected SecondaryStructure secondaryStructure;				// the secondary structure annotation for this pdb object (should never be null)
	protected Scop scop;											// the scop annotation for this pdb object
	protected HashMap<Integer,Double> resser2allrsa;				// internal residue serials to all-atoms rsa
	protected HashMap<Integer,Double> resser2scrsa;					// internal residue serials to SC rsa
	protected HashMap<Integer,Double> resser2consurfhsspscore; 		// internal residue serials to SC rsa
	protected HashMap<Integer,Integer> resser2consurfhsspcolor; 	// internal residue serials to SC rsa
	protected EC ec;												// the ec annotation for this pdb object
	protected CatalSiteSet catalSiteSet;							// the catalytic site annotation for this pdb object
	
	protected String sequence; 		// full sequence as it appears in SEQRES field
	protected String pdbCode;
    // given "external" pdb chain code, i.e. the classic (author's) pdb code ("NULL" if it is blank in original pdb file)	
	protected String pdbChainCode;
    // Our internal chain identifier:
    // - in reading from pdbase or from msdsd it will be set to the internal chain id (asym_id field for pdbase, pchain_id for msdsd)
    // - in reading from pdb file it coincides with pdbChainCode except for "NULL" where we use "A"
	protected String chainCode;
	protected int model;  			// the model serial for NMR structures
	
	// in case we read from pdb file, 2 possibilities:
	// 1) SEQRES field is present: fullLength coincides with sequence length
	// 2) SEQRES not present: fullLength is maximum observed residue (so if residue numbering is correct that only misses possible non-observed residues at end of chain)
	protected int fullLength; // length of full sequence as it appears in SEQRES field 
	protected int obsLength;  // length without unobserved, non standard aas 
	
	protected String db;			// the db from which we have taken the data (if applies)
	protected MySQLConnection conn;
	
	public int checkCSA(String version, boolean online) throws IOException {
		BufferedReader in;
		String inputLine;
		
		if (online) {
			URL csa = new URL("http://www.ebi.ac.uk/thornton-srv/databases/CSA/archive/CSA_"+version.replaceAll("\\.","_")+".dat.gz");
			URLConnection conn = csa.openConnection();
			InputStream inURL = conn.getInputStream();
			OutputStream out = new BufferedOutputStream(new FileOutputStream("CSA_"+version.replaceAll("\\.","_")+".dat.gz"));
			byte[] buffer = new byte[1024];
			int numRead;
			long numWritten = 0;
			while ((numRead = inURL.read(buffer)) != -1) {
				out.write(buffer, 0, numRead);
				numWritten += numRead;
			}
			inURL.close();
			out.close();
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream("CSA_"+version.replaceAll(".","_")+".dat.gz"))));
		} else {
			File csaFile = new File("/project/StruPPi/Databases/CSA/CSA_"+version.replaceAll("\\.","_")+".dat");
			in = new BufferedReader(new FileReader(csaFile));			
		}
		
		int csaMistakes = 0;
		this.catalSiteSet = new CatalSiteSet();
		String curPdbCode, curPdbChainCode, curPdbResSerial, curResType;
		String prevPdbCode = "";
		int curSiteId;
		int prevSiteId = -1;
		CatalyticSite cs = null;
		while ((inputLine = in.readLine()) != null) { 
			String[] fields = inputLine.split(",");
			curPdbCode = fields[0]; 
			curPdbChainCode = (fields[3].equals(""))?"NULL":fields[3];
			if (curPdbCode.equals(pdbCode)) {
				if (curPdbChainCode.equals(pdbChainCode)) {
					curPdbResSerial = fields[4];
					curResType = fields[2];
					curSiteId = Integer.valueOf(fields[1]);
					// only if the pdb residue type is a standard AA
					if (AAinfo.isValidAA(curResType)) {
						// only if the pdb residue type agrees with our residue type
						if (getResTypeFromResSerial(get_resser_from_pdbresser(curPdbResSerial)).equals(curResType)) {
							// each time the site changes except for the first site of a chain,
							// add the site to the set
							if ((curSiteId != prevSiteId) & (prevSiteId != -1)) {
								catalSiteSet.add(cs);
							}
							// each time the site changes or it is the first site of a chain,
							// create a new site
							if ((curSiteId != prevSiteId) | (prevSiteId == -1)) {
								String littEntryPdbCode = fields[7].substring(0,4);
								String littEntryPdbChainCode = fields[7].substring(4).equals("")?"NULL":fields[7].substring(4);
								cs = new CatalyticSite(curSiteId, fields[6], littEntryPdbCode, littEntryPdbChainCode); 
							}
							// add the res to the site
							cs.addRes(get_resser_from_pdbresser(curPdbResSerial), fields[5]);
							prevSiteId = curSiteId;
						} else {
							csaMistakes++;
						}
					}
				}
			} else if (prevPdbCode.equals(pdbCode)) {
				if (prevSiteId != -1) {
					catalSiteSet.add(cs);
				}
				break;
			}
			prevPdbCode = curPdbCode;
		}
		
		in.close();
		if (online) {
			File csaFile = new File("CSA_"+version.replaceAll(".","_")+".dat.gz");
			csaFile.delete();
		}
		
		if ((csaMistakes > 0) | (prevSiteId == -1)) {
			this.catalSiteSet = null;
		}
		
		return csaMistakes;
	}
		
	public void checkEC(boolean online) throws IOException {
		
		this.ec = new EC();	
		ECRegion er = null;
		String startPdbResSer = "", endPdbResSer = "";
		BufferedReader in;
		String inputLine;
		Pattern p = Pattern.compile("^ \\d");
		
		if (online) {
			URL pdb2ecMapping = new URL("http://www.bioinf.org.uk/pdbsprotec/mapping.txt");
			URLConnection p2e= pdb2ecMapping.openConnection();
			in = new BufferedReader(new InputStreamReader(p2e.getInputStream()));
		} else {
			File pdb2ecMapping = new File("/project/StruPPi/Databases/PDBSProtEC/mapping.txt");
			in = new BufferedReader(new FileReader(pdb2ecMapping));
		}
		
        while ((inputLine = in.readLine()) != null) { 
        	Matcher m = p.matcher(inputLine);
        	if (m.find()) {
	        	String curPdbCode = inputLine.substring(0,9).trim();
	        	String curPdbChainCode = (inputLine.charAt(11) == ' ')?"NULL":String.valueOf(inputLine.charAt(11));
	        	if (curPdbCode.equals(pdbCode) && curPdbChainCode.equals(pdbChainCode)) {
	        		startPdbResSer = inputLine.substring(20,26).trim();
	        		endPdbResSer = inputLine.substring(27,33).trim();
	        		String id = inputLine.substring(43).trim();
	        		//System.out.println(curPdbCode+":"+curPdbChainCode+":"+startPdbResSer+"-"+endPdbResSer+":"+ec);
	        		er = new ECRegion(id, startPdbResSer, endPdbResSer, get_resser_from_pdbresser(startPdbResSer), get_resser_from_pdbresser(endPdbResSer));
					ec.add(er);
	        	}
        	}
        }
        
        in.close();
        
	}
	
	public int checkConsurfHssp(boolean online) throws IOException {
		
		BufferedReader in;
		if (online) {
			URL consurfhssp = new URL("http://consurf.tau.ac.il/results/HSSP_ML_"+pdbCode+(pdbChainCode.equals("NULL")?"_":pdbChainCode)+"/pdb"+pdbCode+".gradesPE");
			URLConnection ch = consurfhssp.openConnection();
			in = new BufferedReader(new InputStreamReader(ch.getInputStream()));
		} else {
			File consurfhssp = new File("/project/StruPPi/Databases/ConSurf-HSSP/ConservationGrades/"+pdbCode+(pdbChainCode.equals("NULL")?"_":pdbChainCode)+".grades");
			in = new BufferedReader(new FileReader(consurfhssp));
		}
		String inputLine;
		Pattern p = Pattern.compile("^\\s+\\d+");
		int lineCount = 0;
		int consurfHsspMistakes = 0;
		
		Integer[] ressers = new Integer[resser2restype.size()];
		resser2restype.keySet().toArray(ressers);
		Arrays.sort(ressers);
		
		resser2consurfhsspscore = new HashMap<Integer,Double>();
		resser2consurfhsspcolor = new HashMap<Integer,Integer>();
		
        while ((inputLine = in.readLine()) != null) { 
        	Matcher m = p.matcher(inputLine);
        	if (m.find()) {
        		lineCount++;
        		int resser = ressers[lineCount-1];
        		String[] fields = inputLine.split("\\s+");
    			String pdbresser = fields[3].equals("-")?"-":fields[3].substring(3, fields[3].indexOf(':'));
    			if (fields[2].equals(AAinfo.threeletter2oneletter(getResTypeFromResSerial(resser))) &
    					(pdbresser.equals("-") | pdbresser.equals(get_pdbresser_from_resser(resser)))) {
    				resser2consurfhsspscore.put(resser, Double.valueOf(fields[4]));
    				resser2consurfhsspcolor.put(resser, Integer.valueOf(fields[5]));
    			} else {
    				consurfHsspMistakes++;
    			}
        	}
		}
        in.close();
        
        consurfHsspMistakes += Math.abs(ressers.length - resser2consurfhsspscore.size());
        if (consurfHsspMistakes > 0) {
        	resser2consurfhsspscore.clear();
        	resser2consurfhsspcolor.clear();
        }
        
        return consurfHsspMistakes;
        
    }	
	
	public void runNaccess(String naccessExecutable, String naccessParameters) throws Exception {
		String pdbFileName = pdbCode+chainCode+".pdb";
		dump2pdbfile(pdbFileName, true);
		String line;
		int errorLineCount = 0;
		
		File test = new File(naccessExecutable);
		if(!test.canRead()) throw new IOException("Naccess Executable is not readable");
		
		Process myNaccess = Runtime.getRuntime().exec(naccessExecutable + " " + pdbFileName + " " + naccessParameters);
		BufferedReader naccessOutput = new BufferedReader(new InputStreamReader(myNaccess.getInputStream()));
		BufferedReader naccessError = new BufferedReader(new InputStreamReader(myNaccess.getErrorStream()));
		while((line = naccessOutput.readLine()) != null) {
		}
		while((line = naccessError.readLine()) != null) {
			errorLineCount++;
		}
		naccessOutput.close();
		naccessError.close();
		int exitVal = myNaccess.waitFor();
		if ((exitVal == 1) || (errorLineCount > 0)) {
			throw new Exception("Naccess:Wrong arguments or pdb file format!");
		}
		
		File rsa = new File(pdbCode+chainCode+".rsa");
		if (rsa.exists()) {
			resser2allrsa = new HashMap<Integer,Double>();
			resser2scrsa = new HashMap<Integer,Double>();			
			BufferedReader rsaInput = new BufferedReader(new FileReader(rsa));
			while ((line = rsaInput.readLine()) != null) {
				if (line.startsWith("RES")) {
					int resser = Integer.valueOf(line.substring(9,13).trim());
					double allrsa = Double.valueOf(line.substring(22,28).trim());
					double scrsa = Double.valueOf(line.substring(35,41).trim());
					resser2allrsa.put(resser, allrsa);
					resser2scrsa.put(resser, scrsa);
				}
			}
			rsaInput.close();
		}
		String[] filesToDelete = { ".pdb", ".asa", ".rsa", ".log" };
		for (int i=0; i < filesToDelete.length; i++) {
			File fileToDelete = new File(pdbCode+chainCode+filesToDelete[i]);
			if (fileToDelete.exists()) {
				fileToDelete.delete();
			}
		}
			
	}
	
	public void checkScop(String version, boolean online) throws IOException {
		this.scop = new Scop();	
		ScopRegion sr = null;
		String startPdbResSer = "", endPdbResSer = "";
		BufferedReader in;
		String inputLine;
		
		if (online) {
			URL scop_cla = new URL("http://scop.mrc-lmb.cam.ac.uk/scop/parse/dir.cla.scop.txt_"+version);
			URLConnection sc = scop_cla.openConnection();
			in = new BufferedReader(new InputStreamReader(sc.getInputStream()));
		} else {
			File scop_cla = new File("/project/StruPPi/Databases/SCOP/dir.cla.scop.txt_"+version);
			in = new BufferedReader(new FileReader(scop_cla));
		}
				
        while ((inputLine = in.readLine()) != null) 
			if (inputLine.startsWith(pdbCode,1)) {
				String[] fields = inputLine.split("\\s+");
				String[] regions = fields[2].split(",");
				for (int j=0; j < regions.length; j++) {
					Pattern p = Pattern.compile("^(-)|([a-zA-Z\\d]):(-?\\d+[a-zA-Z]*)-(-?\\d+[a-zA-Z]*)|(-?\\d+[a-zA-Z]*)-(-?\\d+[a-zA-Z]*)|([a-zA-Z\\d]):");
					Matcher m = p.matcher(regions[j]);
					if (m.find()) {
						if (((pdbChainCode.equals("NULL") && ((m.group(1) != null && m.group(1).equals("-")) || m.group(5) != null))) || 
								(m.group(2) != null && m.group(2).equals(pdbChainCode)) || 
								(m.group(7) != null && m.group(7).equals(pdbChainCode))) {
							if (m.group(3) != null) {
								startPdbResSer = m.group(3);
								endPdbResSer = m.group(4);
							} else if (m.group(5) != null) {
								startPdbResSer = m.group(5);
								endPdbResSer = m.group(6);								
							} else {
								startPdbResSer = get_pdbresser_from_resser(Collections.min(resser2pdbresser.keySet()));
								endPdbResSer = get_pdbresser_from_resser(Collections.max(resser2pdbresser.keySet()));
							}
							sr = new ScopRegion(fields[0], fields[3], Integer.parseInt(fields[4]), j, regions.length, startPdbResSer, endPdbResSer, get_resser_from_pdbresser(startPdbResSer), get_resser_from_pdbresser(endPdbResSer));
							scop.add(sr);
						}
					}
				}
				//break;
			}
        in.close();
        
        scop.setVersion(version);
    }	
	
	public void runDssp(String dsspExecutable, String dsspParameters) throws IOException {
		runDssp(dsspExecutable, dsspParameters, SecStrucElement.ReducedState.FOURSTATE, SecStrucElement.ReducedState.FOURSTATE);
	}
	
	/** 
	 * Runs an external DSSP executable and (re)assigns the secondary structure annotation from the parsed output.
	 * Existing secondary structure information will be overwritten.
	 * As of September 2007, a DSSP executable can be downloaded from http://swift.cmbi.ru.nl/gv/dssp/
	 * after filling out a license agreement. For this version, the dsspParamters variable should be
	 * set to "--" (two hyphens).
	 */
	public void runDssp(String dsspExecutable, String dsspParameters, SecStrucElement.ReducedState state4Type, SecStrucElement.ReducedState state4Id) throws IOException {
		String startLine = "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC";
		String line;
		int lineCount = 0;
		char ssType, sheetLabel;
		TreeMap<Integer, Character> ssTypes;
		TreeMap<Integer, Character> sheetLabels;
		int resNum;
		String resNumStr;
		File test = new File(dsspExecutable);
		if(!test.canRead()) throw new IOException("DSSP Executable is not readable");
		Process myDssp = Runtime.getRuntime().exec(dsspExecutable + " " + dsspParameters);
		PrintWriter dsspInput = new PrintWriter(myDssp.getOutputStream());
		BufferedReader dsspOutput = new BufferedReader(new InputStreamReader(myDssp.getInputStream()));
		BufferedReader dsspError = new BufferedReader(new InputStreamReader(myDssp.getErrorStream()));
		writeAtomLines(dsspInput, true);	// pipe atom lines to dssp
		dsspInput.close();
		ssTypes = new TreeMap<Integer,Character>();
		sheetLabels = new TreeMap<Integer,Character>();
		while((line = dsspOutput.readLine()) != null) {
			lineCount++;
			if(line.startsWith(startLine)) {
				//System.out.println("Dssp Output: ");
				break;
			}
		}
		while((line = dsspOutput.readLine()) != null) {
			lineCount++;
			resNumStr = line.substring(5,10).trim();
			ssType = line.charAt(16);			
			sheetLabel = line.charAt(33);
			if (state4Id == SecStrucElement.ReducedState.FOURSTATE && SecStrucElement.getReducedStateTypeFromDsspType(ssType, state4Id) == SecStrucElement.OTHER) {
				sheetLabel = ' ';
			}
			if(!resNumStr.equals("")) {		// this should only happen if dssp inserts a line indicating a chain break
				try {
					resNum = Integer.valueOf(resNumStr);
					ssTypes.put(resNum, ssType);
					sheetLabels.put(resNum, sheetLabel);
				} catch (NumberFormatException e) {
					System.err.println("Error while parsing DSSP output. Expected residue number, found '" + resNumStr + "' in line " + lineCount);
				}
			}
		}
		//for(char c:ssTypes) {System.out.print(c);}; System.out.println(".");
		dsspOutput.close();
		dsspError.close();

		if(ssTypes.size() == 0) {
			throw new IOException("No DSSP output found.");
		}
		
		if(ssTypes.size() != get_length()) {	// compare with number of observed residues
			System.err.println("Error: DSSP output size (" + ssTypes.size() + ") does not match number of observed residues in structure (" + get_length() + ").");
		}

		// assign secondary structure
		this.secondaryStructure = new SecondaryStructure();		// forget the old annotation
		char lastType = SecStrucElement.getReducedStateTypeFromDsspType(ssTypes.get(ssTypes.firstKey()), state4Id);
		int lastResSer = ssTypes.firstKey();
		char lastSheet = sheetLabels.get(lastResSer);
		char thisType, thisSheet, reducedType;
		int start = 0;
		int elementCount = 0;
		SecStrucElement ssElem;
		String ssId;
		for(int resSer:ssTypes.keySet()) {
			thisType = SecStrucElement.getReducedStateTypeFromDsspType(ssTypes.get(resSer), state4Id);
			thisSheet = sheetLabels.get(resSer);
			if(thisType != lastType || thisSheet != lastSheet || resSer > lastResSer+1) {
				// finish previous element, start new one
				elementCount++;
				reducedType = SecStrucElement.getReducedStateTypeFromDsspType(ssTypes.get(lastResSer), state4Type);
				ssId = new Character(lastType).toString() + (lastSheet==' '?"":String.valueOf(lastSheet)) + new Integer(elementCount).toString();
				ssElem = new SecStrucElement(reducedType,start,lastResSer,ssId);
				secondaryStructure.add(ssElem);
				start = resSer;
				lastType = thisType;
				lastSheet = thisSheet;
			}
			lastResSer = resSer;
		}
		// finish last element
		elementCount++;
		reducedType = SecStrucElement.getReducedStateTypeFromDsspType(ssTypes.get(ssTypes.lastKey()), state4Type);
		ssId = new Character(lastType).toString() + (lastSheet==' '?"":String.valueOf(lastSheet)) + new Integer(elementCount).toString();
		ssElem = new SecStrucElement(reducedType, start,ssTypes.lastKey(),ssId);
		secondaryStructure.add(ssElem);
		
		secondaryStructure.setComment("DSSP");
	}

	/** Writes atom lines for this structure to the given output stream */
	private void writeAtomLines(PrintWriter Out, boolean pdbCompatible) {
		TreeMap<Integer,Object[]> lines = new TreeMap<Integer,Object[]>();
		Out.println("HEADER  Dumped from "+db+". pdb code="+pdbCode+", chain='"+chainCode+"'");
		String chainCodeStr = chainCode;
		if (pdbCompatible) {
			chainCodeStr = chainCode.substring(0,1);
		}
		for (String resser_atom:resser_atom2atomserial.keySet()){
			int atomserial = resser_atom2atomserial.get(resser_atom);
			int res_serial = Integer.parseInt(resser_atom.split("_")[0]);
			String atom = resser_atom.split("_")[1];
			String res_type = resser2restype.get(res_serial);
			Point3d coords = atomser2coord.get(atomserial);
			String atomType = atom.substring(0,1);
			Object[] fields = {atomserial, atom, res_type, chainCodeStr, res_serial, coords.x, coords.y, coords.z, atomType};
			lines.put(atomserial, fields);
		}
		for (int atomserial:lines.keySet()){
			// Local.US is necessary, otherwise java prints the doubles locale-dependant (i.e. with ',' for some locales)
			Out.printf(Locale.US,"ATOM  %5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00           %s\n",lines.get(atomserial));
		}
		Out.println("END");
		Out.close();		
	}

	/**
	 * Dumps coordinate data into a file in pdb format (ATOM lines only)
	 * The residue serials dumped are the internal ones.
	 * The chain dumped is the value of the chainCode variable, i.e. our internal chain identifier for Pdb objects
	 * @param outfile
	 * @throws IOException
	 */
	public void dump2pdbfile(String outfile) throws IOException {
		dump2pdbfile(outfile, false);
	}
	
	public void dump2pdbfile(String outfile, boolean pdbCompatible) throws IOException {
		PrintWriter Out = new PrintWriter(new FileOutputStream(outfile));
		writeAtomLines(Out, pdbCompatible);
	}
	
	/**
	 * Dump the full sequence of this Pdb object in fasta file format
	 * @param seqfile
	 * @throws IOException
	 */
	public void dumpseq(String seqfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(seqfile));
		Out.println(">"+pdbCode+"_"+pdbChainCode);
		Out.println(sequence);
		Out.close();
	}
	
	/** 
	 * Returns the number of observed standard residues.
	 * TODO: Refactor method name
	 */
	public int get_length(){
		return obsLength;
	}
	
	/** Returns the number of residues in the sequence of this protein. */
	public int getFullLength() {
		return fullLength;
	}
	
	/**
	 * Returns number of (non-Hydrogen) atoms in the protein
	 * @return
	 */
	public int getNumAtoms() {
		return atomser2atom.size();
	}
	
	/**
	 * Gets a TreeMap with atom serials as keys and their coordinates as values for the given contact type
	 * The contact type can't be a cross contact type, it doesn't make sense here
	 * @param ct
	 * @return
	 */
	private TreeMap<Integer,Point3d> get_coords_for_ct(String ct) {
		TreeMap<Integer,Point3d> coords = new TreeMap<Integer,Point3d>(); 
		for (int resser:resser2restype.keySet()){
			Set<String> atoms = AAinfo.getAtomsForCTAndRes(ct, resser2restype.get(resser));
			for (String atom:atoms){
				if (resser_atom2atomserial.containsKey(resser+"_"+atom)){
					int atomser = resser_atom2atomserial.get(resser+"_"+atom);
					Point3d coord = atomser2coord.get(atomser);
					coords.put(atomser, coord);
				}
				else {
					//System.err.println("Couldn't find "+atom+" atom for resser="+resser+". Continuing without that atom for this resser.");
				}
			}
			// in cts ("ALL","BB") we still miss the OXT, we need to add it now if it is there (it will be there when this resser is the last residue)
			if ((ct.equals("ALL") || ct.equals("BB")) && resser_atom2atomserial.containsKey(resser+"_"+"OXT")){
				int atomser = resser_atom2atomserial.get(resser+"_"+"OXT");
				Point3d coord = atomser2coord.get(atomser);
				coords.put(atomser, coord);
			}
		}
		return coords;
	}

	/**
	 * Gets a TreeMap with residue serials+"_"+atomname as keys and their coordinates as values for the given contact type
	 * The contact type can't be a cross contact type, it doesn't make sense here
	 * Used in rmsd method
	 * @param ct
	 * @return
	 */
	private TreeMap<String,Point3d> get_coords_for_ct_4rmsd(String ct) {
		TreeMap<String,Point3d> coords = new TreeMap<String,Point3d>(); 
		for (int resser:resser2restype.keySet()){
			Set<String> atoms = AAinfo.getAtomsForCTAndRes(ct,resser2restype.get(resser));
			for (String atom:atoms){
				if (resser_atom2atomserial.containsKey(resser+"_"+atom)){
					int atomser = resser_atom2atomserial.get(resser+"_"+atom);
					Point3d coord = atomser2coord.get(atomser);
					coords.put(resser+"_"+atom, coord);
				}
			}
			// in cts ("ALL","BB") we still miss the OXT, we need to add it now if it is there (it will be there when this resser is the last residue)
			if ((ct.equals("ALL") || ct.equals("BB")) && resser_atom2atomserial.containsKey(resser+"_"+"OXT")){
				int atomser = resser_atom2atomserial.get(resser+"_"+"OXT");
				Point3d coord = atomser2coord.get(atomser);
				coords.put(resser+"_"+"OXT", coord);
			}
		}
		return coords;
	}
	
	/**
	 * Returns the distance matrix as a HashMap with Contacts (residue serial pairs) as keys
	 * For multi atom contact types the distance matrix has the minimum distance for each pair of
	 * residues 
	 * AAinfo.isValidSingleAtomCT(ct) can be used to check before calling.
	 * @param ct contact type for which distances are being calculated
	 * @return A map which assings to each edge the corresponding distance
	 * TODO: 
	 */
	public HashMap<Edge, Double> calculate_dist_matrix(String ct){
		HashMap<Edge,Double> distMatrixAtoms = calculate_atom_dist_matrix(ct);

		/**
		 * Helper method which maps atom serials to residue serials.
		 * TODO: Integrate into above method to avoid storing two distance maps in memory
		 */
		HashMap<Edge,Double> distMatrixRes = new HashMap<Edge, Double>();
		for (Edge cont: distMatrixAtoms.keySet()){
			int i_resser = get_resser_from_atomser(cont.i);
			int j_resser = get_resser_from_atomser(cont.j);
			Edge pair = new Edge(i_resser,j_resser);
			if (distMatrixRes.containsKey(pair)) {
				distMatrixRes.put(pair, Math.min(distMatrixRes.get(pair), distMatrixAtoms.get(cont)));
			} else {
				distMatrixRes.put(pair, distMatrixAtoms.get(cont));
			}
		}

		return distMatrixRes;
	}
	
	/**
	 * Returns the distance matrix as a HashMap with Contacts (atom serial pairs) as keys
	 * This method can be used for any contact type
	 * AAinfo.isValidSingleAtomCT(ct) can be used to check before calling.
	 * @param ct contact type for which distances are being calculated
	 * @return A map which assings to each edge the corresponding distance
	 */
	public HashMap<Edge, Double> calculate_atom_dist_matrix(String ct){
		HashMap<Edge,Double> distMatrixAtoms = new HashMap<Edge,Double>();
		if (!ct.contains("/")){
			TreeMap<Integer,Point3d> coords = get_coords_for_ct(ct);
			for (int i_atomser:coords.keySet()){
				for (int j_atomser:coords.keySet()){
					if (j_atomser>i_atomser) {
						Edge pair = new Edge(i_atomser,j_atomser);
						distMatrixAtoms.put(pair, coords.get(i_atomser).distance(coords.get(j_atomser)));
					}
				}
			}
		} else {
			String i_ct = ct.split("/")[0];
			String j_ct = ct.split("/")[1];
			TreeMap<Integer,Point3d> i_coords = get_coords_for_ct(i_ct);
			TreeMap<Integer,Point3d> j_coords = get_coords_for_ct(j_ct);
			for (int i_atomser:i_coords.keySet()){
				for (int j_atomser:j_coords.keySet()){
					if (j_atomser!=i_atomser){
						Edge pair = new Edge(i_atomser,j_atomser);
						distMatrixAtoms.put(pair, i_coords.get(i_atomser).distance(j_coords.get(j_atomser)));
					}
				}
			}
		}

		return distMatrixAtoms;
	}
	
	/**
	 * Get the graph for given contact type and cutoff for this Pdb object.
	 * Returns a Graph object with the contacts
	 * We do geometric hashing for fast contact computation (without needing to calculate full distance matrix)
	 * @param ct
	 * @param cutoff
	 * @return
	 */
	public Graph get_graph(String ct, double cutoff){ 
		TreeMap<Integer,Point3d> i_coords = null;
		TreeMap<Integer,Point3d> j_coords = null;		// only relevant for asymetric edge types
		boolean directed = false;
		if (!ct.contains("/")){
			i_coords = get_coords_for_ct(ct);
			directed = false;
		} else {
			String i_ct = ct.split("/")[0];
			String j_ct = ct.split("/")[1];
			i_coords = get_coords_for_ct(i_ct);
			j_coords = get_coords_for_ct(j_ct);
			directed = true;
		}
		int[] i_atomserials = new  int[i_coords.size()]; // map from matrix indices to atomserials
		int[] j_atomserials = null;
		
		int SCALE=100; // i.e. we use units of hundredths of Amstrongs (thus cutoffs can be specified with a maximum precission of 0.01A)
		
		int boxSize = (int) Math.floor(cutoff*SCALE);
		
		HashMap<Point3i,Box> boxes = new HashMap<Point3i,Box>();
		int i=0;
		for (int i_atomser:i_coords.keySet()){
			//coordinates for atom serial atomser, we will use i as its identifier below
			Point3d coord = i_coords.get(i_atomser);
			int floorX = boxSize*((int)Math.floor(coord.x*SCALE/boxSize));
			int floorY = boxSize*((int)Math.floor(coord.y*SCALE/boxSize));
			int floorZ = boxSize*((int)Math.floor(coord.z*SCALE/boxSize));
			Point3i floor = new Point3i(floorX,floorY,floorZ);
			if (boxes.containsKey(floor)){
				// we put the coords for atom i in its corresponding box (identified by floor)
				boxes.get(floor).put_i_Point(i, coord);
				if (!directed){
					boxes.get(floor).put_j_Point(i, coord);
				}
			} else {
				Box box = new Box(floor);
				box.put_i_Point(i, coord);
				if (!directed){
					box.put_j_Point(i, coord);
				}
				boxes.put(floor,box);
			}
			i_atomserials[i]=i_atomser; //as atomserials in coords were ordered (TreeMap) the new indexing will still be ordered
			i++;
		}
		int j=0;
		if (directed) {
			j_atomserials = new  int[j_coords.size()];
			for (int j_atomser:j_coords.keySet()){
				//coordinates for atom serial atomser, we will use j as its identifier below
				Point3d coord = j_coords.get(j_atomser);
				int floorX = boxSize*((int)Math.floor(coord.x*SCALE/boxSize));
				int floorY = boxSize*((int)Math.floor(coord.y*SCALE/boxSize));
				int floorZ = boxSize*((int)Math.floor(coord.z*SCALE/boxSize));
				Point3i floor = new Point3i(floorX,floorY,floorZ);
				if (boxes.containsKey(floor)){
					// we put the coords for atom j in its corresponding box (identified by floor)
					boxes.get(floor).put_j_Point(j, coord);
				} else {
					Box box = new Box(floor);
					box.put_j_Point(j, coord);
					boxes.put(floor,box);
				}
				j_atomserials[j]=j_atomser; //as atomserials in coords were ordered (TreeMap) the new indexing will still be ordered
				j++;
			}
		} else {
			j_atomserials = i_atomserials;
		}

		
		float[][]distMatrix = new float[i_atomserials.length][j_atomserials.length];
		
		for (Point3i floor:boxes.keySet()){ // for each box
			// distances of points within this box
			boxes.get(floor).getDistancesWithinBox(distMatrix,directed);

			//TODO should iterate only through half of the neighbours here 
			// distances of points from this box to all neighbouring boxes: 26 iterations (26 neighbouring boxes)
			for (int x=floor.x-boxSize;x<=floor.x+boxSize;x+=boxSize){
				for (int y=floor.y-boxSize;y<=floor.y+boxSize;y+=boxSize){
					for (int z=floor.z-boxSize;z<=floor.z+boxSize;z+=boxSize){
						if (!((x==floor.x)&&(y==floor.y)&&(z==floor.z))) { // skip this box
							Point3i neighbor = new Point3i(x,y,z);
							if (boxes.containsKey(neighbor)){
								boxes.get(floor).getDistancesToNeighborBox(boxes.get(neighbor),distMatrix,directed);
							}
						}
					}
				}
			} 
		} 
		
		// getting the contacts (in residue serials) from the atom serials (partial) distance matrix 
		TreeMap<Edge,Integer> weights = new TreeMap<Edge,Integer>(); // we use them to keep the edges together with the weights temporarily
		for (i=0;i<distMatrix.length;i++){
			for (j=0;j<distMatrix[i].length;j++){
				// the condition distMatrix[i][j]!=0.0 takes care of skipping several things: 
				// - diagonal of the matrix in case of undirected
				// - lower half of matrix in case of undirected
				// - cells for which we didn't calculate a distance because the 2 points were not in same or neighbouring boxes (i.e. too far apart)
				if (distMatrix[i][j]!=0.0f && distMatrix[i][j]<=cutoff){
					int i_resser = atomser2resser.get(i_atomserials[i]);
					int j_resser = atomser2resser.get(j_atomserials[j]);
					Edge resser_pair = new Edge(i_resser,j_resser);
					// for multi-atom models (BB, SC, ALL or BB/SC) we need to make sure that we don't have contacts from residue to itself or that we don't have duplicates				
					if (i_resser!=j_resser){ // duplicates are automatically taken care by the TreeMap which doesn't allow duplicates in its keys 
						// as weights we count the number of atom-edges per residue-edge
						if (weights.containsKey(resser_pair)) {
							weights.put(resser_pair, weights.get(resser_pair)+1);
						} else {
							weights.put(resser_pair, 1);
						}
					}
				}

			}
		}
		// finally we put the edges together with their weights into the contacts EdgeSet (right now we just have them in the weights TreeMap):
		EdgeSet contacts = new EdgeSet();
		for (Edge cont: weights.keySet()) {
			contacts.add(new Edge(cont.i,cont.j,weights.get(cont)));
		}

		// creating and returning the graph object
		TreeMap<Integer,String> nodes = new TreeMap<Integer,String>();
		for (int resser:resser2restype.keySet()){
			nodes.put(resser,resser2restype.get(resser));
		}
		Graph graph = new Graph (contacts,nodes,sequence,cutoff,ct,pdbCode,chainCode,pdbChainCode,model, secondaryStructure.copy());
		return graph;
	}
	
	public void calcGridDensity(String ct, double cutoff, Map<Integer, Integer> densityCount) { 
		TreeMap<Integer,Point3d> i_coords = null;
		TreeMap<Integer,Point3d> j_coords = null;		// only relevant for asymmetric edge types
		boolean directed = false;
		if (!ct.contains("/")){
			i_coords = get_coords_for_ct(ct);			// mapping from atom serials to coordinates
			directed = false;
		} else {
			String i_ct = ct.split("/")[0];
			String j_ct = ct.split("/")[1];
			i_coords = get_coords_for_ct(i_ct);
			j_coords = get_coords_for_ct(j_ct);
			directed = true;
		}
		int[] i_atomserials = new  int[i_coords.size()]; // map from matrix indices to atomserials
		int[] j_atomserials = null;
		
		int SCALE=100; // i.e. we use units of hundredths of Angstroms (thus cutoffs can be specified with a maximum precission of 0.01A)
		
		int boxSize = (int) Math.floor(cutoff*SCALE);
		
		HashMap<Point3i,Box> boxes = new HashMap<Point3i,Box>();
		int i=0;
		for (int i_atomser:i_coords.keySet()){
			//coordinates for atom serial atomser, we will use i as its identifier below
			Point3d coord = i_coords.get(i_atomser);
			int floorX = boxSize*((int)Math.floor(coord.x*SCALE/boxSize));
			int floorY = boxSize*((int)Math.floor(coord.y*SCALE/boxSize));
			int floorZ = boxSize*((int)Math.floor(coord.z*SCALE/boxSize));
			Point3i floor = new Point3i(floorX,floorY,floorZ);
			if (boxes.containsKey(floor)){
				// we put the coords for atom i in its corresponding box (identified by floor)
				boxes.get(floor).put_i_Point(i, coord);
				if (!directed){
					boxes.get(floor).put_j_Point(i, coord);
				}
			} else {
				Box box = new Box(floor);
				box.put_i_Point(i, coord);
				if (!directed){
					box.put_j_Point(i, coord);
				}
				boxes.put(floor,box);
			}
			i_atomserials[i]=i_atomser; //as atomserials in coords were ordered (TreeMap) the new indexing will still be ordered
			i++;
		}
		int j=0;
		if (directed) {
			j_atomserials = new  int[j_coords.size()];
			for (int j_atomser:j_coords.keySet()){
				//coordinates for atom serial atomser, we will use j as its identifier below
				Point3d coord = j_coords.get(j_atomser);
				int floorX = boxSize*((int)Math.floor(coord.x*SCALE/boxSize));
				int floorY = boxSize*((int)Math.floor(coord.y*SCALE/boxSize));
				int floorZ = boxSize*((int)Math.floor(coord.z*SCALE/boxSize));
				Point3i floor = new Point3i(floorX,floorY,floorZ);
				if (boxes.containsKey(floor)){
					// we put the coords for atom j in its corresponding box (identified by floor)
					boxes.get(floor).put_j_Point(j, coord);
				} else {
					Box box = new Box(floor);
					box.put_j_Point(j, coord);
					boxes.put(floor,box);
				}
				j_atomserials[j]=j_atomser; //as atomserials in coords were ordered (TreeMap) the new indexing will still be ordered
				j++;
			}
		} else {
			j_atomserials = i_atomserials;
		}
		
		// count density
		for(Point3i floor:boxes.keySet()) {
			//int size = boxes.get(floor).size();
			int size = getNumGridNbs(boxes, floor, boxSize);	// count number of neighbouring grid cells with points in them
			if(densityCount.containsKey(size)) {
				int old = densityCount.get(size);
				densityCount.put(size, ++old);
			} else {
				densityCount.put(size, 1);
			}
		}
		
		
	}
	
	/**
	 * Calculates and returns the difference of the distance maps of this 
	 * structure and another pdb object. On error returns null. Note that 
	 * the distance maps for the single structures have to have the same 
	 * size.
	 * @param contactType1  contact type of this structure
	 * @param pdb2  the second structure
	 * @param contactType2  contact type of the second structure
	 * @return the difference distance map
	 * @see #getDiffDistMap(String, Pdb, String, Alignment, String, String) 
	 */
	public HashMap<Edge,Double> getDiffDistMap(String contactType1, Pdb pdb2, String contactType2) {
		double dist1, dist2, diff;
		HashMap<Edge,Double> otherDistMatrix = pdb2.calculate_dist_matrix(contactType2);
		HashMap<Edge,Double> thisDistMatrix = this.calculate_dist_matrix(contactType1);
		if(thisDistMatrix.size() != otherDistMatrix.size()) {
			System.err.println("Cannot calculate difference distance map. Matrix sizes do not match.");
			return null;
		}
		for(Edge e:otherDistMatrix.keySet()){
			dist1 = otherDistMatrix.get(e);
			if(!thisDistMatrix.containsKey(e)) {
				System.err.println("Error while calculating difference distance map. Entry " + e + " in matrix1 not not found in matrix2.");
				return null;
			}
			dist2 = thisDistMatrix.get(e);
			diff = Math.abs(dist1-dist2);
			otherDistMatrix.put(e, diff);
		}
		return otherDistMatrix;
	}
	
	/**
	 * Calculates the difference distance map of this structure and 
	 * another pdb object given a sequence alignment of the structures. The 
	 * resulting difference distance map may contains non-defined distances. 
	 * This behavior is due to the alignment. If any residue in either 
	 * structures is aligned with a gap one cannot assign a "difference 
	 * distance" to another residue pair.   
	 * @param contactType1  contact type of this structure
	 * @param pdb2  the second structure
	 * @param contactType2  contact type of the second structure
	 * @param ali  sequence alignment of both structures
	 * @param name1  sequence tag of the this structure in the alignment
	 * @param name2  sequence tag og the second structure in the alignment
	 * @return the difference distance map
	 * @see #getDiffDistMap(String, Pdb, String)
	 */
	public HashMap<Edge,Double> getDiffDistMap(String contactType1, Pdb pdb2, String contactType2, Alignment ali, String name1, String name2) {
	    
	    HashMap<Edge,Double> otherDistMatrix = pdb2.calculate_dist_matrix(contactType2);
	    HashMap<Edge,Double> thisDistMatrix = this.calculate_dist_matrix(contactType1);
	    HashMap<Edge,Double> alignedDistMatrix = new HashMap<Edge, Double>(Math.min(this.getFullLength(), pdb2.getFullLength()));
	    int i1,i2,j1,j2;
	    TreeSet<Integer> unobserved1 = new TreeSet<Integer>();
	    TreeSet<Integer> unobserved2 = new TreeSet<Integer>();
	    Edge e1 = new Edge(0,0);
	    Edge e2 = new Edge(0,0);

	    // detect all unobserved residues
	    for(int i = 0; i < ali.getAlignmentLength(); ++i) {
		i1 = ali.al2seq(name1, i);
		i2 = ali.al2seq(name2, i);
		if( i1 != -1 && !hasCoordinates(i1) ) {
		    unobserved1.add(i1);
		}
		if( i2 != -1 && !pdb2.hasCoordinates(i2) ) {
		    unobserved2.add(i2);
		}
	    }

	    // strategy: we always have to look through the alignment to say 
	    // whether a difference distance can be assigned to a pair of 
	    // corresponding residues. To put it differently, for any two 
	    // alignment columns one always has to ensure that both columns 
	    // only contain observed residues (no gaps!), otherwise the one 
	    // cannot obtain a distance in at least one structure as a gap 
	    // indicates "no coordinates available".  
	    
	    for(int i = 0; i < ali.getAlignmentLength()-1; ++i) {
		
		i1 = ali.al2seq(name1, i);
		i2 = ali.al2seq(name2, i);

		// alignment columns must not contain gap characters and both 
		// residues in the current column have to be observed!
		if( i1 == -1 || i2 == -1 || unobserved1.contains(i1) || unobserved2.contains(i2) ) {
		    continue;
		}

		for(int j = i + 1; j < ali.getAlignmentLength(); ++j) {
		  
		    j1 = ali.al2seq(name1, j);
		    j2 = ali.al2seq(name2, j);
		    
		    if( j1 == -1 || j2 == -1 || unobserved1.contains(j1) || unobserved2.contains(j2) ) {
			continue;
		    }
		    
		    // make the edges
		    e1.i = i1;
		    e1.j = j1;
		    e2.i = i2;
		    e2.j = j2;
		    		    
		    alignedDistMatrix.put(new Edge(i+1,j+1),Math.abs(thisDistMatrix.get(e1)-otherDistMatrix.get(e2)));
		}
	    }
	    return alignedDistMatrix;
	}
	// TODO: Version of this where already buffered distance matrices are passed as paremeters

	/** Returns the number of neighbours of this grid cell */
	private int getNumGridNbs(HashMap<Point3i,Box> boxes, Point3i floor, int boxSize) {
		Point3i neighbor;
		int nbs = 0;
		for (int x=floor.x-boxSize;x<=floor.x+boxSize;x+=boxSize){
			for (int y=floor.y-boxSize;y<=floor.y+boxSize;y+=boxSize){
				for (int z=floor.z-boxSize;z<=floor.z+boxSize;z+=boxSize){
						neighbor = new Point3i(x,y,z);
						if (boxes.containsKey(neighbor)) nbs++;
				}
			}
		} 
		// compensate for counting myself as a neighbour
		if(boxes.containsKey(floor)) nbs--;
		return nbs;
	}
	
	/**
	 * Gets the Consurf-HSSP score given an internal residue serial
	 * @param resser
	 * @return
	 */
	public Double getConsurfhsspScoreFromResSerial(int resser){
		if (resser2consurfhsspscore != null) {
			if (resser2consurfhsspscore.get(resser) != null) {
				return resser2consurfhsspscore.get(resser);
			} else {
				return null;
			}		
		} else {
			return null;
		}
	}
	
	/**
	 * Gets the Consurf-HSSP color rsa given an internal residue serial
	 * @param resser
	 * @return
	 */
	public Integer getConsurfhsspColorFromResSerial(int resser){
		if (resser2consurfhsspcolor != null) {
			if (resser2consurfhsspcolor.get(resser) != null) {
				return resser2consurfhsspcolor.get(resser);
			} else {
				return null;
			}		
		} else {
			return null;
		}
	}
	
	/**
	 * Gets the all atoms rsa given an internal residue serial
	 * @param resser
	 * @return
	 */
	public Double getAllRsaFromResSerial(int resser){
		if (resser2allrsa != null) {
			if (resser2allrsa.get(resser) != null) {
				return resser2allrsa.get(resser);
			} else {
				return null;
			}		
		} else {
			return null;
		}
	}
	
	/**
	 * Gets the sc rsa given an internal residue serial
	 * @param resser
	 * @return
	 */
	public Double getScRsaFromResSerial(int resser){
		if (resser2scrsa != null) {
			if (resser2scrsa.get(resser) != null) {
				return resser2scrsa.get(resser);
			} else {
				return null;
			}		
		} else {
			return null;
		}
	}
	
	/**
	 * Gets the internal residue serial (cif) given a pdb residue serial (author assignment)
	 * TODO refactor
	 * @param pdbresser
	 * @return
	 */
	public int get_resser_from_pdbresser (String pdbresser){
		return pdbresser2resser.get(pdbresser);
	}
	
	/**
	 * Gets the pdb residue serial (author assignment) given an internal residue serial (cif)
	 * TODO refactor
	 * @param resser
	 * @return
	 */
	public String get_pdbresser_from_resser (int resser){
		return resser2pdbresser.get(resser);
	}

	/**
	 * Gets the residue serial given an atom serial
	 * TODO refactor
	 * @param atomser
	 * @return
	 */
	public int get_resser_from_atomser(int atomser){
		return atomser2resser.get(atomser);
	}
	
	public String getResTypeFromResSerial(int resser) {
		return resser2restype.get(resser);
	}
	
	/**
	 * Gets the atom serial given the residue serial and atom name
	 * @param resser
	 * @param atom
	 * @return
	 */
	public int getAtomSerFromResSerAndAtom(int resser, String atom) {
		return resser_atom2atomserial.get(resser+"_"+atom);
	}
	
	/**
	 * Checks whether the given residue serial has any associated coordinates.
	 * @param ser  the residue serial
	 * @return true if there is at least one atom with valid coordinates, else false 
	 */
	public boolean hasCoordinates(int resser) {
	    return atomser2resser.values().contains(resser);
	}
		
	/**
	 * Gets the atom coordinates (Point3d object) given the atom serial
	 * @param atomser
	 * @return
	 */
	public Point3d getAtomCoord(int atomser) {
		return this.atomser2coord.get(atomser);
	}
	
	/**
	 * Gets all atom serials in a Set
	 * @return
	 */
	public Set<Integer> getAllAtomSerials() {
		return this.atomser2resser.keySet();
	}
	
	/**
	 * Gets the 4 letter pdb code identifying this structure
	 * @return
	 */
	public String getPdbCode() {
		return this.pdbCode;
	}
	
	/**
	 * Gets the internal chain code (cif)
	 * @return
	 */
	public String getChainCode(){
		return this.chainCode;
	}
	
	/**
	 * Gets the pdb chain code (author assignment)
	 * @return
	 */
	public String getPdbChainCode(){
		return this.pdbChainCode;
	}
	
	/**
	 * Gets the sequence
	 * @return
	 */
	public String getSequence() {
		return sequence;
	}
	
	// csa related methods
	
	/** 
	 * Returns true if csa information is available, false otherwise. 
	 */
	public boolean hasCSA() {
		if (catalSiteSet == null) {
			return false;
		} else if (catalSiteSet.isEmpty()) {
			return false;
		} else {
			return true;
		}
	}
	
	/**
	 * Returns the csa annotation object of this graph.
	 */
	public CatalSiteSet getCSA() {
		return catalSiteSet;
	}
	
	// ec related methods
	
	/** 
	 * Returns true if ec information is available, false otherwise. 
	 */
	public boolean hasEC() {
		if (ec == null) {
			return false;
		} else if (ec.isEmpty()) {
			return false;
		} else {
			return true;
		}
	}
	
	/**
	 * Returns the ec annotation object of this graph.
	 */
	public EC getEC() {
		return ec;
	}
	
	// end of secop related methods
	
	// scop related methods
	
	/** 
	 * Returns true if scop information is available, false otherwise. 
	 */
	public boolean hasScop() {
		if (scop == null) {
			return false;
		} else if (scop.isEmpty()) {
			return false;
		} else {
			return true;
		}
	}
	
	/**
	 * Returns the scop annotation object of this graph.
	 */
	public Scop getScop() {
		return scop;
	}
	
	// end of secop related methods
	
	// secondary structure related methods
	
	/** 
	 * Returns true if secondary structure information is available, false otherwise. 
	 */
	public boolean hasSecondaryStructure() {
		return !this.secondaryStructure.isEmpty();
	}
	
	/**
	 * Returns the secondary structure annotation object of this graph.
	 */
	public SecondaryStructure getSecondaryStructure() {
		return this.secondaryStructure;
	}
	
	// end of secondary structure related methods
	
	/**
	 * Calculates rmsd (on atoms given by ct) of this Pdb object to otherPdb object
	 * Both objects must represent structures with same sequence (save unobserved residues or missing atoms)
	 * 
	 * @param otherPdb
	 * @param ct the contact type (crossed contact types don't make sense here)
	 * @return
	 * @throws ConformationsNotSameSizeError
	 */
	public double rmsd(Pdb otherPdb, String ct) throws ConformationsNotSameSizeError {
		TreeMap<String,Point3d> thiscoords = this.get_coords_for_ct_4rmsd(ct);
		TreeMap<String,Point3d> othercoords = otherPdb.get_coords_for_ct_4rmsd(ct);
		
		// there might be unobserved residues or some missing atoms for a residue
		// here we get the ones that are in common
		ArrayList<String> common = new ArrayList<String>();
		for (String resser_atom:thiscoords.keySet()){
			if (othercoords.containsKey(resser_atom)){
				common.add(resser_atom);
			}
		}
		
		// converting the TreeMaps to Vector3d arrays (input format for calculate_rmsd)
		Vector3d[] conformation1 = new Vector3d[common.size()]; 
		Vector3d[] conformation2 = new Vector3d[common.size()];
		int i = 0;
		for (String resser_atom:common){
			conformation1[i] = new Vector3d(thiscoords.get(resser_atom).x,thiscoords.get(resser_atom).y,thiscoords.get(resser_atom).z);
			conformation2[i] = new Vector3d(othercoords.get(resser_atom).x,othercoords.get(resser_atom).y,othercoords.get(resser_atom).z);
			i++;
		}
				
		// this as well as calculating the rmsd, changes conformation1 and conformation2 to be optimally superposed
		double rmsd = calculate_rmsd(conformation1, conformation2);

//		// printing out individual distances (conformation1 and conformation2 are now optimally superposed)
//		for (i=0;i<conformation1.length;i++){
//			Point3d point1 = new Point3d(conformation1[i].x,conformation1[i].y,conformation1[i].z);
//			Point3d point2 = new Point3d(conformation2[i].x,conformation2[i].y,conformation2[i].z);
//			System.out.println(point1.distance(point2));
//		}
		
		return rmsd;

	}
	
	/**
	 * Calculates the RMSD between two conformations.      
     * conformation1: Vector3d array (matrix of dimensions [N,3])       
     * conformation2: Vector3d array (matrix of dimensions [N,3]) 
     * 
     * Both conformation1 and conformation2 are modified to be optimally superposed
     * 
     * Implementation taken (python) from http://bosco.infogami.com/Root_Mean_Square_Deviation, 
     * then ported to java using Jama matrix package 
     * (page has moved to: http://boscoh.com/protein/rmsd-root-mean-square-deviation)                
	 * @param conformation1
	 * @param conformation2
	 * @return
	 * @throws ConformationsNotSameSizeError
	 */
	private double calculate_rmsd(Vector3d[] conformation1, Vector3d[] conformation2) throws ConformationsNotSameSizeError{
		if (conformation1.length!=conformation2.length) {
			//System.err.println("Conformations not the same size");
			throw new ConformationsNotSameSizeError(
					"Given conformations have different size: conformation1: "+conformation1.length+", conformation2: "+conformation2.length);
		}
		int n_vec = conformation1.length;
		
		// 1st we bring both conformations to the same centre by subtracting their respectives centres
		Vector3d center1 = new Vector3d();
		Vector3d center2 = new Vector3d();
		for (int i=0;i<n_vec;i++){ // summing all vectors in each conformation
			center1.add(conformation1[i]);
			center2.add(conformation2[i]);
		}
		// dividing by n_vec (average)
		center1.scale((double)1/n_vec);
		center2.scale((double)1/n_vec);
		// translating our conformations to the same coordinate system by subtracting centers
		for (Vector3d vec:conformation1){
			vec.sub(center1);
		}
		for (Vector3d vec:conformation2){
			vec.sub(center2);
		}

		//E0: initial sum of squared lengths of both conformations
		double sum1 = 0.0;
		double sum2 = 0.0;
		for (int i=0;i<n_vec;i++){
			sum1 += conformation1[i].lengthSquared();
			sum2 += conformation2[i].lengthSquared();
		}
		double E0 = sum1 + sum2;
		
		// singular value decomposition
		Matrix vecs1 = vector3dAr2matrix(conformation1);
		Matrix vecs2 = vector3dAr2matrix(conformation2);
		
		Matrix correlation_matrix = vecs2.transpose().times(vecs1); //gives a 3x3 matrix

		SingularValueDecomposition svd = correlation_matrix.svd();
		Matrix U = svd.getU();
		Matrix V_trans = svd.getV().transpose(); 
		double[] singularValues = svd.getSingularValues();
		
		boolean is_reflection = false;
		if (U.det()*V_trans.det()<0.0){ 
			is_reflection = true;
		}
		if (is_reflection){
			// reflect along smallest principal axis:
			// we change sign of last coordinate (smallest singular value)
			singularValues[singularValues.length-1]=(-1)*singularValues[singularValues.length-1];  			
		}
		
		// getting sum of singular values
		double sumSV = 0.0;
		for (int i=0;i<singularValues.length;i++){
			sumSV += singularValues[i];
		}
		
		// rmsd square: Kabsch formula
		double rmsd_sq = (E0 - 2.0*sumSV)/((double) n_vec);
		rmsd_sq = Math.max(rmsd_sq, 0.0);

		// finally we modify conformation2 to be aligned to conformation1
		if (is_reflection) { // first we check if we are in is_reflection case: we need to change sign to last row of U
			for (int j=0;j<U.getColumnDimension();j++){
				// we change sign to last row of U
				int lastRow = U.getRowDimension()-1;
				U.set(lastRow, j, (-1)*U.get(lastRow,j));
			}
		}
		Matrix optimal_rotation = U.times(V_trans); 
		Matrix conf2 = vecs2.times(optimal_rotation);
		for (int i=0;i<n_vec;i++){
			conformation2[i].x = conf2.get(i,0);
			conformation2[i].y = conf2.get(i,1);
			conformation2[i].z = conf2.get(i,2);
		}

		return Math.sqrt(rmsd_sq);
	}
	
	/** Gets a Jama.Matrix object from a Vector3d[] (deep copies) */
	private Matrix vector3dAr2matrix(Vector3d[] vecArray) {
		double[][] array = new double[vecArray.length][3];
		for (int i=0;i<vecArray.length;i++){
			vecArray[i].get(array[i]);
		}
		return new Matrix(array);
	}
	
	/**
	 * Write residue info to given db, using our db graph aglappe format, 
	 * i.e. tables: residue_info
	 * @param conn
	 * @param db
	 * @throws SQLException
	 */
	public void writeToDb(MySQLConnection conn, String db) throws SQLException{

		Statement stmt;
		String sql = "";
		
		conn.setSqlMode("NO_UNSIGNED_SUBTRACTION,TRADITIONAL");
		
		for (int resser:resser2restype.keySet()) {
			String resType = AAinfo.threeletter2oneletter(getResTypeFromResSerial(resser));
			String pdbresser = get_pdbresser_from_resser(resser);
			
			String secStructType = null;
			String secStructId = null;
			if (secondaryStructure != null) {
				if (secondaryStructure.getSecStrucElement(resser) != null) {
					secStructType = quote(String.valueOf(secondaryStructure.getSecStrucElement(resser).getType()));
					secStructId = quote(secondaryStructure.getSecStrucElement(resser).getId());
				}
			}
			
			String scopId = null;
			String sccs = null;
			String sunid = null;
			String orderIn = null;
			String domainType = null;
			String domainNumReg = null;
			if (scop != null) {
				if (scop.getScopRegion(resser)!=null) {
					ScopRegion sr = scop.getScopRegion(resser);
					scopId = quote(sr.getSId());
					sccs = quote(sr.getSccs());
					sunid = String.valueOf(sr.getSunid());
					orderIn = String.valueOf(sr.getOrder());
					domainType = quote(String.valueOf(sr.getDomainType()));
					domainNumReg = String.valueOf(sr.getNumRegions());
				}
			}
			
			Double allRsa = getAllRsaFromResSerial(resser);
			Double scRsa = getScRsaFromResSerial(resser);
			
			Double consurfhsspScore = getConsurfhsspScoreFromResSerial(resser);
			Integer consurfhsspColor = getConsurfhsspColorFromResSerial(resser);
			
			String ecId = null;
			if (ec != null) {
				if (ec.getECRegion(resser) != null) {
					ecId = quote(ec.getECNum(resser));
				}
			}
			
			String csaNums = null;
			String csaChemFuncs = null;
			String csaEvids = null;
			if (catalSiteSet != null) {
				if (catalSiteSet.getCatalSite(resser) != null) {
					csaNums = quote(catalSiteSet.getCatalSiteNum(resser));
					csaChemFuncs = quote(catalSiteSet.getCatalSiteChemFunc(resser));
					csaEvids = quote(catalSiteSet.getCatalSiteEvid(resser));
				}
			}
			
			sql = "INSERT IGNORE INTO "+db+".pdb_residue_info (pdb_code, chain_code, pdb_chain_code, res_ser, pdb_res_ser, res_type, sstype, ssid, scop_id, sccs, sunid, order_in, domain_type, domain_num_reg, all_rsa, sc_rsa, consurfhssp_score, consurfhssp_color, ec, csa_site_nums, csa_chem_funcs, csa_evid) " +
			" VALUES ("+quote(pdbCode)+", "+quote(chainCode)+", "+(pdbChainCode.equals("NULL")?quote("-"):quote(pdbChainCode))+","+resser+", "+quote(pdbresser)+", "+quote(resType)+", "+secStructType+", "+secStructId+", "+scopId+", "+sccs+", "+sunid+", "+orderIn+", "+domainType+", "+domainNumReg+", "+allRsa+", "+scRsa+", "+consurfhsspScore+","+consurfhsspColor+","+ecId+","+csaNums+","+csaChemFuncs+","+csaEvids+")";
			//System.out.println(sql);
			stmt = conn.createStatement();
			stmt.executeUpdate(sql);
			stmt.close();
		}			
	}

	/**
	 * Write residue info to given db, using our db graph aglappe format, 
	 * i.e. tables: residue_info
	 * @param conn
	 * @param db
	 * @throws SQLException
	 */
	public void writeToDbFast(MySQLConnection conn, String db) throws SQLException, IOException {

		Statement stmt;
		String sql = "";
		
		conn.setSqlMode("NO_UNSIGNED_SUBTRACTION,TRADITIONAL");
		
		PrintStream resOut = new PrintStream(new FileOutputStream(pdbCode+chainCode+"_residues.txt"));
		
		for (int resser:resser2restype.keySet()) {
			String resType = AAinfo.threeletter2oneletter(getResTypeFromResSerial(resser));
			String pdbresser = get_pdbresser_from_resser(resser);
			
			String secStructType = "\\N";
			String secStructId = "\\N";
			if (secondaryStructure != null) {
				if (secondaryStructure.getSecStrucElement(resser) != null) {
					secStructType = String.valueOf(secondaryStructure.getSecStrucElement(resser).getType());
					secStructId = secondaryStructure.getSecStrucElement(resser).getId();
				}
			}
			
			String scopId = "\\N";
			String sccs = "\\N";
			String sunid = "\\N";
			String orderIn = "\\N";
			String domainType = "\\N";
			String domainNumReg = "\\N";
			if (scop != null) {
				if (scop.getScopRegion(resser)!=null) {
					ScopRegion sr = scop.getScopRegion(resser);
					scopId = sr.getSId();
					sccs = sr.getSccs();
					sunid = String.valueOf(sr.getSunid());
					orderIn = String.valueOf(sr.getOrder());
					domainType = String.valueOf(sr.getDomainType());
					domainNumReg = String.valueOf(sr.getNumRegions());
				}
			}
			
			Double allRsa = getAllRsaFromResSerial(resser);
			Double scRsa = getScRsaFromResSerial(resser);
			
			Double consurfhsspScore = getConsurfhsspScoreFromResSerial(resser);
			Integer consurfhsspColor = getConsurfhsspColorFromResSerial(resser);
			
			String ecId = "\\N";
			if (ec != null) {
				if (ec.getECRegion(resser) != null) {
					ecId = ec.getECNum(resser);
				}
			}
			
			String csaNums = "\\N";
			String csaChemFuncs = "\\N";
			String csaEvids = "\\N";
			if (catalSiteSet != null) {
				if (catalSiteSet.getCatalSite(resser) != null) {
					csaNums = catalSiteSet.getCatalSiteNum(resser);
					csaChemFuncs = catalSiteSet.getCatalSiteChemFunc(resser);
					csaEvids = catalSiteSet.getCatalSiteEvid(resser);
				}
			}
			
			resOut.println(pdbCode+"\t"+chainCode+"\t"+(pdbChainCode.equals("NULL")?"-":pdbChainCode)+"\t"+resser+"\t"+pdbresser+"\t"+resType+"\t"+secStructType+"\t"+secStructId+"\t"+scopId+"\t"+sccs+"\t"+sunid+"\t"+orderIn+"\t"+domainType+"\t"+domainNumReg+"\t"+allRsa+"\t"+scRsa+"\t"+consurfhsspScore+"\t"+consurfhsspColor+"\t"+ecId+"\t"+csaNums+"\t"+csaChemFuncs+"\t"+csaEvids);
			
		}
		resOut.close();
		sql = "LOAD DATA LOCAL INFILE '"+pdbCode+chainCode+"_residues.txt' INTO TABLE "+db+".pdb_residue_info (pdb_code, chain_code, pdb_chain_code, res_ser, pdb_res_ser, res_type, sstype, ssid, scop_id, sccs, sunid, order_in, domain_type, domain_num_reg, all_rsa, sc_rsa, consurfhssp_score, consurfhssp_color, ec, csa_site_nums, csa_chem_funcs, csa_evid);";
		//System.out.println(sql);
		stmt = conn.createStatement();
		stmt.executeUpdate(sql);
		stmt.close();
		File fileToDelete = new File(pdbCode+chainCode+"_residues.txt");
		if (fileToDelete.exists()) {
			fileToDelete.delete();
		}
	}

	private static String quote(String s) {
		return ("'"+s+"'");
	}
}


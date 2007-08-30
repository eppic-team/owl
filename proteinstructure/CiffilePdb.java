package proteinstructure;

import java.io.BufferedReader;
import java.io.File;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Point3d;


/**
 * A single chain pdb protein structure loaded from a mmCIF file 
 * 
 * @author		Jose Duarte
 * Class:		CiffilePdb
 * Package:		proteinstructure
 */
public class CiffilePdb extends Pdb {


	private String ciffile;

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
	private TreeSet<Integer> loopElements; 							// contains list of elements that are of loop type
	private TreeMap<Integer,Interval> loopelements2contentIndex;    // begin and end line index of each loop element
 
	/**
	 * Constructs Pdb object given pdb code and pdb chain code. 
	 * Model will be DEFAULT_MODEL
	 * @param ciffile
	 * @param pdbChainCode
	 * @throws PdbChainCodeNotFoundError
	 * @throws IOException 
	 */
	public CiffilePdb (String ciffile, String pdbChainCode) throws PdbChainCodeNotFoundError, IOException {
		this(ciffile, pdbChainCode, DEFAULT_MODEL);
	}

	/**
	 * Constructs Pdb object given pdb code, pdb chain code, model serial, source db and a MySQLConnection.
	 * The db must be a pdbase database 
	 * @param ciffile
	 * @param pdbChainCode
	 * @param model_serial
	 * @throws PdbChainCodeNotFoundError
	 * @throws IOException 
	 */
	public CiffilePdb (String ciffile, String pdbChainCode, int model_serial) throws PdbChainCodeNotFoundError, IOException {
		this.ciffile = ciffile;
		this.pdbChainCode=pdbChainCode.toUpperCase();	// our convention: chain codes are upper case
		this.model=model_serial;
		
		parseCifFile();	
		
		this.pdbCode = readPdbCode();

		secondaryStructure = new SecondaryStructure();	// create empty secondary structure first to make sure object is not null
		
		parseLoopElements(); // populates resser_atom2atomserial, resser2restype, atomser2coord, atomser2resser, sequence, pdbresser2resser, secondaryStructure
		
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

	private void parseCifFile() throws IOException{
		// data structures to store the parsed fields
		ids2elements = new TreeMap<String, Integer>();
		fields2indices = new TreeMap<String,Integer>();
		fields2values = new TreeMap<String, String>();
		loopElements = new TreeSet<Integer>(); // contains list of elements that are of loop type
		loopelements2contentIndex = new TreeMap<Integer,Interval>();
		TreeMap<String,Integer> fieldsIdx = new TreeMap<String,Integer>(); // this map holds the field index counters for each element id
		
		BufferedReader fcif = new BufferedReader(new FileReader(new File(ciffile)));
		int element = 0;
		String line;
		line = fcif.readLine(); // skip first line
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
				if (!fieldsIdx.containsKey(id)) fieldsIdx.put(id,0);
				Pattern p = Pattern.compile("^"+id+"\\.(\\w+)(?:\\s+(.*))?$");
				Matcher m = p.matcher(line);
				if (m.find()){
					ids2elements.put(id,element);
					String field = id + "." + m.group(1);
					if (!loopElements.contains(element)) { // if not a loop element
						fields2values.put(field, m.group(2)); // 2nd capture group only matches for non-loops where the value of the field is in same line as field name
					} else { // for loop elements we fill the fields2indices TreeMap
						fields2indices.put(field,fieldsIdx.get(id));
					}
					fieldsIdx.put(id,fieldsIdx.get(id)+1); 
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
		return fields2values.get("_entry.id").trim();
	}
	
	private void parseLoopElements() throws IOException, PdbChainCodeNotFoundError {
		resser_atom2atomserial = new HashMap<String,Integer>();
		resser2restype = new HashMap<Integer,String>();
		atomser2coord = new HashMap<Integer,Point3d>();
		atomser2resser = new HashMap<Integer,Integer>();
		pdbresser2resser = new HashMap<String, Integer>();
		sequence = "";
		secondaryStructure = new SecondaryStructure();
		
		ArrayList<String> aalist=AA.aas(); // list of standard 3 letter code aminoacids
		
		String chainCodeStr=pdbChainCode;
		if (pdbChainCode.equals("NULL")) chainCodeStr="A";
		
		ArrayList<String> altLocs = new ArrayList<String>();
		String altLoc = ".";

		Interval intAtomSite = loopelements2contentIndex.get(ids2elements.get(atomSiteId));
		// atom_sites_alt element is optional
		Interval intAtomSitesAlt = null;
		if (ids2elements.containsKey(atomSitesAltId)){
			intAtomSitesAlt = loopelements2contentIndex.get(ids2elements.get(atomSitesAltId));
		}
		Interval intPdbxPoly = loopelements2contentIndex.get(ids2elements.get(pdbxPolySeqId));
		// struct_conf element is optional
		Interval intStructConf = null;
		if (ids2elements.containsKey(structConfId)){
			intStructConf = loopelements2contentIndex.get(ids2elements.get(structConfId));
		}
		// struct_sheet_range element is optional
		Interval intStructSheet = null; 
		if (ids2elements.containsKey(structSheetId)) {
			intStructSheet = loopelements2contentIndex.get(ids2elements.get(structSheetId));
		}
		
		boolean empty = true;
		BufferedReader fcif = new BufferedReader(new FileReader(new File(ciffile)));
		String line;
		int linecount=0;
		while((line = fcif.readLine()) != null ) {
			linecount++; 
			// atom_sites_alt (optional element)
			if (intAtomSitesAlt!=null && linecount>=intAtomSitesAlt.beg && linecount<=intAtomSitesAlt.end){
				int idIdx = fields2indices.get(atomSitesAltId+".id");
				// id=0
				// A ?
				String[] tokens = line.split("\\s+");
				if (!tokens[idIdx].equals(".")) {
					altLocs.add(tokens[idIdx]);
				}
				continue;
			} 
			if (!altLocs.isEmpty()){
				altLoc = Collections.min(altLocs);
			}
			// atom_site
			if (linecount>=intAtomSite.beg && linecount<=intAtomSite.end){ 
				int idIdx = fields2indices.get(atomSiteId+".id");
				int labelAtomIdIdx = fields2indices.get(atomSiteId+".label_atom_id");
				int labelAltIdIdx = fields2indices.get(atomSiteId+".label_alt_id");
				int labelCompIdIdx = fields2indices.get(atomSiteId+".label_comp_id");
				int labelAsymIdIdx = fields2indices.get(atomSiteId+".label_asym_id");
				int labelSeqIdIdx = fields2indices.get(atomSiteId+".label_seq_id");
				int cartnXIdx = fields2indices.get(atomSiteId+".Cartn_x");
				int cartnYIdx = fields2indices.get(atomSiteId+".Cartn_y");
				int cartnZIdx = fields2indices.get(atomSiteId+".Cartn_z");
				int authAsymIdIdx = fields2indices.get(atomSiteId+".auth_asym_id");
				int pdbxPDBModelNumIdx = fields2indices.get(atomSiteId+".pdbx_PDB_model_num");
				// auth_asym_id=22, pdbx_PDB_model_num=24, label_alt_id=4, id=1, label_atom_id=3, label_comp_id=5, label_asym_id=6, label_seq_id=8, Cartn_x=10, Cartn_y=11, Cartn_z=12
				//   0   1    2  3  4   5 6 7 8  9     10    11       12    13    14 151617181920   2122 23 24
				//ATOM   2    C CA  . MET A 1 1  ? 38.591 8.543   15.660  1.00 77.79  ? ? ? ? ? 1  MET A CA  1
				String[] tokens = line.split("\\s+");
				if (tokens[authAsymIdIdx].equals(chainCodeStr) && Integer.parseInt(tokens[pdbxPDBModelNumIdx])==model) { // match our given chain and model 
					empty = false;
					if (tokens[labelAltIdIdx].equals(".") || tokens[labelAltIdIdx].equals(altLoc)) { // don't read lines with something else as "." or altLoc
						int atomserial=Integer.parseInt(tokens[idIdx]); // id
						String atom = tokens[labelAtomIdIdx]; // label_atom_id
						String res_type = tokens[labelCompIdIdx]; // label_comp_id
						chainCode = tokens[labelAsymIdIdx]; // label_asym_id
						int res_serial = Integer.parseInt(tokens[labelSeqIdIdx]); // label_seq_id
						double x = Double.parseDouble(tokens[cartnXIdx]); // Cartn_x
						double y = Double.parseDouble(tokens[cartnYIdx]); // Cartn_y
						double z = Double.parseDouble(tokens[cartnZIdx]); // Cartn_z
						Point3d coords = new Point3d(x,y,z);
						if (aalist.contains(res_type)) {
							atomser2coord.put(atomserial, coords);
							atomser2resser.put(atomserial, res_serial);
							resser2restype.put(res_serial, res_type);
							ArrayList<String> atomlist = aas2atoms.get(res_type);
							atomlist.add("OXT"); // the extra atom OXT is there in the last residue of the chain
							if (atomlist.contains(atom)){
								resser_atom2atomserial.put(res_serial+"_"+atom, atomserial);
							}
						}
					}
				}
				continue;
			}
			// pdbx_poly_seq_scheme
			if (linecount>=intPdbxPoly.beg && linecount<=intPdbxPoly.end){
				int asymIdIdx = fields2indices.get(pdbxPolySeqId+".asym_id");
				int seqIdIdx = fields2indices.get(pdbxPolySeqId+".seq_id");
				int authSeqNumIdx = fields2indices.get(pdbxPolySeqId+".auth_seq_num");
				int pdbInsCodeIdx = fields2indices.get(pdbxPolySeqId+".pdb_ins_code");
				int monIdIdx = fields2indices.get(pdbxPolySeqId+".mon_id");
				// asym_id=0, seq_id=2, auth_seq_num=6, pdb_ins_code=10, mon_id=3 
				// 0 1 2     3 4   5   6     7   8 910
				// A 1 1   ASP 1   1   1   ASP ASP A .
				String[] tokens = line.split("\\s+");
				if (tokens[asymIdIdx].equals(chainCode)) { // we already have chainCode from _atom_site element that should appear before pdbx_poly_seq_scheme in cif file
					int res_serial = Integer.parseInt(tokens[seqIdIdx]); // seq_id
					String pdb_res_serial = tokens[authSeqNumIdx]; // auth_seq_num
					String pdb_ins_code = tokens[pdbInsCodeIdx]; // pdb_ins_code
					String pdb_res_serial_with_icode = pdb_res_serial;
					if (!pdb_ins_code.equals(".")) {
						pdb_res_serial_with_icode=pdb_res_serial+pdb_ins_code;
					}
					String res_type = tokens[monIdIdx]; // mon_id
					// sequence
					if (aalist.contains(res_type)){
		        		sequence+=AA.threeletter2oneletter(res_type);
		        	} else {
		        		sequence+=NONSTANDARD_AA_LETTER;
		        	}
					// pdbresser2resser
					pdbresser2resser.put(pdb_res_serial_with_icode,res_serial);
				}
				continue;
			}
			// struct_conf (optional element), HELIX and TURN secondary structure
			if (intStructConf!=null && linecount>=intStructConf.beg && linecount<=intStructConf.end){
				int idIdx = fields2indices.get(structConfId+".id");
				int begLabelAsymIdIdx = fields2indices.get(structConfId+".beg_label_asym_id");
				int begLabelSeqIdIdx = fields2indices.get(structConfId+".beg_label_seq_id");
				int endLabelSeqIdIdx = fields2indices.get(structConfId+".end_label_seq_id");
				//id=1, beg_label_seq_id=5, end_label_seq_id=9, beg_label_asym_id=4
				//     0       1  2    3 4 5   6   7 8  9 10  111213    1415 16 1718 19
				//HELX_P HELX_P1  1  ASN A 2   ? GLY A 12  ? ASN A 2   GLY A 12  1 ? 11
				String[] tokens = line.split("\\s+");
				if (tokens[begLabelAsymIdIdx].equals(chainCode)) {
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
				String[] tokens = line.split("\\s+");
				if (tokens[begLabelAsymIdIdx].equals(chainCode)){
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
		if (empty) { // no atom data was found for given pdb chain code and model
			throw new PdbChainCodeNotFoundError("Couldn't find _atom_site data for given pdbChainCode: "+pdbChainCode+", model: "+model);
		}
	}
}

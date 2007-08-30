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
	
	private TreeMap<String,Integer> ids2elements;					// map of ids to element serials
	private TreeMap<Integer,ArrayList<String>> elements2fields;		// map of element serials to member fields
	private TreeMap<String,String> fields2values;					// map of fields (id.field) to values (only for non-loop elements)
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
		// fields we will read
		String[] ids = {"_entry","_atom_sites_alt","_atom_site","_pdbx_poly_seq_scheme","_struct_conf","_struct_sheet_range"};
		
		// data structures to store the parsed fields
		ids2elements = new TreeMap<String, Integer>();
		elements2fields = new TreeMap<Integer, ArrayList<String>>();
		fields2values = new TreeMap<String, String>();
		loopElements = new TreeSet<Integer>(); // contains list of elements that are of loop type
		loopelements2contentIndex = new TreeMap<Integer,Interval>();
		
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
				Pattern p = Pattern.compile("^"+id+"\\.(\\w+)(?:\\s+(.*))?$");
				Matcher m = p.matcher(line);
				if (m.find()){
					ids2elements.put(id,element);
					if (!elements2fields.containsKey(element)) {
						ArrayList<String> fields= new ArrayList<String>();
						fields.add(m.group(1));
						elements2fields.put(element,fields);
					} else {
						elements2fields.get(element).add(m.group(1));
					}
					if (!loopElements.contains(element)) { // if not a loop element
						fields2values.put(id+"."+m.group(1), m.group(2));
					}
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
	
		Interval intAtomSite = loopelements2contentIndex.get(ids2elements.get("_atom_site"));
		// atom_sites_alt element is optional
		Interval intAtomSitesAlt = null;
		if (ids2elements.containsKey("_atom_sites_alt")){
			intAtomSitesAlt = loopelements2contentIndex.get(ids2elements.get("_atom_sites_alt"));
		}
		Interval intPdbxPoly = loopelements2contentIndex.get(ids2elements.get("_pdbx_poly_seq_scheme"));
		// struct_conf element is optional
		Interval intStructConf = null;
		if (ids2elements.containsKey("_struct_conf")){
			intStructConf = loopelements2contentIndex.get(ids2elements.get("_struct_conf"));
		}
		// struct_sheet_range element is optional
		Interval intStructSheet = null; 
		if (ids2elements.containsKey("_struct_sheet_range")) {
			intStructSheet = loopelements2contentIndex.get(ids2elements.get("_struct_sheet_range"));
		}
		
		boolean empty = true;
		BufferedReader fcif = new BufferedReader(new FileReader(new File(ciffile)));
		String line;
		int linecount=0;
		while((line = fcif.readLine()) != null ) {
			linecount++; 
			// atom_sites_alt (optional element)
			if (intAtomSitesAlt!=null && linecount>=intAtomSitesAlt.beg && linecount<=intAtomSitesAlt.end){
				String[] tokens = line.split("\\s+");
				if (!tokens[0].equals(".")) {
					altLocs.add(tokens[0]);
				}
				continue;
			} 
			if (!altLocs.isEmpty()){
				altLoc = Collections.min(altLocs);
			}
			// atom_site
			if (linecount>=intAtomSite.beg && linecount<=intAtomSite.end){
				// auth_asym_id=22, pdbx_PDB_model_num=24, label_alt_id=4, id=1, label_atom_id=3, label_comp_id=5, label_asym_id=6, label_seq_id=8, Cartn_x=10, Cartn_y=11, Cartn_z=12
				//   0   1    2  3  4   5 6 7 8  9     10    11       12    13    14 151617181920   2122 23 24
				//ATOM   2    C CA  . MET A 1 1  ? 38.591 8.543   15.660  1.00 77.79  ? ? ? ? ? 1  MET A CA  1
				String[] tokens = line.split("\\s+");
				if (tokens[22].equals(chainCodeStr) && Integer.parseInt(tokens[24])==model) { // match our given chain and model 
					empty = false;
					if (tokens[4].equals(".") || tokens[4].equals(altLoc)) { // don't read lines with something else as "." or altLoc
						int atomserial=Integer.parseInt(tokens[1]); // id
						String atom = tokens[3]; // label_atom_id
						String res_type = tokens[5]; // label_comp_id
						chainCode = tokens[6]; // label_asym_id
						int res_serial = Integer.parseInt(tokens[8]); // label_seq_id
						double x = Double.parseDouble(tokens[10]); // Cartn_x
						double y = Double.parseDouble(tokens[11]); // Cartn_y
						double z = Double.parseDouble(tokens[12]); // Cartn_z
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
				// asym_id=0, seq_id=2, auth_seq_num=6, pdb_ins_code=10, mon_id=3 
				// 0 1 2     3 4   5   6     7   8 910
				// A 1 1   ASP 1   1   1   ASP ASP A .
				String[] tokens = line.split("\\s+");
				if (tokens[0].equals(chainCode)) { // we already have chainCode from _atom_site element that should appear before pdbx_poly_seq_scheme in cif file
					int res_serial = Integer.parseInt(tokens[2]); // seq_id
					String pdb_res_serial = tokens[6]; // auth_seq_num
					String pdb_ins_code = tokens[10]; // pdb_ins_code
					String pdb_res_serial_with_icode = pdb_res_serial;
					if (!pdb_ins_code.equals(".")) {
						pdb_res_serial_with_icode=pdb_res_serial+pdb_ins_code;
					}
					String res_type = tokens[3]; // mon_id
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
				//id=1, beg_label_seq_id=5, end_label_seq_id=9, beg_label_asym_id=4
				//     0       1  2    3 4 5   6   7 8  9 10  111213    1415 16 1718 19
				//HELX_P HELX_P1  1  ASN A 2   ? GLY A 12  ? ASN A 2   GLY A 12  1 ? 11
				String[] tokens = line.split("\\s+");
				if (tokens[4].equals(chainCode)) {
					String id = tokens[1];
					Pattern p = Pattern.compile("^(\\w).+_P(\\d)+$");
					Matcher m = p.matcher(id);
					String ssId="Unknown";
					if (m.find()){
						ssId = m.group(1)+m.group(2); // e.g.: Hnn (helices) or Tnn (turns) 				
					}
					int beg = Integer.parseInt(tokens[5]);
					int end = Integer.parseInt(tokens[9]);
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
				//sheet_id=0, id=1, beg_label_seq_id=4, end_label_seq_id=8, beg_label_asym_id=3
				//0 1   2 3  4 5   6 7  8 910  1112 13  1415 16
				//A 1 ARG A 14 ? LYS A 19 ? ? ARG A 14 LYS A 19
				String[] tokens = line.split("\\s+");
				if (tokens[3].equals(chainCode)){
					String sheetid = tokens[0];
					int id = Integer.parseInt(tokens[1]);
					int beg = Integer.parseInt(tokens[4]);
					int end = Integer.parseInt(tokens[8]);
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

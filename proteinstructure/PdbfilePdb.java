package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Point3d;

/**
 * A single chain PDB protein structure read from a PDB file or CASP TS file
 * We are not tolerant to errors in PDB files, we don't accept PDB files if:
 * - insertion codes are present
 * - residue numbers<=0 present
 * - ATOM lines are not in residue ascending order
 * - SEQRES present and sequence doesn't coincide with ATOM lines sequence
 * - all residues are non-standard for given chain (0 observed residues)
 */
public class PdbfilePdb extends Pdb {
	
	private static final String NULL_chainCode = "A";
	// the regex to use for selecting an atom based on its alt code: either a space or an A are valid
	private static final String ALTCODEREGEX = "[ A]";
	
	private String pdbfile;
	private boolean isCaspTS; // whether we are reading a CASP TS file (true) or a normal PDB file (false)
	private boolean hasSeqRes; // whether we find a non-empty SEQRES field in this pdb file
	
	/**
	 * Constructs an empty Pdb object given a pdbfile name
	 * Data will be loaded from pdb file upon call of load(pdbChainCode, modelSerial) 
	 * @param pdbfile
	 */
	public PdbfilePdb (String pdbfile) {
		this.pdbfile = pdbfile;
		this.isCaspTS = false; //we assume by default this is not a CASP TS file, when reading we set it to true if we detect CASP TS headers
		this.hasSeqRes = false; // we assume that there is no SEQRES in this pdb file, if SEQRES is found when reading we set it to true
		
		// we initialise the secondary structure to empty, if no sec structure info is found then it remains empty
		this.secondaryStructure = new SecondaryStructure("");		
	}

	/**
	 * Loads PDB data (coordinates, sequence, etc.) from the PDB file
	 * for given pdbChainCode and modelSerial
	 * @param pdbChainCode
	 * @param modelSerial
	 * @throws PdbLoadError
	 */
	public void load(String pdbChainCode, int modelSerial) throws PdbLoadError {
		try {
			initialiseResidues();
			this.model=modelSerial;
			this.pdbChainCode=pdbChainCode;			// NOTE! pdb chain codes are case sensitive!
			// we set chainCode to pdbChainCode except for case Pdb.NULL_CHAIN_CODE where we use "A"
			this.chainCode=pdbChainCode;
			if (pdbChainCode.equals(Pdb.NULL_CHAIN_CODE)) this.chainCode=NULL_chainCode;

			this.sequence = ""; // we initialize to blank so we can append to the string in read_pdb_data_from_file
			parse();
			
			secondaryStructure.setSequence(this.sequence);
			if(!secondaryStructure.isEmpty()) {
				secondaryStructure.setComment("Author");
				this.initialiseResiduesSecStruct();
			}
			
			// when reading from pdb file we have no information of residue numbers or author's (original) pdb residue number, so we fill the mapping with the residue numbers we know
			this.resser2pdbresser = new TreeMap<Integer, String>();
			this.pdbresser2resser = new TreeMap<String, Integer>();
			for (int resser:getAllSortedResSerials()){
				resser2pdbresser.put(resser, String.valueOf(resser));
				pdbresser2resser.put(String.valueOf(resser), resser);
			}
			this.initialiseMaps();
			dataLoaded = true;
			
		} catch (PdbfileFormatError e) {
			throw new PdbLoadError(e);
		} catch (IOException e) {
			throw new PdbLoadError(e);
		} catch (PdbChainCodeNotFoundError e) {
			throw new PdbLoadError(e);
		}
	}
	
	/**
	 * Returns all PDB chain codes present in the PDB file
	 * @return array with all pdb chain codes
	 */
	public String[] getChains() throws PdbLoadError {
		TreeSet<String> chains = new TreeSet<String>();
		try {
			BufferedReader fpdb = new BufferedReader(new FileReader(new File(pdbfile)));
			String  line;
			while ((line=fpdb.readLine())!=null) {
				if (line.startsWith("ATOM")) {
					String chain = line.substring(21, 22);
					if (chain.equals(" ")) chain=Pdb.NULL_CHAIN_CODE;
					chains.add(chain);
				}
			}
			fpdb.close();
		} catch (IOException e) {
			throw new PdbLoadError(e);
		}
		
		if (chains.isEmpty()) return null;
		
		String[] chainsArray = new String[chains.size()];
		chains.toArray(chainsArray);
		return chainsArray;
	}
	
	/**
	 * Returns all model serials present in the PDB file
	 * @return array with all model serials
	 */
	public Integer[] getModels() throws PdbLoadError {
		TreeSet<Integer> models = new TreeSet<Integer>();
		try {
			BufferedReader fpdb = new BufferedReader(new FileReader(new File(pdbfile)));
			String  line;
			while ((line=fpdb.readLine())!=null) {
				// The model serial numbers should occur in columns 11-14 (official PDB format spec)
				// Here we are less strict: we allow for the numbers to appear in any column after the MODEL keyword (with any number of spaces in between)
				// We mainly need this because of CASP TS file which don't follow the PDB format 100% 
				Pattern p = Pattern.compile("^MODEL\\s+(\\d+)");
				Matcher m = p.matcher(line);
				if (m.find()){
					int model = Integer.parseInt(m.group(1));
					models.add(model);
				}
			}
			fpdb.close();
		} catch (IOException e) {
			throw new PdbLoadError(e);
		} catch (NumberFormatException e) {
			throw new PdbLoadError("Wrong format for MODEL lines!");
		}
		
		if (models.isEmpty()) models.add(DEFAULT_MODEL);//return null;		
		Integer[] modelsArray = new Integer[models.size()];
		models.toArray(modelsArray);
		return modelsArray;
	}
	
	/**
	 * To read the pdb data (atom coordinates, residue serials, atom serials) from file.
	 * <code>chainCode</code> gets set to same as <code>pdbChainCode</code>, except if input chain code 
	 * is Pdb.NULL_CHAIN_CODE then chainCode will be <code>NULL_chainCode</code>
	 * <code>pdbCode</code> gets set to the one parsed in HEADER or to <code>NO_PDB_CODE</code> 
	 * if not found
	 * The sequence is either read from SEQRES if present or from the residues read from ATOM 
	 * lines using <code>Pdb.UNKNOWN_UNOBSERVED_RES_LETTER</code> as gaps.
	 * We are not tolerant to errors in PDB files: PdbfileFormatErrors will be thrown if:
	 * - insertion codes are present
	 * - residue numbers<=0 present
	 * - ATOM lines are not in residue ascending order
	 * - SEQRES present and sequence doesn't coincide with ATOM lines sequence
	 * - all residues are non-standard for given chain (0 observed residues)
	 * We never renumber the residues: we take them as they are and check all of the above
	 * @param pdbfile the PDB file name
	 * @throws IOException
	 * @throws PdbfileFormatError if file is empty, if file is a CASP TS file 
	 * and no TARGET line found, if sequences from ATOM lines and SEQRES do not coincide, 
	 * if an insertion code is found, if a residue number<=0 found, if ATOM lines are not 
	 * in residue ascending order, if 0 observed residues are found for given chain
	 * @throws PdbChainCodeNotFoundError if no ATOM lines are found for given 
	 * pdbChainCode and model
	 */
	private void parse() throws IOException, PdbfileFormatError, PdbChainCodeNotFoundError {
		Pattern p;
		Matcher m;
		boolean empty = true; // controls whether we don't find any atom line for given pdbChainCode and model
		// we set chainCodeStr (for regex) to pdbChainCode except for case Pdb.NULL_CHAIN_CODE where we use " " (Pdb.NULL_CHAIN_CODE is a blank chain code in pdb files)
		String chainCodeStr=pdbChainCode;
		if (pdbChainCode.equals(Pdb.NULL_CHAIN_CODE)) chainCodeStr=" ";
		
		int lastResSerial = 0; // we store the last residue serial read from the ATOM lines to check for correct ascending order
		int lastAtomSerial = -1; // same for atoms
		boolean atomAtOriginSeen = false; // if we've read at least 1 atom at the origin (0,0,0) it is set to true
		int thismodel=DEFAULT_MODEL; // we initialise to DEFAULT_MODEL, in case file doesn't have MODEL lines 
		BufferedReader fpdb = new BufferedReader(new FileReader(new File(pdbfile)));
		int linecount=0;
		String line;
		LinkedList<String> parentList = null;

		while((line = fpdb.readLine()) != null ) {
			linecount++;
			if (linecount==1) {
				// HEADER
				p = Pattern.compile("^HEADER");
				m = p.matcher(line);
				if (m.find()){
					Pattern ph = Pattern.compile("^HEADER.{56}(\\d\\w{3})");
					Matcher mh = ph.matcher(line);
					if (mh.find()) {
						pdbCode=mh.group(1).toLowerCase();
					}
				} else { // header not found
					// check whether this is a Casp prediction file
					p = Pattern.compile("^PFRMAT\\s+TS");
					m = p.matcher(line);
					if(m.find()) {
						// ok, it is
						isCaspTS = true; 
						pdbCode = NO_PDB_CODE;
						// we try to read the TARGET from the next line, if there's no TARGET line appearing this is not respecting the format: exception
						if((line = fpdb.readLine()) != null ) {
							linecount++;
							p = Pattern.compile("^TARGET\\s+[Tt](\\d+)");
							m = p.matcher(line);
							if (m.find()) {
								this.targetNum = Integer.parseInt(m.group(1));
							} else {
								throw new PdbfileFormatError("The CASP TS file "+pdbfile+" does not have a TARGET line");
							}
						} else {
							throw new PdbfileFormatError("The CASP TS file "+pdbfile+" is empty after the PFRMAT line");
						}
					}
				}
			}
			// SEQRES
			//SEQRES   1 A  348  VAL ASN ILE LYS THR ASN PRO PHE LYS ALA VAL SER PHE
			p = Pattern.compile("^SEQRES.{5}"+chainCodeStr);
			m = p.matcher(line);
			if (m.find()){
				for (int i=19;i<=67;i+=4) {
					if (!line.substring(i, i+3).equals("   ")) {
						if (AAinfo.isValidAA(line.substring(i, i+3))) { // for non-standard aas
							sequence+= AAinfo.threeletter2oneletter(line.substring(i, i+3));
						} else {
							sequence+= AAinfo.NONSTANDARD_AA_ONE_LETTER;
						}
					}
				}
				if (!sequence.equals("")) {// if SEQRES was not empty then we have a sequence
					hasSeqRes = true;
				}
			}
			// SECONDARY STRUCTURE
			// helix
			//HELIX    1   1 LYS A   17  LEU A   26  1
			//							helix ser				beg res ser					end res ser
			p = Pattern.compile("^HELIX..(...).{9}"+chainCodeStr+".(....).{6}"+chainCodeStr+".(....)");
			m = p.matcher(line);
			if (m.find()){
				int serial = Integer.valueOf(m.group(1).trim());
				int beg = Integer.valueOf(m.group(2).trim());
				int end = Integer.valueOf(m.group(3).trim());
				String ssId = new Character(SecStrucElement.HELIX).toString()+serial;
				SecStrucElement ssElem = new SecStrucElement(SecStrucElement.HELIX,beg,end,ssId);
				secondaryStructure.add(ssElem);
			}
			// sheet
			//SHEET    2   A 5 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107
			//                       strand ser sheet id			 beg res ser                 end res ser
			p = Pattern.compile("^SHEET..(...).(...).{7}"+chainCodeStr+"(....).{6}"+chainCodeStr+"(....)");
			m = p.matcher(line);
			if (m.find()){
				int strandSerial = Integer.valueOf(m.group(1).trim());
				String sheetId = m.group(2).trim();
				int beg = Integer.valueOf(m.group(3).trim());
				int end = Integer.valueOf(m.group(4).trim());
				String ssId = new Character(SecStrucElement.STRAND).toString()+sheetId+strandSerial;
				SecStrucElement ssElem = new SecStrucElement(SecStrucElement.STRAND,beg,end,ssId);
				secondaryStructure.add(ssElem);
			}
			// we've stored the sec structure info in the strands2begEnd and sheets2strands maps.
			// the assignment to resser2secstruct is done when we reach the ATOM lines, see below
			// turn
			//TURN     1 S1A GLY A  16  GLN A  18     SURFACE
			//							turn ser				beg res ser					end res ser
			p = Pattern.compile("^TURN...(...).{9}"+chainCodeStr+"(....).{6}"+chainCodeStr+"(....)");
			m = p.matcher(line);
			if (m.find()){
				int serial = Integer.valueOf(m.group(1).trim());
				int beg = Integer.valueOf(m.group(2).trim());
				int end = Integer.valueOf(m.group(3).trim());
				String ssId = new Character(SecStrucElement.TURN).toString()+serial;
				SecStrucElement ssElem = new SecStrucElement(SecStrucElement.TURN,beg,end,ssId);
				secondaryStructure.add(ssElem);
			}			
			// MODEL
			// The model serial numbers should occur in columns 11-14 (official PDB format spec)
			// Here we are less strict: we allow for the numbers to appear in any column after the MODEL keyword (with any number of spaces in between)
			// We mainly need this because of CASP TS file which don't follow the PDB format 100% 
			p = Pattern.compile("^MODEL\\s+(\\d+)");
			m = p.matcher(line);
			if (m.find()){
				thismodel=Integer.parseInt(m.group(1));
			}
			if (thismodel!=model) continue; // we skip reading of atom lines if we are not in the desired model
			
			// PARENT (optional for Casp TS files)
			p = Pattern.compile("^PARENT\\s+.+");
			m = p.matcher(line);
			if(m.find()) {
				parentList = new LinkedList<String>();
				Pattern p2 = Pattern.compile("(\\d\\w\\w\\w_?\\w?)");
				m = p2.matcher(line);
				//System.out.printf("| ");
				while(m.find()) {
					parentList.add(parseParent(m.group()));
					//System.out.printf("%s ", parentList.getLast());
				}
				//System.out.printf("(%d) ", parentList.size());
				this.caspParents = new String[0];
				this.caspParents = parentList.toArray(this.caspParents);
			}
			// ATOM
			p = Pattern.compile("^ATOM");
			m = p.matcher(line);
			if (m.find()){
				try {
					//                                 serial    atom       altcode     restype      chain 	   res_ser icode    x     y     z
					Pattern pl = Pattern.compile("^.{6}(.....).{2}(...)"+ALTCODEREGEX+"(...).{1}"+chainCodeStr+"(.{4})(.).{3}(.{8})(.{8})(.{8})",
												Pattern.CASE_INSENSITIVE);
					Matcher ml = pl.matcher(line);
					if (ml.find()) {
						empty=false;
						int atomserial=Integer.parseInt(ml.group(1).trim());
						String atom = ml.group(2).trim();
						String res_type = ml.group(3).trim();
						int res_serial = Integer.parseInt(ml.group(4).trim());
						if (res_serial<1) {
							throw new PdbfileFormatError("A residue serial <=0 was found in the ATOM lines");
						}
						if (res_serial<lastResSerial) {
							throw new PdbfileFormatError("Residue serials do not occur in ascending order in ATOM lines of the PDB file "+pdbfile);
						}
						if(atomserial <= lastAtomSerial) {
							throw new PdbfileFormatError("Atom serials do not occur in ascending order in PDB file " + pdbfile + "(atom=" + atomserial + ")");
						}
						lastResSerial = res_serial;
						lastAtomSerial = atomserial;
						String iCode = ml.group(5);
						if (!iCode.equals(" ")) {
							throw new PdbfileFormatError("PDB file "+pdbfile+" contains insertion codes. Please use cif file instead.");
						}
						double x = Double.parseDouble(ml.group(6).trim());
						double y = Double.parseDouble(ml.group(7).trim());
						double z = Double.parseDouble(ml.group(8).trim());
						Point3d coords = new Point3d(x,y,z);
	
						if (isCaspTS && coords.equals(new Point3d(0.0,0.0,0.0))) {
							// in CASP TS (0,0,0) coordinates are considered unobserved (see http://predictioncenter.org/casp7/doc/casp7-format.html)
							if (!atomAtOriginSeen) {
								// first atom at origin we see is valid, we set the flag that we've seen it to true and later it will be read
								atomAtOriginSeen = true;
							} else {
								// more than 1 atom at origin: we don't want to read it: we skip it by continuing to next line
								continue;
							}
						} 
						
						if (AAinfo.isValidAA(res_type)) {
							if (!this.containsResidue(res_serial)) {
								this.addResidue(new Residue(AminoAcid.getByThreeLetterCode(res_type),res_serial,this));
							}
							if (AAinfo.isValidAtomWithOXT(res_type,atom)){
								Residue residue = this.getResidue(res_serial);
								residue.addAtom(new Atom(atomserial, atom, coords, residue));
							}
						}
	
					}
				} catch(NumberFormatException e) {
					throw new PdbfileFormatError("Wrong number format: " + e.getMessage());
				}
				
			}
		}
		fpdb.close();
		if (empty) {
			throw new PdbChainCodeNotFoundError("Couldn't find any ATOM line for given pdbChainCode: "+pdbChainCode+", model: "+model);
		}

		// we check also that resser2restype is not empty: happens when all residues in chain are non-standard
		if (this.get_length()==0) {
			throw new PdbfileFormatError("No standard aminoacids found for given chain in ATOM lines");
		}
		
		if (!hasSeqRes){ // no SEQRES could be read
			
			// 1) check whether first observed residue number is <100, if so we warn because it's likely to be an error
			int minResSer = this.getMinObsResSerial(); 
			if (minResSer>=100) {
				System.err.println("Warning! PDB file "+pdbfile+" starts with residue number "+minResSer+ ". "+
						"This is likely to be an error in the residue numbering of this PDB file. " +
						"A gap of size "+(minResSer-1)+" will be inserted at the beginning of the sequence.");
			}
			
			// 2) we take the sequence from the resser2restype Map, 
			// we assume the numbering is correct: we introduce non-observed gaps whenever we don't have the resser in resser2restype, 
			// that only misses non-observed gaps at end of chain which we can't know about
			sequence = "";
			for (int resser=1;resser<=getMaxObsResSerial();resser++) {
				if (this.containsResidue(resser)) {
					sequence += this.getResidue(resser).getAaType().getOneLetterCode();
				} else {
					sequence += AAinfo.UNKNOWN_UNOBSERVED_RES_ONE_LETTER;
				}
			}
			// 3) we set fullLength
			fullLength = sequence.length(); 
			
		} else { // we could read the sequence from SEQRES
			// 1) we check that the sequences from ATOM lines and SEQRES coincide (except for unobserved residues)
			for (int resser:getAllSortedResSerials()) {
				// before checking the sequence matching we need to be sure that the resser is not out of the range of the SEQRES sequence length 
				if (resser>sequence.length()) {
					throw new PdbfileFormatError("Residue serial "+resser+" from ATOM line is bigger than SEQRES sequence length. Incorrect residue numbering in PDB file "+pdbfile);
				}
				if (sequence.charAt(resser-1)!=this.getResidue(resser).getAaType().getOneLetterCode()) {
					throw new PdbfileFormatError("Sequences from ATOM lines and SEQRES do not match for position "+resser+". Incorrect residue numbering in PDB file "+pdbfile);
				}
			}
			// 2) we set fullLength
			fullLength = sequence.length();
		}
	}
	
	/**
	 * Reformats a parent record in one of the various styles found in Casp files to our standard form
	 * @param s the input parent string
	 * @return a parent string in our convention
	 */
	private String parseParent(String s) {
		String chain = "";
		String pdb = s.substring(0,4).toLowerCase();
		if(s.length() > 4 && !s.substring(s.length()-1).equals("_")) {
			chain = s.substring(s.length()-1).toUpperCase();
		}
		return pdb+chain;
		
	}
	
	/**
	 * Sets the sequence for this PdbasePdb object, overriding any current sequence information.
	 * If the observed residues (those having 3d coordinates) do not match the new sequence,
	 * an exception will be thrown.
	 * @param seq the new sequence
	 * @throws PdbLoadError if the given sequence does not match observed sequence from ATOM lines
	 */
	public void setSequence(String seq) throws PdbLoadError {
		// we check that the sequences from ATOM lines and the new sequence coincide (except for unobserved residues)
		for (int resser:getAllSortedResSerials()) {
			if (seq.charAt(resser-1)!=this.getResidue(resser).getAaType().getOneLetterCode()) {
				throw new PdbLoadError("Given sequence does not match observed sequence from ATOM lines for position "+resser+".");
			}
		}
		this.sequence = seq;
		fullLength = sequence.length();
	}
	
	/**
	 * Gets the PDB file name from which this Pdb is read from.
	 * @return
	 */
	public String getPdbFileName() {
		return pdbfile;
	}
	
}

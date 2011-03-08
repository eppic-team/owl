package owl.core.structure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Point3d;

import owl.core.sequence.alignment.PairwiseSequenceAlignment;
import owl.core.sequence.alignment.PairwiseSequenceAlignment.PairwiseSequenceAlignmentException;
import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;
import owl.core.util.FileFormatException;

/**
 * A single chain PDB protein structure read from a PDB file or CASP TS file
 * The sequences will be read from both SEQRES and ATOM lines and aligned if
 * they don't match. The internal residue numbers will then match the SEQRES.
 * If only ATOM lines are present the ATOM line sequence will be used as SEQRES. 
 * Exceptions will be thrown when:
 * - all residues are non-standard for a given chain (0 observed residues)
 */
public class PdbfilePdb extends Pdb {
	
	private static final long serialVersionUID = 1L;

	// the values for the parameters of the alignment have been fine tuned, checking that it works with as many PDB files as possible
	private static final float	GAP_OPEN_SCORE =	0.2f; // default 10f
	private static final float	GAP_EXTEND_SCORE =	0.1f; // default 0.5f
	private static final String ALI_SCORING_MATRIX = "IDENTITY"; //so that we force matching of identities only
	
	private static final String NULL_chainCode = "A";
	
	private String pdbfile;
	private boolean isCaspTS; // whether we are reading a CASP TS file (true) or a normal PDB file (false)
	private boolean hasSeqRes; // whether we find a non-empty SEQRES field in this pdb file
	
	private boolean hasAltCodes;
	
	private ArrayList<Residue> tmpResiduesList; // the temp list to store the residues read from ATOM lines (need to keep them in order as they appear in file)
	private HashMap<String,Residue> tmpResiduesMap; // the temp map of all seen pdb residue serials to residues (need it for fast searches)
	
	
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
	 * @throws PdbLoadException
	 */
	public void load(String pdbChainCode, int modelSerial) throws PdbLoadException {
		try {
			hasAltCodes = false;
			initialiseResidues();
			this.model=modelSerial;
			this.pdbChainCode=pdbChainCode;			// NOTE! pdb chain codes are case sensitive!
			// we set chainCode to pdbChainCode except for case Pdb.NULL_CHAIN_CODE where we use "A"
			this.chainCode=pdbChainCode;
			if (pdbChainCode.equals(Pdb.NULL_CHAIN_CODE)) this.chainCode=NULL_chainCode;

			parse();
			
			secondaryStructure.setSequence(this.sequence);
			if(!secondaryStructure.isEmpty()) {
				secondaryStructure.setComment("Author");
				this.initialiseResiduesSecStruct();
			}			

			this.initialiseMaps(); //TODO revise, is the remapping of pdbresser2resser in this necessary?
			dataLoaded = true;
			
			// so that the GC releases memory (hopefully) for the tmp residue objects
			tmpResiduesList = null;
			tmpResiduesMap = null;
			
		} catch (FileFormatException e) {
			throw new PdbLoadException(e);
		} catch (IOException e) {
			throw new PdbLoadException(e);
		} 
	}
	
	/**
	 * Returns all PDB chain codes present in the PDB file
	 * @return array with all pdb chain codes or null if no chains found
	 */
	public String[] getChains() throws PdbLoadException {
		TreeSet<String> chains = new TreeSet<String>();
		try {
			BufferedReader fpdb = new BufferedReader(new FileReader(new File(pdbfile)));
			String  line;
			while ((line=fpdb.readLine())!=null) {
				if (line.startsWith("ATOM")) {
					if (line.length()>22) {
						String chain = line.substring(21, 22);
						if (chain.equals(" ")) chain=Pdb.NULL_CHAIN_CODE;
						chains.add(chain);						 
					}
				}
			}
			fpdb.close();
		} catch (IOException e) {
			throw new PdbLoadException(e);
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
	public Integer[] getModels() throws PdbLoadException {
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
			throw new PdbLoadException(e);
		} catch (NumberFormatException e) {
			throw new PdbLoadException("Wrong format for MODEL lines!");
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
	 * lines. If the alignment given by the residue numbers does not match, then we realign
	 * and reassign internal residue numbers.
	 * @throws IOException
	 * @throws FileFormatException if file is empty, if file is a CASP TS file 
	 * and no TARGET line found
	 * @throws PdbLoadException if no ATOM lines are found for given 
	 * pdbChainCode and model or the space group found is not recognised
	 */
	private void parse() throws IOException, FileFormatException, PdbLoadException { 
		tmpResiduesList = new ArrayList<Residue>();
		tmpResiduesMap = new HashMap<String,Residue>();
		Pattern p;
		Matcher m;
		// we set chainCodeStr (for regex) to pdbChainCode except for case Pdb.NULL_CHAIN_CODE where we use " " (Pdb.NULL_CHAIN_CODE is a blank chain code in pdb files)
		String chainCodeStr=pdbChainCode;
		if (pdbChainCode.equals(Pdb.NULL_CHAIN_CODE)) chainCodeStr=" ";
		this.sequence = ""; // we will put here the sequence we find (either from SEQRES or ATOM lines)
		int lastAtomSerial = -1;
		int totalInsCodesFound = 0;
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
								fpdb.close();
								throw new FileFormatException("The CASP TS file "+pdbfile+" does not have a TARGET line");
							}
						} else {
							fpdb.close();
							throw new FileFormatException("The CASP TS file "+pdbfile+" is empty after the PFRMAT line");
						}
					}
				}
			}
			// EXPDTA
			if (line.startsWith("EXPDTA")) {
				 String exp = line.substring(6, line.length()).trim();
				 expMethod = exp.split(";\\s")[0]; // in some (strange) cases there are several exp methods, we simply take first, e.g. 2krl
			}
			// REMARK 3 (for resolution)
			if (line.startsWith("REMARK   3   RESOLUTION RANGE HIGH")){
				Pattern pR = Pattern.compile("^REMARK   3   RESOLUTION RANGE HIGH \\(ANGSTROMS\\) :\\s+(\\d+\\.\\d+).*");
				Matcher mR = pR.matcher(line);
				if (mR.matches()) {
					resolution = Double.parseDouble(mR.group(1));
				}
			}
			// REMARK 3 (for R free)
			if (line.startsWith("REMARK   3   FREE R VALUE")) {
				Pattern pR = Pattern.compile("^REMARK   3   FREE R VALUE\\s+(?:\\(NO CUTOFF\\))?\\s+:\\s+(\\d\\.\\d+).*");
				Matcher mR = pR.matcher(line);
				if (mR.matches()) {
					rFree = Double.parseDouble(mR.group(1));
				}				
			}
			// REMARK 200 (for R merge)
			if (line.startsWith("REMARK 200  R MERGE                    (I)")) {
				Pattern pR = Pattern.compile("^REMARK 200  R MERGE                    \\(I\\)\\s+:\\s+(\\d\\.\\d+).*");
				Matcher mR = pR.matcher(line);
				if (mR.matches()) {
					if (this.rSym == -1) {
						this.rSym = Double.parseDouble(mR.group(1));
					}
				}								
			}
			// REMARK 200 (for R sym)
			// if both Rsym/Rmerge are present, we don't compare them but take the Rsym value to be 
			// the right one (there's not much consensus in the field as to what's the 
			// right thing to do anyway!)			
			if (line.startsWith("REMARK 200  R SYM                      (I)")) {
				Pattern pR = Pattern.compile("^REMARK 200  R SYM                      \\(I\\)\\s+:\\s+(\\d\\.\\d+).*");
				Matcher mR = pR.matcher(line);
				if (mR.matches()) {
					this.rSym = Double.parseDouble(mR.group(1));
				}												
			}
			// CRYST1
			if (line.startsWith("CRYST1")) {
				double a = Double.parseDouble(line.substring( 6, 15));
				double b = Double.parseDouble(line.substring(15, 24));
				double c = Double.parseDouble(line.substring(24, 33));
				double alpha = Double.parseDouble(line.substring(33, 40));
				double beta = Double.parseDouble(line.substring(40, 47));
				double gamma = Double.parseDouble(line.substring(47, 54));
				int endsg = 66;
				if (line.length()<66) {
					endsg = line.length(); // some programs like phenix don't write the z number and thus can have no trailing spaces so that the line is shorter
				}
				String sg = line.substring(55,endsg).trim();
				crystalCell = new CrystalCell(a, b, c, alpha, beta, gamma);
				spaceGroup = SymoplibParser.getSpaceGroup(sg);
				if (spaceGroup==null) {
					throw new PdbLoadException("The space group found '"+sg+"' is not recognised as a standard space group");
				}
			}
			// SEQRES
			//SEQRES   1 A  348  VAL ASN ILE LYS THR ASN PRO PHE LYS ALA VAL SER PHE
			p = Pattern.compile("^SEQRES.{5}"+chainCodeStr);
			m = p.matcher(line);
			if (m.find()){
				for (int i=19;i<=67;i+=4) {
					// most pdb files have blank spaces up to 80 characters, but some don't.
					// because of that we need to check that (in the last line of SEQRES) the line is long enough
					// else we'd get an out of bounds error
					if (line.length()>=i+3) {  
						if (!line.substring(i, i+3).equals("   ")) {
							if (AminoAcid.isStandardAA(line.substring(i, i+3))) { // for non-standard aas
								sequence+= AminoAcid.three2one(line.substring(i, i+3));
							} else {
								sequence+= AminoAcid.XXX.getOneLetterCode();
							}
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
					if (line.length()<54) {
						// the least we admit is a PDB file with coordinates up to z
						fpdb.close();
						throw new FileFormatException("ATOM/HETATM line is too short to contain the minimum fields required. PDB file "+pdbfile+" at line "+linecount);
					}
					if (line.substring(21, 22).matches(chainCodeStr)) {
						int atomserial=Integer.parseInt(line.substring(6,11).trim());
						String atom = line.substring(12,16).trim();
						String res_type = line.substring(17,20).trim();
						String pdbResSerial = line.substring(22,27).trim();
						if(atomserial <= lastAtomSerial) {
							fpdb.close();
							throw new FileFormatException("Atom serials do not occur in ascending order in PDB file " + pdbfile + "(atom=" + atomserial + ")");
						}
						lastAtomSerial = atomserial;
						String altCode = line.substring(16, 17);
						if (!altCode.equals(" ")) hasAltCodes = true;
						double x = Double.parseDouble(line.substring(30,38).trim());
						double y = Double.parseDouble(line.substring(38,46).trim());
						double z = Double.parseDouble(line.substring(46,54).trim());
						Point3d coords = new Point3d(x,y,z);
						double occupancy = Atom.DEFAULT_OCCUPANCY;
						if (line.length()>=60)
							occupancy = Double.parseDouble(line.substring(54,60).trim());
						double bfactor = Atom.DEFAULT_B_FACTOR;
						if (line.length()>=66)
							bfactor = Double.parseDouble(line.substring(60,66).trim());
						String element = null;
						if (line.length()>=78)
							element = line.substring(76,78).trim();
	
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
						
						if (AminoAcid.isStandardAA(res_type)) {
							int resSerial = 0;
							if (Character.isDigit(pdbResSerial.charAt(pdbResSerial.length()-1))) {
								resSerial = Integer.parseInt(pdbResSerial);
							} else {
								resSerial = Integer.parseInt(pdbResSerial.substring(0,pdbResSerial.length()-1));
								if (!tmpResiduesMap.containsKey(pdbResSerial)) {
									totalInsCodesFound++; // if it is the first time seen we increase counter
									// the strategy of summing the totalInsCodesFound to the parsed number without ins code
									// only works when pdbSerial is same as last with ins code e.g. 27 -> 27A (27+1) -> 27B (27+2), 
									// doesn't work for 0A (0+1) -> 1B (1+2)
									// in the latter case we simply use the wrong numbering and rely on realigning later
								}
							}

							if (!tmpResiduesMap.containsKey(pdbResSerial)) {
								Residue residue = new Residue(AminoAcid.getByThreeLetterCode(res_type),resSerial+totalInsCodesFound,this); 
								tmpResiduesList.add(residue);
								residue.setPdbSerial(pdbResSerial);
								tmpResiduesMap.put(pdbResSerial,residue);
							}
							if (AminoAcid.isValidAtomWithOXT(res_type,atom)){
								Residue residue = tmpResiduesMap.get(pdbResSerial);
								// for alt codes we take either blanks or As
								// this is slightly different from what we do in CifFile or Pdbase where we take either blanks or first (alphabetically) letter
								// for some entries like 2imf there can be discrepancies (in that one there's no As but only B,C)
								if (altCode.equals(" ") || altCode.equals("A")) {
									residue.addAtom(new Atom(atomserial, atom, element, coords, residue, occupancy, bfactor));
								}
							}
						}

					}
				} catch(NumberFormatException e) {
					fpdb.close();
					throw new FileFormatException("Wrong number format in PDB file "+pdbfile+" at line "+linecount+". Error: " + e.getMessage());
				}
				
			}
		}
		fpdb.close();
		// we check that there was at least one observed residue for the chain
		if (tmpResiduesList.size()==0) {
			throw new PdbLoadException("Couldn't find any ATOM line for given pdbChainCode: "+pdbChainCode+", model: "+model);
		}
		
		checkForEmptyResidues();
		checkSeqResMatching();
	}
	
	/**
	 * We need to check for residues that are in the ATOM lines but for which we 
	 * did not find any valid atoms.
	 * This happens when there are alt codes in the file but not used for all residues 
	 * in the file, e.g. in 2heu several alt codes are present (A,B,C,D) but for some 
	 * residues there's no A or no D atoms. This is not standard practice.
	 * Our approach to alt codes is to take either blanks or the first one encountered in the file (usually 'A')
	 */
	private void checkForEmptyResidues() {
		ArrayList<Residue> emptyResidues = new ArrayList<Residue>();
		for (Residue residue:tmpResiduesList) {
			if (residue.getNumAtoms()==0) {
				emptyResidues.add(residue);
			}
		}
		for (Residue emptyRes:emptyResidues) {
			tmpResiduesList.remove(emptyRes);
			tmpResiduesMap.remove(emptyRes.getPdbSerial());
		}
			
	}
	
	private void checkSeqResMatching() throws PdbLoadException {

		if (!hasSeqRes){ // no SEQRES could be read
			
			// we take the sequence from the residues Map
			// and renumber
			sequence = "";
			int newResSer = 1;
			for (Residue residue:tmpResiduesList) {
				sequence += residue.getAaType().getOneLetterCode();
				residue.setSerial(newResSer);
				newResSer++;
			}
			
		} else { // we could read the sequence from SEQRES
			boolean aligned = true;
			if (tmpResiduesList.size()>sequence.length()) {
				throw new PdbLoadException("The sequence from ATOM lines is longer than the SEQRES sequence.");
			}
			// we check that the sequences from ATOM lines and SEQRES coincide (except for unobserved residues)
			for (Residue residue:tmpResiduesList) {
				int resser = residue.getSerial();
				if (resser<1 || resser>sequence.length()) {
					aligned = false;
					break;
				}
				if (residue.getAaType().getOneLetterCode()!=sequence.charAt(resser-1)){
					aligned = false;
					break;
				}				
			}
			if (!aligned) {
				if (!checkIfShifted()) {
					reAlignSeqRes();
				}
			} 
		}
		// now that we have realigned and renumbered we fill the final Pdb object
		for (Residue residue:tmpResiduesList){
			this.addResidue(residue);
		}
	
		// finally we initialise the pdb 2 resser maps
		this.resser2pdbresser = new TreeMap<Integer, String>();
		this.pdbresser2resser = new TreeMap<String, Integer>();
		for (int i=0;i<sequence.length();i++) {
			if (this.containsResidue(i+1)) {
				Residue residue = getResidue(i+1);
				resser2pdbresser.put(residue.getSerial(),residue.getPdbSerial());
				pdbresser2resser.put(residue.getPdbSerial(),residue.getSerial());
			} 
		}

	}
	
	private void reAlignSeqRes() throws PdbLoadException {
		String obsSequence = "";
		for (Residue residue:tmpResiduesList) {
			obsSequence += residue.getAaType().getOneLetterCode();
		}
		try {
			PairwiseSequenceAlignment psa = new PairwiseSequenceAlignment(sequence, obsSequence, "SEQRES", "ATOM",
					GAP_OPEN_SCORE,GAP_EXTEND_SCORE,ALI_SCORING_MATRIX);
			//psa.printAlignment();
			// 1st check: all positions of ATOM sequence are identical to a position in SEQRES
			if (psa.getIdentity()!=obsSequence.length()) {
				throw new PdbLoadException("The ATOM lines sequence does not align with identities for all its positions.\n"+psa.getAlignmentString());
			}
			// 2nd check: there's no gaps on the SEQRES side
			if (psa.getLength()!=psa.getAlignedSequences()[0].length()) {
				throw new PdbLoadException("Alignment of ATOM to SEQRES sequences contains gap on the SEQRES sequence.");
			}

			// now we can trust it's a full match: renumber
			for (int i=0;i<tmpResiduesList.size();i++) {
				int posInSeqRes = psa.getMapping2To1(i);
				Residue res = tmpResiduesList.get(i);
				res.setSerial(posInSeqRes+1);
			}


		} catch (PairwiseSequenceAlignmentException e) {
			throw new PdbLoadException("Could not create alignment of SEQRES and ATOM lines sequences to realign them");
		}
		
	}
	
	private boolean checkIfShifted() {
		// strategy is to take first window residues and see if they match to seqres, 
		// from that we know the offset and check the rest of sequence
		// If the whole thing match then we renumber according to the found offset
		// and fill the residues with addResidue
		int window = 5;
		String seqPattern = "";
		int i = 0;
		int lastResSerial = -1;
		while (tmpResiduesList.size()>i && i<window) { 	// if ATOM sequence smaller than window we have to take only the first 1, 2, 3, 4 depending on size
			Residue currentObsRes = tmpResiduesList.get(i);
			if (lastResSerial!=-1 && (currentObsRes.getSerial()-lastResSerial)>1) {
				// if there are gaps in the first window residues we consider them to be non-standard (hoping there are no inner gaps in the first window residues!)
				// if unobserved residues are present in the first window residues (e.g. 2pvb) then this doesn't find the right shifting and we have to rely on alignment
				for (int j=1;j<(currentObsRes.getSerial()-lastResSerial);j++) {
					seqPattern+=AminoAcid.XXX.getOneLetterCode();
				}
			} 
			seqPattern+=currentObsRes.getAaType().getOneLetterCode();
			lastResSerial = currentObsRes.getSerial();
			i++;
		}
		Pattern p = Pattern.compile(seqPattern);
		Matcher m = p.matcher(this.sequence);
		int offSet = -1;
		if (m.find()) {
			offSet = m.start();
		} else {
			return false; // couldn't find it, will continue with the alignment
		}
		int shift = tmpResiduesList.get(0).getSerial()-offSet;
		boolean aligned = true;
		for (Residue residue:tmpResiduesList) {
			int seqIndex = residue.getSerial()-shift;
			if (seqIndex<0 || seqIndex>=this.sequence.length()) {
				aligned = false;
				break;
			}
			if (residue.getAaType().getOneLetterCode()!=this.sequence.charAt(seqIndex)) {
				aligned = false;
				break;
			}
		}
		// if the above doesn't fail for any position, it means there was a shift, we renumber accordingly
		if (aligned) {
			for (Residue residue:tmpResiduesList) {
				residue.setSerial(residue.getSerial()-shift+1);
			}
		}
		
		return aligned;
		
		
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
	 * @throws PdbLoadException if the given sequence does not match observed sequence from ATOM lines
	 */
	public void setSequence(String seq) throws PdbLoadException {
		// we check that the sequences from ATOM lines and the new sequence coincide (except for unobserved residues)
		for (int resser:getAllSortedResSerials()) {
			if (seq.charAt(resser-1)!=this.getResidue(resser).getAaType().getOneLetterCode()) {
				throw new PdbLoadException("Given sequence does not match observed sequence from ATOM lines for position "+resser+".");
			}
		}
		this.sequence = seq;
	}
	
	/**
	 * Gets the PDB file name from which this Pdb is read from.
	 * @return
	 */
	public String getPdbFileName() {
		return pdbfile;
	}
	
	/**
	 * True if the PDB file contains alt codes for this chain, false otherwise.
	 * @return
	 */
	public boolean hasAltCodes() {
		return hasAltCodes;
	}
}

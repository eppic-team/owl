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
 * A PDB file format parser, reads PDB files or CASP TS files
 * The sequences will be read from both SEQRES and ATOM lines and aligned if
 * they don't match. The internal residue numbers will then match the SEQRES.
 * If only ATOM lines are present the ATOM line sequence will be used as SEQRES. 
 * Exceptions will be thrown when:
 * - a chain contains 0 observed residues
 */
public class PdbfileParser {
	
	private static boolean DEBUG = false;
	
	// the values for the parameters of the alignment have been fine tuned, checking that it works with as many PDB files as possible
	private static final float	GAP_OPEN_SCORE =	0.2f; // default 10f
	private static final float	GAP_EXTEND_SCORE =	0.1f; // default 0.5f
	private static final String ALI_SCORING_MATRIX = "IDENTITY"; //so that we force matching of identities only

	private static final Pattern INNER_GAPS_REGEX = Pattern.compile("\\w-+\\w");
	
	public static final String NULL_chainCode = "A";
	
	private String pdbfile;
	private boolean isCaspTS; // whether we are reading a CASP TS file (true) or a normal PDB file (false)
	private boolean hasSeqRes; // whether we find a non-empty SEQRES field in this pdb file
	
	private HashMap<String,ArrayList<Residue>> tmpResiduesLists; // a map of pdb chain codes to the temp lists to store the residues read from ATOM lines (need to keep them in order as they appear in file)
	private HashMap<String,HashMap<String,Residue>> tmpResiduesMaps; // a map of pdb chain codes to the temp maps of all seen pdb residue serials to residues (need it for fast searches)
	
	private String pdbCode;
	private String title;
	private String expMethod;
	private CrystalCell crystalCell;
	private SpaceGroup spaceGroup;
	private double resolution;
	private double rFree;
	private double rSym; 
	
	private int caspTargetNum;
	private String[] caspParents;
	
	private int model;
	
	private HashMap<String,String> sequences;
	
	private String[] chainsArray;
	private Integer[] modelsArray;
	
	private boolean missingSeqResPadding;
	
	/**
	 * Constructs a pdb file parser object given a pdbfile name
	 * Use readChain(pdbChainCode, modelSerial) to generate PdbChain objects. 
	 * @param pdbfile
	 */
	public PdbfileParser (String pdbfile) {
		this.pdbfile = pdbfile;
		this.isCaspTS = false; //we assume by default this is not a CASP TS file, when reading we set it to true if we detect CASP TS headers
		
		this.pdbCode = PdbAsymUnit.NO_PDB_CODE;
		this.resolution = -1;
		this.rFree = -1;
		this.rSym = -1;
		this.missingSeqResPadding = true;
	}

	/**
	 * Constructs a pdb file parser object given a pdbfile name
	 * Use readChain(pdbChainCode, modelSerial) to generate PdbChain objects. 
	 * @param pdbfile
	 * @param missingSeqResPadding if true in a PDB file without SEQRES, the sequence will be padded with X 
	 * for any missing residue (deduced from residue numbering) in ATOM line. If false no padding will be
	 */
	public PdbfileParser (String pdbfile, boolean missingSeqResPadding) {
		this.pdbfile = pdbfile;
		this.isCaspTS = false; //we assume by default this is not a CASP TS file, when reading we set it to true if we detect CASP TS headers
		
		this.pdbCode = PdbAsymUnit.NO_PDB_CODE;
		this.resolution = -1;
		this.rFree = -1;
		this.rSym = -1;
		this.missingSeqResPadding = missingSeqResPadding;

	}
	
	/**
	 * Reads PDB data (coordinates, sequence, etc.) from the PDB file
	 * for given pdbChainCode and modelSerial
	 * @param pdbAsymUnit
	 * @param modelSerial
	 * @throws PdbLoadException
	 */
	public void readChains(PdbAsymUnit pdbAsymUnit, int modelSerial) throws PdbLoadException {
		
		try {
			this.model=modelSerial;

			parse(pdbAsymUnit);
			
			for (PdbChain pdb:pdbAsymUnit.getAllChains()) {
				pdb.initialiseMaps();
			}
			
			// so that the GC releases memory (hopefully) for the tmp residue objects
			tmpResiduesLists = null;
			tmpResiduesMaps = null;
			
		} catch (FileFormatException e) {
			throw new PdbLoadException(e);
		} catch (IOException e) {
			throw new PdbLoadException(e);
		} 
//		return pdb;
	}
	
	/**
	 * Returns all alphabetically sorted PDB chain codes present in the PDB file by quickly parsing the chain 
	 * info from the ATOM lines of the PDB file.
	 * It caches the result so that next time called no parsing has to be done.
	 * @return array with all pdb chain codes or null if no chains found
	 */
	public String[] getChains() throws PdbLoadException {
		if (chainsArray!=null) {
			return chainsArray;
		}
		TreeSet<String> chains = new TreeSet<String>();
		try {
			BufferedReader fpdb = new BufferedReader(new FileReader(new File(pdbfile)));
			String  line;
			while ((line=fpdb.readLine())!=null) {
				if (line.startsWith("ATOM")) {
					if (line.length()>22) {
						String chain = line.substring(21, 22);
						if (chain.equals(" ")) chain=PdbAsymUnit.NULL_CHAIN_CODE;
						chains.add(chain);						 
					}
				}
			}
			fpdb.close();
		} catch (IOException e) {
			throw new PdbLoadException(e);
		}
		
		if (chains.isEmpty()) return null;
		
		chainsArray = new String[chains.size()];
		chains.toArray(chainsArray);
		return chainsArray;
	}
	
	/**
	 * Returns all model serials present in the PDB file by parsing the MODEL lines
	 * of the PDB file.
	 * It caches the result so that next time called no parsing has to be done. 
	 * @return array with all model serials
	 */
	public Integer[] getModels() throws PdbLoadException {
		if (modelsArray!=null) {
			return modelsArray;
		}
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
		
		if (models.isEmpty()) models.add(PdbAsymUnit.DEFAULT_MODEL);//return null;		
		modelsArray = new Integer[models.size()];
		models.toArray(modelsArray);
		return modelsArray;
	}
	
	protected String getPdbCode() {
		return pdbCode;
	}
	
	protected String getTitle(){
		return title;
	}
	
	protected String getExpMethod() {
		return expMethod;
	}
	
	protected CrystalCell getCrystalCell() {
		return crystalCell;
	}
	
	protected SpaceGroup getSpaceGroup() {
		return spaceGroup;
	}
	
	protected double getResolution() {
		return resolution;
	}
	
	protected double getRfree() {
		return rFree;
	}
	
	protected double getRsym() {
		return rSym;
	}
	
	/**
	 * To read the pdb data (atom coordinates, residue serials, atom serials) from file.
	 * <code>chainCode</code> gets set to same as <code>pdbChainCode</code>, except if input chain code 
	 * is PdbChain.NULL_CHAIN_CODE then chainCode will be <code>NULL_chainCode</code>
	 * <code>pdbCode</code> gets set to the one parsed in HEADER or to <code>NO_PDB_CODE</code> 
	 * if not found
	 * The sequence is either read from SEQRES if present or from the residues read from ATOM 
	 * lines. If the alignment given by the residue numbers does not match, then we realign
	 * and reassign internal residue numbers.
	 * @param pdbAsymUnit
	 * @throws IOException
	 * @throws FileFormatException if file is empty, if file is a CASP TS file 
	 * and no TARGET line found
	 * @throws PdbLoadException if no ATOM lines are found for given 
	 * pdbChainCode and model or the space group found is not recognised
	 */
	private void parse(PdbAsymUnit pdbAsymUnit) throws IOException, FileFormatException, PdbLoadException {
		AtomLineList atomLines = new AtomLineList();
		sequences = new HashMap<String,String>();
		//HashMap<String,SecondaryStructure> secStructures = new HashMap<String,SecondaryStructure>();
		ArrayList<SecStructureLine> secStructureLines = new ArrayList<SecStructureLine>();
		tmpResiduesLists = new HashMap<String, ArrayList<Residue>>();
		tmpResiduesMaps = new HashMap<String, HashMap<String,Residue>>();
		Pattern p;
		Matcher m;
		this.title = "";
		boolean outOfPolyChain = false; // true when out of poly chain (after TER record seen), false again when an ATOM record found
		boolean atomAtOriginSeen = false; // if we've read at least 1 atom at the origin (0,0,0) it is set to true
		boolean terRecordSeen = false;
		int thismodel=PdbAsymUnit.DEFAULT_MODEL; // we initialise to DEFAULT_MODEL, in case file doesn't have MODEL lines 
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
						pdbCode = PdbAsymUnit.NO_PDB_CODE;
						// we try to read the TARGET from the next line, if there's no TARGET line appearing this is not respecting the format: exception
						if((line = fpdb.readLine()) != null ) {
							linecount++;
							p = Pattern.compile("^TARGET\\s+[Tt](\\d+)");
							m = p.matcher(line);
							if (m.find()) {
								caspTargetNum = Integer.parseInt(m.group(1));
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
			// TITLE
			if (line.startsWith("TITLE")) {
				if (!line.substring(8,10).equals("  ")) {
					char lastChar = title.charAt(title.length()-1);
					if (lastChar!='-') { 
						this.title += " ";
					}
				}
				this.title += line.substring(10,line.length()).trim();
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
			if (line.startsWith("SEQRES")){
				String chain = line.substring(11,12);
				if (!sequences.containsKey(chain)) {
					sequences.put(chain,"");
				}
				String sequence = sequences.get(chain);
				for (int i=19;i<=67;i+=4) {
					// most pdb files have blank spaces up to 80 characters, but some don't.
					// because of that we need to check that (in the last line of SEQRES) the line is long enough
					// else we'd get an out of bounds error
					if (line.length()>=i+3) {  
						if (!line.substring(i, i+3).equals("   ")) {
							if (AminoAcid.isStandardAA(line.substring(i, i+3))) { // for non-standard aas
								sequence+= AminoAcid.three2one(line.substring(i, i+3));
							} else if (Nucleotide.isStandardNuc(line.substring(i, i+3).trim())) {
								sequence+=Nucleotide.getByCode(line.substring(i, i+3).trim()).getOneLetterCode();
							} else {
								sequence+=AminoAcid.XXX.getOneLetterCode();
							}
						}
					}
				}
				sequences.put(chain,sequence);
				// if SEQRES was not empty then we have sequences
				if (!sequences.isEmpty()) {
					hasSeqRes = true;
				}
			}
			// SECONDARY STRUCTURE
			// helix
			//HELIX    1   1 LYS A   17  LEU A   26  1
			if (line.startsWith("HELIX")){
				String begChain = line.substring(19,20);
				String endChain = line.substring(31,32);
				int serial = Integer.valueOf(line.substring(7,10).trim());
				String beg = line.substring(21,26).trim();
				String end = line.substring(33,38).trim();
				secStructureLines.add(new SecStructureLine("HELIX", begChain, endChain, beg, end, serial, ""));
			}
			// sheet
			//SHEET    2   A 5 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107
			if (line.startsWith("SHEET")){
				String begChain = line.substring(21,22);
				String endChain = line.substring(32,33);
				int strandSerial = Integer.valueOf(line.substring(7,10).trim());
				String sheetId = line.substring(11,14).trim();
				String beg = line.substring(22,27).trim();
				String end = line.substring(33,38).trim();
				secStructureLines.add(new SecStructureLine("SHEET", begChain, endChain, beg, end, strandSerial, sheetId));
			}
			
			//TURN     1 S1A GLY A  16  GLN A  18     SURFACE
			//							turn ser				beg res ser					end res ser
			// turn has been deprecated (see PDB file format ver 3.20)
			
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
				caspParents = new String[0];
				caspParents = parentList.toArray(caspParents);
			}
			if (line.startsWith("TER ")) {
				outOfPolyChain = true;
				terRecordSeen = true;
			}
			// ATOM
			p = Pattern.compile("^(ATOM|HETATM)");
			m = p.matcher(line);
			if (m.find()){
				try {
					if (line.length()<54) {
						// the least we admit is a PDB file with coordinates up to z
						fpdb.close();
						throw new FileFormatException("ATOM/HETATM line is too short to contain the minimum fields required. PDB file "+pdbfile+" at line "+linecount);
					}
					if (!line.substring(17,20).trim().equals("HOH")) {
						if (m.group(1).equals("ATOM")) outOfPolyChain = false;
						int atomserial=Integer.parseInt(line.substring(6,11).trim());
						String atom = line.substring(12,16).trim();
						String res_type = line.substring(17,20).trim();
						String pdbChainCode = line.substring(21, 22);
						String pdbResSerialField = line.substring(22,27).trim();
						int pdbResSerial = 0;
						String insCode = ".";
						if (Character.isDigit(pdbResSerialField.charAt(pdbResSerialField.length()-1))) {
							pdbResSerial = Integer.parseInt(pdbResSerialField);
						} else {
							pdbResSerial = Integer.parseInt(pdbResSerialField.substring(0,pdbResSerialField.length()-1));
							insCode = pdbResSerialField.substring(pdbResSerialField.length()-1);
						}
						
						String altCode = line.substring(16, 17);
						if (altCode.equals(" ")) altCode=".";
						double x = Double.parseDouble(line.substring(30,38).trim());
						double y = Double.parseDouble(line.substring(38,46).trim());
						double z = Double.parseDouble(line.substring(46,54).trim());
						Point3d coords = new Point3d(x,y,z);
						double occupancy = Atom.DEFAULT_OCCUPANCY;
						if (line.length()>=60) {
							String occupancyStr = line.substring(54,60).trim();
							if (!occupancyStr.equals("")) {
								occupancy = Double.parseDouble(occupancyStr);
							}
						}
						double bfactor = Atom.DEFAULT_B_FACTOR;
						if (line.length()>=66) {
							String bfactorStr = line.substring(60,66).trim();
							if (!bfactorStr.equals("")) {
								bfactor = Double.parseDouble(bfactorStr);
							}
						}
						String element = null;
						if (line.length()>=78) {
							element = line.substring(76,78).trim();
							if (element.equals("") || Character.isDigit(element.charAt(0)) || (element.length()==2 && Character.isDigit(element.charAt(1)))) 
								element = null;
						}
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
						
						atomLines.addAtomLine(
							new AtomLine(null, altCode, atomserial, atom, element, res_type, 0, pdbResSerial, insCode, coords, occupancy, bfactor, 
									pdbChainCode, outOfPolyChain, m.group(1).equals("HETATM")));
					}
				} catch(NumberFormatException e) {
					fpdb.close();
					throw new FileFormatException("Wrong number format in PDB file "+pdbfile+" at line "+linecount+". Error: " + e.getMessage());
				}
				
			}
		}
		fpdb.close();
		// we check that there was at least one observed residue
		if (atomLines.isEmpty()) {
			throw new PdbLoadException("Couldn't find any ATOM/HETATM line for model: "+model);
		}
		
		//if (!terRecordSeen) System.err.println("Warning: TER records are missing in PDB file, chain assignments can be unreliable");
		
		atomLines.sortIntoChains(terRecordSeen); // this assigns the labelAsymIds (chainCodes) and the isNonPoly fields and finds out the pdbchaincode2chaincode mapping
		pdbAsymUnit.setPdbchaincode2chaincode(atomLines.getPdbChainCode2chainCode());
		String altLoc = atomLines.getAtomAltLoc();

		for (ArrayList<AtomLine> group:atomLines.getAtomLineGroups().values()) {
			AtomLine firstAtomLine = group.get(0);
			
			PdbChain pdb = new PdbChain();
			pdb.setChainCode(firstAtomLine.labelAsymId);
			pdb.setPdbChainCode(firstAtomLine.authAsymId);
			pdb.setParent(pdbAsymUnit);

			if (firstAtomLine.isNonPoly) {
				pdb.setIsNonPolyChain(true);
				pdbAsymUnit.setNonPolyChain(firstAtomLine.labelAsymId,pdb);
			} else {
				pdb.setIsNonPolyChain(false);
				pdbAsymUnit.setPolyChain(firstAtomLine.labelAsymId,pdb);
				tmpResiduesLists.put(firstAtomLine.authAsymId, new ArrayList<Residue>());
				tmpResiduesMaps.put(firstAtomLine.authAsymId, new HashMap<String,Residue>());
			}
			
			for (AtomLine atomLine:group) {	
				
				// we read only the alt locs we want
				if (altLoc!=null && !atomLine.labelAltId.equals(altLoc) && !atomLine.labelAltId.equals(".")) continue;

				// creating the residues
				if (atomLine.isNonPoly==false) {
					// for polymer chains we don't add them yet to the chains but to the temp lists/maps
					HashMap<String,Residue> tmpResiduesMap = tmpResiduesMaps.get(atomLine.authAsymId);
					ArrayList<Residue> tmpResiduesList = tmpResiduesLists.get(atomLine.authAsymId);

					if (!tmpResiduesMap.containsKey(atomLine.getPdbResSerialWithInsCode())) {
						Residue residue = null;
						int resSerial = atomLine.pdbResSerial+atomLines.getNumInsCodeForChain(atomLine.authAsymId);
						if (AminoAcid.isStandardAA(atomLine.res_type)) {
							residue = new AaResidue(AminoAcid.getByThreeLetterCode(atomLine.res_type), resSerial, pdbAsymUnit.getChain(atomLine.authAsymId));
						} else if (Nucleotide.isStandardNuc(atomLine.res_type)) {
							residue = new NucResidue(Nucleotide.getByCode(atomLine.res_type),resSerial,pdbAsymUnit.getChain(atomLine.authAsymId));
						} else {
							residue = new HetResidue(atomLine.res_type, resSerial, pdbAsymUnit.getChain(atomLine.authAsymId));
						}
						tmpResiduesList.add(residue);
						residue.setPdbSerial(atomLine.getPdbResSerialWithInsCode());
						Residue oldVal = tmpResiduesMap.put(atomLine.getPdbResSerialWithInsCode(),residue);
						if (oldVal!=null) throw new PdbLoadException("Duplicate residue number for residue "+oldVal);

					}
					Residue residue = tmpResiduesMap.get(atomLine.getPdbResSerialWithInsCode());
					residue.addAtom(new Atom(atomLine.atomserial, atomLine.atom, atomLine.element, atomLine.coords, residue, atomLine.occupancy, atomLine.bfactor));


				} else {
					// for non-poly chains we add already residues to the chains (there's no realignment to do) 
					pdb = pdbAsymUnit.getChainForChainCode(atomLine.labelAsymId);
					if (!pdb.containsResidue(atomLine.pdbResSerial)) {
						Residue residue = null;
						if (AminoAcid.isStandardAA(atomLine.res_type)) {
							residue = new AaResidue(AminoAcid.getByThreeLetterCode(atomLine.res_type), atomLine.pdbResSerial, pdb);
						} else if (Nucleotide.isStandardNuc(atomLine.res_type)) {
							residue = new NucResidue(Nucleotide.getByCode(atomLine.res_type),atomLine.pdbResSerial,pdb);
						} else { 
							residue = new HetResidue(atomLine.res_type,atomLine.pdbResSerial,pdb);
						}
						residue.setPdbSerial(String.valueOf(atomLine.pdbResSerial));
						pdb.addResidue(residue);
					}
					Residue residue = pdb.getResidue(atomLine.pdbResSerial);
					residue.addAtom(new Atom(atomLine.atomserial, atomLine.atom, atomLine.element, atomLine.coords, residue, atomLine.occupancy, atomLine.bfactor));
				}
			}
		}

		if (altLoc!=null) {
			for (PdbChain pdb:pdbAsymUnit.getAllChains()) {
				pdb.setHasAltCodes(true);
			}
		}
		
		for (String pdbChainCode:pdbAsymUnit.getPdbChainCodes()){
			ArrayList<Residue> tmpResiduesList = tmpResiduesLists.get(pdbChainCode);
			HashMap<String,Residue> tmpResiduesMap = tmpResiduesMaps.get(pdbChainCode);

			checkForEmptyResidues(tmpResiduesList,tmpResiduesMap);
			checkSeqResMatching(tmpResiduesList,tmpResiduesMap,pdbAsymUnit,pdbChainCode);

			PdbChain pdb = pdbAsymUnit.getChain(pdbChainCode);
			
			boolean protein = false;
			boolean nucleotide = false;
			
			for (Residue residue:tmpResiduesList){
				if (residue instanceof AaResidue) protein = true;
				if (residue instanceof NucResidue) nucleotide = true;
				pdb.addResidue(residue);
			}
			
			// now that we have realigned and renumbered we fill the final PdbChain object
			if (protein && nucleotide) {
				throw new PdbLoadException("Mix of protein and nucleotide sequences in chain");
			}
			if (!protein && !nucleotide) {
				// in case where not a single std aa or nucleotide are found then we have a chain of Xs, we assume a protein
				protein = true;
			}
			
			pdb.setSequence(sequences.get(pdbChainCode), protein);
			
			
			// finally we initialise the pdb 2 resser maps
			TreeMap<Integer,String> resser2pdbresser = new TreeMap<Integer, String>();
			TreeMap<String,Integer> pdbresser2resser = new TreeMap<String, Integer>();
			for (int i=0;i<sequences.get(pdbChainCode).length();i++) {
				if (pdb.containsResidue(i+1)) {
					Residue residue = pdb.getResidue(i+1);
					resser2pdbresser.put(residue.getSerial(),residue.getPdbSerial());
					pdbresser2resser.put(residue.getPdbSerial(),residue.getSerial());
				} 
			}
			pdb.setResser2pdbresserMap(resser2pdbresser);
			pdb.setPdbresser2resserMap(pdbresser2resser);

			// we set casp members
			pdb.setCaspParents(caspParents);
			pdb.setTargetNum(caspTargetNum);
			
			// we initialise the secondary structure elements to blank and then we assign them if ss lines were read
			SecondaryStructure secondaryStructure = new SecondaryStructure("");
			if (hasSeqRes) {
				secondaryStructure.setSequence(sequences.get(pdbChainCode));
			}
			pdb.setSecondaryStructure(secondaryStructure);
		}
		
		// finally we assign the secondary structure from the parsed sec structure lines
		int generatedSerial = 1; // we have to generate our own serials because some times those in PDB files are wrong (whilst they have been fixed in CIF!)
		String lastType = null;
		String lastId = null;
		for (SecStructureLine ssline:secStructureLines) {
			if (!ssline.begChain.equals(ssline.endChain)) 
				throw new PdbLoadException(ssline.type+" element beg and end chain id differ for ss element with serial "+ssline.serial);
			SecondaryStructure secStructure = pdbAsymUnit.getChain(ssline.begChain).getSecondaryStructure();
			String ssId = "";
			char secStructElem = 0;
			int beg = pdbAsymUnit.getChain(ssline.begChain).getResSerFromPdbResSer(String.valueOf(ssline.begPdbChainCode));
			int end = pdbAsymUnit.getChain(ssline.begChain).getResSerFromPdbResSer(String.valueOf(ssline.endPdbChainCode));
			if ((lastType!=null && !lastType.equals(ssline.type)) || (lastId!=null && !lastId.equals(ssline.id))) {
				generatedSerial = 1;
			}
			if (ssline.type.equals("HELIX")) {
				ssId = new Character(SecStrucElement.HELIX).toString()+generatedSerial;//ssline.serial;
				secStructElem = SecStrucElement.HELIX;
			} else if (ssline.type.equals("SHEET")) {
				ssId = new Character(SecStrucElement.STRAND).toString()+ssline.id+generatedSerial;//ssline.serial;
				secStructElem = SecStrucElement.STRAND;
			}
			secStructure.add(new SecStrucElement(secStructElem,beg,end,ssId));
			generatedSerial++;
			lastType = ssline.type;
			lastId = ssline.id;
		}
		for (PdbChain pdb:pdbAsymUnit.getPolyChains()) {
			SecondaryStructure secStructure = pdb.getSecondaryStructure();
			if (!secStructure.isEmpty()) {
				secStructure.setComment("Author");
			}
			pdb.initialiseResiduesSecStruct();
		}

	

	}
	
	/**
	 * We need to check for residues that are in the ATOM lines but for which we 
	 * did not find any valid atoms.
	 * This happens when there are alt codes in the file but not used for all residues 
	 * in the file, e.g. in 2heu several alt codes are present (A,B,C,D) but for some 
	 * residues there's no A or no D atoms. This is not standard practice.
	 * Our approach to alt codes is to take either blanks or 'A'
	 */
	private void checkForEmptyResidues(ArrayList<Residue> tmpResiduesList, HashMap<String,Residue> tmpResiduesMap) {
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
	
	private void checkSeqResMatching(ArrayList<Residue> tmpResiduesList, HashMap<String,Residue> tmpResiduesMap, PdbAsymUnit pdbAsymUnit, String pdbChainCode) throws PdbLoadException {

		if (!hasSeqRes){ // no SEQRES could be read
			boolean canUseResidueNumberingAsIs = true;
			int lastResSerial = 0;
			for (Residue residue:tmpResiduesList) {
				// negative or 0 residue numbers
				if (residue.getSerial()<=0) {
					canUseResidueNumberingAsIs = false;
					break;
				}
				// insertion codes
				if (!Character.isDigit(residue.getPdbSerial().charAt(residue.getPdbSerial().length()-1))) {
					canUseResidueNumberingAsIs = false;
					break;
				}
				// residue serials in ascending order
				if (residue.getSerial()<lastResSerial) {
					canUseResidueNumberingAsIs = false;
					break;					
				}
				lastResSerial = residue.getSerial();
			}
			
			if (canUseResidueNumberingAsIs && missingSeqResPadding) { // we only do padding if padding mode was specified
				// we take the residue serials to be a valid numbering and fill the unknown gaps with Xs
				String sequence = "";
				for (int resser=1;resser<=tmpResiduesList.get(tmpResiduesList.size()-1).getSerial();resser++) {
					if (tmpResiduesMap.containsKey(String.valueOf(resser))) {
						sequence += tmpResiduesMap.get(String.valueOf(resser)).getShortCode();
					} else {
						sequence += AminoAcid.XXX.getOneLetterCode();
					}
				}
				sequences.put(pdbChainCode, sequence);
			} else {
				// we take the sequence from the residues Map
				// and renumber
				String sequence = "";
				int newResSer = 1;
				for (Residue residue:tmpResiduesList) {
					sequence += residue.getShortCode();
					residue.setSerial(newResSer);
					newResSer++;
				}
				sequences.put(pdbChainCode, sequence);
			}
			
		} else { // we could read the sequence from SEQRES
			boolean aligned = true;
			
			if (tmpResiduesList.size()>sequences.get(pdbChainCode).length()) {
				throw new PdbLoadException("The sequence from ATOM lines is longer than the SEQRES sequence.");
			}
			// we check that the sequences from ATOM lines and SEQRES coincide (except for unobserved residues)
			for (Residue residue:tmpResiduesList) {
				int resser = residue.getSerial();
				if (resser<1 || resser>sequences.get(pdbChainCode).length()) {
					aligned = false;
					break;
				}
				if (residue.getShortCode()!=sequences.get(pdbChainCode).charAt(resser-1)){
					aligned = false;
					break;
				}				
			}
			if (!aligned) {
				if (!checkIfShifted(tmpResiduesList,sequences.get(pdbChainCode))) {
					reAlignSeqRes(tmpResiduesList,sequences.get(pdbChainCode));
				}
			} 
		}
	}
	
	private void reAlignSeqRes(ArrayList<Residue> tmpResiduesList, String sequence) throws PdbLoadException {
		String obsSequence = "";
		for (Residue residue:tmpResiduesList) {
			obsSequence += residue.getShortCode();
		}
		try {
			PairwiseSequenceAlignment psa = new PairwiseSequenceAlignment(sequence, obsSequence, "SEQRES", "ATOM",
					GAP_OPEN_SCORE,GAP_EXTEND_SCORE,ALI_SCORING_MATRIX);
			if (DEBUG) {
				System.out.println("Realigning ATOM to SEQRES");
				psa.printAlignment();
			}
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

			// the 3D contiguity check we can only do for amino acids and not for nucleotides
			if (isProtein(tmpResiduesList)) {
				check3Dcontiguity(psa,tmpResiduesList,sequence);
			}
			
		} catch (PairwiseSequenceAlignmentException e) {
			throw new PdbLoadException("Could not create alignment of SEQRES and ATOM lines sequences to realign them");
		}
		

	}
	
	private boolean isProtein(ArrayList<Residue> tmpResiduesList) {
		for (Residue res:tmpResiduesList){
			if (res instanceof NucResidue) return false;
		}
		return true;
	}
	
	private boolean checkIfShifted(ArrayList<Residue> tmpResiduesList, String sequence) {
		if (DEBUG) System.out.println("Checking ATOM residue numbering is shifted");
		// strategy is to take first window residues and see if they match to seqres, 
		// from that we know the offset and check the rest of sequence
		// If the whole thing match then we renumber according to the found offset
		// and fill the residues with addResidue
		int window = 5;
		String seqPattern = "";
		int i = 0;
		while (tmpResiduesList.size()>i && i<window) { 	// if ATOM sequence smaller than window we have to take only the first 1, 2, 3, 4 depending on size
			Residue currentObsRes = tmpResiduesList.get(i);
			seqPattern+=currentObsRes.getShortCode();
			i++;
		}
		Pattern p = Pattern.compile(seqPattern);
		Matcher m = p.matcher(sequence);
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
			if (seqIndex<0 || seqIndex>=sequence.length()) {
				aligned = false;
				break;
			}
			if (residue.getShortCode()!=sequence.charAt(seqIndex)) {
				aligned = false;
				break;
			}
		}
		// if the above doesn't fail for any position, it means there was a shift, we renumber accordingly
		if (aligned) {			
			if (DEBUG) System.out.println("Shift of "+shift+" residues found, renumbering");
			for (Residue residue:tmpResiduesList) {
				residue.setSerial(residue.getSerial()-shift+1);
			}
		}
		
		return aligned;
		
		
	}
	
	private void check3Dcontiguity(PairwiseSequenceAlignment psa,ArrayList<Residue> tmpResiduesList, String sequence) {
		if (DEBUG) System.out.println("Checking 3D contiguity");
		// There are still ambiguities even if the alignment matches fully, e.g. 2ofz:
		//SEQRES             1 MGSDKIHHHHHHNTASWFTALTQHGKEELRFPRGQGVPINTNSGPDDQIG     50
        //                       |||||    |||||||||||||||||||||||||||||||||||||||
        //ATOM               1 --SDKIH----HNTASWFTALTQHGKEELRFPRGQGVPINTNSGPDDQIG     44
		// the H at both sides of the gap could be swapped to the other side (it's an alignment with the same score)
		// in these cases we want to check for 3D contiguity to determine if the H in left belongs to right or
		// the H in right belongs to left
		// Letters can be different e.g.:
		//SEQRES             1 MGSDKISGHHHGNTASWFTALTQHGKEELRFPRGQGVPINTNSGPDDQIG     50
        //                       |||||    |||||||||||||||||||||||||||||||||||||||
        //ATOM               1 --SDKIS----GNTASWFTALTQHGKEELRFPRGQGVPINTNSGPDDQIG     44
		// in here the only possibility is that the right G can be swapped to the left
		
		// from our tests this method fixes entries: 2ofz and 1dki (without it one residue doesn't match to CIF)
		
		//psa.printAlignment();
		String alignedAtomSeq = psa.getAlignedSequences()[1];
		Matcher m = INNER_GAPS_REGEX.matcher(alignedAtomSeq);
		while (m.find()) { // for each gap
			//String gap = m.group(1);
			
			int leftSeqresIndex = m.start();
			int rightSeqresIndex = m.end()-1;
			
			int leftAtomIndex = psa.getMapping1To2(leftSeqresIndex);
			int rightAtomIndex = psa.getMapping1To2(rightSeqresIndex);

			//System.err.println("gap: "+gap);
			//System.err.println("aln: "+sequence.substring(m.start(),m.end()));
			
			Residue leftResidue = (AaResidue)tmpResiduesList.get(leftAtomIndex);
			Residue rightResidue = (AaResidue)tmpResiduesList.get(rightAtomIndex);
			if (!(leftResidue instanceof AaResidue) || !(rightResidue instanceof AaResidue)) {
				// one of them is not an amino acid, there's not much we can do
				return;
			}

			AaResidue leftRes = (AaResidue) leftResidue;
			AaResidue rightRes = (AaResidue) rightResidue;
			
			char letterLeftSwapRes = sequence.charAt(leftSeqresIndex+1);
			char letterRightSwapRes = sequence.charAt(rightSeqresIndex-1);
			
			// we only want to check if there is a possibility of swapping the residues flanking the gap (one of the 2 cases above)
			if (leftRes.getShortCode()==letterRightSwapRes || rightRes.getShortCode()==letterLeftSwapRes) {
				if (leftRes.isContiguous(rightRes)) { // if not there's nothing to do, they are correctly aligned
					if (DEBUG) System.out.println("Possible 3D discontiguity found between residues "+leftRes+" and "+rightRes+". Will try to fix.");
					AaResidue leftResMin1 = null;
					if (leftAtomIndex-1>0) {
						leftResMin1 = (AaResidue)tmpResiduesList.get(leftAtomIndex-1);
					}					
					AaResidue rightResPlus1 = null;
					if (rightAtomIndex+1<tmpResiduesList.size()) {
						rightResPlus1 = (AaResidue)tmpResiduesList.get(rightAtomIndex+1);
					}
					if (leftResMin1!=null && leftResMin1.isContiguous(leftRes)) {
						//System.err.println("right to left");
						if (rightRes.getAaType().getOneLetterCode()==letterLeftSwapRes) {
							//System.err.println("letters match! it's a right to left");
							rightRes.setSerial(leftSeqresIndex+1+1); 
						}
					}
					if (rightResPlus1!=null && rightRes.isContiguous(rightResPlus1)) {
						//System.err.println("left to right");
						if (leftRes.getAaType().getOneLetterCode()==letterRightSwapRes) {
							//System.err.println("letters match! it's a left to right");
							leftRes.setSerial(rightSeqresIndex-1+1);
						}

					}

				}

			}
			
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
	 * Gets the PDB file name from which this PdbChain is read from.
	 * @return
	 */
	public String getPdbFileName() {
		return pdbfile;
	}
	
}

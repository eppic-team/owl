package owl.core.structure;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.ArrayList;
import java.util.TreeSet;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import owl.core.features.Feature;
import owl.core.features.FeatureType;
import owl.core.features.HasFeatures;
import owl.core.features.InvalidFeatureCoordinatesException;
import owl.core.features.OverlappingFeatureException;
import owl.core.sequence.alignment.MultipleSequenceAlignment;
import owl.core.structure.features.CatalSiteSet;
import owl.core.structure.features.EC;
import owl.core.structure.features.Scop;
import owl.core.structure.features.ScopRegion;
import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;
import owl.core.structure.graphs.AIGEdge;
import owl.core.structure.graphs.AIGNode;
import owl.core.structure.graphs.AIGraph;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.Box;
import owl.core.util.Interval;
import owl.core.util.IntervalSet;
import owl.core.util.MySQLConnection;

import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

import Jama.Matrix;
import Jama.SingularValueDecomposition;


import java.sql.SQLException;
import java.sql.Statement;

/**
 * A single chain PDB protein structure
 * 
 */
public class Pdb implements HasFeatures {

	/*------------------------------------  constants ---------------------------------------------*/
	public static final int    DEFAULT_MODEL   = 1;			// default model serial (NMR structures)
	public static final String DEFAULT_CHAIN   = "A";		// used when constructing models, this will be the chain assigned to them
	public static final String NULL_CHAIN_CODE = "NULL";	// to specify the NULL (blank in pdb file) chain code. 
															// Should be safe now to change the value of this constant from "NULL" to something else,
															// all hard coded "NULL" strings have been checked now (Jose svn rev. 609)
	public static final String NO_PDB_CODE       = "";		// to specify no pdb code
	public static final String NO_PDB_CHAIN_CODE = "";		// to specify no pdb chain code
	public static final String NO_CHAIN_CODE     = "";		// to specify no internal chain code
	public static final String DEFAULT_CASP_TS_CHAINCODE = " "; // Casp TS format allows only empty chain codes
	
	private static final double MIN_AVRG_NUM_ATOMS_RES = 6.5;	// the cutoff to consider that the average number of atoms per residue 
																// corresponds to that of an all atoms protein. See isAllAtom()
		
	/*-------------------------------------  members ---------------------------------------------*/

	// atom/residue data
	private TreeMap<Integer, Residue> residues;			// residue serials to residue object references (only observed residues)
	private TreeMap<Integer, Atom>    atomser2atom;		// atom serials to Atom object references, we keep this is as a separate map to speed up searches
	protected TreeMap<Integer,String> resser2pdbresser; // internal residue serials to pdb (author) residue serials (can include insertion codes so they are strings)
	protected TreeMap<String,Integer> pdbresser2resser; // pdb (author) residue serials (can include insertion codes so they are strings) to internal residue serials
	protected int fullLength; 							// length of full sequence as it appears in SEQRES field 
														// to get observed length use method getObsLength()
	
	// sequence features (annotations)
	protected SecondaryStructure secondaryStructure;	// the secondary structure annotation for this pdb object (should never be null)
	protected Scop scop;								// the scop annotation for this pdb object
	protected EC ec;									// the ec annotation for this pdb object
	protected CatalSiteSet catalSiteSet;				// the catalytic site annotation for this pdb object
	
	private Map<FeatureType, Collection<Feature>> features; // all other features. Eventually all features above should be implemented by using  
															// this object and the HasFeature interface
	
	// identifiers
	protected String sequence; 			// full sequence as it appears in SEQRES field
	protected String pdbCode;			// the 4 letters PDB code. By convention we always use lower case. 
	protected String pdbChainCode;		// Given "external" pdb chain code, i.e. the classic (author's) pdb code 
										// If it is blank in original PDB file then it is: Pdb.NULL_CHAIN_CODE
	protected String chainCode;			// Our internal chain identifier:
										// - in reading from pdbase/cif file it will be set to the cif chain id (asym_id field)
										// - in reading from PDB file it coincides with pdbChainCode except for 
										//   Pdb.NULL_CHAIN_CODE where we use "A"
	protected int model;  				// the model serial for NMR structures
	protected String title;				// the title of the structure (e.g. from the PDB)
	protected String sid;				// the scop id if Pdb has been restricted (restrictToScopDomain)
	
	// optional fields for structures based on casp predictions
	protected int targetNum;
	protected int caspModelNum;
	protected String caspAuthorStr;
	protected String caspMethodStr;
	protected int groupNum;
	protected String[] caspParents;		// optional list of parents used for modelling, may be null

	// flags
	protected boolean dataLoaded;		// true if this object has been loaded with pdb data, false when is empty
	private boolean hasASA; 			// true if naccess has been run and ASA values assigned
	private boolean hasBfactors; 		// true if atom b-factors have been assigned

	/*----------------------------------  constructors -----------------------------------------------*/

	/**
	 * Constructs an empty Pdb with no sequence, no residues and no atoms 
	 */
	public Pdb() {
		this.chainCode = DEFAULT_CHAIN;
		this.pdbChainCode = DEFAULT_CHAIN;
		this.model = DEFAULT_MODEL;
		this.pdbCode = NO_PDB_CODE;
		
		this.dataLoaded = false;
		this.hasASA = false;
		this.hasBfactors = false;
		
		this.initialiseResidues();
	}
	
	/**
	 * Constructs a Pdb for given sequence with empty residues (residues with no atoms)
	 * @param sequence
	 * @throws IllegalArgumentException if sequence contains invalid characters
	 */
	public Pdb(String sequence) {
		this.chainCode = DEFAULT_CHAIN;
		this.pdbChainCode = DEFAULT_CHAIN;
		this.model = DEFAULT_MODEL;
		this.pdbCode = NO_PDB_CODE;
		
		this.dataLoaded = true;
		this.hasASA = false;
		this.hasBfactors = false;
		
		this.sequence = sequence;
		this.fullLength = sequence.length();
		this.initialiseResidues();
		this.resser2pdbresser = new TreeMap<Integer, String>();
		this.pdbresser2resser = new TreeMap<String, Integer>();

		for (int i=0;i<sequence.length();i++) {
			int resser = i+1;
			char one = sequence.charAt(i);
			if (!AminoAcid.isStandardAA(one)) {
				throw new IllegalArgumentException("Given input sequence to construct a pdb model contains an invalid aminoacid "+one);
			}
			this.addResidue(new Residue(AminoAcid.getByOneLetterCode(one),resser,this));
			resser2pdbresser.put(resser, String.valueOf(resser));
			pdbresser2resser.put(String.valueOf(resser), resser);
		}
		
		this.initialiseMaps();

	}
	
	/**
	 * Constructs a Pdb for given sequence setting coordinates of the given atom type to 
	 * the given coordinates coords
	 * @param sequence
	 * @param coords the array of all coordinates, must be ordered as the sequence and be
	 * of the same size 
	 * @param atom the atom code for which to set the coordinates
	 * @throws IllegalArgumentException if sequence contains invalid characters or if array
	 * of coordinates is of different length as sequence
	 */
	public Pdb(String sequence, Vector3d[] coords, String atom) {
		this.chainCode = DEFAULT_CHAIN;
		this.pdbChainCode = DEFAULT_CHAIN;
		this.model = DEFAULT_MODEL;
		this.pdbCode = NO_PDB_CODE;
		
		this.dataLoaded = true;
		this.hasASA = false;
		this.hasBfactors = false;
		
		if (coords.length!=sequence.length()) {
			throw new IllegalArgumentException("Array of coordinates is not of same length as given sequence");
		}
		
		this.sequence = sequence;
		this.fullLength = sequence.length();
		this.initialiseResidues();
		this.resser2pdbresser = new TreeMap<Integer, String>();
		this.pdbresser2resser = new TreeMap<String, Integer>();

		for (int i=0;i<sequence.length();i++) {
			int resser = i+1;
			char one = sequence.charAt(i);
			if (!AminoAcid.isStandardAA(one)) {
				throw new IllegalArgumentException("Given input sequence to construct a pdb model contains an invalid aminoacid "+one);
			}
			this.addResidue(new Residue(AminoAcid.getByOneLetterCode(one),resser,this));
			resser2pdbresser.put(resser, String.valueOf(resser));
			pdbresser2resser.put(String.valueOf(resser), resser);
		}
		
		int i = 0;
		for (Residue residue:this.residues.values()) {
			residue.addAtom(new Atom(i+1, atom, new Point3d(coords[i]), residue));
			i++;
		}
		this.initialiseMaps();
	}
	
	/*-------------------------   methods to be overriden by subclasses  -----------------------------*/
	// these methods should really be abstract, but we can't as Pdb is not an abstract class anymore
	// to work around that we implement them here just throwing an error if they are called at all
	// because in fact it's an error to call them!
	
	public void load(String pdbChainCode, int model) throws PdbLoadError {
		throw new PdbLoadError("Fatal error. This method shouldn't be called. This is a bug!");
	}
	
	public String[] getChains() throws PdbLoadError {
		throw new PdbLoadError("Fatal error. This method shouldn't be called. This is a bug!");
	}
	
	public Integer[] getModels() throws PdbLoadError {
		throw new PdbLoadError("Fatal error. This method shouldn't be called. This is a bug!");
	}

	// this method doesn't actually call load(String,int) above but rather the overridden ones in subclasses
	/**
	 * Load (from file/db) the given pdbChainCode and the default model {@link #DEFAULT_MODEL}
	 * @param pdbChainCode
	 * @throws PdbLoadError
	 */
	public void load(String pdbChainCode) throws PdbLoadError {
		load(pdbChainCode, DEFAULT_MODEL);
	}
	
	/*---------------------------------  public methods ----------------------------------------------*/

	/**
	 * Returns true if this Pdb has been loaded with pdb data (i.e. when 
	 * load(pdbChainCode) has been called), false if it is empty
	 */
	public boolean isDataLoaded() {
		return dataLoaded;
	}
	
	/**
	 * Initialises the residues to an empty map. 
	 */
	protected void initialiseResidues() {
		residues = new TreeMap<Integer, Residue>();
	}
	
	/**
	 * Adds a residue to this Pdb
	 * @param residue
	 */
	public void addResidue(Residue residue) {
		residues.put(residue.getSerial(), residue);
	}
	
	/**
	 * Gets a Residue given its serial
	 * @param resSerial
	 * @return the reference to the Residue object or null if residue serial doesn't exist
	 * for instance because it's not observed
	 */
	public Residue getResidue(int resSerial) {
		
		Residue res = residues.get(resSerial);
		return res;
	}
	
	/**
	 * Tells whether residue of given residue number is a (observed) residue in this Pdb
	 * instance. See also {@link #hasCoordinates(int)}
	 * @param resSerial
	 * @return
	 */
	public boolean containsResidue(int resSerial) {
		return residues.containsKey(resSerial);
	}
	
	/**
	 * Gets all residues of specified type in a List
	 * @param aa the amino acid type (see AminoAcid class)
	 * @return
	 */
	public ArrayList<Residue> getResiduesOfType(AminoAcid aa) {
		ArrayList<Residue> list = new ArrayList<Residue>();
		for (Residue residue:residues.values()) {
			if (residue.getAaType().equals(aa)) {
				list.add(residue);
			}
		}
		return list;
	}
	
	/**
	 * Gets a new Map with residue serials to Residues that contain only the atoms for the 
	 * given contact type and given interval set. The Residue objects are new, but the Atom 
	 * objects to which they point to are the same old references.
	 * @param ct the contact type
	 * @param intervSet only residues of this intervals will be considered, if null then
	 * all residues taken
	 * @return
	 */
	private TreeMap<Integer, Residue> getReducedResidues(String ct, IntervalSet intervSet) {
		TreeMap<Integer,Residue> reducedResidues = new TreeMap<Integer, Residue>();
		
		if (intervSet!=null) {
			for (Interval interv:intervSet) {
				for (int resser=interv.beg;resser<=interv.end;resser++) {
					Residue residue = getResidue(resser);
					if (residue==null) 
						throw new IllegalArgumentException("Invalid interval specified, residue "+resser+" is not part of this Pdb");
					reducedResidues.put(resser,residue.getReducedResidue(ct));
				}
			}
		} else { // we take all observed residues
			for (Residue residue:residues.values()) {
				reducedResidues.put(residue.getSerial(),residue.getReducedResidue(ct));
			}
		}
		return reducedResidues;
	}
	
	/**
	 * Gets an Atom given its serial
	 * @param atomSerial
	 * @return
	 */
	public Atom getAtom(int atomSerial) {
		return atomser2atom.get(atomSerial);
	}
	
	/**
	 * Populates the atomser2atom map (from the Residues and Atoms objects)
	 * and the pdbResSerials of the Residue objects from resser2pdbresser map
	 */
	protected void initialiseMaps() {
		atomser2atom = new TreeMap<Integer, Atom>();
		for (Residue residue:residues.values()) {
			for (Atom atom:residue.getAtoms()) {
				atomser2atom.put(atom.getSerial(), atom);
			}
			if (resser2pdbresser.containsKey(residue.getSerial())) {
				residue.setPdbSerial(getPdbResSerFromResSer(residue.getSerial()));
			}
		}
		
	}
	
	/**
	 * Sets the secondary structure members of the Residue objects from the SecondaryStructure object
	 */
	protected void initialiseResiduesSecStruct() {
		Iterator<SecStrucElement> it = this.getSecondaryStructure().iterator();
		while (it.hasNext()) {
			SecStrucElement ssElem = it.next();
			for (int resser=ssElem.getInterval().beg;resser<=ssElem.getInterval().end;resser++) {
				if (this.containsResidue(resser)) {
					this.getResidue(resser).setSsElem(ssElem);
				} else {
					// we don't warn because this does happen really often!
					//System.err.println("Warning: the residue serial "+resser+" can't be assigned with the secondary structure element "+
					//		ssElem.getId()+" with interval "+ssElem.getInterval()+" because it's not present in this Pdb (e.g. it's not observed)");
				}
			}			
		}
	}
	
	/**
	 * Calculates compactness coefficient as defined by:
	 * Micheal H. Zehfus and George D. Rose (1986) Compact Units in Proteins. Biochemistry 25: 5759-5765.
	 * DOI: 10.1021/bi00367a062
	 * @param surface
	 * @param volume
	 * @return
	 */
	public double calcCompactnessCoefficient(double surface, double volume) {
		return (surface/Math.pow(36*Math.PI*Math.pow(volume,2), 1.0/3.0));
	}
	
	/**
	 * Returns the per-residue relative solvent accessible surface areas (SASA) as 
	 * calculated by NACCESS.
	 * Returns null if SASA has not previously been calculated with {@link runner.NaccessRunner}.
	 * @return a map from residue serial to SASA value (null if not calculated yet)
	 */
	public HashMap<Integer, Double> getSurfaceAccessibilities() {
		if (hasASA()) {
			HashMap<Integer, Double> resser2allrsa = new HashMap<Integer, Double>();
			for (Residue residue:residues.values()) {
				resser2allrsa.put(residue.getSerial(), residue.getRsa());
			}
			return resser2allrsa;
		} else {
			return null;
		}
	}

	/**
	 * Returns the per-residue side-chain relative solvent accessible surface areas (SASA)
	 * as calculated by NACCESS.
	 * Returns null if SASA has not previously been calculated with {@link runner.NaccessRunner}.
	 * @return a map from residue serial to SASA value (null if not calculated yet)
	 */
	public HashMap<Integer, Double> getSideChainSurfaceAccessibilities() {
		if (hasASA()) {
			HashMap<Integer, Double> resser2scrsa = new HashMap<Integer, Double>();
			for (Residue residue:residues.values()) {
				resser2scrsa.put(residue.getSerial(), residue.getScRsa());
			}
			return resser2scrsa;
		} else {
			return null;
		}
	}

	/**
	 * Assigns b-factor values to the atoms of this structure. If structure is written to pdb file,
	 * these values will appear in the b-factor column. Currently, if no b-factors are assigned, the
	 * default value 0 will be written.
	 * @param bfactorsPerAtom a map of atom serials to b-factors
	 */
	public void setBFactorsPerAtom(HashMap<Integer, Double> bfactorsPerAtom) {
		for(int atomser:bfactorsPerAtom.keySet()) {
			getAtom(atomser).setBfactor(bfactorsPerAtom.get(atomser));
		}
		hasBfactors = true;
	}
	
	/**
	 * Calculates for each atom in this structure the deviation to the corresponding atom
	 * in the reference structure and returns a map from atom serials to differences in Angstrom.
	 * The corresponding atom is the atom with the same residue serial and the same PDB residue code
	 * (e.g. CA, CB, N), ignoring residues which are mutated between this and the reference structure.
	 * If no corresponding atom is found, the value Atom.DEFAULT_B_FACTOR will be used. To give useful
	 * results the two structures need to have the same residue numbering (save mutations) and have to
	 * be properly superimposed. This holds for example for two consecutive states in an MD trajectory
	 * or for different FoldX mutant models based in the same original structure.
	 * The returned values can be assigned to the b-factor column using this.setBFactorsPerAtom()
	 * and can then be visualized in PyMol with the command 'spectrum b'. 
	 * @param referencePdb the reference structures to compare to
	 * @return the per atom distances between this and the reference structure
	 */
	public HashMap<Integer, Double> getPerAtomDistances(Pdb referencePdb) {
		HashMap<Integer, Double> distances = new HashMap<Integer, Double>();
		
		// first set all distances to the default value, so that the map contains the full set of atoms
		for(int atomSer: this.getAllAtomSerials()) {
				distances.put(atomSer,Atom.DEFAULT_B_FACTOR);
		}
		
		// then set the real distance whereever possible
		for(int resSer: this.getAllSortedResSerials()) {
			Residue r = this.getResidue(resSer);
			AminoAcid thisType = r.getAaType();
			if(referencePdb.hasCoordinates(resSer)) {
				AminoAcid refType = referencePdb.getResidue(resSer).getAaType();
				if(thisType == refType) {
					for(Atom a1:r.getAtoms()) {
						int atomSer = a1.getSerial();
						if(referencePdb.hasCoordinates(resSer, a1.getCode())) {
							Atom a2 = referencePdb.getAtom(referencePdb.getAtomSerFromResSerAndAtom(resSer, a1.getCode()));
							// a2 should never be null
							double dist = Math.abs(a1.getCoords().distance(a2.getCoords()));
							distances.put(atomSer,dist);
						}
					}
				} else {
					System.err.println("Skipping mismatching residue at position " + resSer + " ("+ thisType + "!=" + refType + ")");
				}
			} else {
				System.err.println("Skipping " + thisType.getThreeLetterCode() + " at position " + resSer + ". Not found in reference structure.");
			}
		}
		return distances;
	}
	
	/**
	 * Assigns b-factor values to all atoms for the given residues.
	 * @param bfactorsPerResidue a map of residue serials to b-factors
	 * @throws IllegalArgumentException when residue serials given in input are not
	 * present in this Pdb instance
	 */
	public void setBFactorsPerResidue(HashMap<Integer, Double> bfactorsPerResidue) {
		for(int resser:bfactorsPerResidue.keySet()) {
			if (this.containsResidue(resser)) {
				Residue residue = this.getResidue(resser);
				double bfactor = bfactorsPerResidue.get(resser);
				for(Atom atom:residue.getAtoms()) {
					atom.setBfactor(bfactor);
				}
			} else {
				System.err.println("Warning! Can't assign bfactor for residue serial "+resser+", it is not present in this Pdb instance, pdbCode "+pdbCode+", pdbChainCode "+pdbChainCode);
			}
		}
		hasBfactors = true;
	}
	
	/**
	 * Writes to given PrintWriter the PDB file format HEADER line
	 * @param Out
	 */
	public void writePDBFileHeader(PrintStream Out) {
		String source = "";
		if (this instanceof CiffilePdb) {
			source = ((CiffilePdb) this).getCifFile().getAbsolutePath();
		} else if (this instanceof PdbfilePdb) {
			source = ((PdbfilePdb) this).getPdbFileName();
		} else if (this instanceof PdbasePdb){
			source = ((PdbasePdb) this).getDb();
		} else {
			source = "model";
		}
		Out.println("HEADER  Source: "+source+". "+pdbCode+", chain='"+chainCode+"', model="+model);		
	}
	
	/**
	 * Write CASP TS file headers complying with the CASP 8 specification.
	 * Note that the CASP target number, CASP model number, CASP author
	 * and CASP method will be written from the internally 
	 * set values (targetNum, caspModelNum, caspAuthorStr, caspMethodStr) 
	 * so they must be set before trying to write them out.
	 * If the method string has new line characters it will be splitted into
	 * several METHOD lines
	 * @param Out
	 * @param parents PDB entries in which this homology prediction is based on or
	 * null if not applicable i.e. if this is an ab-initio prediction
	 */
	private void writeCaspTSHeader(PrintStream Out, String[] parents) {
		Out.println("PFRMAT TS");
		Out.printf("TARGET T%04d\n",targetNum);
		if(caspAuthorStr != null) Out.printf("AUTHOR %s\n", caspAuthorStr);
		if(caspMethodStr != null) {
			String[] caspMethodLines = caspMethodStr.split("\n");
			for (String caspMethodLine: caspMethodLines) {
				Out.printf("METHOD %s\n", caspMethodLine);
			}
		}
		Out.println("MODEL "+caspModelNum);
		String parentStr = "";
		if (parents==null) {
			parentStr = "N/A";
		} else {
			for (int i=0;i<parents.length;i++) {
				Pattern p = Pattern.compile("(\\d\\w\\w\\w)(\\w)");
				Matcher m = p.matcher(parents[i]);
				if (m.matches()) {
					parentStr += m.group(1)+"_"+m.group(2) + " ";
				} else {
					System.err.println("Parent string "+parents[i]+" to be written to file CASP TS file not in the expected format 1abcA ");
					parentStr += parents[i] + " ";
				}
				
			}
		}
		Out.println("PARENT "+parentStr);
	}
	
	/** 
	 * Writes atom lines for this structure to the given output stream
	 * @param Out
	 */
	public void writeAtomLines(PrintStream Out) {
		String chainCodeStr = chainCode;
		if (chainCode.length()>1) {
			System.err.println("Warning! Chain code with more than 1 character ("+chainCode+"), only first character will be written to ATOM lines");
			chainCodeStr = chainCode.substring(0,1);
		}
		for (int atomser:this.getAllAtomSerials()) {
			Atom atom = getAtom(atomser);
			int resser = atom.getParentResSerial();
			String atomCode = atom.getCode();
			String atomType = atom.getType().getSymbol();
			String res = atom.getParentResidue().getAaType().getThreeLetterCode();
			Point3d coords = atom.getCoords();
			double occupancy = atom.getOccupancy();
			double bFactor = atom.getBfactor();
			
			String printfStr = "ATOM  %5d  %-3s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %s\n";
			if (atomCode.length()==4) { // some hydrogens have a 4 letter code and it is not aligned to column 14 but to 13 instead 
				printfStr = "ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %s\n";
			}
			// Local.US is necessary, otherwise java prints the doubles locale-dependant (i.e. with ',' for some locales)
			Out.printf(Locale.US, printfStr,
					atomser, atomCode, res, chainCodeStr, resser, coords.x, coords.y, coords.z, occupancy, bFactor, atomType);
		}
	}

	/**
	 * Dumps coordinate data into a file in PDB format (ATOM lines only)
	 * The residue serials written are the internal ones.
	 * The chain dumped is the value of the chainCode variable, i.e. our internal
	 * chain identifier for Pdb objects 
	 * @param outfile
	 * @throws FileNotFoundException
	 */
	public void writeToPDBFile(String outfile) throws FileNotFoundException {
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		writePDBFileHeader(Out);
		writeAtomLines(Out);
		Out.println("END");
		Out.close();
	}

	/**
	 * Writes coordinates to given File in CASP TS format.
     * Note that the CASP target number, CASP model number, CASP author, 
     * CASP method and CASP parents will be written from the internally set values
     * (targetNum, caspModelNum, caspAuthorStr, caspMethodStr, caspParents) so they
     * must be set before trying to write them out.
	 * @param outFile
	 * @throws FileNotFoundException
	 */
	public void writeToCaspTSFile(File outFile) throws FileNotFoundException {
		PrintStream Out = new PrintStream(new FileOutputStream(outFile));
		writeCaspTSHeader(Out, this.caspParents);
		String oldChainCode = this.chainCode;
		this.chainCode = DEFAULT_CASP_TS_CHAINCODE;
		writeAtomLines(Out);
		this.chainCode = oldChainCode;
		Out.println("TER"); // note that CASP TS requires a TER field at the end of each model
		Out.println("END");
		Out.close();
	}
	
	/**
	 * Dump the full sequence of this Pdb object in FASTA file format
	 * The FASTA tag is written as the concatenation of pdbCode and pdbChainCode
	 * @param seqfile
	 * @throws IOException if file can't be written
	 */
	public void writeSeqToFasta(String seqfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(seqfile));
		int len = 80;
		Out.println(">"+pdbCode+pdbChainCode);
		for(int i=0; i<sequence.length(); i+=len) {
			Out.println(sequence.substring(i, Math.min(i+len,sequence.length())));
		}		
		Out.close();
	}

	/** 
	 * Returns the number of observed standard residues.
	 * @return number of observed standard residues
	 */
	public int getObsLength(){
		return residues.size();
	}

	/** 
	 * Returns the number of residues in the sequence of this protein.
	 * @return number of residues in the full sequence
	 */
	public int getFullLength() {
		return fullLength;
	}

	/**
	 * Returns number of atoms in the protein, including Hydrogens if they are present
	 * @return number of atoms
	 */
	public int getNumAtoms() {
		return atomser2atom.size();
	}
	
	/**
	 * Returns number of heavy (non-Hydrogen) atoms in the protein
	 * @return number of (non-Hydrogen) atoms
	 */
	public int getNumHeavyAtoms() {
		int numAtoms = 0;
		for (Residue residue: residues.values()) {
			numAtoms+=residue.getNumHeavyAtoms();
		}
		return numAtoms;
	}

	/**
	 * Returns the number of chain breaks in this structure, i.e. how many unobserved 
	 * gaps are there in the chain (excluding unobserved regions in the terminals)
	 * @return
	 */
	public int getNumChainBreaks() {
		int numChainBreaks = 0;
		int lastResser = 0;
		for (int resser:getAllSortedResSerials()) {
			if (lastResser!=0 && resser-lastResser>1) 
				numChainBreaks++;
			lastResser = resser;
		}
		return numChainBreaks;
	}
	
	/**
	 * Gets a TreeMap with atom serials as keys and their coordinates as values for the 
	 * given contact type.
	 * The contact type can't be a cross contact type, it doesn't make sense here
	 * @param ct
	 * @return
	 */
	private TreeMap<Integer,Point3d> getCoordsForCt(String ct) {
		TreeMap<Integer,Point3d> coords = new TreeMap<Integer,Point3d>();
		for (Residue residue:residues.values()) {
			Set<String> atomCodes = AAinfo.getAtomsForCTAndRes(ct, residue.getAaType().getThreeLetterCode());
			for (String atomCode:atomCodes){
				if (residue.containsAtom(atomCode)) {
					coords.put(residue.getAtom(atomCode).getSerial(),residue.getAtom(atomCode).getCoords());
				} 
			}
			// in cts ("ALL","BB") we still miss the OXT, we need to add it now if it is there (it will be there when this resser is the last residue)
			if ((ct.equals("ALL") || ct.equals("BB")) &&  residue.containsOXT()) { 
				coords.put(residue.getAtom("OXT").getSerial(), residue.getAtom("OXT").getCoords());
			}
		}
		return coords;
	}

	/**
	 * Returns the distance matrix (one distance per pair of residues) as a Jama Matrix object. 
	 * For multi atom contact types the distance matrix has the minimum distance for each pair of
	 * residues. 
	 * The indices of the matrix will be 0 to get_length()-1 and will be in the same order 
	 * as the residue serials
	 * @param ct contact type for which distances are being calculated
	 * @return
	 */
	public Matrix calcDistMatrixJamaFormat(String ct) {
		HashMap<Integer,Integer> resser2ind = new HashMap<Integer, Integer>();
		int i = 0;
		for (int resser:this.getAllSortedResSerials()) {
			resser2ind.put(resser,i);
			i++;
		}
		HashMap<Pair<Integer>, Double> distHM = this.calcDistMatrix("Ca");
		Matrix matrix = new Matrix(this.getObsLength(), this.getObsLength());
		for (Pair<Integer> pair: distHM.keySet()) {
			double currentElem = distHM.get(pair);
			matrix.set(resser2ind.get(pair.getFirst()), resser2ind.get(pair.getSecond()), currentElem);
			matrix.set(resser2ind.get(pair.getSecond()), resser2ind.get(pair.getFirst()), currentElem);
		}
		return matrix;
	}
	
	/**
	 * Returns the distance matrix as a HashMap with residue serial pairs as keys
	 * For multi atom contact types the distance matrix has the minimum distance for each pair of
	 * residues 
	 * AAinfo.isValidSingleAtomCT(ct) can be used to check before calling.
	 * @param ct contact type for which distances are being calculated
	 * @return a map which assigns to each edge the corresponding distance 
	 */
	public HashMap<Pair<Integer>, Double> calcDistMatrix(String ct){
		HashMap<Pair<Integer>,Double> distMatrixAtoms = calcAtomDistMatrix(ct);

		 // mapping atom serials to residue serials
		 // TODO: we could integrate this with the code in calculate_atom_dist_matrix to avoid storing two distance maps in memory
		HashMap<Pair<Integer>,Double> distMatrixRes = new HashMap<Pair<Integer>, Double>();
		for (Pair<Integer> cont: distMatrixAtoms.keySet()){
			int i_resser = getResSerFromAtomSer(cont.getFirst());
			int j_resser = getResSerFromAtomSer(cont.getSecond());
			Pair<Integer> edge = new Pair<Integer>(i_resser,j_resser);
			if (distMatrixRes.containsKey(edge)) {
				distMatrixRes.put(edge, Math.min(distMatrixRes.get(edge), distMatrixAtoms.get(cont)));
			} else {
				distMatrixRes.put(edge, distMatrixAtoms.get(cont));
			}
		}

		return distMatrixRes;
	}
	
	/**
	 * Returns the distance matrix as a HashMap with atom serial pairs as keys
	 * This method can be used for any contact type 
	 * AAinfo.isValidSingleAtomCT(ct) can be used to check before calling.
	 * @param ct contact type for which distances are being calculated
	 * @return a map which assings to each atom pair the corresponding distance
	 */
	public HashMap<Pair<Integer>, Double> calcAtomDistMatrix(String ct){
		HashMap<Pair<Integer>,Double> distMatrixAtoms = new HashMap<Pair<Integer>,Double>();
		if (!ct.contains("/")){
			TreeMap<Integer,Point3d> coords = getCoordsForCt(ct);
			for (int i_atomser:coords.keySet()){
				for (int j_atomser:coords.keySet()){
					if (j_atomser>i_atomser) {
						Pair<Integer> pair = new Pair<Integer>(i_atomser,j_atomser);
						distMatrixAtoms.put(pair, coords.get(i_atomser).distance(coords.get(j_atomser)));
					}
				}
			}
		} else {
			String i_ct = ct.split("/")[0];
			String j_ct = ct.split("/")[1];
			TreeMap<Integer,Point3d> i_coords = getCoordsForCt(i_ct);
			TreeMap<Integer,Point3d> j_coords = getCoordsForCt(j_ct);
			for (int i_atomser:i_coords.keySet()){
				for (int j_atomser:j_coords.keySet()){
					if (j_atomser!=i_atomser){
						Pair<Integer> pair = new Pair<Integer>(i_atomser,j_atomser);
						distMatrixAtoms.put(pair, i_coords.get(i_atomser).distance(j_coords.get(j_atomser)));
					}
				}
			}
		}

		return distMatrixAtoms;
	}
	
	/**
	 * Calculates the radius of gyration of this Pdb 
	 * (defined as half of the maximum distance between any 2 CA atoms)
	 * @return
	 */
	public double calcRadGyration() {
		//TODO this is a very raw implementation o(n^2): should optimise it if that's really needed
		return Collections.max(this.calcAtomDistMatrix("Ca").values())/2;
	}
	
	/**
	 * Get the graph for given contact type and cutoff for this Pdb object.
	 * Returns a Graph object with the contacts
	 * We do geometric hashing for fast contact computation (without needing to calculate full distance matrix)
	 * @param ct
	 * @param cutoff
	 * @return
	 */
	private AIGraph getAIGraph(String ct, double cutoff){ 
		TreeMap<Integer,Point3d> i_coords = null;
		TreeMap<Integer,Point3d> j_coords = null;		// only relevant for asymetric edge types
		boolean crossed = false;
		if (!ct.contains("/")){
			i_coords = getCoordsForCt(ct);
			crossed = false;
		} else {
			String i_ct = ct.split("/")[0];
			String j_ct = ct.split("/")[1];
			i_coords = getCoordsForCt(i_ct);
			j_coords = getCoordsForCt(j_ct);
			crossed = true;
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
				if (!crossed){
					boxes.get(floor).put_j_Point(i, coord);
				}
			} else {
				Box box = new Box(floor);
				box.put_i_Point(i, coord);
				if (!crossed){
					box.put_j_Point(i, coord);
				}
				boxes.put(floor,box);
			}
			i_atomserials[i]=i_atomser; //as atomserials in coords were ordered (TreeMap) the new indexing will still be ordered
			i++;
		}
		int j=0;
		if (crossed) {
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
			boxes.get(floor).getDistancesWithinBox(distMatrix,crossed);

			//TODO should iterate only through half of the neighbours here 
			// distances of points from this box to all neighbouring boxes: 26 iterations (26 neighbouring boxes)
			for (int x=floor.x-boxSize;x<=floor.x+boxSize;x+=boxSize){
				for (int y=floor.y-boxSize;y<=floor.y+boxSize;y+=boxSize){
					for (int z=floor.z-boxSize;z<=floor.z+boxSize;z+=boxSize){
						if (!((x==floor.x)&&(y==floor.y)&&(z==floor.z))) { // skip this box
							Point3i neighbor = new Point3i(x,y,z);
							if (boxes.containsKey(neighbor)){
								boxes.get(floor).getDistancesToNeighborBox(boxes.get(neighbor),distMatrix,crossed);
							}
						}
					}
				}
			} 
		} 

		// creating the AIGraph
		AIGraph graph = new AIGraph();
		TreeMap<Integer,RIGNode> rignodemap = new TreeMap<Integer,RIGNode>();
		TreeSet<Integer> atomSerials = new TreeSet<Integer>();
		atomSerials.addAll(i_coords.keySet());
		if (j_coords!=null){
			atomSerials.addAll(j_coords.keySet());
		}
		// adding the AIGNodes (including parent RIGNode references)
		SecondaryStructure secondaryStructureCopy = secondaryStructure.copy();
		for (int atomSer:atomSerials) {
			int resser = getResSerFromAtomSer(atomSer);
			SecStrucElement sselem = secondaryStructureCopy.getSecStrucElement(resser);
			if (!rignodemap.containsKey(resser)) {
				// NOTE!: we are passing references to the SecStrucElement objects! they point to the same objects as secondaryStructureCopy 
				RIGNode resNode = new RIGNode(resser, getResidue(resser).getAaType().getThreeLetterCode(), sselem);
				rignodemap.put(resser,resNode);
			}
			AIGNode atomNode = new AIGNode(atomSer,getAtom(atomSer).getCode(),rignodemap.get(resser));
			graph.addVertex(atomNode); // this also adds the atomNode to the serials2nodes map
		}
		
		graph.setSecondaryStructure(secondaryStructureCopy);
		graph.setCutoff(cutoff);
		graph.setSequence(sequence);
		graph.setPdbCode(pdbCode);
		graph.setChainCode(chainCode);
		graph.setPdbChainCode(pdbChainCode);
		graph.setModel(model);
		graph.setSid(sid);
		graph.setTargetNum(targetNum);
		graph.setGroupNum(groupNum);
		graph.setCaspModelNum(caspModelNum);
		graph.setMethodStr(caspMethodStr);
		graph.setAuthorStr(caspAuthorStr);
		graph.setParents(caspParents);
		graph.setCrossed(crossed);
		
		// populating the AIGraph with AIGEdges 
		for (i=0;i<distMatrix.length;i++){
			for (j=0;j<distMatrix[i].length;j++){
				// the condition distMatrix[i][j]!=0.0 takes care of skipping several things: 
				// - diagonal of the matrix in case of non-crossed
				// - lower half of matrix in case of non-crossed
				// - cells for which we didn't calculate a distance because the 2 points were not in same or neighbouring boxes (i.e. too far apart)
				if (distMatrix[i][j]!=0.0f && distMatrix[i][j]<=cutoff){
					if (!crossed) {
						graph.addEdge(new AIGEdge(distMatrix[i][j]), graph.getNodeFromSerial(i_atomserials[i]), graph.getNodeFromSerial(j_atomserials[j]), EdgeType.UNDIRECTED);
					}
					// This condition is to take care of crossed contact types that have overlapping sets of atoms: 
					//   the matrix would contain both i,j and j,i but that's only 1 edge in the AIGraph
					//TODO if our AIGraph didn't allow parallel edges, this extra check woudn't be necessary
					else if (!graph.containsEdgeIJ(i_atomserials[i], j_atomserials[j])) {
						graph.addEdge(new AIGEdge(distMatrix[i][j]), graph.getNodeFromSerial(i_atomserials[i]), graph.getNodeFromSerial(j_atomserials[j]), EdgeType.UNDIRECTED);
					}
				}

			}
		}

		return graph;
	}

	/**
	 * Returns a RIGraph for given contact type, cutoff and directionality
	 * @param ct  the contact type
	 * @param cutoff  the distance cutoff
	 * @param directed  true if we want a directed graph, false for undirected
	 * @return
	 */
	public RIGraph getRIGraph(String ct, double cutoff, boolean directed) {
		//NOTE: At the moment we don't allow directed graphs for overlapping contact types e.g. directed ALL/BB
		//      because the code wouldn't be able to cope correctly with them
		// To lift this restriction, one possibility would be to make AIGraph directed and
		//- remove the check for adding parallel atomic edges if crossed
		//- add atomic edges in both directions if !crossed
		//- make sure in getRIGraph for undirected graphs not to double-count atomic edges
		if (directed && AAinfo.isOverlapping(ct)) {
			throw new IllegalArgumentException("Contact type "+ct+" is overlapping. Generating directed graphs for it is unsupported");
		}
		
		String[] cts = ct.split("\\+");		
		AIGraph atomGraph = getAIGraph(cts[0], cutoff);
		for(int i=1; i< cts.length; i++) {
			atomGraph.addGraph(getAIGraph(cts[i], cutoff));
		}
		RIGraph graph = atomGraph.getRIGraph(directed);
		
		graph.setContactType(ct);
		graph.setCutoff(cutoff);

		return graph;
	}
	
	/**
	 * Returns a RIGraph for given contact type and cutoff
	 * Crossed contact types (i.e. those containing a "/") will be considered always as 
	 * directed. If one wants an undirected graph for a crossed contact type then method
	 * {@link #getRIGraph(String, double, boolean)} should be used instead.
	 * @param ct  the contact type
	 * @param cutoff  the distance cutoff
	 * @return
	 */
	public RIGraph getRIGraph(String ct, double cutoff) {
		// TODO eventually we should use the ContactType class as parameter so we could encapsulate all properties of a contact type in there
		boolean crossed = false;
		if (ct.contains("/")) {
			crossed = true;
		}		
		return getRIGraph(ct, cutoff, crossed);
	}
	
	/**
	 * Returns an all atom graph in a AIGraph object
	 * @param cutoff  the distance cutoff
	 * @return
	 */
	public AIGraph getAllAtomGraph(double cutoff) {
		return this.getAIGraph("ALL", cutoff);
	}
	
	public void calcGridDensity(String ct, double cutoff, Map<Integer, Integer> densityCount) { 
		TreeMap<Integer,Point3d> i_coords = null;
		TreeMap<Integer,Point3d> j_coords = null;		// only relevant for asymmetric edge types
		boolean directed = false;
		if (!ct.contains("/")){
			i_coords = getCoordsForCt(ct);			// mapping from atom serials to coordinates
			directed = false;
		} else {
			String i_ct = ct.split("/")[0];
			String j_ct = ct.split("/")[1];
			i_coords = getCoordsForCt(i_ct);
			j_coords = getCoordsForCt(j_ct);
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
	 */
	public HashMap<Pair<Integer>,Double> getDiffDistMap(String contactType1, Pdb pdb2, String contactType2, MultipleSequenceAlignment ali, String name1, String name2) {

		HashMap<Pair<Integer>,Double> otherDistMatrix = pdb2.calcDistMatrix(contactType2);
		HashMap<Pair<Integer>,Double> thisDistMatrix = this.calcDistMatrix(contactType1);
		HashMap<Pair<Integer>,Double> alignedDistMatrix = new HashMap<Pair<Integer>, Double>(Math.min(this.getFullLength(), pdb2.getFullLength()));
		int i1,i2,j1,j2;
		TreeSet<Integer> unobserved1 = new TreeSet<Integer>();
		TreeSet<Integer> unobserved2 = new TreeSet<Integer>();

		// detect all unobserved residues
		for(int i = 1; i <= ali.getAlignmentLength(); ++i) {
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

		for(int i = 1; i <= ali.getAlignmentLength()-1; ++i) {

			i1 = ali.al2seq(name1, i);
			i2 = ali.al2seq(name2, i);

			// alignment columns must not contain gap characters and both 
			// residues in the current column have to be observed!
			if( i1 == -1 || i2 == -1 || unobserved1.contains(i1) || unobserved2.contains(i2) ) {
				continue;
			}

			for(int j = i + 1; j <= ali.getAlignmentLength(); ++j) {

				j1 = ali.al2seq(name1, j);
				j2 = ali.al2seq(name2, j);

				if( j1 == -1 || j2 == -1 || unobserved1.contains(j1) || unobserved2.contains(j2) ) {
					continue;
				}

				// make the edges
				Pair<Integer> e1 = new Pair<Integer>(i1,j1);
				Pair<Integer> e2 = new Pair<Integer>(i2,j2);

				alignedDistMatrix.put(new Pair<Integer>(i,j),Math.abs(thisDistMatrix.get(e1)-otherDistMatrix.get(e2)));
			}
		}
		return alignedDistMatrix;
	}
	// TODO: Version of this where already buffered distance matrices are passed as paremeters

	/** 
	 * Returns the number of neighbours of this grid cell
	 * @param boxes
	 * @param floor
	 * @param boxSize
	 * @return 
	 */
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
	 * @throws NullPointerException if residue serial not present in this Pdb instance
	 */
	public Double getConsurfhsspScoreFromResSerial(int resser){
		return getResidue(resser).getConsurfScore();
	}

	/**
	 * Gets the Consurf-HSSP color rsa given an internal residue serial
	 * @param resser
	 * @throws NullPointerException if residue serial not present in this Pdb instance 
	 * @return
	 */
	public Integer getConsurfhsspColorFromResSerial(int resser){
		return getResidue(resser).getConsurfColor();
	}

	/**
	 * Gets the all atoms rsa given an internal residue serial
	 * @param resser
	 * @return the rsa or null if rsa has not been calculated yet or the residue number cannot be found
	 * @throws NullPointerException if residue serial not present in this Pdb instance 
	 */
	public Double getAllRsaFromResSerial(int resser){
		if (hasASA()) {
			return getResidue(resser).getRsa();
		} else {
			return null;
		}
	}

	/**
	 * Gets the sc rsa given an internal residue serial
	 * @param resser
	 * @return
	 * @throws NullPointerException if residue serial not present in this Pdb instance 
	 */
	public Double getScRsaFromResSerial(int resser){
		if (hasASA()) {
			return getResidue(resser).getScRsa();
		} else {
			return null;
		}
	}

	/**
	 * Gets the internal residue serial (cif) given a pdb residue serial (author assignment)
	 * @param pdbresser
	 * @return
	 */
	public int getResSerFromPdbResSer (String pdbresser){
		return pdbresser2resser.get(pdbresser);
	}

	/**
	 * Gets the pdb residue serial (author assignment) given an internal residue serial (cif)
	 * @param resser
	 * @return 
	 */
	public String getPdbResSerFromResSer (int resser){
		return this.resser2pdbresser.get(resser);
	}

	/**
	 * Gets the residue serial given an atom serial
	 * @param atomser
	 * @return
	 */
	public int getResSerFromAtomSer(int atomser){
		return atomser2atom.get(atomser).getParentResSerial();
	}

	/**
	 * Convenience method to get residue type (3-letter code) from residue serial.
	 * It is also possible to use {@link #getResidue(int)} and from that get any
	 * other residue properties.
	 * @param resser
	 * @return
	 * @throws NullPointerException if residue serial not present in this Pdb instance
	 * or if residue serial is unobserved. Use {@link #hasCoordinates(int)} to check 
	 */
	public String getResTypeFromResSerial(int resser) {
		return getResidue(resser).getAaType().getThreeLetterCode();
	}

	/**
	 * Gets the atom serial given the residue serial and atom code.
	 * The caller of this method needs to check whether the resser, atom and combination of the two exists
	 * using {@link #hasCoordinates(int, String)}. Otherwise, if this function
	 * is called and no atom serial exists, a null pointer exception will be thrown.
	 * @param resser
	 * @param atomCode
	 * @return the atom serial
	 * @throws NullPointerException if residue serial not present in this Pdb 
	 * instance or if no atom of given type exists for given residue
	 */
	public int getAtomSerFromResSerAndAtom(int resser, String atomCode) {
		return getResidue(resser).getAtom(atomCode).getSerial();			
	}

	/**
	 * Returns the set of all atom serials for observed atoms of the given residue serial.
	 * @param resser the residue serial
	 * @return an ordered set of serials of the (observed) atoms in this residue
	 * @throws NullPointerException if given residue serial not present in this Pdb
	 * instance
	 */
	public Set<Integer> getAtomSersFromResSer(int resser) {
		TreeSet<Integer> atomSers = new TreeSet<Integer>();
		for (Atom atom:this.getResidue(resser).getAtoms()) {
			atomSers.add(atom.getSerial());
		}
		return atomSers;
	}
	
	/**
	 * Checks whether this Pdb is an all-atom one (disregarding Hydrogens), i.e. is not a 
	 * CA only or BB only or has no other major group of atoms missing.
	 * Even PDB structures with all atoms can still have missing atoms for some residues, 
	 * here what we check is that the average number of (non-Hydrogen) atoms per residue is above the 
	 * threshold {@value #MIN_AVRG_NUM_ATOMS_RES} . This threshold has been obtained from statistics of a set of non-redundant 
	 * PDB structures.
	 * 
	 * @return true if above average atoms per residue threshold, false otherwise
	 */
	public boolean isAllAtom() {
		if (((double)this.getNumHeavyAtoms()/(double)this.getObsLength())<MIN_AVRG_NUM_ATOMS_RES) {
			return false;
		} else {
			return true;
		}
	}
	
	/**
	 * Checks whether the given residue serial has any associated coordinates
	 * for its atoms.
	 * @param resser  the residue serial
	 * @return true if there is at least one atom with coordinates, else false 
	 */
	public boolean hasCoordinates(int resser) {
		if (this.containsResidue(resser)) {
			if (this.getResidue(resser).getNumAtoms()>0) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}

	/**
	 * Checks whether the given atom type for the given residue serial has any associated coordinates
	 * @param resser the residue serial
	 * @param atomCode  
	 * @return true if there is coordinates for given resser-atom, else false 
	 */
	public boolean hasCoordinates(int resser, String atomCode) {
		if (!this.containsResidue(resser)) 
			return false;
		return getResidue(resser).containsAtom(atomCode);
	}
	
	/**
	 * Returns true if this Pdb has been restricted to a specific SCOP domain 
	 * @return
	 */
	public boolean isRestrictedToScopDomain() {
		return sid!=null;
	}
	
	/**
	 * Returns the sid of this Pdb 
	 * It is set when restrictToScopDomain is run
	 */
	public String getSid() {
		return sid;
	}
	
	/**
	 * @return the title of this pdb object (may be null).
	 */
	public String getTitle() {
		return this.title;
	}
	
	/**
	 * Sets the title for this pdb object (may be null).
	 * @param title the new title
	 */
	public void setTitle(String title) {
		this.title = title;
	}
	
	/**
	 * Gets the atom coordinates (Point3d object) given the atom serial
	 * @param atomser
	 * @return
	 */
	public Point3d getAtomCoord(int atomser) {
		return getAtom(atomser).getCoords();
	}

	/**
	 * Gets the atom coordinates (a Point3d) given the residue serial and atom code 
	 * (standard PDB atom name, e.g. CA, N, C, O, CB, ...) 
	 * @param resser
	 * @param atomCode
	 * @return
	 */
	public Point3d getAtomCoord(int resser, String atomCode) {
		return getResidue(resser).getAtom(atomCode).getCoords();
	}
	
	/**
	 * Gets the atom code (standard PDB atom name, e.g. CA, N, C, O, CB, ...) given the 
	 * atom serial
	 * @param atomser
	 * @return
	 */
	public String getAtomNameFromAtomSer(int atomser) {
		return this.getAtom(atomser).getCode();
	}
	
	/**
	 * Returns a Set of all ordered atom serials (only observed atoms)
	 * The order is according to atom serials (not necessarily ordered by residue
	 * although they almost always are)
	 * @return
	 */
	public Set<Integer> getAllAtomSerials() {
		return this.atomser2atom.keySet();
	}

	/**
	 * Returns a Set of all ordered residue serials (only observed residues)
	 * @return
	 */
	public Set<Integer> getAllSortedResSerials() {
		return residues.keySet();
	}
	
	/**
	 * Returns the lowest observed residue serial 
	 * @return
	 */
	public int getMinObsResSerial() {
		return residues.firstEntry().getKey();
	}
	
	/**
	 * Returns the highest observed residue serial 
	 * @return
	 */
	public int getMaxObsResSerial() {
		return residues.lastEntry().getKey();
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
	 * Gets the model number of this Pdb 
	 * @return
	 */
	public int getModel() {
		return this.model;
	}
	
	/**
	 * Gets the CASP target number
	 * @return
	 */
	public int getTargetNum() {
		return this.targetNum;
	}
	
	/**
	 * Sets the CASP target number
	 * @param targetNum
	 */
	public void setTargetNum(int targetNum) {
		this.targetNum = targetNum;
	}
	
	/**
	 * Gets the CASP model number
	 * @return
	 */
	public int getCaspModelNum() {
		return this.caspModelNum;
	}

	/**
	 * Sets the CASP model number
	 * @param caspModelNum
	 */
	public void setCaspModelNum(int caspModelNum) {
		this.caspModelNum = caspModelNum;
	}
	
	/**
	 * Sets the author string required for casp submissions, set to null to supress output.
	 * @param authorStr the author string to be set
	 */
	public void setCaspAuthorStr(String authorStr) {
		this.caspAuthorStr = authorStr;
	}

	/**
	 * Sets the method string required for casp submissions, set to null to supress output.
	 * @param methodStr the method string to be set
	 */
	public void setCaspMethodStr(String methodStr) {
		this.caspMethodStr = methodStr;
	}
	
	/**
	 * Sets the parents record which Casp models may optionally have to specify parents used in modelling the strucure.
	 * The value is used when writing to Casp TS files or by the getParents() method.
	 * @param parents an array of parent strings used in modelling this structure.
	 */
	public void setParents(String[] parents) {
		this.caspParents = parents;
	}	
	
	/**
	 * Returns the list of parents which Casp models may optionally have.
	 * The value will be set when reading from Casp TS files or by the setParents method.
	 * @return an array of parent strings or null if no parents have been specified.
	 */
	public String[] getParents() {
		return this.caspParents;
	}
	
	/**
	 * Gets the sequence
	 * @return
	 */
	public String getSequence() {
		return sequence;
	}
	
	/**
	 * Gets the observed sequence, i.e. the sequence as it appears in the ATOM 
	 * lines of the PDB file (non-standard aas are not in this sequence even if 
	 * they have coordinates)
	 * @return
	 */
	public String getObsSequence() {
		String obsSequence = "";
		for (Residue residue:residues.values()) {
			obsSequence += residue.getAaType().getOneLetterCode();
		}
		return obsSequence;
	}
	
	/**
	 * True if this Pdb has the sequence field set to not blank 
	 * @return
	 */
	public boolean hasSequence() {
		if (sequence==null) return false;
		return !sequence.equals("");
	}
	
	/**
	 * Returns true if this Pdb instance has been assigned atom b-factors 
	 * @return
	 */
	public boolean hasBfactors() {
		return hasBfactors;
	}
	
	/**
	 * Returns true if naccess has been run and thus ASA values are present
	 * in this Pdb instance
	 * @return
	 */
	public boolean hasASA() {
		return hasASA;
	}

	/**
	 * Sets the hasASA flag. To be used when running Naccess and setting the ASA values 
	 * of residues.
	 * @param hasASA
	 */
	public void setHasASA(boolean hasASA) {
		this.hasASA = hasASA;
	}
	
	/** 
	 * Returns true if csa information is available, false otherwise. 
	 * @return
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
	 * Returns the csa annotation object of this Pdb.
	 * @return
	 */
	public CatalSiteSet getCSA() {
		return catalSiteSet;
	}

	/** 
	 * Returns true if ec information is available, false otherwise. 
	 * @return
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
	 * @return
	 */
	public EC getEC() {
		return ec;
	}

	/** 
	 * Returns true if scop information is available, false otherwise.
	 * @return 
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
	 * @return
	 */
	public Scop getScop() {
		return scop;
	}

	// secondary structure related methods

	/** 
	 * Returns true if secondary structure information is available, false otherwise.
	 * @return 
	 */
	public boolean hasSecondaryStructure() {
		return !this.secondaryStructure.isEmpty();
	}

	/**
	 * Returns the secondary structure annotation object of this Pdb.
	 * @return
	 */
	public SecondaryStructure getSecondaryStructure() {
		return this.secondaryStructure;
	}
	
	/**
	 * Sets the secondary structure annotation for this Pdb, overwriting the existing one. 
	 * @param secondaryStructure
	 */
	public void setSecondaryStructure(SecondaryStructure secondaryStructure) {
		this.secondaryStructure = secondaryStructure;
		initialiseResiduesSecStruct();
	}

	// end of secondary structure related methods
	
	/**
	 * Sets the catalitic site set annotation object of this Pdb
	 * @param catalSiteSet
	 */
	public void setCatalSiteSet(CatalSiteSet catalSiteSet) {
		this.catalSiteSet = catalSiteSet;
	}

	/**
	 * Sets the scop annotation object of this Pdb
	 * @param scop
	 */
	public void setScop(Scop scop) {
		this.scop = scop;
	}
	
	/**
	 * Sets the EC annotation object of this Pdb
	 * @param ec
	 */
	public void setEC(EC ec) {
		this.ec = ec;
	}
	
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
		return rmsd(otherPdb, ct, null);
	}
	
	/**
	 * Calculates rmsd (on atoms given by ct) of this Pdb object to otherPdb object
	 * restricted only to the given set of intervals
	 * Both objects must represent structures with same sequence (save unobserved residues or missing atoms)
	 * 
	 * @param otherPdb
	 * @param ct the contact type (crossed contact types don't make sense here)
	 * @param intervSet an interval set of residues for which the rmsd will be calculated, if 
	 * null then all residues will be considered
	 * @return
	 * @throws ConformationsNotSameSizeError
	 */
	public double rmsd(Pdb otherPdb, String ct, IntervalSet intervSet) throws ConformationsNotSameSizeError {
		TreeMap<Integer, Residue> thisResidues = this.getReducedResidues(ct,intervSet);
		TreeMap<Integer, Residue> otherResidues = otherPdb.getReducedResidues(ct,intervSet);

		ArrayList<Vector3d> conf1AL = new ArrayList<Vector3d>();
		ArrayList<Vector3d> conf2AL = new ArrayList<Vector3d>();	
		// there might be unobserved residues or some missing atoms for a residue
		// here we get the ones that are in common
		for (int resser:thisResidues.keySet()) {
			Residue thisRes = thisResidues.get(resser);
			if (otherResidues.containsKey(resser)) {
				Residue otherRes = otherResidues.get(resser);
				for (Atom atom:thisRes.getAtoms()) {
					if (otherRes.containsAtom(atom.getCode())) {
						conf1AL.add(new Vector3d(atom.getCoords()));
						conf2AL.add(new Vector3d(otherRes.getAtom(atom.getCode()).getCoords()));
					}
				}
			}
		}

		// converting the ArrayLists to arrays
		Vector3d[] conformation1 = new Vector3d[conf1AL.size()]; 
		Vector3d[] conformation2 = new Vector3d[conf2AL.size()];
		conf1AL.toArray(conformation1);
		conf2AL.toArray(conformation2);
		
		// this as well as calculating the rmsd, changes conformation1 and conformation2 to be optimally superposed
		double rmsd = calculate_rmsd(conformation1, conformation2);

//		// printing out individual distances (conformation1 and conformation2 are now optimally superposed)
//		for (i=0;i<conformation1.length;i++){
//		Point3d point1 = new Point3d(conformation1[i].x,conformation1[i].y,conformation1[i].z);
//		Point3d point2 = new Point3d(conformation2[i].x,conformation2[i].y,conformation2[i].z);
//		System.out.println(point1.distance(point2));
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
	public static double calculate_rmsd(Vector3d[] conformation1, Vector3d[] conformation2) throws ConformationsNotSameSizeError{
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
	private static Matrix vector3dAr2matrix(Vector3d[] vecArray) {
		double[][] array = new double[vecArray.length][3];
		for (int i=0;i<vecArray.length;i++){
			vecArray[i].get(array[i]);
		}
		return new Matrix(array);
	}

	/**
	 * Write residue info to given db, using our db graph OWL format, 
	 * i.e. tables: residue_info
	 * @param conn
	 * @param db
	 * @throws SQLException
	 */
	public void writeToDb(MySQLConnection conn, String db) throws SQLException{

		Statement stmt;
		String sql = "";

		conn.setSqlMode("NO_UNSIGNED_SUBTRACTION,TRADITIONAL");

		for (int resser:getAllSortedResSerials()) {
			String resType = String.valueOf(this.getResidue(resser).getAaType().getOneLetterCode());
			String pdbresser = getPdbResSerFromResSer(resser);

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
			" VALUES ("+quote(pdbCode)+", "+quote(chainCode)+", "+(pdbChainCode.equals(NULL_CHAIN_CODE)?quote("-"):quote(pdbChainCode))+","+resser+", "+quote(pdbresser)+", "+quote(resType)+", "+secStructType+", "+secStructId+", "+scopId+", "+sccs+", "+sunid+", "+orderIn+", "+domainType+", "+domainNumReg+", "+allRsa+", "+scRsa+", "+consurfhsspScore+","+consurfhsspColor+","+ecId+","+csaNums+","+csaChemFuncs+","+csaEvids+")";
			//System.out.println(sql);
			stmt = conn.createStatement();
			stmt.executeUpdate(sql);
			stmt.close();
		}			
	}

	/**
	 * Write residue info to given db, using our db graph OWL format, 
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
		
		for (int resser:getAllSortedResSerials()) {
			String resType = String.valueOf(this.getResidue(resser).getAaType().getOneLetterCode());
			String pdbresser = getPdbResSerFromResSer(resser);
			
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
			
			resOut.println(pdbCode+"\t"+chainCode+"\t"+(pdbChainCode.equals(NULL_CHAIN_CODE)?"-":pdbChainCode)+"\t"+resser+"\t"+pdbresser+"\t"+resType+"\t"+secStructType+"\t"+secStructId+"\t"+scopId+"\t"+sccs+"\t"+sunid+"\t"+orderIn+"\t"+domainType+"\t"+domainNumReg+"\t"+allRsa+"\t"+scRsa+"\t"+consurfhsspScore+"\t"+consurfhsspColor+"\t"+ecId+"\t"+csaNums+"\t"+csaChemFuncs+"\t"+csaEvids);
			
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
	
	/**
	 * Restricts thisPdb object to residues that belong to the given sunid
	 * Can only be used if SCOP annotation is loaded. Check it with {@link #hasScop()}
	 * @param sunid
	 */
	public void restrictToScopDomain (int sunid) {
		Vector<ScopRegion> scopRegions = this.scop.getScopRegions(sunid);
		if (scopRegions.size()!=0) {
			this.sid = scopRegions.get(0).getSId();
			if (scopRegions.get(0).getDomainType() != ScopRegion.DomainType.WHOLECHAIN) {
				restrictToScopRegions(scopRegions);
				
				Iterator<ScopRegion> allScopRegions = this.scop.getIterator();
				while (allScopRegions.hasNext()) {
					ScopRegion scopRegion = allScopRegions.next();
					if (!scopRegion.getSId().equals(sid)) {
						allScopRegions.remove();
					}					
				}
			}
		}
	}
	
	/**
	 * Restricts thisPdb object to residues that belong to the given sid
	 * Can only be used if SCOP annotation is loaded. Check it with {@link #hasScop()} 
	 * @param sid
	 */
	public void restrictToScopDomain (String sid) {

		Vector<ScopRegion> scopRegions = this.scop.getScopRegions(sid);
		if (scopRegions.size()!=0) {
			this.sid = sid;
			if (scopRegions.get(0).getDomainType() != ScopRegion.DomainType.WHOLECHAIN) {
				restrictToScopRegions(scopRegions);
				
				Iterator<ScopRegion> allScopRegions = this.scop.getIterator();
				while (allScopRegions.hasNext()) {
					ScopRegion scopRegion = allScopRegions.next();
					if (!scopRegion.getSId().equals(sid)) {
						allScopRegions.remove();
					}					
				}
			}
		}
	}
	
	/**
	 * Restricts thisPdb object to residues that belong to the given sids
	 * Can only be used after calling checkScop() 
	 * @param sid
	 */
	public void restrictToScopDomains (String[] sids) {
		
		TreeSet<String> keepSids = new TreeSet<String>();
		
		Vector<ScopRegion> scopRegions = new Vector<ScopRegion>();
		for (String sid: sids) {
			if (!keepSids.contains(sid)) {
				Vector<ScopRegion> curScopRegions = this.scop.getScopRegions(sid);
				if (curScopRegions.size() != 0) {
					scopRegions.addAll(curScopRegions);
					keepSids.add(sid);
				}
			}
		}
		
		if (scopRegions.size()!=0) {
			Iterator<String> it = keepSids.iterator();
			this.sid = it.next();
			while (it.hasNext()) {
				this.sid += ","+it.next(); 
			}
			
			if (scopRegions.get(0).getDomainType() != ScopRegion.DomainType.WHOLECHAIN) {
				restrictToScopRegions(scopRegions);
				
				Iterator<ScopRegion> allScopRegions = this.scop.getIterator();
				while (allScopRegions.hasNext()) {
					ScopRegion scopRegion = allScopRegions.next();
					if (!keepSids.contains(scopRegion.getSId())) {
						allScopRegions.remove();
					}					
				}
			}
		}
	}
	
	/**
	 * Restricts thisPdb object to residues that belong to the given sunids
	 * Can only be used after calling checkScop() 
	 * @param sid
	 */
	public void restrictToScopDomains (int[] sunids) {
		
		TreeSet<String> keepSids = new TreeSet<String>();
		TreeSet<Integer> keepSunids = new TreeSet<Integer>();
		
		Vector<ScopRegion> scopRegions = new Vector<ScopRegion>();
		for (int sunid: sunids) {
			if (!keepSunids.contains(sunid)) {
				Vector<ScopRegion> curScopRegions = this.scop.getScopRegions(sunid);
				if (curScopRegions.size() != 0) {
					scopRegions.addAll(curScopRegions);
					keepSunids.add(sunid);
					keepSids.add(curScopRegions.get(0).getSId());
				}
			}
		}
		
		if (scopRegions.size()!=0) {
			Iterator<String> it = keepSids.iterator();
			this.sid = it.next();
			while (it.hasNext()) {
				this.sid += ","+it.next(); 
			}
			
			if (scopRegions.get(0).getDomainType() != ScopRegion.DomainType.WHOLECHAIN) {
				restrictToScopRegions(scopRegions);
				
				Iterator<ScopRegion> allScopRegions = this.scop.getIterator();
				while (allScopRegions.hasNext()) {
					ScopRegion scopRegion = allScopRegions.next();
					if (!keepSids.contains(scopRegion.getSId())) {
						allScopRegions.remove();
					}					
				}
			}
		}
	}
	
	/**
	 * Restricts this Pdb object to residues within the given ScopRegions 
	 * @param scopRegions
	 */
	private void restrictToScopRegions (Vector<ScopRegion> scopRegions) {
		IntervalSet intervSet = new IntervalSet();
		Iterator<ScopRegion> it = scopRegions.iterator();
		while(it.hasNext()) {
			ScopRegion scopRegion = it.next();
			intervSet.add(scopRegion.getInterval());
		}
		restrictToIntervalSet(intervSet);
	}
	
	/**
	 * Restricts this Pdb object to residues within the given IntervalSet
	 * @param intervSet a set of internal residue serials
	 */
	public void restrictToIntervalSet(IntervalSet intervSet) {
		
		// getting list of the residue serials to keep
		TreeSet<Integer> resSersToKeep = intervSet.getIntegerSet();

		// removing residues
		Iterator<Residue> resIt = residues.values().iterator();
		while (resIt.hasNext()) {
			int resser = resIt.next().getSerial();
			if (!resSersToKeep.contains(resser)) {
				resIt.remove();
				if (catalSiteSet != null) {
					catalSiteSet.removeCatalSiteRes(resser);
				}
			}
		}
		// reinitialise the maps
		initialiseMaps();
			
		// setting sequence to scop sequence and obsLength and fullLength respectively
		Iterator<Interval> regionsToKeep = intervSet.iterator();
		String newSequence = "";
		while (regionsToKeep.hasNext()) {
			Interval region = regionsToKeep.next();
			newSequence += sequence.substring((region.beg-1),region.end);
		}
		sequence = newSequence;
		fullLength = sequence.length();
		
	}
	
	/**
	 * Mirror this Pdb structure by inverting through the origin.
	 */
	public void mirror() {
		for (int atomserial:getAllAtomSerials()){
			Point3d coords = getAtomCoord(atomserial);
			coords.x *= -1;
			coords.y *= -1;
			coords.z *= -1;
		}
	}

	/**
	 * Rotates or translates this structure.
	 * @param m the rotation/translation matrix
	 */
	public void transform(Matrix4d m) {
		for(int atomserial:getAllAtomSerials()) {
			Point3d coords = getAtomCoord(atomserial);
			m.transform(coords);
		}
	}
	
	/**
	 * Transform the coordinates of this structure translating them to the given center
	 * and rotating them so that the given axis aligns with the z-axis
	 * @param center
	 * @param axis
	 */
	public void transformToCenterAndAxis(Point3d center, Vector3d axis) {
		// finding the rotation matrix to align z axis to the given inertia axis
		Vector3d r = new Vector3d();
		Vector3d k = new Vector3d(0,0,1);
		r.cross(axis, k); // this is the axis of rotation
		double alpha = axis.angle(k); // this is the angle to rotate
		AxisAngle4d axisAngle = new AxisAngle4d(r, alpha);
		// note that the matrix needs to be initialised to the unit matrix otherwise setRotation() doesn't work properly
		Matrix4d rot = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1); 
		rot.setRotation(axisAngle);
		for(int atomserial:getAllAtomSerials()) {
			Point3d coords = getAtomCoord(atomserial);
			// translate to new origin
			coords.sub(center);
			// rotate so that z axis is the given axis
			rot.transform(coords);
		}
	}
	
	/**
	 * Rotates this structure around rotAxis with the given rotAngle 
	 * @param rotAxis the vector around which the rotation will be performed
	 * @param rotAngle the rotation angle in radians
	 */
	public void rotate(Vector3d rotAxis, double rotAngle) {
		AxisAngle4d axisAngle = new AxisAngle4d(rotAxis, rotAngle);
		// note that the matrix needs to be initialised to the unit matrix otherwise setRotation() doesn't work properly
		Matrix4d rot = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1); 
		rot.setRotation(axisAngle);
		for(int atomserial:getAllAtomSerials()) {
			Point3d coords = getAtomCoord(atomserial);
			// rotate
			rot.transform(coords);			
		}			
	}

	
	/**
	 * Moves this structure such that the center of mass is at the origin using all atoms
	 */
	public void moveToOrigin() {
		Vector3d sumVector = new Vector3d();
		int numVectors = 0;
		for(int atomserial:getAllAtomSerials()) {
			Point3d coords = getAtomCoord(atomserial);
			sumVector.add(coords);
			numVectors++;
		}
		sumVector.scale(1.0/numVectors);
		//System.out.println(sumVector);
		
		for(int atomserial:getAllAtomSerials()) {
			Point3d coords = getAtomCoord(atomserial);
			coords.sub(sumVector);
		}		
	}
	
	/**
	 * Moves this structure such that the center of mass of the given subset of residues is at the origin using only CA atoms
	 */
	public void moveToOrigin(TreeSet<Integer> residues) {
		Vector3d sumVector = new Vector3d();
		int numVectors = 0;
		for(int resser:residues) {
			Point3d coords = getAtomCoord(resser, "CA");
			sumVector.add(coords);
			numVectors++;			
		}

		sumVector.scale(1.0/numVectors);
		//System.out.println(sumVector);
		
		for(int atomserial:getAllAtomSerials()) {
			Point3d coords = getAtomCoord(atomserial);
			coords.sub(sumVector);
		}		
	}
	
	/**
	 * Gets the phi angle in degrees for given residue serial
	 * @param i
	 * @return the phi angle or NaN if there are no coordinates for given i or i-1
	 */
	public double getPhiAngle(int i) { 
		if (!hasCoordinates(i-1, "C") || !hasCoordinates(i, "N") || !hasCoordinates(i, "CA") || !hasCoordinates(i, "C")) {
			return Double.NaN;
		}
		Point3d Ciminus1 = getAtomCoord(i-1, "C");
		Point3d Ni = getAtomCoord(i, "N");
		Point3d CAi = getAtomCoord(i, "CA");
		Point3d Ci = getAtomCoord(i, "C");
		return getTorsionAngle(Ciminus1, Ni, CAi, Ci);	
	}
	
	/**
	 * Gets the psi angle in degrees for the given residue serial
	 * @param i
	 * @return the psi angle or NaN if there are no coordinates for given i or i+1
	 */
	public double getPsiAngle(int i) {
		if (!hasCoordinates(i, "N") || !hasCoordinates(i, "CA") || !hasCoordinates(i, "C") || !hasCoordinates(i+1,"N")) {
			return Double.NaN;
		}
		Point3d Ni = getAtomCoord(i, "N");
		Point3d CAi = getAtomCoord(i, "CA");
		Point3d Ci = getAtomCoord(i, "C");
		Point3d Niplus1 = getAtomCoord(i+1, "N");
		return getTorsionAngle(Ni, CAi, Ci, Niplus1);
	}
	
	/**
	 * Gets the omega angle in degrees for the given residue serial
	 * @param i
	 * @return the omega angle or NaN if there are no coordinates for given i or i+1
	 */
	public double getOmegaAngle(int i) {
		if (!hasCoordinates(i, "CA") || !hasCoordinates(i, "C") || !hasCoordinates(i+1, "N") || !hasCoordinates(i+1,"CA")) {
			return Double.NaN;
		}
		Point3d CAi = getAtomCoord(i, "CA");
		Point3d Ci = getAtomCoord(i, "C");
		Point3d Niplus1 = getAtomCoord(i+1, "N");
		Point3d CAiplus1 = getAtomCoord(i+1, "CA");
		return getTorsionAngle(CAi, Ci, Niplus1, CAiplus1);
	}
	
	/**
	 * Gets the torsion angle in degrees for the 4 given atoms 
	 * See http://en.wikipedia.org/wiki/Dihedral_angle for how to calculate it.
	 * @param atom1
	 * @param atom2
	 * @param atom3
	 * @param atom4
	 * @return
	 */
	private double getTorsionAngle(Point3d atom1, Point3d atom2, Point3d atom3, Point3d atom4) {
		Vector3d b1 = new Vector3d();
		Vector3d b2 = new Vector3d();
		Vector3d b3 = new Vector3d();
		Vector3d b1xb2 = new Vector3d();
		Vector3d b2xb3 = new Vector3d();
		Vector3d b1scalelengthb2 = new Vector3d();
		b1.sub(atom2, atom1);
		b2.sub(atom3, atom2);
		b3.sub(atom4, atom3);
		b1xb2.cross(b1, b2);
		b2xb3.cross(b2, b3);
		b1scalelengthb2.scale(b2.length(), b1);
		return Math.toDegrees(Math.atan2(b1scalelengthb2.dot(b2xb3), b1xb2.dot(b2xb3)));		
	}
	
	/**
	 * Gets all phi/psi angles in degrees for this Pdb structure 
	 * @return residue serials as keys, values arrays of 2 doubles: first phi angle, second psi angle 
	 */
	public TreeMap<Integer, double[]> getAllPhiPsi() {
		TreeMap<Integer, double[]> phipsi = new TreeMap<Integer, double[]>();
		for (int resser: getAllSortedResSerials()) {
			double[] angles = {getPhiAngle(resser), getPsiAngle(resser)};
			phipsi.put(resser, angles);
		}
		return phipsi;
	}
	
	/**
	 * Prints in 2 columns all phi, psi angles for this protein structure.
	 * Use to generate a Ramachandran plot.
	 */
	public void printAllPhiPsi() {
		TreeMap<Integer, double[]> phipsi = getAllPhiPsi();
		for (int resser:phipsi.keySet()) {
			System.out.printf("%7.2f %7.2f\n",phipsi.get(resser)[0],phipsi.get(resser)[1]);
		}
	}
	
	/**
	 * Returns the number of unobserved residues
	 * @return
	 */
	public int countUnobserved() {
		return getFullLength()-getObsLength();
	}
	
	public TreeMap<Integer,Character> getUnobservedResidues() {
		TreeMap<Integer,Character> unobserved = new TreeMap<Integer,Character>();

		// detect all unobserved residues
		for(int i = 1; i <= fullLength; ++i) {
			if(!residues.containsKey(i)) {
				unobserved.put(i,sequence.charAt(i-1));
			}
		}
		return unobserved;
		
	}

	
	/*---------------------------- static methods ---------------------------*/
	
	/**
	 * Loads a pdb structure where arg can be a pdbcode+chaincode or a pdb file name.
	 * If something goes wrong, prints an error message and exits.
	 * @param arg a pdbcode+chaincode (e.g. 1tdrB) or a pdb file name
	 * @return the newly created pdb object
	 */
	public static Pdb readStructureOrExit(String arg) {
		return readFromFileOrPdbCode(arg, true, true);
	}
	
	/**
	 * Loads a pdb structure where arg can be a pdbcode+chaincode or a pdb file name.
	 * If something goes wrong, prints an error message and returns null;
	 * @param arg a pdbcode+chaincode (e.g. 1tdrB) or a pdb file name
	 * @return the newly created pdb object or null
	 */	
	public static Pdb readStructureOrNull(String arg) {
		return readFromFileOrPdbCode(arg, false, true);		
	}
	
	public static Pdb readStructureOrNull(String arg, String chain) {
		return readFromFileOrPdbCode(arg, chain, true, true);		
	}	

	/**
	 * Loads a pdb structure given a pdbcode+chaincode or a pdb file name and chain code.
	 * Common exceptions are caught internally. The behvaiour in case of an error is
	 * specified by the parameters <code>exit</code> and <code>silent</code>.
	 * Parameter chain code is only used if first parameter is a file, otherwise ignored.
	 * @param arg a pdbcode+chaincode (e.g. 1tdrB) or a pdb file name
	 * @param chain a chain code or null (first chain code in file)
	 * @param exit if true, system.exit(1) is called on error
	 * @param silent if false, error messages will be printed
	 * @return the structure object
	 */
	public static Pdb readFromFileOrPdbCode(String arg, String chain, boolean exit, boolean silent) {
		Pdb pdb = null;
		
		// check if argument is a filename
		File inFile = new File(arg);
		if(inFile.canRead()) {
			if(!silent) System.out.println("Reading file " + arg);
			pdb = new PdbfilePdb(arg);
			try {
				if(chain == null) {
					String[] chains = pdb.getChains();
					if(!silent) System.out.println("Loading chain " + chains[0]);
					pdb.load(chains[0]);
				} else {
					if(!silent) System.out.println("Loading chain " + chain);
					pdb.load(chain);
				}
			} catch (PdbLoadError e) {
				if(!silent) System.err.println("Error loading file " + arg + ":" + e.getMessage());
			}
		} else {
			// check if argument is a pdb code
			if(arg.length() < 4 || arg.length() > 5) {
				if(!silent) System.err.println(arg + "is neither a valid file name nor a valid pdb(+chain) code");
				if(exit) System.exit(1);
			} else {
				String pdbCode = arg.substring(0,4);
				String chainCode = arg.substring(4,5);
				try {
					if(!silent) System.out.println("Loading pdb code " + pdbCode);
					pdb = new PdbasePdb(pdbCode);
					if(chainCode.length() == 0) {
						try {
							chainCode = pdb.getChains()[0];
						} catch (PdbLoadError e) {
							if(!silent) System.err.println("Error loading pdb structure:" + e.getMessage());
							if(exit) System.exit(1);
						}
					}
					try {
						if(!silent) System.out.println("Loading chain " + chainCode);
						pdb.load(chainCode);
					} catch (PdbLoadError e) {
						if(!silent) System.err.println("Error loading pdb structure:" + e.getMessage());
						if(exit) System.exit(1);
					}
				} catch (SQLException e) {
					if(!silent) System.err.println("Database error: " + e.getMessage());
					if(exit) System.exit(1);
				} catch (PdbCodeNotFoundError e) {
					if(!silent) System.err.println("Pdb code " + pdbCode + " not found in database.");
					if(exit) System.exit(1);
				}
			}
		}
		return pdb;	
	}
	
	
	/**
	 * Loads a pdb structure given a pdbcode+chaincode or a pdb file name.
	 * Common exceptions are caught internally. The behvaiour in case of an error is
	 * specified by the parameters <code>exit</code> and <code>silent</code>.
	 * @param arg a pdbcode+chaincode (e.g. 1tdrB) or a pdb file name
	 * @param exit if true, system.exit(1) is called on error
	 * @param silent if false, error messages will be printed
	 * @return the structure object
	 */
	public static Pdb readFromFileOrPdbCode(String arg, boolean exit, boolean silent) {
		return readFromFileOrPdbCode(arg, null, exit, silent);
	}
	
	/*------------------------ HasFeature interface implementation -----------------------*/
	
	@Override
	public boolean addFeature(Feature feature) throws InvalidFeatureCoordinatesException,
														OverlappingFeatureException {

		boolean result = false;
		FeatureType ft = feature.getType();	// will throw a null pointer exception if feature == null
		if(this.features == null) {
			this.features = new HashMap<FeatureType, Collection<Feature>>();
		}
		Collection<Feature> fc = this.features.get(ft);
		if(fc==null) {
			fc = new LinkedList<Feature>();
			features.put(ft, fc);
		}
		IntervalSet intervSet = feature.getIntervalSet();
		for (Interval interv:intervSet){
			if (interv.beg<1 || interv.end>this.getFullLength())
				throw new InvalidFeatureCoordinatesException("Feature being added "+feature.getDescription()+" of type "+feature.getType()+" contains invalid coordinates for this Pdb.\n"
						+"Interval: "+interv+". Max residue serial for this Pdb: "+this.getFullLength());
		}
		for (Feature f:fc) {
			if (intervSet.overlaps(f.getIntervalSet())) {
				throw new OverlappingFeatureException("Feature being added "+feature.getDescription()+" of type "+feature.getType()+" overlaps existing feature "+f.getDescription()+" of type "+f.getType()+"\n" +
						"New interval set: "+intervSet+". Existing interval set: "+f.getIntervalSet());
			}
		}
		result = fc.add(feature);
		return result;	
	}

	@Override
	public Collection<FeatureType> getFeatureTypes() {
		return features.keySet();
	}

	@Override
	public Collection<Feature> getFeatures() {
		Collection<Feature> allfeatures = new LinkedList<Feature>();
		for (Collection<Feature> coll:features.values()) {
			allfeatures.addAll(coll);
		}
		return allfeatures;
	}

	@Override
	public Collection<Feature> getFeaturesForPositon(int position) {
		Collection<Feature> result = new LinkedList<Feature>(); 
		for(Feature f:this.getFeatures()) {
			if(f.getIntervalSet().getIntegerSet().contains(position)) result.add(f);
		}
		return result;		
	}

	@Override
	public Collection<Feature> getFeaturesOfType(FeatureType featureType) {
		return features.get(featureType);
	}

	@Override
	public Collection<Feature> getFeaturesOfTypeForPosition(FeatureType featureType, int position) {
		Collection<Feature> result = new LinkedList<Feature>(); 
		if(this.features.get(featureType) != null) {
			for(Feature f:this.features.get(featureType)) {
				if(f.getIntervalSet().getIntegerSet().contains(position)) result.add(f);
			}
		}
		return result;
	}
}


package owl.core.structure;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.net.URL;
import java.net.URLConnection;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import owl.core.util.BoundingBox;
import owl.core.util.FileFormatException;
import owl.core.util.FileTypeGuesser;
import owl.core.util.MySQLConnection;

/**
 * A protein crystal's asymmetric unit, i.e. a PDB entry.
 *  
 * @author duarte_j
 *
 */
public class PdbAsymUnit implements Serializable { //, Iterable<PdbChain>
	
	/*------------------------------------  constants ---------------------------------------------*/
	
	private static final long serialVersionUID = 1L;

	public static final int    DEFAULT_MODEL   = 1;			// default model serial (NMR structures)
	
	// to specify the NULL (blank in pdb file) chain code
	// NOTE: the value of this constant used to be "NULL"
	public static final String NULL_CHAIN_CODE = "A";	
	
	public static final String NO_PDB_CODE       = "";		// to specify no pdb code


	
	protected static final Matrix4d IDENTITY_TRANSFORM = new Matrix4d(1,0,0,0,
																	0,1,0,0,
																	0,0,1,0,
																	0,0,0,1);
	/*------------------------------------  members -----------------------------------------------*/
	private String pdbCode;
	private int model;
	private String title;
	
	/**
	 * Polymer chains (peptide/nucleotide): chainCodes to PdbChains
	 */
	private TreeMap<String, PdbChain> chains;		
	
	/**
	 * Non-polymer chains (ligands and other het residues): chainCodes to PdbChains
	 */
	private TreeMap<String, PdbChain> nonPolyChains;
	
	/**
	 * A map of pdbChainCodes to chainCodes of polymer chains. 
	 */
	private TreeMap<String, String> pdbchaincode2chaincode;
	
	/**
	 * A map of all PDB chain codes to the representative pdbChainCode
	 * (the first chain alphabetically from the group of identical chains)
	 */
	private HashMap<String, String> chain2repChain;
	
	/**
	 * A map of representative PDB chain code to the members of the 
	 * group of identical chains
	 */
	private HashMap<String,List<String>> repChain2members; 
	
	/**
	 * Experimental method: x-ray crystallograpy, NMR, etc.
	 */
	private String expMethod;
	
	// crystallographic data
	
	/**
	 * The parameters of the crystal cell (in CRYST1 field in PDB files)
	 */
	private CrystalCell crystalCell;
	
	/**
	 * The space group (in CRYST1 field in PDB files)
	 */
	private SpaceGroup spaceGroup;
	
	// quality data (for crystal structures), all 3 parameters initialised to -1 (for non-crystal structures they just don't apply)
	private double resolution;					
	private double rFree;
	private double rSym;
	private BoundingBox bounds; // cached bounds for speed up of interface calculations, filled in getAllAtoms()
	private boolean boundsProtOnly; // to track whether the cached bounds are for protein chains only (true) or for any chain poly/non-poly
	private Point3d centroid;  // cached centroid
	
	/**
	 * Transformation matrix used to generate this asym unit (includes 
	 * the possible translation from the original cell)
	 */
	private Matrix4d transform; 
	
	/**
	 * An identifier of the space group transformation used (0 is identity i.e. original asymmetric unit) 
	 * does not count the translations: 2 equivalent asym units of 2 different cells will have the same identifier)
	 * i.e. it is unique within the unit cell but equivalent units of different crystal cells will have same id
	 * goes from 1 to m (m=number of symmetry operations of the space group)
	 */
	private int transformId;    


	
	/*------------------------------------  constructors -----------------------------------------------*/

	/**
	 * Constructs an empty PdbAsymUnit with default meta-data and no chains.
	 */
	public PdbAsymUnit() {
		this.expMethod = null;
		this.resolution = -1;
		this.rFree = -1;
		this.rSym = -1;
		this.pdbCode = NO_PDB_CODE;
		this.model = DEFAULT_MODEL;
		this.title = null;
		
		this.chains = new TreeMap<String,PdbChain>();
		this.nonPolyChains = new TreeMap<String,PdbChain>();
		this.transform = IDENTITY_TRANSFORM;
		this.transformId = 0;

	}
	
	/**
	 * Constructs a PdbAsymUnit by reading model {@value #DEFAULT_MODEL} for all chains 
	 * from given pdbSourceFile. 
	 * The pdbSourceFile can be in PDB or mmCIF format.
	 * @param pdbSourceFile
	 * @throws IOException if IO problems while reading PDB data from file
	 * @throws FileFormatException if given file is not in PDB or mmCIF format
	 * @throws PdbLoadException if problems while loading the PDB data from file
	 */
	public PdbAsymUnit(File pdbSourceFile) throws IOException, FileFormatException, PdbLoadException {
		this(pdbSourceFile, DEFAULT_MODEL);
	}
	
	/**
	 * Constructs a PdbAsymUnit by reading given model for all chains 
	 * from given pdbSourceFile. 
	 * The pdbSourceFile can be in PDB or mmCIF format. 
	 * @param pdbSourceFile
	 * @param model
	 * @throws IOException if IO problems while reading PDB data from file
	 * @throws FileFormatException if given file is not in PDB or mmCIF format
	 * @throws PdbLoadException if problems while loading the PDB data from file
	 */
	public PdbAsymUnit(File pdbSourceFile, int model) throws IOException, FileFormatException, PdbLoadException {
		this.expMethod = null;
		this.resolution = -1;
		this.rFree = -1;
		this.rSym = -1;
		this.pdbCode = NO_PDB_CODE;
		this.model = model;
		this.title = null;
		this.transform = IDENTITY_TRANSFORM;
		this.transformId = 0;
		this.chains = new TreeMap<String, PdbChain>();
		this.nonPolyChains = new TreeMap<String,PdbChain>();
		int type = FileTypeGuesser.guessFileType(pdbSourceFile);
		if (type==FileTypeGuesser.PDB_FILE || type ==FileTypeGuesser.RAW_PDB_FILE || type==FileTypeGuesser.CASP_TS_FILE) {
			loadFromPdbFile(pdbSourceFile, true);
		} else if (type==FileTypeGuesser.CIF_FILE) {
			loadFromCifFile(pdbSourceFile);
		} else {
			throw new FileFormatException("The given file does not seem to be neither a PDB file nor a mmCIF file");
		}
		
	}
	
	/**
	 * Constructs a PdbAsymUnit by reading given model for all chains 
	 * from given pdbSourceFile. 
	 * The pdbSourceFile can be in PDB or mmCIF format. 
	 * @param pdbSourceFile
	 * @param model
	 * @param missingSeqResPadding if true in a PDB file without SEQRES, the sequence will be padded with X 
	 * for any missing residue (deduced from residue numbering) in ATOM line. If false no padding will be
	 * done if no SEQRES found 
	 * @throws IOException if IO problems while reading PDB data from file
	 * @throws FileFormatException if given file is not in PDB or mmCIF format
	 * @throws PdbLoadException if problems while loading the PDB data from file
	 */
	public PdbAsymUnit(File pdbSourceFile, int model, boolean missingSeqResPadding) throws IOException, FileFormatException, PdbLoadException {
		this.expMethod = null;
		this.resolution = -1;
		this.rFree = -1;
		this.rSym = -1;
		this.pdbCode = NO_PDB_CODE;
		this.model = model;
		this.title = null;
		this.transform = IDENTITY_TRANSFORM;
		this.transformId = 0;
		this.chains = new TreeMap<String, PdbChain>();
		this.nonPolyChains = new TreeMap<String,PdbChain>();
		int type = FileTypeGuesser.guessFileType(pdbSourceFile);
		if (type==FileTypeGuesser.PDB_FILE || type ==FileTypeGuesser.RAW_PDB_FILE || type==FileTypeGuesser.CASP_TS_FILE) {
			loadFromPdbFile(pdbSourceFile, missingSeqResPadding);
		} else if (type==FileTypeGuesser.CIF_FILE) {
			loadFromCifFile(pdbSourceFile);
		} else {
			throw new FileFormatException("The given file does not seem to be neither a PDB file nor a mmCIF file");
		}
		
	}
	
	/**
	 * Constructs a PdbAsymUnit by reading model {@value #DEFAULT_MODEL} for all chains 
	 * from given pdbase database. 
	 * @param pdbCode
	 * @param conn
	 * @param dbName
	 * @throws PdbLoadException
	 * @throws PdbCodeNotFoundException
	 */
	public PdbAsymUnit(String pdbCode, MySQLConnection conn, String dbName) throws PdbLoadException, PdbCodeNotFoundException {
		this(pdbCode, DEFAULT_MODEL, conn, dbName);
	}
	
	/**
	 * Constructs a PdbAsymUnit by reading given model for all chains 
	 * from given pdbase database. 
	 * @param pdbCode
	 * @param model
	 * @param conn
	 * @param dbName
	 * @throws PdbLoadException
	 * @throws PdbCodeNotFoundException
	 */
	public PdbAsymUnit(String pdbCode, int model, MySQLConnection conn, String dbName) throws PdbLoadException, PdbCodeNotFoundException {
		this.pdbCode = pdbCode.toLowerCase();
		this.model = model;

		this.expMethod = null;
		this.resolution = -1;
		this.rFree = -1;
		this.rSym = -1;
		this.title = null;
		this.transform = IDENTITY_TRANSFORM;
		this.transformId = 0;
		
		this.chains = new TreeMap<String, PdbChain>();
		this.nonPolyChains = new TreeMap<String,PdbChain>();

		loadFromPdbase(conn, dbName);
	}
	
	/*----------------------------------  loader methods ----------------------------------------*/
	
	private void loadFromPdbFile(File pdbFile, boolean missingSeqResPadding) throws PdbLoadException {
		PdbfileParser parser = new PdbfileParser(pdbFile.getAbsolutePath(),missingSeqResPadding);
		
		parser.readChains(this,model);
		
		this.pdbCode = parser.getPdbCode();
		this.expMethod = parser.getExpMethod();
		this.title = parser.getTitle();
		this.crystalCell = parser.getCrystalCell();
		this.spaceGroup = parser.getSpaceGroup();
		this.resolution = parser.getResolution();
		this.rFree = parser.getRfree();
		this.rSym = parser.getRsym();
		
		// check for those very weird few PDB entries that are in a non-standard frame (e.g. 1bbb, 1bab)
		Matrix4d scaleMatrix = parser.getScaleMatrix();
		if (isXrayDiffraction() && scaleMatrix!=null && this.crystalCell!=null && !this.crystalCell.checkScaleMatrix(scaleMatrix)) {
			throw new PdbLoadException("This PDB entry is in a non-standard crystal frame. Please tell the PDB about it.");
		}

	}
	
	private void loadFromCifFile(File cifFile) throws PdbLoadException, IOException, FileFormatException {
		CiffileParser parser = new CiffileParser(cifFile);
		this.pdbCode = parser.readPdbCode();
		this.expMethod = parser.readExpMethod();
		this.title = parser.readTitle();
		this.crystalCell = parser.readCrystalCell();
		this.spaceGroup = parser.readSpaceGroup();
		double[] qParams = parser.readQparams();
		this.resolution = qParams[0];
		this.rFree = qParams[1];
		this.rSym = qParams[2];
		
		parser.readChains(this, model);
		for (PdbChain chain:getAllChains()) {
			chain.setParent(this);
		}

		// check for those very weird few PDB entries that are in a non-standard frame (e.g. 1bbb, 1bab)
		Matrix4d scaleMatrix = parser.readScaleMatrix();
		if (isXrayDiffraction() && scaleMatrix!=null && this.crystalCell!=null && !this.crystalCell.checkScaleMatrix(scaleMatrix)) {
			throw new PdbLoadException("This PDB entry is in a non-standard crystal frame. Please tell the PDB about it.");
		}
		
		parser.closeFile();		
	}
	
	private void loadFromPdbase(MySQLConnection conn, String dbName) throws PdbLoadException, PdbCodeNotFoundException {
		try {
			PdbaseParser parser = new PdbaseParser(pdbCode,dbName,conn);
			this.title = parser.readTitle();
			this.expMethod = parser.readExpMethod();
			this.crystalCell = parser.readCrystalCell();
			this.spaceGroup = parser.readSpaceGroup();
			double[] qParams = parser.readQparams();
			this.resolution = qParams[0];
			this.rFree = qParams[1];
			this.rSym = qParams[2];

			parser.readChains(this, model);
			for (PdbChain chain:getAllChains()) {
				chain.setParent(this);
			}
				
			// check for those very weird few PDB entries that are in a non-standard frame (e.g. 1bbb,1bab)
			Matrix4d scaleMatrix = parser.readScaleMatrix();
			if (isXrayDiffraction() && scaleMatrix!=null && this.crystalCell!=null && !this.crystalCell.checkScaleMatrix(scaleMatrix)) {
				throw new PdbLoadException("This PDB entry is in a non-standard crystal frame. Please tell the PDB about it.");
			}
			
		} catch(SQLException e) {
			throw new PdbLoadException(e);
		}
	}
	
	/*----------------------------------  public methods ----------------------------------------*/
	
	/**
	 * Gets a polymer chain given its PDB chain code or null if no such chain exists.
	 * @return
	 */
	public PdbChain getChain(String pdbChainCode) {
		String chainCode = getChainCodeForPdbChainCode(pdbChainCode);
		if (chainCode==null) {
			return null;
		} 
		return chains.get(chainCode);
	}
	
	/**
	 * Gets a chain (polymer/non-polymer) given its CIF chain code or null if no such chain exists.
	 * @param chainCode
	 * @return
	 */
	public PdbChain getChainForChainCode(String chainCode) {
		if (chains.containsKey(chainCode)) return chains.get(chainCode);
		if (nonPolyChains.containsKey(chainCode)) return nonPolyChains.get(chainCode);
		return null;
	}

	/**
	 * Returns the first polymer chain (in ascending alphabetical order of the PDB chain codes)
	 * of this PDB entry.
	 * @return
	 */
	public PdbChain getFirstChain() {
		String chainCode = pdbchaincode2chaincode.firstEntry().getValue();
		return chains.get(chainCode);
	}
	
	/**
	 * Returns a Collection of all polymer (peptide/nucleotide) chains
	 * @return
	 */
	public Collection<PdbChain> getPolyChains() {
		return chains.values();
	}
	
	/**
	 * Returns a Collection of all protein chains (no nucleotide chains)
	 * @return
	 */
	protected Collection<PdbChain> getProtChains() {
		Collection<PdbChain> protChains = new ArrayList<PdbChain>();
		for (PdbChain chain:chains.values()){
			if (chain.getSequence().isProtein()) protChains.add(chain);
		}
		return protChains;
	}
	
	/**
	 * Returns a Collection of all non-polymer (ligands and het residues) chains
	 * @return
	 */
	public Collection<PdbChain> getNonPolyChains() {
		return nonPolyChains.values();
	}
	
	/**
	 * Returns a Collection of all polymer and non-polymer chains of this asym unit
	 * @return
	 */
	public Collection<PdbChain> getAllChains() {
		Collection<PdbChain> all = new ArrayList<PdbChain>();
		all.addAll(getPolyChains());
		all.addAll(getNonPolyChains());
		return all;
	}
	
	/**
	 * Returns a sorted set of all CIF chain codes in this entry
	 * @return
	 */
	public Set<String> getChainCodes() {
		Set<String> all = new TreeSet<String>();
		all.addAll(chains.keySet());
		all.addAll(nonPolyChains.keySet());
		return all;
	}
	
	/**
	 * Returns a sorted set of all CIF chain codes of polymer chains (protein/nucleic acid) in this entry
	 * @return
	 */
	@SuppressWarnings("unused")
	private Set<String> getPolyChainCodes() {
		Set<String> all = new TreeSet<String>();
		all.addAll(chains.keySet());
		return all;		
	}
	
	/**
	 * Returns a sorted set of all CIF chain codes of protein chains in this entry (no nucleic acid chains)
	 * @return
	 */
	protected Set<String> getProtChainCodes() {
		Set<String> all = new TreeSet<String>();
		for (PdbChain polyChain:chains.values()) {
			if (polyChain.getSequence().isProtein()) all.add(polyChain.getChainCode());
		}
		return all;				
	}
	
	/**
	 * Returns a sorted set of all PDB chain codes of polymer chains in this entry
	 * @return
	 */
	public Set<String> getPdbChainCodes() {
		Set<String> all = new TreeSet<String>();
		for (PdbChain chain:getPolyChains()) {
			all.add(chain.getPdbChainCode());	
		}
		return all;		
	}
	
	/**
	 * Returns the total number of chains present in this asym unit. Includes both 
	 * polymer (protein/nucleotide) and non-polymer (het residues only) chains
	 * @return
	 */
	public int getNumChains() {
		return chains.size()+nonPolyChains.size();
	}
	
	/**
	 * Returns the total number of polymer chains present in this asym unit.
	 * @return
	 */
	public int getNumPolyChains() {
		return chains.size();
	}
	
	/**
	 * Returns the total number of protein polymer chains present in this asym unit.
	 * @return
	 */
	public int getNumProtChains() {
		int count = 0;
		for (PdbChain chain:chains.values()) {
			if (chain.getSequence().isProtein()) count++;
		}
		return count;
	}
	
	/**
	 * Returns the total number of non-polymer chains present in this asym unit.
	 * @return
	 */
	public int getNumNonPolyChains() {
		return nonPolyChains.size();
	}
	
	public String getPdbCode() {
		return pdbCode;
	}
	
	public void setPdbCode(String pdbCode) {
		this.pdbCode = pdbCode;
	}
	
	/**
	 * Returns this PDB's title (may be null)
	 * @return
	 */
	public String getTitle() {
		return title;
	}
	
	public void setTitle(String title) {
		this.title = title;
	}
	
	public String getExpMethod() {
		return this.expMethod;
	}
	
	public void setExpMethod(String expMethod) {
		this.expMethod = expMethod;
	}
	
	/**
	 * Returns true if the experimental method of this PdbAsymUnit
	 * is a crystallographic one, i.e. one of :
	 * X-RAY DIFFRACTION, NEUTRON DIFFRACTION or ELECTRON CRYSTALLOGRAPHY
	 * If the experimental method is not set (null) then also true 
	 * is returned, i.e. we consider the default PDB to be crystallographic
	 * @return
	 */
	public boolean isCrystallographicExpMethod() {

		return (getExpMethod()==null ||
				getExpMethod().equals("X-RAY DIFFRACTION") || 
				getExpMethod().equals("NEUTRON DIFFRACTION") || 
				getExpMethod().equals("ELECTRON CRYSTALLOGRAPHY"));
	}
	
	public int getModel() {
		return model;
	}
	
	public void setModel(int model) {
		this.model = model;
	}
	
	public CrystalCell getCrystalCell() {
		return crystalCell;
	}
	
	public void setCrystalCell(CrystalCell crystalCell){
		this.crystalCell = crystalCell;
	}
	
	public SpaceGroup getSpaceGroup() {
		return spaceGroup;
	}
	
	public void setSpaceGroup(SpaceGroup spaceGroup) {
		this.spaceGroup = spaceGroup;
	}
	
	/**
	 * Returns the resolution or -1 if does not apply.
	 * @return
	 */
	public double getResolution() {
		return this.resolution;
	}
	
	public void setResolution(double resolution) {
		this.resolution = resolution;
	}
	
	/**
	 * Returns the R free value or -1 if does not apply or not available
	 * @return
	 */
	public double getRfree() {
		return this.rFree;
	}
	
	public void setRfree(double rFree) {
		this.rFree = rFree;
	}
	
	/**
	 * Returns the R sym/R merge value or -1 if does not apply or not available.
	 * @return
	 */
	public double getRsym() {
		return this.rSym;
	}
	
	public void setRsym(double rSym) {
		this.rSym = rSym;
	}
	
	/**
	 * Returns true if the experimental method for this entry is X-RAY DIFFRACTION
	 * False if experimental method not set or different from x-ray diffraction
	 * @return
	 */
	public boolean isXrayDiffraction() {
		if (getExpMethod()==null) return false;
		return getExpMethod().equals("X-RAY DIFFRACTION");
	}
	
	/**
	 * Adds the given polymer chain with given CIF chain code to this asym unit.
	 * @param chainCode
	 * @param chain
	 */
	public void setPolyChain(String chainCode, PdbChain chain) {
		chains.put(chainCode, chain);
	}
	
	/**
	 * Adds the given non-polymer chain with gicen CIF chain code to this asym unit. 
	 * @param chainCode
	 * @param chain
	 */
	public void setNonPolyChain(String chainCode, PdbChain chain) {
		nonPolyChains.put(chainCode,chain);
	}
	
	/**
	 * Tells whether a polymer chain of this PdbAsymUnit contains a (observed) residue (standard amino acid,
	 * het residue or nucleotide) with the given residue serial and given CIF chain code
	 * @param resSerial
	 * @param chainCode
	 * @return
	 */
	public boolean containsResidue(int resSerial, String chainCode) {
		if (!containsChainCode(chainCode)) {
			return false;
		}
		return this.getChainForChainCode(chainCode).containsResidue(resSerial);
	}
	
	/**
	 * Tells whether given polymer PdbChain object is contained in this PdbAsymUnit.
	 * Note that this depends on the implementation of PdbChain.equals(). If not 
	 * implemented then the references are compared.
	 * @param pdb
	 * @return
	 */
	public boolean containsChain(PdbChain pdb) {
		return this.chains.values().contains(pdb);
	}
	
	/**
	 * Tells whether a polymer chain with given PDB chain code is contained in this PdbAsymUnit.
	 * @param pdbChainCode
	 * @return
	 */
	public boolean containsPdbChainCode(String pdbChainCode) {
		String chainCode = getChainCodeForPdbChainCode(pdbChainCode);
		if (chainCode==null) return false;
		return this.chains.containsKey(chainCode);
	}

	/**
	 * Tells whether a polymer or non-polymer chain with given CIF chain code is contained in this PdbAsymUnit.
	 * @param pdbChainCode
	 * @return
	 */
	public boolean containsChainCode(String chainCode) {
		if (this.chains.containsKey(chainCode)) return true;
		if (this.nonPolyChains.containsKey(chainCode)) return true;
		return false;
	}

	public void setPdbchaincode2chaincode(TreeMap<String,String> pdbchaincode2chaincode) {
		this.pdbchaincode2chaincode = pdbchaincode2chaincode;	
	}
	
	/**
	 * Gets the CIF chain code corresponding to the given PDB chain code of a polymer (protein/nucleotide) chain.
	 * If given PDB chain code does not correspond to one of a polymer chain, null will be returned.
	 * @param pdbChainCode the PDB chain code
	 * @return
	 */
	public String getChainCodeForPdbChainCode(String pdbChainCode) {
		return this.pdbchaincode2chaincode.get(pdbChainCode);
	}
	
	/**
	 * Returns the Residue for the given residue serial and CIF chain code
	 * @param resSerial
	 * @param chainCode
	 * @return
	 * @throws the Residue or null if the chainCode/resSerial combination not contained in this PdbAsymUnit
	 */
	public Residue getResidue(int resSerial, String chainCode) {
		if (chainCode==null) return null;
		if (!containsChainCode(chainCode)) {
			return null;
		}
		return this.getChainForChainCode(chainCode).getResidue(resSerial);
	}
	
	/**
	 * Returns the centroid (average coordinate) of all atoms in this structure 
	 * (i.e. center of mass disregarding atom masses)
	 * Atoms of all chains (polymer/non-polymer) and of all kinds of residues (standard aminoacids, 
	 * nucleotides, hets) are considered
	 * The centroid valued is cached and reused upon subsequent calls
	 * @return
	 */
	public Point3d getCentroid() {
		if (centroid!=null) return centroid;
		Vector3d sumVector = new Vector3d();
		int numAtoms = 0;
		for (PdbChain chain:getAllChains()) {
			for (Residue residue:chain) {
				for (Atom atom:residue) {
					sumVector.add(atom.getCoords());
					numAtoms++;
				}
			}
		}
		sumVector.scale(1.0/numAtoms);
		centroid = new Point3d(sumVector);
		return centroid;
	}
	
	/**
	 * Returns the separation in the three crystal axes (unit cell units) of this PdbAsymUnit's 
	 * centroid with respect to the given one's centroid
	 * @param pdb
	 * @return
	 */
	public Point3d getCrystalSeparation(PdbAsymUnit pdb) {
		Point3d thisCoM  = this.getCentroid();
		Point3d otherCoM = pdb.getCentroid();
		crystalCell.getCrystalFromOrthCoords(thisCoM);
		crystalCell.getCrystalFromOrthCoords(otherCoM);
		double asep = otherCoM.x-thisCoM.x;
		double bsep = otherCoM.y-thisCoM.y;
		double csep = otherCoM.z-thisCoM.z;
		return new Point3d(asep,bsep,csep);
	}
	
	/**
	 * Transforms (rotation+translation) this structure in place as indicated by the given matrix. 
	 * @param m
	 */
	public void transform(Matrix4d m) {
		this.bounds = null; // cached bounds must be reset whenever we transform the coordinates
		this.centroid = null;
		for (PdbChain pdb:getAllChains()) {
			pdb.transform(m);
		}
	}
	
	/**
	 * Mirror this structure by inverting through the origin.
	 */
	public void mirror() {
		this.bounds = null; // cached bounds must be reset whenever we transform the coordinates
		this.centroid = null;
		for (PdbChain chain:getAllChains()){
			chain.mirror();
		}
	}
	
	private List<PdbAsymUnit> getSymRelatedObjects() {
		List<PdbAsymUnit> syms = new ArrayList<PdbAsymUnit>();
		int i = 1; // we start at 1 because we want to skip the identity
		for (Matrix4d m:this.getTransformations()) {
			PdbAsymUnit sym = this.copy();
			sym.bounds = null;
			sym.centroid = null;
			for (PdbChain chain:sym.getAllChains()) {
				// note that the call to transform in PdbChain also resets the cached bounds in the chain
				chain.transform(m);
			}
			sym.setTransform((Matrix4d)this.spaceGroup.getTransformation(i).clone());
			sym.setTransformId(i);
			syms.add(sym);
			i++;
		}
		return syms;
	}
	
	/**
	 * Generates all symmetry-related objects from this asym unit and returns the whole
	 * unit cell (this asymmetric unit plus the symmetry-related objects)
	 * @return
	 */
	public PdbUnitCell getUnitCell() {
		PdbUnitCell cell = new PdbUnitCell();
		List<PdbAsymUnit> syms = this.getSymRelatedObjects();
		cell.addUnit(this);
		for (PdbAsymUnit sym:syms) {
			cell.addUnit(sym);
		}
		return cell;
	}
	
	/**
	 * Translates this PdbAsymUnit to the given unit cell (direction).
	 * e.g. doCrystalTranslation(new Vector3d(1,1,1)) will translate this PdbAsymUnit to 
	 * crystal cell (1,1,1), considering always this PdbAsymUnit's cell to be (0,0,0)
	 * @param direction
	 */
	public void doCrystalTranslation(Vector3d direction) {
		// in order to speed up things we translate bounds instead of resetting them
		// so no need for recalculating them will exist after
		if (bounds!=null || centroid!=null) {
			Matrix4d m = getCrystalCell().getTransform(direction);
			if (bounds!=null) this.bounds.translate(new Vector3d(m.m03,m.m13,m.m23)); 
			if (centroid!=null) this.centroid.add(new Vector3d(m.m03,m.m13,m.m23));
		}
		
		for (PdbChain pdb:getAllChains()) {
			// note that the bounds in the chain are also translated by this call
			pdb.doCrystalTranslation(direction);
		}		

		transform.m03 = transform.m03+direction.x;
		transform.m13 = transform.m13+direction.y;
		transform.m23 = transform.m23+direction.z;
		// note that transformId doesn't change here: 
		// that's the whole point of having such an id: to identify equivalent crystal symmetry units  
	}
	
	/**
	 * Deep copies this PdbAsymUnit
	 * @return
	 */
	public PdbAsymUnit copy() {
		PdbAsymUnit newAsym = new PdbAsymUnit();
		newAsym.pdbCode = pdbCode;
		newAsym.model = model;
		newAsym.title = title;
		newAsym.crystalCell = crystalCell;
		newAsym.spaceGroup = spaceGroup;
		newAsym.expMethod = expMethod;
		newAsym.resolution = resolution;
		newAsym.rFree = rFree;
		newAsym.rSym = rSym;
		
		for (PdbChain chain:getPolyChains()){
			newAsym.chains.put(chain.getChainCode(), chain.copy(newAsym));
		}
		for (PdbChain chain:getNonPolyChains()){
			newAsym.nonPolyChains.put(chain.getChainCode(), chain.copy(newAsym));
		}

		newAsym.pdbchaincode2chaincode = new TreeMap<String, String>();
		newAsym.pdbchaincode2chaincode.putAll(this.pdbchaincode2chaincode);
		newAsym.chain2repChain = null;
		newAsym.repChain2members = null;
		if (this.bounds!=null) {
			newAsym.bounds = new BoundingBox(bounds.xmin, bounds.xmax, bounds.ymin, bounds.ymax, bounds.zmin, bounds.zmax);
		}
		newAsym.boundsProtOnly = this.boundsProtOnly;
		
		if (this.centroid!=null) {
			newAsym.centroid = new Point3d(this.centroid);
		}
		
		newAsym.setTransform(new Matrix4d(this.transform));
		newAsym.setTransformId(this.getTransformId());
		return newAsym;
	}
	
	/**
	 * Writes PDB data out in PDB file format.
	 * The chain codes written are the CIF chain codes. If they are longer than one character
	 * then only the first one is used and a warning issued.
	 * @param outFile
	 * @throws FileNotFoundException
	 */
	public void writeToPdbFile(File outFile) throws FileNotFoundException {
		PrintStream ps = new PrintStream(outFile);
		writePDBFileHeader(ps);
		for (PdbChain chain:getPolyChains()) {
			chain.writeSeqresRecord(ps, chain.getChainCode());
		}
		for (PdbChain chain:getAllChains()) {
			chain.writeAtomLines(ps);
		}
		ps.println("END");
		ps.close();
	}
	
	/**
	 * Writes to given PrintWriter the PDB file format HEADER line
	 * @param out
	 */
	public void writePDBFileHeader(PrintStream out) {
		out.printf("HEADER%56s",pdbCode);
		out.println();
	}
	

	
	/**
	 * Gets all symmetry transformation operators corresponding to this Pdb's space group 
	 * (except for the identity) expressed in the orthonormal basis. Using PDB's axes 
	 * convention (NCODE=1).
	 * @return
	 */	
	public List<Matrix4d> getTransformations() {
		List<Matrix4d> transfs = new ArrayList<Matrix4d>();
		for (int i=1;i<this.getSpaceGroup().getNumOperators();i++) {
			transfs.add(this.crystalCell.transfToOrthonormal(this.getSpaceGroup().getTransformation(i)));
		}
		return transfs;
	}
	
	/**
	 * Gets the transformation used to generate the asymmetric unit, if no transformation used on
	 * it yet, this will return the identity matrix.
	 * @return
	 */
	public Matrix4d getTransform() {
		return transform;
	}
	
	/**
	 * Gets the translation component of the transformation used to generate the asymmetric unit
	 * @return
	 */
	protected Vector3d getTranslation() {
		return new Vector3d(transform.m03,transform.m13,transform.m23);
	}
	
	/**
	 * Sets the transformation used to generate the asymmetric unit.
	 * @param transform
	 */
	public void setTransform(Matrix4d transform) {
		this.transform = transform;
	}

	/**
	 * Returns the transform id, which identifies which space group symmetry operator was
	 * used to generate this asym unit, 2 equivalent units related by a crystal translation
	 * will have the same id. 
	 * @param id
	 */
	public int getTransformId() {
		return transformId;
	}
	
	/**
	 * Sets the transform id, which identifies which space group symmetry operator was
	 * used to generate this asym unit, 2 equivalent units related by a crystal translation
	 * will have the same id. 
	 * @param id
	 */
	public void setTransformId(int id) {
		this.transformId = id;
	}

	/**
	 * Returns a sorted (decreasing area) list of all interfaces that this PdbAsymUnit has upon 
	 * generation of all crystal symmetry objects. An interface is defined as any pair of chains 
	 * that contact, i.e. for which there is at least a pair of atoms (one from each chain) within 
	 * the given cutoff distance.
	 * The interface areas and BSAs are calculated with either our implementation of the rolling
	 * ball algorithm (naccessExe set to null) or the external NACCESS program (naccessExe must 
	 * be passed)
	 * @param cutoff the distance cutoff for 2 chains to be considered in contact
	 * @param naccessExe the NACCESS executable if null our rolling ball algorithm implementation
	 * will be used
	 * @param nSpherePoints
	 * @param nThreads
	 * @param hetAtoms whether to consider HETATOMs in surface area calculations or not
	 * @param nonPoly if true interfaces will be calculated for non-polymer chains and 
	 * protein/nucleic acid polymers, if false only interfaces between protein polymer chains calculated
	 * @return
	 * @throws IOException when problems when running NACCESS (if NACCESS used)
	 */
	public ChainInterfaceList getAllInterfaces(double cutoff, File naccessExe, int nSpherePoints, int nThreads, boolean hetAtoms, boolean nonPoly) throws IOException {
		InterfacesFinder interfFinder = new InterfacesFinder(this);
		return interfFinder.getAllInterfaces(cutoff, naccessExe, nSpherePoints, nThreads, hetAtoms, nonPoly);
	}

	/**
	 * Calculate the Accessible Surface Areas using our implementation of the 
	 * rolling ball algorithm. Sets both the Atoms' and Residues' asa members.
	 * See Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms. 
	 * Lysozyme and Insulin." JMB (1973) 79:351-371.
	 * @param nSpherePoints
	 * @param nThreads
	 * @param hetAtoms whether to consider hetAtoms or not for area calculations 
	 */
	public void calcASAs(int nSpherePoints, int nThreads, boolean hetAtoms) {
		
		int numAtoms = this.getNumAtoms();
		if (!hetAtoms) {
			numAtoms = this.getNumNonHetAtoms();
		}
		Atom[] atoms = new Atom[numAtoms];

		int i = 0;
		for (PdbChain pdb:getAllChains()) {
			pdb.setAtomRadii();
			for (Residue residue:pdb) {
				if (!hetAtoms && (residue instanceof HetResidue)) continue;
				for (Atom atom:residue) {
					atoms[i] = atom;
					i++;
				}
			}
		}
		
		double[] asas = Asa.calculateAsa(atoms, Asa.DEFAULT_PROBE_SIZE, nSpherePoints, nThreads);
		for (i=0;i<atoms.length;i++){
			atoms[i].setAsa(asas[i]);
		}

		// and finally sums per residue
		for (PdbChain pdb:getAllChains()) {
			for (Residue residue: pdb) {
				double tot = 0;
				for (Atom atom:residue) {
					tot+=atom.getAsa();
				}
				residue.setAsa(tot);
				if (residue instanceof AaResidue) {
					AaResidue aares = (AaResidue) residue;
					residue.setRsa(tot/aares.getAaType().getAsaInExtTripept());
				}
			}
		}
	}

	/**
	 * Calculate the Buried Surface Areas by calculating the ASAs of this PdbAsymUnit 
	 * as a complex and then using the (previously calculated) isolated chain ASA values
	 * to compute the difference and set the bsa members of Atoms and Residues
	 * @param nSpherePoints
	 * @param nThreads
	 * @param hetAtoms whether to consider hetAtoms or not for area calculations
	 */
	public void calcBSAs(int nSpherePoints, int nThreads, boolean hetAtoms) {
		
		int numAtoms = this.getNumAtoms();
		if (!hetAtoms) {
			numAtoms = this.getNumNonHetAtoms();
		}
		Atom[] atoms = new Atom[numAtoms];
		
		int i = 0;
		for (PdbChain pdb:getAllChains()) {
			pdb.setAtomRadii();
			for (Residue residue:pdb) {
				if (!hetAtoms && (residue instanceof HetResidue)) continue;
				for (Atom atom:residue) {
					atoms[i] = atom;
					i++;
				}
			}
		}
		
		double[] asas = Asa.calculateAsa(atoms, Asa.DEFAULT_PROBE_SIZE, nSpherePoints, nThreads);
		for (i=0;i<atoms.length;i++){
			atoms[i].setBsa(atoms[i].getAsa()-asas[i]);
		}
		// and finally sums per residue
		for (PdbChain pdb:getAllChains()) {
			for (Residue residue: pdb) {
				double tot = 0;
				for (Atom atom:residue) {
					tot+=atom.getBsa();
				}
				residue.setBsa(tot);
			}
		}		
	}	

	protected BoundingBox getBoundingBox(boolean protOnly) {
		if (bounds!=null && boundsProtOnly==protOnly) {
			return bounds;
		}
		int numChains = getNumChains();
		if (protOnly) numChains = getNumProtChains();
		
		BoundingBox[] boxes = new BoundingBox[numChains];
		Collection<PdbChain> cs = null;
		if (!protOnly) cs = getAllChains();
		else cs = getProtChains();
		
		int i = 0;
		for (PdbChain chain:cs) {
			boxes[i] = chain.getBoundingBox();
			i++;
		}
		bounds = new BoundingBox(boxes);
		boundsProtOnly = protOnly;
		return bounds;
	}
	
	/**
	 * Returns true if the two bounding boxes of this AsymUnit and other AsymUnit
	 * are not within the given cutoff, false if they are within the given cutoff
	 * @param other
	 * @param cutoff
	 * @param protOnly whether to consider protein chains only or all (poly and non-poly)
	 * @return
	 */
	public boolean areNotOverlapping(PdbAsymUnit other, double cutoff, boolean protOnly) {
		BoundingBox thisbb = this.getBoundingBox(protOnly);
		BoundingBox otherbb = other.getBoundingBox(protOnly);
		return !thisbb.overlaps(otherbb, cutoff);
	}
	
	/**
	 * Returns number of atoms in the protein, including Hydrogens if they are present
	 * Includes all chains (polymer/non-polymer) and all residues (standard aminoacids, 
	 * non-standard aminoacids, nucleotides and hets) 
	 * @return number of atoms
	 */
	public int getNumAtoms() {
		int tot = 0;
		for (PdbChain pdb:getPolyChains()) {
			tot+=pdb.getNumAtoms();
		}
		for (PdbChain pdb:getNonPolyChains()) {
			tot+=pdb.getNumAtoms();
		}

		return tot;
	}
	
	/**
	 * Returns number of atoms in the protein, including Hydrogens if they are present
	 * Includes all chains (polymer/non-polymer) but only standard amino acids and nucleotides
	 * (no het residues)
	 * @return number of atoms
	 */
	public int getNumNonHetAtoms() {
		int tot = 0;
		for (PdbChain pdb:getPolyChains()) {
			tot+=pdb.getNumNonHetAtoms();
		}
		for (PdbChain pdb:getNonPolyChains()) {
			tot+=pdb.getNumNonHetAtoms();
		}

		return tot;		
	}
	
	/**
	 * Initialises the map of chain codes to representative chain codes.
	 * @see #getRepChain(String)
	 */
	private void initialiseRepChainsMaps() {
		// map of sequences to list of PDB chain codes
		Map<String, List<String>> uniqSequences = new HashMap<String, List<String>>();
		// finding the entities (groups of identical chains)
		for (PdbChain pdb:getPolyChains()) {
			
			if (uniqSequences.containsKey(pdb.getSequence().getSeq())) {
				uniqSequences.get(pdb.getSequence().getSeq()).add(pdb.getPdbChainCode());
			} else {
				List<String> list = new ArrayList<String>();
				list.add(pdb.getPdbChainCode());
				uniqSequences.put(pdb.getSequence().getSeq(),list);
			}		
		}
		
		chain2repChain = new HashMap<String,String>();
		repChain2members = new HashMap<String,List<String>>();
		for (List<String> entity:uniqSequences.values()) {
			for (int i=0;i<entity.size();i++) {
				chain2repChain.put(entity.get(i),entity.get(0));
			}
			repChain2members.put(entity.get(0), entity);
		}
		
	}
	
	/**
	 * Returns the representative chain's PDB chain code for the given PDB chain 
	 * code, i.e. the first chain (alphabetically) that represents a set of 
	 * sequence-identical chains.
	 * A PDB entry is composed of several chains, some of them can be identical in sequence
	 * (an entity). For each of those groups (entities) we define a representative chain as (the 
	 * first PDB chain code alphabetically). 
	 * That's the chain returned here, given any PDB chain code.
	 * @param pdbChainCode the chain for which we want to get the representative chain
	 * @return
	 * @see {@link #getSeqIdenticalGroup(String)}
	 */
	public String getRepChain(String pdbChainCode) {
		if (chain2repChain==null) {
			initialiseRepChainsMaps();
		}
		return chain2repChain.get(pdbChainCode);
	}
	
	/**
	 * Returns the map containing as keys all PDB chain codes and values their corresponding 
	 * representative chains (the first chain alphabetically from the group of identical chains)
	 * @return
	 * @see #getRepChain(String)
	 */
	public HashMap<String,String> getChain2repChainMap() {
		if (chain2repChain==null) {
			initialiseRepChainsMaps();
		}
		return chain2repChain;
	}
	
	/**
	 * For a given PDB chain code of a representative chain returns all the PDB chain codes
	 * of the sequence-identical chains that it represents.
	 * @param repPdbChainCode
	 * @return
	 * @see {@link #getRepChain(String)}
	 */
	public List<String> getSeqIdenticalGroup(String repPdbChainCode) {
		if (repChain2members==null) {
			initialiseRepChainsMaps();
		}
		return repChain2members.get(repPdbChainCode);
	}
	
	/**
	 * Returns a sorted list of the representative chains' PDB chain codes. 
	 * @return
	 * @see {@link #getRepChain(String)} and {@link #getSeqIdenticalGroup(String)}
	 */
	public List<String> getAllRepChains() {
		if (repChain2members==null) {
			initialiseRepChainsMaps();
		}
		List<String> reps = new ArrayList<String>();
		for (String rep:repChain2members.keySet()) {
			reps.add(rep);
		}
		Collections.sort(reps);
		return reps;
	}
	
	/**
	 * Returns a string containing the given representative PDB chain code and the
	 * PDB chain codes of all sequence-identical chains that it represents. 
	 * The string is formatted like: A (B,C,D,E)  
	 * @param repPdbChainCode
	 * @return
	 */
	public String getSeqIdenticalGroupString(String repPdbChainCode) {
		if (repChain2members==null) {
			initialiseRepChainsMaps();
		}
		
		List<String> members =  repChain2members.get(repPdbChainCode);
		String str = repPdbChainCode;
		
		if (members.size()>1) str+=" (";
		
		for (int i=0;i<members.size();i++) {
			if (!members.get(i).equals(repPdbChainCode)) {
				if (i==members.size()-1) {
					str+= members.get(i)+")";
				}else {
					str+= members.get(i)+",";
				}
			}
		}
		return str;
	}
	
	/**
	 * Calculates molecular weight of this PDB asym unit for all chains (polymer or non-polymer)
	 * and all atoms of all residues
	 * @return
	 */
	public double getMass() {
		double mass = 0;
		for (PdbChain chain:getAllChains()) {
			mass+=chain.getMass();
		}
		return mass;
	}
	
	/**
	 * Removes all Hydrogen atoms from this structure
	 */
	public void removeHatoms() {
		for (PdbChain chain:this.getAllChains()) {
			chain.removeHatoms();
		}
	}
	
	/**
	 * Returns true if at least one Hydrogen atom is present in a standard amino-acid 
	 * of at least one protein chain of this PDB structure
	 * @see #removeHatoms()
	 * @return
	 */
	public boolean hasHydrogens() {
		for (PdbChain chain:this.getPolyChains()) {
			if (chain.hasHydrogens()) return true;
		}
		return false;
	}
	
	/*--------------------------------------- static methods -----------------------------------------*/
	
	/**
	 * Grabs and unzips a cif file from either the online PDB ftp or a local directory 
	 * containing zipped cif files. The file is written to the given cifFile. 
	 * @param localCifDir
	 * @param pdbFtpCifUrl
	 * @param pdbCode
	 * @param cifFile the file where the data will be written to
	 * @param online if true taken from online ftp, otherwise from local dir
	 * @throws IOException
	 */
	public static void grabCifFile(String localCifDir, String pdbFtpCifUrl, String pdbCode, File cifFile, boolean online) throws IOException {
		pdbCode = pdbCode.toLowerCase();
		String gzCifFileName = pdbCode+".cif.gz";
		File gzCifFile = null;
		if (!online) {	
			gzCifFile = new File(localCifDir,gzCifFileName);
		} else {
			gzCifFile = File.createTempFile(pdbCode+".", ".cif.gz");
			System.out.println("Downloading cif file "+gzCifFileName+" from ftp...");
			// getting gzipped cif file from ftp
			URL url = new URL(pdbFtpCifUrl+gzCifFileName);
			URLConnection urlc = url.openConnection();
			InputStream is = urlc.getInputStream();
			FileOutputStream os = new FileOutputStream(gzCifFile);
			int b;
			while ( (b=is.read())!=-1) {
				os.write(b);
			}
			is.close();
			os.close();
			gzCifFile.deleteOnExit();
		} 

		// unzipping file
		GZIPInputStream zis = new GZIPInputStream(new FileInputStream(gzCifFile));
		FileOutputStream os = new FileOutputStream(cifFile);
		int b;
		while ( (b=zis.read())!=-1) {
			os.write(b);
		}
		zis.close();
		os.close();
	}




}

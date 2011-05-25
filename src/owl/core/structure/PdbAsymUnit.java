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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import edu.uci.ics.jung.graph.util.EdgeType;

import owl.core.structure.graphs.AICGEdge;
import owl.core.structure.graphs.AICGraph;
import owl.core.util.BoundingBox;
import owl.core.util.FileFormatException;
import owl.core.util.FileTypeGuesser;
import owl.core.util.Grid;
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


	
	private static final Matrix4d IDENTITY_TRANSFORM = new Matrix4d(1,0,0,0,
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
			loadFromPdbFile(pdbSourceFile);
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
	
	private void loadFromPdbFile(File pdbFile) throws PdbLoadException {
		PdbfileParser parser = new PdbfileParser(pdbFile.getAbsolutePath());
		
//		String[] chainCodes = parser.getChains();
//		for (int i=0;i<chainCodes.length;i++) {
//			PdbChain chain = parser.readChain(chainCodes[i],model);
//			chain.setParent(this);
//			chains.put(chainCodes[i],chain);
//		}
		
		parser.readChains(this,model);
		
		this.pdbCode = parser.getPdbCode();
		this.expMethod = parser.getExpMethod();
		this.title = parser.getTitle();
		this.crystalCell = parser.getCrystalCell();
		this.spaceGroup = parser.getSpaceGroup();
		this.resolution = parser.getResolution();
		this.rFree = parser.getRfree();
		this.rSym = parser.getRsym();
		

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
	 * Returs a sorted set of all CIF chain codes of polymer chains in this entry
	 * @return
	 */
	private Set<String> getPolyChainCodes() {
		Set<String> all = new TreeSet<String>();
		all.addAll(chains.keySet());
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
	 * Returns the average coordinate of all atoms in this structure (disregarding atom masses)
	 * Atoms of all chains (polymer/non-polymer) and of all kinds of residues (standard aminoacids, 
	 * nucleotides, hets) are considered
	 * @return
	 */
	public Point3d getCenterOfMass() {
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
		return new Point3d(sumVector);
	}
	
	/**
	 * Returns the separation in the three crystal axes (unit cell units) of this PdbAsymUnit's 
	 * centre of mass with respect to the given one's centre of mass.
	 * @param pdb
	 * @return
	 */
	public Point3d getCrystalSeparation(PdbAsymUnit pdb) {
		Point3d thisCoM  = this.getCenterOfMass();
		Point3d otherCoM = pdb.getCenterOfMass();
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
		for (PdbChain pdb:getAllChains()) {
			pdb.transform(m);
		}
	}
	
	/**
	 * Mirror this structure by inverting through the origin.
	 */
	public void mirror() {
		for (PdbChain chain:getAllChains()){
			chain.mirror();
		}
	}
	
	private List<PdbAsymUnit> getSymRelatedObjects() {
		List<PdbAsymUnit> syms = new ArrayList<PdbAsymUnit>();
		int i = 1; // we start at 1 because we want to skip the identity
		for (Matrix4d m:this.getTransformations()) {
			PdbAsymUnit sym = this.copy();
			for (PdbChain chain:sym.getAllChains()) {
				chain.transform(m);
			}
			// the transformed object might end up in another cell but we want it in the original cell
			// that's why we check it now and translate if it wasn't
			Point3d sep3d = this.getCrystalSeparation(sym);
			Point3i sep = new Point3i((int)Math.round(sep3d.x),(int)Math.round(sep3d.y),(int)Math.round(sep3d.z));
			if (!sep.equals(new Point3i(0,0,0))) {
				// we don't use here doCrystalTranslation method because we don't want sym's transf member to be reset
				for (PdbChain pdb:sym.getAllChains()) {
					pdb.doCrystalTranslation(new Vector3d(-sep.x,-sep.y,-sep.z));
				}	
			}
			sym.setTransform(this.spaceGroup.getTransformation(i));
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
		this.bounds = null; //bounds must be reset whenever the coordinates are transformed
		for (PdbChain pdb:getAllChains()) {
			pdb.doCrystalTranslation(direction);
		}		

		transform.m03 = transform.m03+direction.x;
		transform.m13 = transform.m13+direction.y;
		transform.m23 = transform.m23+direction.z;
		// note that transformId doesn't change here. That's the whole point of having such an id: to identify equivalent crystal symmtry units  
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
		newAsym.bounds = null;
		
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
	 * Returns a sorted (decreasing area) list of all interfaces (any 2 atoms under cutoff) 
	 * that this chain has upon generation of all crystal symmetry objects. 
	 * The interface areas and BSAs are calculated with either our implementation of the rolling
	 * ball algorithm (naccessExe set to null) or the external NACCESS program (naccessExe must 
	 * be passed).
	 * No debugging output is produced.
	 * @param cutoff the distance cutoff for 2 chains to be considered in contact
	 * @param naccessExe the NACCESS executable if null our rolling ball algorithm implementation
	 * will be used
	 * @param nSpherePoints
	 * @param nThreads
	 * @param hetAtoms whether to consider HETATOMs in surface area calculations or not
	 * @return
	 * @throws IOException when problems when running NACCESS (if NACCESS used)
	 */
	public ChainInterfaceList getAllInterfaces(double cutoff, File naccessExe, int nSpherePoints, int nThreads, boolean hetAtoms) throws IOException {
		return getAllInterfaces(cutoff, naccessExe, nSpherePoints, nThreads, hetAtoms, false);
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
	 * @param debug set to true to produce some debugging output (run times of each part of the calculation)
	 * @return
	 * @throws IOException when problems when running NACCESS (if NACCESS used)
	 */
	public ChainInterfaceList getAllInterfaces(double cutoff, File naccessExe, int nSpherePoints, int nThreads, boolean hetAtoms, boolean debug) throws IOException {	
		// TODO also take care that for longer cutoffs or for very small angles and small molecules one might need to go to the 2nd neighbour
		// TODO pathological cases, 3hz3: one needs to go to the 2nd neighbour
		
		// the set takes care of eliminating duplicates, comparison is based on the equals() 
		// and hashCode() of ChainInterface and that in turn on that of AICGraph and Atom
		Set<ChainInterface> set = new HashSet<ChainInterface>();

		// 0. generate complete unit cell
		PdbUnitCell cell = null;
		if (this.crystalCell!=null) {
			cell = this.getUnitCell();
		}
		
		long start = -1; 
		long end = -1;
		int trialCount = 0, countSkipped = 0, duplicatesCount1=0, duplicatesCount2=0, duplicatesCount3=0;
		if (debug) {
			trialCount = 0;
			start= System.currentTimeMillis();
			System.out.println("Interfaces within asymmetric unit");
		}
		// 1. interfaces within unit cell
		// 1.1 within asymmetric unit
		for (String iChainCode:this.getPolyChainCodes()) { //getChainCodes
			for (String jChainCode:this.getPolyChainCodes()) { // getChainCodes
				if (iChainCode.compareTo(jChainCode)<=0) continue;
				if (debug) {
					System.out.print(".");
					trialCount++;
				}
				PdbChain chaini = this.getChainForChainCode(iChainCode);
				PdbChain chainj = this.getChainForChainCode(jChainCode);
				AICGraph graph = chaini.getAICGraph(chainj, cutoff);
				if (graph.getEdgeCount()>0) {
					if (debug) System.out.print("x");
					// because of the bsas are values of the residues of each chain we need to make a copy so that each interface has independent residues
					PdbChain chainiCopy = chaini.copy(this);
					PdbChain chainjCopy = chainj.copy(this);
					if (!set.add(new ChainInterface(chainiCopy,chainjCopy,graph,IDENTITY_TRANSFORM,IDENTITY_TRANSFORM))) {
						duplicatesCount1++;
					}
				}											
			}
		}
		if (debug) {
			end = System.currentTimeMillis();
			System.out.println("\n"+trialCount+" trials done. Time "+(end-start)/1000+"s");
		}

		
		if (debug) {
			trialCount = 0;
			start= System.currentTimeMillis();
			System.out.println("Interfaces within the rest of the unit cell");
		}
		if (cell!=null) { // for NMR structures or structures with no crystal data we can't do more than AU contacts
			// 1.2 between the original asymmetric unit and the others resulting from applying the symmetry transformations
			for (int j=0;j<cell.getNumAsymUnits();j++) {
				PdbAsymUnit jAsym = cell.getAsymUnit(j);
				if (jAsym==this) continue; // we want to compare this to all others but not to itself
				for (PdbChain chaini:this.getPolyChains()) { // getAllChains
					for (PdbChain chainj:jAsym.getPolyChains()) { // getAllChains
						if (debug) {
							System.out.print(".");
							trialCount++;
						}
						AICGraph graph = chaini.getAICGraph(chainj, cutoff);
						if (graph.getEdgeCount()>0) {
							if (debug) System.out.print("x");
							// because of the bsas are values of the residues of each chain we need to make a copy so that each interface has independent residues
							PdbChain chainiCopy = chaini.copy(this);
							PdbChain chainjCopy = chainj.copy(jAsym);
							ChainInterface interf = new ChainInterface(chainiCopy,chainjCopy,graph,this.getTransform(),jAsym.getTransform()); 
							if (!set.add(interf)) {
								duplicatesCount2++;
							}
						}													
					}
				}

			}
			if (debug) {
				end = System.currentTimeMillis();
				System.out.println("\n"+trialCount+" trials done. Time "+(end-start)/1000+"s");
			}

			if (debug) {
				trialCount = 0;
				start= System.currentTimeMillis();
				int trials = this.getNumChains()*cell.getNumAsymUnits()*this.getNumChains()*26;
				System.out.println("Interfaces between the original asym unit and the 26 neighbouring whole unit cells ("+trials+")");
			}
			// 2. interfaces between original asymmetric unit and 26 neighbouring whole unit cells
			for (int i=-1;i<=1;i++) {
				for (int j=-1;j<=1;j++) {
					for (int k=-1;k<=1;k++) {
						if (i==0 && j==0 && k==0) continue; // that would be the identity translation, we calculate that before
						PdbUnitCell translated = cell.copy();
						Vector3d trans = new Vector3d(i,j,k);
						translated.doCrystalTranslation(trans);

						for (PdbAsymUnit jAsym:translated.getAllAsymUnits()) {
							Point3d sep = this.getCrystalSeparation(jAsym);
							if (Math.abs(sep.x)>1.1 || Math.abs(sep.y)>1.1 || Math.abs(sep.z)>1.1) {
								if (debug) {
									//System.out.println("\nskipping:");
									//System.out.printf("(%2d,%2d,%2d) - %2d : %5.2f,%5.2f,%5.2f (%2d,%2d,%2d)\n",i,j,k,jAsym.getTransformId(),
									//		sep.x,sep.y,sep.z,
									//		(int)Math.round(sep.x),(int)Math.round(sep.y),(int)Math.round(sep.z));
									countSkipped++;
								}
								continue;
							}
							for (PdbChain chainj:jAsym.getPolyChains()) {
								//try {
								//	chainj.writeToPDBFile("/home/duarte_j/"+pdbCode+"."+i+"."+j+"."+k+"."+jAsym.getTransformId()+".pdb");
								//} catch (FileNotFoundException e) {
								//	e.printStackTrace();
								//}

								for (PdbChain chaini:this.getPolyChains()) { // we only have to compare the original asymmetric unit to every full cell around
									if (debug) {
										System.out.print(".");
										trialCount++;
									}
									AICGraph graph = chaini.getAICGraph(chainj, cutoff);
									if (graph.getEdgeCount()>0) {
										if (debug) System.out.print("x");
										// because of the bsas are values of the residues of each chain we need to make a copy so that each interface has independent residues
										PdbChain chainiCopy = chaini.copy(this);
										PdbChain chainjCopy = chainj.copy(jAsym);
										ChainInterface interf = new ChainInterface(chainiCopy,chainjCopy,graph,this.getTransform(),jAsym.getTransform());
										if (!set.add(interf)){
											duplicatesCount3++;
										}
									}							
								}
							}
						}
					}
				}
			}
			if (debug) {
				end = System.currentTimeMillis();
				System.out.println("\n"+trialCount+" trials done ("+
						countSkipped+" branches skipped). Total "+(trialCount+countSkipped*this.getNumChains()*this.getNumChains())+
						" trials. Time "+(end-start)/1000+"s");
				System.out.println("Duplicates: "+duplicatesCount1+" "+duplicatesCount2+" "+duplicatesCount3);
				System.out.println("Found "+set.size()+" interfaces.");
			}
		}
		// bsa calculation 
		// NOTE in principle it is more efficient to calculate asas only once per isolated chain
		// BUT! surprisingly the rolling ball algorithm gives slightly different values for same molecule in different 
		// orientations! (can't really understand why!). Both NACCESS and our own implementation behave like that.
		// That's why we calculate always for the 2 separate members of interface and the complex, otherwise 
		// we get (not very big but annoying) discrepancies and also things like negative (small) bsa values
		if (debug) {
			start= System.currentTimeMillis();
			System.out.println("Calculating ASAs with "+nThreads+" threads and "+nSpherePoints+" sphere points");
		}
		for (ChainInterface interf:set) {
			//System.out.print(".");
			if (naccessExe!=null) {
				interf.calcSurfAccessNaccess(naccessExe,hetAtoms);
			} else {
				interf.calcSurfAccess(nSpherePoints, nThreads,hetAtoms);
			}
		}
		if (debug) {
			end = System.currentTimeMillis();
			System.out.println("\nDone. Time "+(end-start)/1000+"s");
		}
		
		// now that we have the areas we can put them into a list and sort them
		ChainInterfaceList.AsaCalcMethod asaCalcMethod = ChainInterfaceList.AsaCalcMethod.INTERNAL;
		if (naccessExe!=null) {
			asaCalcMethod = ChainInterfaceList.AsaCalcMethod.NACCESS;
		}
		ChainInterfaceList list = new ChainInterfaceList(asaCalcMethod);
		if (asaCalcMethod == ChainInterfaceList.AsaCalcMethod.INTERNAL) {
			list.setAsaCalcAccuracyParam(nSpherePoints);
		}
		for (ChainInterface interf:set) {
			list.addInterface(interf);
		}
		list.sort(); // this sorts the returned list and assigns ids to the ChainInterface members
		
		return list;
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
	
	private Atom[] getAllAtoms() {
 
		Atom[][] allAtoms = new Atom[this.getNumChains()][];
		BoundingBox[] boxes = new BoundingBox[getNumChains()];
		int numThisAtoms = 0;
		int i = 0;
		for (PdbChain chain:getAllChains()) {
			allAtoms[i] = chain.getAllAtoms();
			numThisAtoms += allAtoms[i].length;
			boxes[i] = chain.getBoundingBox();
			i++;
		}
		Atom[] thisAtoms = new Atom[numThisAtoms];
		int k = 0;
		for (i =0;i<this.getNumChains();i++) {
			for (int j=0;j<allAtoms[i].length;j++) {
				thisAtoms[k] = allAtoms[i][j];
				k++;
			}
		}
		if (bounds==null) {
			bounds = new BoundingBox(boxes); 
		}
		return thisAtoms;		
	}
	
	protected BoundingBox getBoundingBox() {
		if (bounds!=null) {
			return bounds;
		}
		BoundingBox[] boxes = new BoundingBox[getNumChains()];
		int i = 0;
		for (PdbChain chain:getAllChains()) {
			boxes[i] = chain.getBoundingBox();
			i++;
		}
		bounds = new BoundingBox(boxes);
		return bounds;
	}
	
	public AICGraph getAIAUGraph(PdbAsymUnit other, double cutoff) {
		Atom[] thisAtoms = this.getAllAtoms();
		Atom[] otherAtoms = other.getAllAtoms();
		
		Grid grid = new Grid(cutoff);
		grid.addAtoms(thisAtoms,this.bounds,otherAtoms,other.bounds);
		
		AICGraph graph = new AICGraph();

		float[][] distMatrix = grid.getDistMatrix(true);
		
		for (int i=0;i<distMatrix.length;i++){ 
			for (int j=0;j<distMatrix[i].length;j++){
				// the condition distMatrix[i][j]!=0.0 takes care of skipping cells for which we 
				// didn't calculate a distance because the 2 points were not in same or neighbouring boxes (i.e. too far apart)
				if (distMatrix[i][j]!=0.0f && distMatrix[i][j]<=cutoff){
					graph.addEdge(new AICGEdge(distMatrix[i][j]), thisAtoms[i], otherAtoms[j], EdgeType.UNDIRECTED);
				}

			}
		}
		return graph;
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

//	@Override
//	public Iterator<PdbChain> iterator() {
//		return chains.values().iterator();
//	}
	


}

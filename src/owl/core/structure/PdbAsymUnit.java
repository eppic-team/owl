package owl.core.structure;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import owl.core.structure.graphs.AICGraph;
import owl.core.util.FileFormatError;
import owl.core.util.FileTypeGuesser;
import owl.core.util.MySQLConnection;

/**
 * A protein crystal's asymmetric unit, i.e. a PDB entry.
 * Note that due to historical reasons this is called PdbAsymUnit and the single chain
 * object is called Pdb. 
 * 
 * TODO must refactor at some point, it would make a lot more sense to call this Pdb
 *  
 * @author duarte_j
 *
 */
public class PdbAsymUnit {
	
	private static final String DEFAULT_TRANSFORM = "X,Y,Z";
	
	private String pdbCode;
	private int model;
	private String title;
	
	private TreeMap<String, Pdb> chains;		// pdbChainCodes to Pdbs
	
	private CrystalCell crystalCell;
	private SpaceGroup spaceGroup;
	
	private String transform;
	
	/**
	 * Constructs a new PdbAsymUnit with no chains and all meta-data information passed.
	 * @param pdbCode
	 * @param model
	 * @param title
	 * @param crystalCell
	 * @param spaceGroup
	 */
	public PdbAsymUnit(String pdbCode, int model, String title, CrystalCell crystalCell, SpaceGroup spaceGroup) {
		this.pdbCode = pdbCode;
		this.model = model;
		this.title = title;
		this.crystalCell = crystalCell;
		this.spaceGroup = spaceGroup;
		this.chains = new TreeMap<String, Pdb>();
		this.transform = DEFAULT_TRANSFORM;
	}
	
	public PdbAsymUnit(File pdbSourceFile) throws IOException, FileFormatError, PdbLoadError {
		this.transform = DEFAULT_TRANSFORM;
		chains = new TreeMap<String, Pdb>();
		int type = FileTypeGuesser.guessFileType(pdbSourceFile);
		if (type==FileTypeGuesser.PDB_FILE) {
			loadFromPdbFile(pdbSourceFile);
		} else if (type==FileTypeGuesser.CIF_FILE) {
			loadFromCifFile(pdbSourceFile);
		} else {
			throw new FileFormatError("The given file does not seem to be neither a PDB file nor a mmCIF file");
		}
		
	}
	
	public PdbAsymUnit(String pdbCode, MySQLConnection conn, String dbName) throws PdbLoadError, PdbCodeNotFoundError {
		this.transform = DEFAULT_TRANSFORM;
		chains = new TreeMap<String, Pdb>();
		loadFromPdbase(pdbCode, conn, dbName);
	}
	
	private void loadFromPdbFile(File pdbFile) throws PdbLoadError {
		Pdb pdb = new PdbfilePdb(pdbFile.getAbsolutePath());
		String[] chainCodes = pdb.getChains();
		for (int i=0;i<chainCodes.length;i++) {
			Pdb chain = new PdbfilePdb(pdbFile.getAbsolutePath());
			chain.load(chainCodes[i]);
			chains.put(chainCodes[i],chain);
			if (i==0) {
				this.pdbCode = chain.getPdbCode();
				this.title = chain.getTitle();
				this.model = chain.getModel();
				this.crystalCell = chain.getCrystalCell();
				this.spaceGroup = chain.getSpaceGroup();
				
			}

		}
		
	}
	
	private void loadFromCifFile(File cifFile) throws PdbLoadError {
		Pdb pdb = new CiffilePdb(cifFile);
		String[] chainCodes = pdb.getChains();
		for (int i=0;i<chainCodes.length;i++) {
			Pdb chain = new CiffilePdb(cifFile);
			chain.load(chainCodes[i]);
			chains.put(chainCodes[i],chain);
			if (i==0) {
				this.pdbCode = chain.getPdbCode();
				this.title = chain.getTitle();
				this.model = chain.getModel();
				this.crystalCell = chain.getCrystalCell();
				this.spaceGroup = chain.getSpaceGroup();
				
			}

		}
		
	}
	
	private void loadFromPdbase(String pdbCode, MySQLConnection conn, String dbName) throws PdbLoadError, PdbCodeNotFoundError {
		try {
			Pdb pdb = new PdbasePdb(pdbCode,dbName,conn);
			String[] chainCodes = pdb.getChains();
			for (int i=0;i<chainCodes.length;i++) {
				Pdb chain = new PdbasePdb(pdbCode,dbName,conn);
				chain.load(chainCodes[i]);
				chains.put(chainCodes[i],chain);
				if (i==0) {
					this.pdbCode = chain.getPdbCode();
					this.title = chain.getTitle();
					this.model = chain.getModel();
					this.crystalCell = chain.getCrystalCell();
					this.spaceGroup = chain.getSpaceGroup();
				}
			}
		} catch(SQLException e) {
			throw new PdbLoadError(e);
		}
	}
	
	public Pdb getChain(String pdbChainCode) {
		return chains.get(pdbChainCode);
	}
	
	public Collection<Pdb> getAllChains() {
		return chains.values();
	}
	
	public Set<String> getPdbChainCodes() {
		return chains.keySet();
	}
	
	public int getNumChains() {
		return chains.size();
	}
	
	public String getPdbCode() {
		return pdbCode;
	}
	
	public String getTitle() {
		return title;
	}
	
	public int getModel() {
		return model;
	}
	
	public CrystalCell getCrystalCell() {
		return crystalCell;
	}
	
	public SpaceGroup getSpaceGroup() {
		return spaceGroup;
	}
	
	public void setChain(String pdbChainCode, Pdb chain) {
		chains.put(pdbChainCode, chain);
	}
	
	/**
	 * Tells whether this PdbAsymUnit contains the given residue serial for the given pdbChainCode
	 * @param resSerial
	 * @param pdbChainCode
	 * @return
	 */
	public boolean containsResidue(int resSerial, String pdbChainCode) {
		if (!this.chains.containsKey(pdbChainCode)) {
			return false;
		}
		return this.getChain(pdbChainCode).containsResidue(resSerial);
	}
	
	/**
	 * Tells whether given Pdb object is contained in this PdbAsymUnit.
	 * At the moment the comparison will be done based on reference, if the Pdb.equals() is 
	 * implemented then that will change.
	 * @param pdb
	 * @return
	 */
	public boolean containsChain(Pdb pdb) {
		return this.chains.values().contains(pdb);
	}
	
	/**
	 * Returns the Residue for the given residue serial and pdbChainCode
	 * @param resSerial
	 * @param pdbChainCode
	 * @return
	 * @throws the Residue or null if the pdbChaiCode/resSerial combination not contained in this PdbAsymUnit
	 */
	public Residue getResidue(int resSerial, String pdbChainCode) {
		if (!this.chains.containsKey(pdbChainCode)) {
			return null;
		}
		return this.getChain(pdbChainCode).getResidue(resSerial);
	}
	
	public Point3d getCenterOfMass() {
		Vector3d sumVector = new Vector3d();
		int numAtoms = 0;
		for (Pdb chain:this.chains.values()) {
			for(int atomserial:chain.getAllAtomSerials()) {
				Point3d coords = chain.getAtomCoord(atomserial);
				sumVector.add(coords);
				numAtoms++;
			}
		}
		sumVector.scale(1.0/numAtoms);
		return new Point3d(sumVector);
	}
	
	/**
	 * Returns an integer triplet indicating the crystal coordinates of this PdbAsymUnit with respect to the given one.
	 * e.g. if given PdbAsymUnit is a translation of this to adjacent cell (-1,0,0) then the triplet (-1,0,0) is returned  
	 * @param pdb
	 * @return
	 */
	public Point3i getCrystalSeparation(PdbAsymUnit pdb) {
		Point3d thisCoM  = this.getCenterOfMass();
		Point3d otherCoM = pdb.getCenterOfMass();
		crystalCell.getCrystalFromOrthCoords(thisCoM);
		crystalCell.getCrystalFromOrthCoords(otherCoM);
		double asep = otherCoM.x-thisCoM.x;
		double bsep = otherCoM.y-thisCoM.y;
		double csep = otherCoM.z-thisCoM.z;
		//System.out.printf("a: %5.2f ",asep);
		//System.out.printf("b: %5.2f ",bsep);
		//System.out.printf("c: %5.2f \n",csep);
		return new Point3i((int)asep, (int)bsep, (int)csep);
	}

	public void transform(Matrix4d m) {
		for (Pdb pdb:this.chains.values()) {
			pdb.transform(m);
		}
	}
	
	public void transform(Matrix4d m, PdbAsymUnit pdb) {

		for (String pdbChainCode:getPdbChainCodes()) {
			Pdb newChain = this.getChain(pdbChainCode).copy();
			this.getChain(pdbChainCode).transform(m, newChain);
			pdb.setChain(pdbChainCode, newChain);
		}
	}
	
	public List<PdbAsymUnit> getSymRelatedObjects() {
		List<PdbAsymUnit> syms = new ArrayList<PdbAsymUnit>();
		int i = 1; // we start at 1 because we want to skip the identity
		for (Matrix4d m:this.getTransformations()) {
			PdbAsymUnit sym = new PdbAsymUnit(this.pdbCode, this.model, this.title, this.crystalCell, this.spaceGroup);
			for (String pdbChainCode:getPdbChainCodes()) {
				Pdb newChain = this.getChain(pdbChainCode).copy();
				newChain.transform(m);
				sym.setChain(pdbChainCode, newChain);
			}
			// the transformed object might end up in another cell but we want it in the original cell
			// that's why we check it now and translate if it wasn't
			Point3i sep = this.getCrystalSeparation(sym);
			if (!sep.equals(new Point3i(0,0,0))) {
				sym.doCrystalTranslation(new Vector3d(-sep.x,-sep.y,-sep.z));
			}
			sym.setTransform(this.spaceGroup.getTransfAlgebraic(i));
			syms.add(sym);
			i++;
		}
		return syms;
	}
	
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
		for (Pdb pdb:this.chains.values()) {
			pdb.doCrystalTranslation(direction);
		}		
		
	}
	
	public PdbAsymUnit copy() {
		PdbAsymUnit newAsym = new PdbAsymUnit(this.pdbCode, this.model, this.title, this.crystalCell, this.spaceGroup);
		for (String pdbChainCode:getPdbChainCodes()){
			newAsym.setChain(pdbChainCode, this.getChain(pdbChainCode).copy());
		}
		return newAsym;
	}
	
	public void writeToPdbFile(File outFile) throws FileNotFoundException {
		PrintStream ps = new PrintStream(outFile);
		for (Pdb chain:getAllChains()) {
			chain.writeAtomLines(ps);
		}
		ps.close();
	}
	
	/**
	 * Gets all symmetry transformation operators corresponding to this Pdb's space group 
	 * (except for the identity) expressed in the orthonormal basis. Using PDB's axes 
	 * convention (NCODE=1).
	 * @return
	 */	
	public List<Matrix4d> getTransformations() {
		return this.chains.firstEntry().getValue().getTransformations();
	}
	
	public String getTransform() {
		return transform;
	}
	
	public void setTransform(String transform) {
		this.transform = transform;
	}
	
	/**
	 * Returns a list of all interfaces (any 2 atoms under cutoff) that this chain has 
	 * upon generation of all crystal symmetry objects. 
	 * @param cutoff the distance cutoff for 2 chains to be considered in contact 
	 * @return
	 */
	public List<ChainInterface> getAllInterfaces(double cutoff) {
		
		// TODO also take care that for longer cutoffs or for very small angles and small molecules one might need to go to the 2nd neighbour
		// TODO pathological cases, 3hz3: one needs to go to the 2nd neighbour
		
		List<ChainInterface> list = new ArrayList<ChainInterface>();

		// 0. generate complete unit cell
		PdbUnitCell cell = this.getUnitCell();
		
		// 1. interfaces within unit cell
		// 1.1 within asymmetric unit
		for (String iChainCode:this.getPdbChainCodes()) {
			for (String jChainCode:this.getPdbChainCodes()) {
				if (iChainCode.compareTo(jChainCode)<=0) continue;
				//System.out.print(".");
				Pdb chaini = this.getChain(iChainCode);
				Pdb chainj = this.getChain(jChainCode);
				AICGraph graph = chaini.getAICGraph(chainj, "ALL", cutoff);
				if (graph.getEdgeCount()>0) {
					list.add(new ChainInterface(chaini,chainj,graph,"X,Y,Z","X,Y,Z"));
				}											
			}
			//System.out.println();
		}
		//System.out.println();
		
		// 1.2 between the original asymmetric unit and the others resulting in applying the symmetry transformations
		for (int j=0;j<cell.getNumAsymUnits();j++) {
			PdbAsymUnit jAsym = cell.getAsymUnit(j);
			if (jAsym==this) continue; // we want to compare this to all others but not to itself
			//System.out.print(".");
			for (Pdb chaini:this.getAllChains()) {
				for (Pdb chainj:jAsym.getAllChains()) {
					AICGraph graph = chaini.getAICGraph(chainj, "ALL", cutoff);
					if (graph.getEdgeCount()>0) {
						list.add(new ChainInterface(chaini,chainj,graph,this.getTransform(),jAsym.getTransform()));
					}													
				}
			}
			
		}
		//System.out.println();
		
		// 2. interfaces between original asymmetric unit and 26 neighbouring whole unit cells, we only need to check half of them (13)
		// to choose the half we simply take the left half of the tree (i<=0, if i==0 then j<=0, if i==0 & j==0 then k<0) (PISA does it in the same way) 
		for (int i=-1;i<=0;i++) {
			for (int j=-1;j<=1;j++) {
				if (i==0 && j>0) continue;
				for (int k=-1;k<=1;k++) {
					if (i==0 && j==0 && k>=0) continue;
					PdbUnitCell translated = cell.copy();
					translated.doCrystalTranslation(new Vector3d(i,j,k));
					//System.out.print(".");
					for (Pdb chaini:this.getAllChains()) { // we only have to compare the original asymmetric unit to every full cell around
						for (Pdb chainj:translated.getAllChains()) {
							AICGraph graph = chaini.getAICGraph(chainj, "ALL", cutoff);
							if (graph.getEdgeCount()>0) {
								list.add(new ChainInterface(chaini,chainj,graph,"X,Y,Z",String.format("X%+d,Y%+d,Z%+d",i,j,k)));
							}							
						}
					}
				}
				//System.out.println();
			}
			//System.out.println();
		}
 
		return list;
	}
	
}

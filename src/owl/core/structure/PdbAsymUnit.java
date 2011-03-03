package owl.core.structure;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
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
import java.util.zip.GZIPInputStream;

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
	
	private static final Matrix4d IDENTITY_TRANSFORM = new Matrix4d(1,0,0,0,
																	0,1,0,0,
																	0,0,1,0,
																	0,0,0,1);
	
	private String pdbCode;
	private int model;
	private String title;
	
	private TreeMap<String, Pdb> chains;		// pdbChainCodes to Pdbs
	private HashMap<String, String> chain2repChain; // a map of all pdb chain codes to the representative pdb chain code (the first chain alphabetically from the group of identical chains)
	private HashMap<String,List<String>> repChain2members; // a map of representative pdb chain code to the members of the group of identical chains
	
	private CrystalCell crystalCell;
	private SpaceGroup spaceGroup;
	
	private Matrix4d transform; // transformation matrix used to generate this asym unit (includes the possible translation from the original cell)
	
	private int transformId;    // and identifier of the space group transformation used (0 is identity i.e. original asymmetric unit) 
								// (does not count the translations: 2 equivalent asym units of 2 different cells will have the same identifier)
								// i.e. it is unique within the unit cell but equivalent units of different crystal cells will have same id
								// goes from 1 to m (m=number of symmetry operations of the space group)
	
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
		this.transform = IDENTITY_TRANSFORM;
		this.transformId = 0;
	}
	
	public PdbAsymUnit(File pdbSourceFile) throws IOException, FileFormatError, PdbLoadError {
		this.transform = IDENTITY_TRANSFORM;
		this.transformId = 0;
		chains = new TreeMap<String, Pdb>();
		int type = FileTypeGuesser.guessFileType(pdbSourceFile);
		if (type==FileTypeGuesser.PDB_FILE || type ==FileTypeGuesser.RAW_PDB_FILE) {
			loadFromPdbFile(pdbSourceFile);
		} else if (type==FileTypeGuesser.CIF_FILE) {
			loadFromCifFile(pdbSourceFile);
		} else {
			throw new FileFormatError("The given file does not seem to be neither a PDB file nor a mmCIF file");
		}
		
	}
	
	public PdbAsymUnit(String pdbCode, MySQLConnection conn, String dbName) throws PdbLoadError, PdbCodeNotFoundException {
		this.transform = IDENTITY_TRANSFORM;
		this.transformId = 0;
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
	
	private void loadFromPdbase(String pdbCode, MySQLConnection conn, String dbName) throws PdbLoadError, PdbCodeNotFoundException {
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
			Point3d sep3d = this.getCrystalSeparation(sym);
			Point3i sep = new Point3i((int)Math.round(sep3d.x),(int)Math.round(sep3d.y),(int)Math.round(sep3d.z));
			if (!sep.equals(new Point3i(0,0,0))) {
				// we don't use here doCrystalTranslation method because we don't want sym's transf member to be reset
				for (Pdb pdb:sym.chains.values()) {
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
		transform.m03 = transform.m03+direction.x;
		transform.m13 = transform.m13+direction.y;
		transform.m23 = transform.m23+direction.z;
		// note that transformId doesn't change here. That's the whole point of having such an id: to identify equivalent crystal symmtry units  
	}
	
	public PdbAsymUnit copy() {
		PdbAsymUnit newAsym = new PdbAsymUnit(this.pdbCode, this.model, this.title, this.crystalCell, this.spaceGroup);
		for (String pdbChainCode:getPdbChainCodes()){
			newAsym.setChain(pdbChainCode, this.getChain(pdbChainCode).copy());
		}
		newAsym.setTransform(new Matrix4d(this.transform));
		newAsym.setTransformId(this.getTransformId());
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
	 * be passed)
	 * @param cutoff the distance cutoff for 2 chains to be considered in contact
	 * @param naccessExe the NACCESS executable if null our rolling ball algorithm implementation
	 * will be used
	 * @param nSpherePoints
	 * @param nThreads
	 * @return
	 * @throws IOException when problems when running NACCESS (if NACCESS used)
	 */
	public ChainInterfaceList getAllInterfaces(double cutoff, File naccessExe, int nSpherePoints, int nThreads) throws IOException {	
		// TODO also take care that for longer cutoffs or for very small angles and small molecules one might need to go to the 2nd neighbour
		// TODO pathological cases, 3hz3: one needs to go to the 2nd neighbour
		
		// the set takes care of eliminating duplicates, comparison is based on the equals() 
		// and hashCode() of ChainInterface and that in turn on that of AICGraph and Atom
		Set<ChainInterface> set = new HashSet<ChainInterface>();

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
					// because of the bsas are values of the residues of each chain we need to make a copy so that each interface has independent residues
					Pdb chainiCopy = chaini.copy();
					Pdb chainjCopy = chainj.copy();
					set.add(new ChainInterface(chainiCopy,chainjCopy,graph,IDENTITY_TRANSFORM,IDENTITY_TRANSFORM));
				}											
			}
		}

		// 1.2 between the original asymmetric unit and the others resulting in applying the symmetry transformations
		for (int j=0;j<cell.getNumAsymUnits();j++) {
			PdbAsymUnit jAsym = cell.getAsymUnit(j);
			if (jAsym==this) continue; // we want to compare this to all others but not to itself
			for (Pdb chaini:this.getAllChains()) {
				for (Pdb chainj:jAsym.getAllChains()) {
					AICGraph graph = chaini.getAICGraph(chainj, "ALL", cutoff);
					if (graph.getEdgeCount()>0) {
						// because of the bsas are values of the residues of each chain we need to make a copy so that each interface has independent residues
						Pdb chainiCopy = chaini.copy();
						Pdb chainjCopy = chainj.copy();
						ChainInterface interf = new ChainInterface(chainiCopy,chainjCopy,graph,this.getTransform(),jAsym.getTransform()); 
						set.add(interf);
					}													
				}
			}
			
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
							//System.out.println("skipping:");
							//System.out.printf("(%2d,%2d,%2d) - %2d : %5.2f,%5.2f,%5.2f (%2d,%2d,%2d)\n",i,j,k,jAsym.getTransformId(),
							//		sep.x,sep.y,sep.z,
							//		(int)Math.round(sep.x),(int)Math.round(sep.y),(int)Math.round(sep.z));
							continue;
						}
						for (Pdb chainj:jAsym.getAllChains()) {
							//try {
							//	chainj.writeToPDBFile("/home/duarte_j/"+pdbCode+"."+i+"."+j+"."+k+"."+jAsym.getTransformId()+".pdb");
							//} catch (FileNotFoundException e) {
							//	e.printStackTrace();
							//}

							for (Pdb chaini:this.getAllChains()) { // we only have to compare the original asymmetric unit to every full cell around
								AICGraph graph = chaini.getAICGraph(chainj, "ALL", cutoff);
								if (graph.getEdgeCount()>0) {
									// because of the bsas are values of the residues of each chain we need to make a copy so that each interface has independent residues
									Pdb chainiCopy = chaini.copy();
									Pdb chainjCopy = chainj.copy();
									ChainInterface interf = new ChainInterface(chainiCopy,chainjCopy,graph,this.getTransform(),jAsym.getTransform());
									set.add(interf);
								}							
							}
						}
					}
				}
			}
		}
		
		// bsa calculation 
		// NOTE in principle it is more efficient to calculate asas only once per isolated chain
		// BUT! surprisingly the rolling ball algorithm gives slightly different values for same molecule in different 
		// orientations! (can't really understand why!). Both NACCESS and our own implementation behave like that.
		// That's why we calculate always for the 2 separate members of interface and the complex, otherwise 
		// we get (not very big but annoying) discrepancies and also things like negative (small) bsa values
		//long start = System.currentTimeMillis();
		for (ChainInterface interf:set) {
			//System.out.print(".");
			if (naccessExe!=null) {
				interf.calcSurfAccessNaccess(naccessExe);
			} else {
				interf.calcSurfAccess(nSpherePoints, nThreads);
			}
		}
		//long end = System.currentTimeMillis();
		//System.out.println("ASA computation time: "+(end-start)/1000+" s");
		
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
	 */
	public void calcASAs(int nSpherePoints, int nThreads) {
		Atom[] atoms = new Atom[this.getNumAtoms()];
		
		int i = 0;
		for (Pdb pdb:this.getAllChains()) {
			pdb.setAtomRadii();
			for (int atomser: pdb.getAllAtomSerials()) {
				atoms[i] = pdb.getAtom(atomser);
				i++;
			}
		}
		
		double[] asas = Asa.calculateAsa(atoms, Asa.DEFAULT_PROBE_SIZE, nSpherePoints, nThreads);
		for (i=0;i<atoms.length;i++){
			atoms[i].setAsa(asas[i]);
		}

		// and finally sums per residue
		for (Pdb pdb:this.getAllChains()) {
			for (Residue residue: pdb.getResidues().values()) {
				double tot = 0;
				for (Atom atom:residue.getAtoms()) {
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
	 */
	public void calcBSAs(int nSpherePoints, int nThreads) {
		Atom[] atoms = new Atom[this.getNumAtoms()];
		
		int i = 0;
		for (Pdb pdb:this.getAllChains()) {
			pdb.setAtomRadii();
			for (int atomser: pdb.getAllAtomSerials()) {
				atoms[i] = pdb.getAtom(atomser);
				i++;
			}
		}
		
		double[] asas = Asa.calculateAsa(atoms, Asa.DEFAULT_PROBE_SIZE, nSpherePoints, nThreads);
		for (i=0;i<atoms.length;i++){
			atoms[i].setBsa(atoms[i].getAsa()-asas[i]);
		}
		// and finally sums per residue
		for (Pdb pdb:this.getAllChains()) {
			for (Residue residue: pdb.getResidues().values()) {
				double tot = 0;
				for (Atom atom:residue.getAtoms()) {
					tot+=atom.getBsa();
				}
				residue.setBsa(tot);
			}
		}		
	}
	
	public int getNumAtoms() {
		int tot = 0;
		for (Pdb pdb:this.getAllChains()) {
			tot+=pdb.getNumAtoms();
		}
		return tot;
	}
	
	/**
	 * Initialises the map of chain codes to representative chain codes.
	 * @see #getRepChain(String)
	 */
	private void initialiseRepChainsMaps() {
		// map of sequences to list of chain codes
		Map<String, List<String>> uniqSequences = new HashMap<String, List<String>>();
		// finding the entities (groups of identical chains)
		for (String chain:this.chains.keySet()) {
			Pdb pdb = getChain(chain);
			if (uniqSequences.containsKey(pdb.getSequence())) {
				uniqSequences.get(pdb.getSequence()).add(chain);
			} else {
				List<String> list = new ArrayList<String>();
				list.add(chain);
				uniqSequences.put(pdb.getSequence(),list);
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
	 * A PDB entrity is composed of several chains, some of them can be identical in sequence
	 * (an entity). For each of those groups (entities) we define a representative chain (the 
	 * first PDB chain code alphabetically). 
	 * That's the chain returned here, given a particular chain.
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

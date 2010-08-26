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
	
	private String pdbCode;
	private int model;
	private String title;
	
	private TreeMap<String, Pdb> chains;		// pdbChainCodes to Pdbs
	
	private CrystalCell crystalCell;
	private SpaceGroup spaceGroup;
	
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
	}
	
	public PdbAsymUnit(File pdbSourceFile) throws IOException, FileFormatError, PdbLoadError {
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
		for (Matrix4d m:this.getTransformations()) {
			PdbAsymUnit sym = new PdbAsymUnit(this.pdbCode, this.model, this.title, this.crystalCell, this.spaceGroup);
			for (String pdbChainCode:getPdbChainCodes()) {
				Pdb newChain = this.getChain(pdbChainCode).copy();
				newChain.transform(m);
				sym.setChain(pdbChainCode, newChain);
			}
			syms.add(sym);
		}
		return syms;
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
	
	/**
	 * Returns a list of all interfaces (any 2 atoms under cutoff) that this chain has 
	 * upon generation of all crystal symmetry objects. 
	 * @param cutoff the distance cutoff for 2 chains to be considered in contact 
	 * @return
	 */
	public List<ChainInterface> getAllInterfaces(double cutoff) {
		
		// the strategy here should be first generate the unit cell from the asymmetric unit with the space group's symmetry operator
		// once you have that it's like having a P1, you need 2 things: check interfaces within unit cell and with 26 (actually just half) 
		// of the neighbouring cells
		// also take care that for longer cutoffs or for very small angles and small molecules one might need to go to the 2nd neighbour
		
		//Set<ChainInterface> list = new HashSet<ChainInterface>();
		List<ChainInterface> list = new ArrayList<ChainInterface>();
		// 0. interfaces of original asymmetric unit with the space-group symmetry generated objects
		List<PdbAsymUnit> syms = getSymRelatedObjects();
		int t = 1;
		for (PdbAsymUnit sym:syms) {
			String transf = this.getSpaceGroup().getTransfAlgebraic(t);
			for (int i=-1;i<=1;i++) {
				for (int j=-1;j<=1;j++) {
					for (int k=-1;k<=1;k++) {
						PdbAsymUnit symtrans = sym.copy();
						symtrans.doCrystalTranslation(new Vector3d(i,j,k));
						for (Pdb chaini:this.getAllChains()) {
							System.out.print(".");
							for (Pdb chainj:symtrans.getAllChains()) {
								//if (!chaini.isCrystalNeighbor(chainj)) {
								//	System.out.println("Not neighbors: "+i+","+j+","+k+" "+chaini.getPdbChainCode()+" "+chainj.getPdbChainCode());
								//	continue;
								//}
								AICGraph graph = chaini.getAICGraph(chainj, "ALL", cutoff);
								if (graph.getEdgeCount()>0) {
									list.add(new ChainInterface(chaini,chainj,graph,transf));
								}							
							}
						}
						System.out.println();
					}
				}
			}
			t++;
		}
		
		// 1. interfaces within unit cell
		for (String iChain:this.chains.keySet()) {
			for (String jChain:this.chains.keySet()) {
				// we skip half of the pairs
				if (iChain.compareTo(jChain)<=0) continue; 
				Pdb chaini = chains.get(iChain);
				Pdb chainj = chains.get(jChain);
				AICGraph graph = chaini.getAICGraph(chainj, "ALL", cutoff);
				if (graph.getEdgeCount()>0) {
					list.add(new ChainInterface(chaini,chainj,graph,"X,Y,Z"));
				}							
			}
		}		
		// 2. interfaces between unit cell and 26 neighbouring cells, we only need to check half of them (13)
		// to choose the half we simply take the right half of the tree (i>=0, if i==0 then j>=0, if i==0 & j==0 then k>0) 
		for (int i=0;i<=1;i++) {
			for (int j=-1;j<=1;j++) {
				if (i==0 && j<0) continue;
				for (int k=-1;k<=1;k++) {
					if (i==0 && j==0 && k<=0) continue;
					PdbAsymUnit translated = this.copy();
					translated.doCrystalTranslation(new Vector3d(i,j,k));
					//try {
					//	translated.writeToPdbFile(new File("/home/duarte_j/"+pdbCode+"."+i+"."+j+"."+k+".pdb"));
					//} catch (IOException e) {
					//	e.printStackTrace();
					//}
					System.out.print(".");
					for (Pdb chaini:this.getAllChains()) {
						for (Pdb chainj:translated.getAllChains()) {
							AICGraph graph = chaini.getAICGraph(chainj, "ALL", cutoff);
							if (graph.getEdgeCount()>0) {
								list.add(new ChainInterface(chaini,chainj,graph,String.format("X%+d,Y%+d,Z%+d",i,j,k)));
							}							
						}
					}
				}
				System.out.println();
			}
			System.out.println();
		}
 
		return list;
	}
	
}

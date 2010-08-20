package owl.core.structure;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.Collection;
import java.util.Set;
import java.util.TreeMap;

import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

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
}

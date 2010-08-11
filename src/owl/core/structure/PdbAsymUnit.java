package owl.core.structure;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.TreeMap;

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
	
}

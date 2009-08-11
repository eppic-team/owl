package proteinstructure;

import java.util.HashMap;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

/**
 * A Pdb structure to be constructed from scratch (by passing coordinates)
 * 
 * @author duarte
 *
 */
public class ModelPdb extends Pdb {
	
	private static final String DEFAULT_CHAIN = "A";
	private static final int DEFAULT_MODEL = 1;
	

	/**
	 * Constructs a ModelPdb for given sequence setting all coordinates of all atoms to (0,0,0)
	 * @param sequence
	 * @throws IllegalArgumentException if sequence contains invalid characters
	 */
	public ModelPdb(String sequence) {
		this.chainCode = DEFAULT_CHAIN;
		this.pdbChainCode = DEFAULT_CHAIN;
		this.model = DEFAULT_MODEL;
		this.pdbCode = NO_PDB_CODE;
		
		this.dataLoaded = true;
		
		this.setSequence(sequence);
		this.setResidueData();
		this.setAtomData();
	}
	
	/**
	 * Constructs a ModelPdb for given sequence setting coordinates of the given atom type to 
	 * the given coordinates coords
	 * @param sequence
	 * @param coords the array of all coordinates, must be ordered as the residue serials and of 
	 * the same size as this Pdb's observed size
	 * @param atom the atom type for which to set the coordinates
	 * @throws IllegalArgumentException if sequence contains invalid characters 
	 */
	public ModelPdb(String sequence, Vector3d[] coords, String atom) {
		this.chainCode = DEFAULT_CHAIN;
		this.pdbChainCode = DEFAULT_CHAIN;
		this.model = DEFAULT_MODEL;
		this.pdbCode = NO_PDB_CODE;
		
		this.dataLoaded = true;
		this.setSequence(sequence);
		this.setResidueData();
		this.setAtomData();
		this.setAtomsCoords(coords, atom);
	}
	
	
	@Override
	public String[] getChains() throws PdbLoadError {
		// this doesn't make much sense here, it makes only sense for subclasses that load from some external source
		// but anyway better than returning a null we return the one chain
		String[] chains = {this.pdbChainCode};
		return chains;
	}

	@Override
	public Integer[] getModels() throws PdbLoadError {
		// this doesn't make much sense here, it makes only sense for subclasses that load from some external source
		// but anyway better than returning a null we return the one model
		Integer[] models = {this.model};
		return models;
	}

	@Override
	public void load(String pdbChainCode, int modelSerial) throws PdbLoadError {
		System.err.println("Warning: load(pdbChainCode,modelSerial) was called for a ModelPdb, this shouldn't happen!");
		// do nothing
	}
	
	/**
	 * Initializes sequence, fullLength and obsLength fields 
	 * @param sequence
	 */
	private void setSequence(String sequence) {
		this.sequence = sequence;
		this.fullLength = sequence.length();
		this.obsLength = this.fullLength;		
	}

	/**
	 * Initializes all residue data of this Pdb object
	 */
	private void setResidueData() {
		resser2restype = new HashMap<Integer, String>();
		pdbresser2resser = new HashMap<String, Integer>();
		resser2pdbresser = new HashMap<Integer, String>();
		for (int i=0;i<sequence.length();i++) {
			int resser = i+1;
			String aa = AAinfo.oneletter2threeletter(Character.toString(sequence.charAt(i)));
			if (!AAinfo.isValidAA(aa)) {
				throw new IllegalArgumentException("Given input sequence to construct a pdb model contains an invalid aminoacid "+aa);
			}
			// initialization of resser2restype
			resser2restype.put(resser, aa);	
			// initialization of resser2pdbresser and pdbresser2resser
			pdbresser2resser.put(String.valueOf(resser),resser);
			resser2pdbresser.put(resser,String.valueOf(resser));			
		}
	}
	
	/**
	 * Initializes all atom data for this Pdb data.
	 * All (non-H) atoms are added and assigned coordinate (0,0,0)
	 */
	private void setAtomData() {
		resser_atom2atomserial = new HashMap<String, Integer>();
		atomser2resser = new HashMap<Integer, Integer>();
		atomser2atom = new HashMap<Integer, String>();
		atomser2coord = new HashMap<Integer, Point3d>(); // we only initialise the Map but don't assing any coordinates	

		int atomSerial = 1; 
		for (int resser:getAllSortedResSerials()) {
			String aa = getResTypeFromResSerial(resser);
			for (String atom:AAinfo.getAtomsForCTAndRes("ALL", aa)) {
				// initialization of resser_atom2atomserial
				resser_atom2atomserial.put(resser+"_"+atom,atomSerial);
				// initialization of atomser2resser
				atomser2resser.put(atomSerial, resser);
				// initialization of atomser2atom
				atomser2atom.put(atomSerial, atom);
				// atomser2coord
				atomser2coord.put(atomSerial, new Point3d(0,0,0));
				atomSerial++;
			}
		}
		//atomser2bfactor;				// should we initialize it to 0 too?		
	}
	
	/**
	 * Resets coordinates of all atoms to (0,0,0)
	 */
	private void resetConformation() {
		for (int atomser:getAllAtomSerials()){
			this.atomser2coord.put(atomser,new Point3d(0,0,0));
		}
	}
	
	/**
	 * Sets the coordinates for given atom and residue serial
	 * @param resser
	 * @param atom
	 * @param coord
	 */
	private void setAtomCoord(int resser, String atom, Point3d coord) {
		this.atomser2coord.put(getAtomSerFromResSerAndAtom(resser, atom),coord);
	}
	
	/**
	 * Resets all coordinates and sets the coordinates of all atoms of given atom type 
	 * to given coords.
	 * @param coords the array of all coordinates, must be ordered as the residue serials and of 
	 * the same size as this Pdb's observed size
	 * @param atom the atom type for which to set the coordinates
	 */
	private void setAtomsCoords(Vector3d[] coords, String atom) {
		if (coords.length!=this.get_length()) {
			throw new IllegalArgumentException("The given vector of coordinates has different size than the number of observed residues in the original pdb entry");
		}
		resetConformation();
		int i = 0;
		for (int resser:getAllSortedResSerials()){
			setAtomCoord(resser, atom, new Point3d(coords[i]));
			i++;
		}
	}

	/**
	 * Removes all atoms that have (0,0,0) coordinates
	 */
	public void removeAtomsWith0Coords() {
		for (int atomser:getAllAtomSerials()) {
			if (getAtomCoord(atomser).epsilonEquals(new Point3d(0,0,0), 0.0000001)) {
				int resser = get_resser_from_atomser(atomser);
				String atom = getAtomNameFromAtomSer(atomser);
				resser_atom2atomserial.put(resser+"_"+atom,atomser);
				atomser2resser.remove(atomser);
				atomser2atom.remove(atomser);
				atomser2coord.remove(atomser);
			}
		}
	}
	
}

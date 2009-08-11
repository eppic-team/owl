package proteinstructure;

import java.util.Collection;
import java.util.Set;
import java.util.TreeMap;

import javax.vecmath.Point3d;

/**
 * A class representing a residue as part of a PDB protein structure 
 * @author duarte
 *
 */
public class Residue {

	
	private AminoAcid aaType;
	private int serial;
	private String pdbSerial;
	
	private Pdb parentPdb;
	
	private SecStrucElement ssElem;

	// following variables are Double/Integer objects instead of primitives
	// to be able to keep the old behaviour of using nulls when something
	// is missing
	private Double consurfScore;
	private Integer consurfColor;
	private Double rsa;
	private Double scRsa;
	
	private TreeMap<String, Atom> atoms; // atom codes to atoms
	
	/**
	 * Constructs a new Residue given its type, serial and parentPdb. The Residue will 
	 * have an empty list of Atoms until they are added using {@link #addAtom(Atom)} 
	 * @param aaType
	 * @param serial
	 * @param parentPdb
	 */
	public Residue(AminoAcid aaType, int serial, Pdb parentPdb) {
		atoms = new TreeMap<String, Atom>();
		this.aaType = aaType;
		this.serial = serial;
		this.parentPdb = parentPdb;
	}

	/**
	 * Adds an atom to this Residue. If the atom code for given Atom is already present
	 * in this Residue a warning will be printed. No duplicate atom codes are allowed in 
	 * one Residue.
	 * @param atom
	 */
	public void addAtom(Atom atom) {
		if (this.containsAtom(atom.getCode())) {
			System.err.println("Warning: atom "+atom+" being added to residue "+this+" is already present in this residue.");
		} else {
			atoms.put(atom.getCode(), atom);
		}
	}
	
	/**
	 * Returns the Atom object given its atomCode or if no such atom for this residue returns null.
	 * @param atomCode
	 * @return the Atom object or null if no such atom exists in this residue
	 */
	public Atom getAtom(String atomCode) {
		if (this.containsAtom(atomCode)) {
			return atoms.get(atomCode);
		} else {
			// we used to have this warning in the Pdb class, copied here to keep it where it should be
			//System.err.println("Couldn't find "+atomCode+" atom for resser="+this.getSerial()+" in protein "+pdbCode+" and chain "+chainCode+". Continuing without that atom for this resser.");
			return null;
		}
	}
	
	/**
	 * Returns true if this Residue instance has the given atomCode
	 * @param atomCode
	 * @return
	 */
	public boolean containsAtom(String atomCode) {
		return atoms.containsKey(atomCode);
	}

	/**
	 * Returns true if this Residue contains an OXT atom (terminal Oxygen atom 
	 * present in the last residue of a chain)
	 * @return
	 */
	public boolean containsOXT() {
		return containsAtom("OXT");
	}
	
	/**
	 * Returns the number of atoms that this Residue contains
	 * @return
	 */
	public int getNumAtoms() {
		return atoms.size();
	}
	
	public AminoAcid getAaType() {
		return aaType;
	}

	public void setAaType(AminoAcid aaType) {
		this.aaType = aaType;
	}

	public int getSerial() {
		return serial;
	}

	public void setSerial(int serial) {
		this.serial = serial;
	}

	public String getPdbSerial() {
		return pdbSerial;
	}

	public void setPdbSerial(String pdbSerial) {
		this.pdbSerial = pdbSerial;
	}

	public String getChainCode() {
		return parentPdb.getChainCode();
	}

	public Pdb getParentPdb() {
		return parentPdb;
	}

	public void setParentPdb(Pdb pdb) {
		this.parentPdb = pdb;
	}

	public SecStrucElement getSsElem() {
		return ssElem;
	}

	public void setSsElem(SecStrucElement ssElem) {
		this.ssElem = ssElem;
	}

	public Double getConsurfScore() {
		return consurfScore;
	}

	public void setConsurfScore(Double consurfScore) {
		this.consurfScore = consurfScore;
	}

	public Integer getConsurfColor() {
		return consurfColor;
	}

	public void setConsurfColor(Integer consurfColor) {
		this.consurfColor = consurfColor;
	}

	public Double getRsa() {
		return rsa;
	}

	public void setRsa(Double rsa) {
		this.rsa = rsa;
	}

	public Double getScRsa() {
		return scRsa;
	}

	public void setScRsa(Double scrsa) {
		this.scRsa = scrsa;
	}

	/**
	 * Get Collection of Atoms belonging to this Residue sorted by atom codes
	 * @return
	 */
	public Collection<Atom> getAtoms() {
		return atoms.values();
	}
	
	/**
	 * Gets a new Residue object containing only the atoms for the given contact type
	 * (the new atom instances are references to the atom instances of this Residue)
	 * @param ct
	 * @return
	 */
	public Residue getReducedResidue(String ct) {
		Residue reducedResidue = new Residue(this.getAaType(), this.getSerial(), this.getParentPdb());
		Set<String> atomCodes = AAinfo.getAtomsForCTAndRes(ct, this.getAaType().getThreeLetterCode());
		for (String atomCode:atomCodes) {
			if (this.containsAtom(atomCode)) {
				reducedResidue.addAtom(this.getAtom(atomCode));
			}
		}
		return reducedResidue;
	}

	/**
	 * Returns the centre of mass of the heavy side chains atoms present in this Residue,
	 * if no atoms in the sice chain, then the CA coordinates are returned or null if there's
	 * no CA.
	 * @return
	 */
	public Point3d getScCentroid() {
		Residue scOnly = getReducedResidue("SC");
		if (scOnly.getNumAtoms()==0) {
			if (this.containsAtom("CA")) {
				return this.getAtom("CA").getCoords();
			} else {
				return null;
			}
		}
		Point3d centroid = new Point3d();
		double massSum = 0;
		for (Atom atom:scOnly.getAtoms()) {
			Point3d coord = new Point3d(atom.getCoords());
			coord.scale(atom.getType().getAtomicMass());
			centroid.add(coord);
			massSum+=atom.getType().getAtomicMass();
		}
		centroid.scale(1.0/massSum);
		
		return centroid;
	}
	
	public String toString() {
		return this.getChainCode()+serial+aaType.getOneLetterCode();
	}
	
}

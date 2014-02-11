package owl.core.structure;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

import owl.core.structure.features.SecStrucElement;

/**
 * A class representing a standard amino acid residue as part of a PDB protein structure
 *  
 * @author duarte
 *
 */
public class AaResidue implements Residue {

	private static final long serialVersionUID = 1L;

	public enum Chirality {
		L( 1,"L","L-form"),
		D(-1,"D","D-form"),
		U(10,"U","undetermined"),
		C( 0,"C","coplanar");
		private int number;
		private String abbrev;
		private String name;
		private Chirality(int number, String abbrev, String name) {
			this.number = number;
			this.abbrev = abbrev;
			this.name = name;
		}
		public int getNumber() {
			return number;
		}
		public String getAbbrev() {
			return abbrev;
		}
		public String getName() {
			return name;
		}
	}
	
	private AminoAcid aaType;
	private int serial;
	private String pdbSerial;
	
	private PdbChain parent;
	
	private SecStrucElement ssElem;

	private double rsa;   // relative all-atoms accessible surface area
	private double scRsa; // relative side chain-atoms accessible surface area
	
	private double asa;   // all-atoms accessible surface area (absolute)
	private double bsa;   // all-atoms accessible surface area upon burial (this can only be calculated by running naccess for 2 partners and a complex, see ChainInterface class)
	
	private TreeMap<String, Atom> atoms; // atom codes to atoms
	
	/**
	 * Constructs an empty Residue. Use setters to add content to it.
	 */
	public AaResidue() {
		atoms = new TreeMap<String, Atom>();
	}
	
	/**
	 * Constructs a new Residue given its type, serial and parent. The Residue will 
	 * have an empty list of Atoms until they are added using {@link #addAtom(Atom)} 
	 * @param aaType
	 * @param serial
	 * @param parent
	 */
	public AaResidue(AminoAcid aaType, int serial, PdbChain parentPdb) {
		atoms = new TreeMap<String, Atom>();
		this.aaType = aaType;
		this.serial = serial;
		this.parent = parentPdb;
	}

	@Override
	public void addAtom(Atom atom) {
		if (this.containsAtom(atom.getCode())) {
			System.err.println("Warning: atom "+atom+" being added to residue "+this+" is already present in this residue.");
		} else {
			atoms.put(atom.getCode(), atom);
		}
	}
	
	@Override
	public Atom getAtom(String atomCode) {
		if (this.containsAtom(atomCode)) {
			return atoms.get(atomCode);
		} else {
			// we used to have this warning in the PdbChain class, copied here to keep it where it should be
			//System.err.println("Couldn't find "+atomCode+" atom for resser="+this.getSerial()+" in protein "+pdbCode+" and chain "+chainCode+". Continuing without that atom for this resser.");
			return null;
		}
	}
	
	@Override
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
	
	@Override
	public int getNumAtoms() {
		return atoms.size();
	}
	
	@Override
	public int getNumHeavyAtoms() {
		int numAtoms = 0;
		for (Atom atom:atoms.values()) {
			if (atom.getType()!=AtomType.H) {
				numAtoms++;
			}
		}
		return numAtoms;
	}
	
	public AminoAcid getAaType() {
		return aaType;
	}

	public void setAaType(AminoAcid aaType) {
		this.aaType = aaType;
	}

	@Override
	public int getSerial() {
		return serial;
	}

	@Override
	public void setSerial(int serial) {
		this.serial = serial;
	}

	@Override
	public String getPdbSerial() {
		return pdbSerial;
	}

	@Override
	public void setPdbSerial(String pdbSerial) {
		this.pdbSerial = pdbSerial;
	}

	public String getChainCode() {
		return parent.getChainCode();
	}
	
	public String getPdbChainCode() {
		return parent.getPdbChainCode();
	}

	@Override
	public SecStrucElement getSsElem() {
		return ssElem;
	}

	public void setSsElem(SecStrucElement ssElem) {
		this.ssElem = ssElem;
	}

	@Override
	public double getRsa() {
		return rsa;
	}

	@Override
	public void setRsa(double rsa) {
		this.rsa = rsa;
	}
	
	@Override
	public double getScRsa() {
		return scRsa;
	}

	@Override
	public void setScRsa(double scrsa) {
		this.scRsa = scrsa;
	}

	@Override
	public double getAsa() {
		return asa;
	}

	@Override
	public void setAsa(double asa) {
		this.asa = asa;
	}

	@Override
	public double getBsa() {
		return bsa;
	}
	
	@Override
	public void setBsa(double bsa) {
		this.bsa = bsa;
	}

	@Override
	public double getBsaToAsaRatio() {
		return (double)this.bsa/ (double)this.asa;
	}
	
	@Override
	public Collection<Atom> getAtoms() {
		return atoms.values();
	}
	
	/**
	 * Get Collection of Atoms belonging to this Residue sorted by atom codes
	 * @return
	 */
	public TreeMap<String,Atom> getAtomsMap() {
		return atoms;
	}
	
	/**
	 * Gets a new Residue object containing only the atoms for the given contact type
	 * (the new atom instances are references to the atom instances of this Residue)
	 * @param ct
	 * @return
	 */
	public AaResidue getReducedResidue(String ct) {
		AaResidue reducedResidue = new AaResidue(this.getAaType(), this.getSerial(), this.getParent());
		Set<String> atomCodes = ContactType.getAtomsForCTAndRes(ct, this.getAaType().getThreeLetterCode());
		for (String atomCode:atomCodes) {
			if (this.containsAtom(atomCode)) {
				reducedResidue.addAtom(this.getAtom(atomCode));
			}
		}
		// in cts ("ALL","BB") we still miss the OXT, we need to add it now if it is there (it will be there when this resser is the last residue)
		if ((ct.equals("ALL") || ct.equals("BB")) &&  this.containsOXT()) { 
			reducedResidue.addAtom(this.getAtom("OXT"));
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
		AaResidue scOnly = getReducedResidue("SC");
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
	
	/**
	 * Finds out wheter this AaResidue is of the L-form or the D-form or has 
	 * no chiral center, either because it is a Gly or because one of the 4 atoms 
	 * around the CA is missing: C, CB, HA or N 
	 * See Chapter "The Mathematics of Chirality", Distance Geometry and Molecular Conformation,
	 * G.M. Crippen, T.M. Havel
	 * @return
	 */
	public Chirality getChirality() {
		if (!containsAtom("C") || !containsAtom("CB") || !containsAtom("HA") || !containsAtom("N")) {
			return Chirality.U;
		}
		Point3d h = getAtom("HA").getCoords();
		Point3d c = getAtom("C").getCoords();
		Point3d r = getAtom("CB").getCoords();
		Point3d n = getAtom("N").getCoords();
		
		// see equation 2.1 of book mentioned above
		Matrix4d oriVolMat = new Matrix4d(1,   1,   1,   1, 
										h.x, c.x, r.x, n.x, 
										h.y, c.y, r.y, n.y,
										h.z, c.z, r.z, n.z);
		double oriVol = oriVolMat.determinant();
		if (oriVol>0) {
			return Chirality.L;
		} else if (oriVol<0) {
			return Chirality.D;
		} else {
			// this shouldn't ever happen because oriVol can't really be 0 exactly
			// anyway we keep this just in case
			return Chirality.C;
		}
	}
	
	/**
	 * Tells whether this residue is contiguous in the chain to given i+1
	 * residue, i.e. whether they form a peptide bond.
	 * It tests whether this residue's C atom is below 1.4A of iPlus1Residue's N atom.
	 * @param iPlus1Residue
	 * @return true if they are contiguous, false otherwise or if this residue does 
	 * not have a C atom or iPlus1residue does not have a N atom
	 */
	public boolean isContiguous(AaResidue iPlus1Residue) {
		if (!this.containsAtom("C") || !iPlus1Residue.containsAtom("N")) {
			return false;
		}
		if (this.getAtom("C").getCoords().distance(iPlus1Residue.getAtom("N").getCoords())<1.4) {
			return true;
		}
		return false;
	}
	
	public String toString() {
		return this.getChainCode()+serial+aaType.getOneLetterCode();
	}

	@Override
	public Iterator<Atom> iterator() {
		return atoms.values().iterator();
	}

	@Override
	public AaResidue copy(PdbChain parentPdb) {
		AaResidue newResidue = new AaResidue(this.aaType, this.serial, parentPdb);
		newResidue.pdbSerial = this.pdbSerial;
		newResidue.rsa = this.rsa;
		newResidue.asa = this.asa;
		newResidue.bsa = this.bsa;
		newResidue.scRsa = this.scRsa;
		if (parentPdb.getSecondaryStructure()!=null) {
			newResidue.ssElem = parentPdb.getSecondaryStructure().getSecStrucElement(serial);
		}
		newResidue.atoms = new TreeMap<String, Atom>();
		for (String code:this.atoms.keySet()) {
			newResidue.atoms.put(code,this.atoms.get(code).copy(newResidue));
		}
 
		return newResidue;
	}

	@Override
	public String getLongCode() {
		return aaType.getThreeLetterCode();
	}

	@Override
	public char getShortCode() {
		return aaType.getOneLetterCode();
	}
	
	@Override
	public boolean isPeptideLinked() {
		return true;
	}

	@Override
	public PdbChain getParent() {
		return parent;
	}
	
	@Override
	public void removeHatoms() {
		List<String> toRemove = new ArrayList<String>();
		for (Atom atom:this.atoms.values()) {
			if (atom.getType()==AtomType.H) {
				toRemove.add(atom.getCode());
			}
		}
		for (String atomCode:toRemove) {
			atoms.remove(atomCode);
		}
	}
}

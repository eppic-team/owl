package owl.core.structure;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import owl.core.structure.features.SecStrucElement;

/**
 * A HET residue in a PDB structure i.e. residues that come from HETATM records in PDB:
 * ligands, non-standard amino acids, ions etc.
 * 
 * @author duarte_j
 *
 */
public class HetResidue implements Residue {

	private static final long serialVersionUID = 1L;
	
	// a (non-comprehensive) list of het residues that are within protein chains (peptide-linked)
	// some parents: 
	// MSE:MET NLE:LEU SMC:CYS AIB:ALA ABA:ALA
	private static final String[] BACKBONE_HET_RESIDUES = 
		{"ACE", "NH2", "SUI", "PYR", "GL3", "MSE", "SNN", "CRO", "AKZ", "GLK", "LLP", "NLE", "SMC", "AIB", "ABA", "L2O"};
	
	private String mol3letterCode;
	
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
	public HetResidue() {
		atoms = new TreeMap<String, Atom>();
	}
	
	/**
	 * Constructs a new Residue given its HET molecule 3-letter code, serial and parentPdb. The Residue will 
	 * have an empty list of Atoms until they are added using {@link #addAtom(Atom)}  
	 * @param mol3letterCode
	 * @param serial
	 * @param parentPdb
	 */
	public HetResidue(String mol3letterCode, int serial, PdbChain parent) {
		atoms = new TreeMap<String, Atom>();
		this.mol3letterCode = mol3letterCode;
		this.serial = serial;
		this.parent = parent;
	}

	@Override
	public Iterator<Atom> iterator() {
		return atoms.values().iterator();
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

	@Override
	public String getLongCode() {
		return mol3letterCode;
	}
	
	@Override
	public char getShortCode() {
		return AminoAcid.XXX.getOneLetterCode();
	}

	@Override
	public boolean isPeptideLinked() {
		
		for (String backboneRes:BACKBONE_HET_RESIDUES) {
			if (mol3letterCode.equals(backboneRes)) {
				return true;
			}
		}
		
		// this works for non-standard residues but not for other peptide linked residues like ACE or SUI
		return (this.containsAtom("C") &&
				this.containsAtom("N") &&
				//this.containsAtom("O") && // in principle this atom is not strictly necessary, e.g. in GL3 where it is a S
				this.containsAtom("CA"));
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

	@Override
	public PdbChain getParent() {
		return parent;
	}

	@Override
	public Collection<Atom> getAtoms() {
		return atoms.values();
	}

	@Override
	public Residue copy(PdbChain parent) {
		HetResidue newResidue = new HetResidue(this.mol3letterCode, this.serial, parent);
		newResidue.pdbSerial = this.pdbSerial;
		newResidue.asa = this.asa;
		newResidue.bsa = this.bsa;
		newResidue.atoms = new TreeMap<String, Atom>();
		for (String code:this.atoms.keySet()) {
			newResidue.atoms.put(code,this.atoms.get(code).copy(newResidue));
		}
 
		return newResidue;
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
	public void setAtomRadii() {
		for (Atom atom:this) {
			atom.setRadius(AtomRadii.getRadius(this.mol3letterCode, atom));
		}
	}
	
	public String toString() {
		return this.parent.getChainCode()+serial+getLongCode();
	}
	
	@Override
	public SecStrucElement getSsElem() {
		return ssElem;
	}

	/**
	 * Sets the secondary structure element for this HetResidue
	 * @param ssElem
	 * @throws IllegalArgumentException if this residue is not part of a protein chain
	 */
	public void setSsElem(SecStrucElement ssElem) {
		if (!this.parent.getSequence().isProtein()) {
			throw new IllegalArgumentException("Residue "+this+" is not part of a protein chain, can't add a secondary structure element to it.");
		}
		this.ssElem = ssElem;
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

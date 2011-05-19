package owl.core.structure;

import java.util.Collection;
import java.util.Iterator;
import java.util.TreeMap;

import owl.core.structure.features.SecStrucElement;

public class NucResidue implements Residue {

	private static final long serialVersionUID = 1L;
	
	private Nucleotide nucType;
	private int serial;
	private String pdbSerial;
	
	private PdbChain parent;
	
	private double rsa;   // relative all-atoms accessible surface area
	private double scRsa; // relative side chain-atoms accessible surface area
	private double asa;   // all-atoms accessible surface area (absolute)
	private double bsa;   // all-atoms accessible surface area upon burial (this can only be calculated by running naccess for 2 partners and a complex, see ChainInterface class)
	
	private TreeMap<String, Atom> atoms; // atom codes to atoms
	
	/**
	 * Constructs an empty Residue. Use setters to add content to it.
	 */
	public NucResidue() {
		atoms = new TreeMap<String, Atom>();
	}
	
	/**
	 * Constructs a new Residue given its Nucleotide type, serial and parentPdb. The Residue will 
	 * have an empty list of Atoms until they are added using {@link #addAtom(Atom)}  
	 * @param nucType
	 * @param serial
	 * @param parentPdb
	 */
	public NucResidue(Nucleotide nucType, int serial, PdbChain parentPdb) {
		atoms = new TreeMap<String, Atom>();
		this.nucType = nucType;
		this.serial = serial;
		this.parent = parentPdb;
	}
	
	public Nucleotide getNucType(){
		return this.nucType;
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
		return this.nucType.getTwoLetterCode();
	}

	@Override
	public char getShortCode() {
		return this.nucType.getOneLetterCode();
	}
	
	@Override
	public boolean isPeptideLinked() {
		return false;
	}

	@Override
	public int getSerial() {
		return this.serial;
	}

	@Override
	public void setSerial(int serial) {
		this.serial = serial;
	}

	@Override
	public String getPdbSerial() {
		return this.pdbSerial;
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
		NucResidue newResidue = new NucResidue(this.nucType, this.serial, parent);
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
			atom.setRadius(AtomRadii.getRadius(this.nucType, atom));
		}
	}
	
	@Override
	public SecStrucElement getSsElem() {
		return null;
	}
	
	public String toString() {
		return this.parent.getChainCode()+serial+nucType.getOneLetterCode();
	}

}

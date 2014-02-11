package owl.core.structure;

import java.io.Serializable;
import java.util.Collection;

import owl.core.runners.NaccessRunner;
import owl.core.structure.features.SecStrucElement;

/**
 * A residue as a part of a protein structure: a group of atoms that form a distinct
 * chemical entity.
 * 
 * @author duarte_j
 *
 */
public interface Residue extends Iterable<Atom>, Serializable {

	/**
	 * Adds an atom to this Residue. If the atom code for given Atom is already present
	 * in this Residue a warning will be printed and no atom added. 
	 * No duplicate atom codes are allowed in one Residue.
	 * @param atom
	 */
	public void addAtom(Atom atom);
	
	/**
	 * Returns the Atom object given its atomCode or if no such atom for this residue returns null.
	 * @param atomCode
	 * @return the Atom object or null if no such atom exists in this residue
	 */
	public Atom getAtom(String atomCode);
	
	/**
	 * Returns true if this Residue instance has the given atomCode
	 * @param atomCode
	 * @return
	 */
	public boolean containsAtom(String atomCode);
	
	/**
	 * Returns the number of atoms in this Residue, including Hydrogens if they are 
	 * present
	 * @return
	 */
	public int getNumAtoms();
	
	/**
	 * Returns the number of heavy (non-Hydrogen) atoms in this Residue
	 * @return
	 */
	public int getNumHeavyAtoms();
	
	/**
	 * Returns this residue's 3 letter code (PDB standard) identifying the type of molecule.
	 * For nucleotides it is only 2 letters rather than 3.
	 * See the PDB chemical component dictionary ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz 
	 * for all possible codes. Nice web interfaces to it are:
	 * - PDBeChem (http://www.ebi.ac.uk/msd-srv/msdchem/cgi-bin/cgi.pl)
	 * - PDB ligand expo browser (http://ligand-expo.rcsb.org) 
	 * @return
	 */
	public String getLongCode();
	
	/**
	 * Returns this residue's 1 letter code (PDB standard) identifying the type of molecule.
	 * If HetResidue then an X is returned.
	 * @return
	 */
	public char getShortCode();
	
	/**
	 * Tells whether this Residue is peptide linked to the rest of the chain.
	 * For nucleotides it always returns false, for standard amino acid always true and
	 * for HET residues it tries to guess whether it is peptide linked by using 
	 * a (non-comprehensive) list of known peptide linked compounds or from the
	 * presence of backbone atoms C, N and CA. 
	 * For the case of HET residues it is safer to check whether the parent PdbChain 
	 * is a protein chain (all HETs that are not part of a chain should be assigned to 
	 * their independent non-polymer chains).
	 * @return
	 */
	public boolean isPeptideLinked();
	
	/**
	 * Gets this Residue's serial (corresponding to SEQRES sequence)
	 * @return
	 */
	public int getSerial();

	/**
	 * Sets this Residue's serial (corresponding to SEQRES sequence)
	 * @param serial
	 */
	public void setSerial(int serial);
	
	/**
	 * Gets this Residue's PDB serial (can have insertion codes and be negative or 0)
	 * @return
	 */
	public String getPdbSerial();

	/**
	 * Sets this Residue's PDB serial (can have insertion codes and be negative or 0)
	 * @param pdbSerial
	 */
	public void setPdbSerial(String pdbSerial);
	
	/**
	 * Gets this Residue's parent PdbChain
	 * @return
	 */
	public PdbChain getParent();
	
	/**
	 * Get Collection of Atoms belonging to this Residue sorted by atom codes
	 * @return
	 */
	public Collection<Atom> getAtoms();
	
	/**
	 * Returns a deep copy of this Residue
	 * @return
	 */
	public Residue copy(PdbChain parent);
	
	/**
	 * Returns the absolute all-atoms accessible surface area (square Angstroms) 
	 * @see AsaCalculator
	 * @see NaccessRunner
	 * @return
	 */
	public double getAsa();
	
	/**
	 * Sets the absolute all-atoms acessible surface area (square Angstroms)
	 * @param asa
	 */
	public void setAsa(double asa);
	
	/**
	 * Returns the absolute all-atoms accessible surface area upon burial (square Angstroms)
	 * This can be calculated by computing asa for 2 partners separately and the complex of both.
	 * @see ChainInterface
	 * @see AsaCalculator
	 * @see NaccessRunner
	 * @return
	 */
	public double getBsa();
	
	/**
	 * Sets the absolute all-atoms accessible surface area upon burial (square Angstroms)
	 * @param bsa
	 */
	public void setBsa(double bsa);
	
	/**
	 * Returns the bsa to asa ratio (bsa/asa)
	 * @return
	 */
	public double getBsaToAsaRatio();
	
	/**
	 * Returns the relative all-atoms relative accessible surface area as calculated by NACCESS 
	 * @return the rsa value or null if NACCESS has not been run
	 */
	public double getRsa();

	/**
	 * Sets the all-atoms relative accessible surface area
	 * @param rsa
	 */
	public void setRsa(double rsa);
	
	/**
	 * Returns the relative sidechain-atoms accessible surface area as calculated by NACCESS 
	 * @return the sc-rsa value or null if NACCESS has not been run
	 */
	public double getScRsa();

	/**
	 * Sets the relative sidechain-atoms accessible surface area
	 * @param scrsa
	 */
	public void setScRsa(double scrsa);
	
	/**
	 * Returns the secondary structure element to which this Residue belongs. 
	 * It returns always null for HetResidues not peptide linked and for NucResidues
	 * @return
	 */
	public SecStrucElement getSsElem();
	
	/**
	 * Removes all Hydrogen atoms from this Residue
	 */
	public void removeHatoms();
}

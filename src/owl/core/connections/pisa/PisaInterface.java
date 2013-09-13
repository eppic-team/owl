package owl.core.connections.pisa;

import java.io.PrintStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.ChainInterface;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.Residue;

public class PisaInterface implements Comparable<PisaInterface> {

	private int id;						// note this is call NN in the PISA web tables
	private int type;					// note this is call id in the PISA web tables
	private double interfaceArea;
	private double solvEnergy;
	private double solvEnergyPvalue;
	
	private PisaMolecule firstMolecule;
	private PisaMolecule secondMolecule;
	
	// our own fields
	private int protprotId; 		// the serial number of this interface if only protein-protein interfaces are counted
	
	public PisaInterface() {
		
	}

	public void printTabular(PrintStream ps) {
		ps.print("# ");
		ps.printf("%d\t%d\t%9.2f\t%5.2f\t%4.2f\n",this.getId(),this.getType(),this.getInterfaceArea(),this.getSolvEnergy(),this.getSolvEnergyPvalue());
		ps.print("## ");
		this.getFirstMolecule().printTabular(ps);
		ps.print("## ");
		this.getSecondMolecule().printTabular(ps);
	}	
	
	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}

	/**
	 * @param id the id to set
	 */
	public void setId(int id) {
		this.id = id;
	}

	/**
	 * @return the type
	 */
	public int getType() {
		return type;
	}

	/**
	 * @param type the type to set
	 */
	public void setType(int type) {
		this.type = type;
	}

	/**
	 * @return the interfaceArea
	 */
	public double getInterfaceArea() {
		return interfaceArea;
	}

	/**
	 * @param interfaceArea the interfaceArea to set
	 */
	public void setInterfaceArea(double interfaceArea) {
		this.interfaceArea = interfaceArea;
	}

	/**
	 * @return the solvEnergy
	 */
	public double getSolvEnergy() {
		return solvEnergy;
	}

	/**
	 * @param solvEnergy the solvEnergy to set
	 */
	public void setSolvEnergy(double solvEnergy) {
		this.solvEnergy = solvEnergy;
	}

	/**
	 * @return the solvEnergyPvalue
	 */
	public double getSolvEnergyPvalue() {
		return solvEnergyPvalue;
	}

	/**
	 * @param solvEnergyPvalue the solvEnergyPvalue to set
	 */
	public void setSolvEnergyPvalue(double solvEnergyPvalue) {
		this.solvEnergyPvalue = solvEnergyPvalue;
	}

	/**
	 * @return the firstMolecule
	 */
	public PisaMolecule getFirstMolecule() {
		return firstMolecule;
	}

	/**
	 * @param firstMolecule the firstMolecule to set
	 */
	public void setFirstMolecule(PisaMolecule firstMolecule) {
		this.firstMolecule = firstMolecule;
	}

	/**
	 * @return the secondMolecule
	 */
	public PisaMolecule getSecondMolecule() {
		return secondMolecule;
	}

	/**
	 * @param secondMolecule the secondMolecule to set
	 */
	public void setSecondMolecule(PisaMolecule secondMolecule) {
		this.secondMolecule = secondMolecule;
	}
	
	/**
	 * Returns true if both members of this interface are proteins.
	 * @return
	 */
	public boolean isProtein() {
		return (firstMolecule.isProtein() && secondMolecule.isProtein());
	}

	/**
	 * Return the serial number of this interface if only protein-protein interfaces are
	 * counted. 1 would be the top prot-prot interface in the PISA list (usually biggest prot-prot interface) 
	 * @return
	 */
	public int getProtProtId() {
		return protprotId;
	}
	
	/**
	 * Sets the protprotId
	 * @param protprotId
	 */
	public void setProtProtId(int protprotId) {
		this.protprotId = protprotId;
	}
	
	@Override
	public int compareTo(PisaInterface o) {
		return new Double(this.getInterfaceArea()).compareTo(o.getInterfaceArea());
	}
	
	/**
	 * Converts a PisaInterface into our own ChainInterface which contains the full coordinates
	 * of the PDB entry with properly transformed chains.
	 * The chains read from the given pdb are deep-copied and then transformed so that they stay
	 * independent from the input pdb. 
	 * @param pdb the externally read PDB entry (cif file, pdb file, pdbase)
	 * @return
	 */
	public ChainInterface convertToChainInterface(PdbAsymUnit pdb) {
		ChainInterface interf = new ChainInterface();
		interf.setInterfaceArea(this.interfaceArea);
		interf.setId(this.id);
		interf.setFirstTransf(firstMolecule.getTransf());
		interf.setFirstTransfOrth(firstMolecule.getTransfOrth());
		interf.setSecondTransf(secondMolecule.getTransf());
		interf.setSecondTransfOrth(secondMolecule.getTransfOrth());
		interf.setName(firstMolecule.getChainId()+"+"+secondMolecule.getChainId());
		
		PdbChain pdb1 = findChainForPisaMolecule(this.firstMolecule, pdb).copy(pdb);
		pdb1.transform(firstMolecule.getTransfOrth());
		// this might be confusing: the setAsaAndBsas methods sets the asa/bsa values of the passed pdb object from the PisaMolecule
		firstMolecule.setAsaAndBsas(pdb1);
		interf.setFirstMolecule(pdb1);
		PdbChain pdb2 = findChainForPisaMolecule(this.secondMolecule, pdb).copy(pdb);
		pdb2.transform(secondMolecule.getTransfOrth());
		// this might be confusing: the setAsaAndBsas methods sets the asa/bsa values of the passed pdb object from the PisaMolecule
		secondMolecule.setAsaAndBsas(pdb2);
		interf.setSecondMolecule(pdb2);		
		
		return interf;
	}
	
	/**
	 * For a given PisaMolecule and a corresponding PdbAsymUnit it finds what is the 
	 * PdbChain corresponding to the PisaMolecule: if is a protein/nucleotide chain it is straight 
	 * forward from the PDB chain code, but if it is a non-polymer chain then it has to find the corresponding
	 * chain by matching the 3 parts of the PISA identifier (PDB chain code, residue type and PDB residue serial). 
	 * @param molecule
	 * @param pdb the chain or null if nothing found (a warning is printed as well)
	 * @return
	 */
	private PdbChain findChainForPisaMolecule(PisaMolecule molecule, PdbAsymUnit pdb) {
		if (molecule.isProtein()) {
			return pdb.getChain(molecule.getChainId());
		}
		String pisaNonPolyChainId = molecule.getChainId();
		Pattern p = Pattern.compile("^\\[(\\w+)\\](\\w):(\\d+)$");
		Matcher m = p.matcher(pisaNonPolyChainId);
		String resCode = null;
		String chain = null;
		String pdbResSerial = null;
		if (m.matches()){
			resCode = m.group(1).trim();
			chain = m.group(2);
			pdbResSerial = m.group(3);
		}
		PdbChain selectedChain = null; // here we put the chain that we think matches the pisa non-poly chain identifier 
		for (PdbChain nonPolyChain: pdb.getNonPolyChains()){
			if (nonPolyChain.getPdbChainCode().equals(chain)) {
				for (Residue res:nonPolyChain) {
					if (res.getLongCode().equals(resCode) && res.getPdbSerial().equals(pdbResSerial)) {
						selectedChain = nonPolyChain;
						break;
					}
				}
			}
		}
		if (selectedChain==null) System.err.println("Warning! couldn't find a corresponding non-polymer chain for PISA identifier "+pisaNonPolyChainId);
		return selectedChain;
	}
}

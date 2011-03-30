package owl.core.connections.pisa;

import java.io.PrintStream;

import owl.core.structure.ChainInterface;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;

public class PisaInterface implements Comparable<PisaInterface> {

	private int id;						// note this is call NN in the PISA web tables
	private int type;					// note this is call id in the PISA web tables
	private double interfaceArea;
	private double solvEnergy;
	private double solvEnergyPvalue;
	
	private PisaMolecule firstMolecule;
	private PisaMolecule secondMolecule;
	
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
		interf.setScore(this.solvEnergy);
		interf.setInterfaceArea(this.interfaceArea);
		interf.setId(this.id);
		interf.setFirstTransf(firstMolecule.getTransf());
		interf.setFirstTransfOrth(firstMolecule.getTransfOrth());
		interf.setSecondTransf(secondMolecule.getTransf());
		interf.setSecondTransfOrth(secondMolecule.getTransfOrth());
		interf.setFirstMolType(firstMolecule.getMolClass());
		interf.setSecondMolType(secondMolecule.getMolClass());
		interf.setName(firstMolecule.getChainId()+"+"+secondMolecule.getChainId());
		
		if (firstMolecule.isProtein()) {
			PdbChain pdb1 = pdb.getChain(this.firstMolecule.getChainId()).copy(pdb);
			pdb1.transform(firstMolecule.getTransfOrth());
			// this might be confusing: the setAsaAndBsas methods sets the asa/bsa values of the passed pdb object from the PisaMolecule
			firstMolecule.setAsaAndBsas(pdb1);
			interf.setFirstMolecule(pdb1);
		}
		if (secondMolecule.isProtein()) {
			PdbChain pdb2 = pdb.getChain(secondMolecule.getChainId()).copy(pdb);
			pdb2.transform(secondMolecule.getTransfOrth());
			// this might be confusing: the setAsaAndBsas methods sets the asa/bsa values of the passed pdb object from the PisaMolecule
			secondMolecule.setAsaAndBsas(pdb2);
			interf.setSecondMolecule(pdb2);
		}
		
		return interf;
	}
	
	
}

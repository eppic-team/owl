package owl.core.connections.pisa;

import java.io.PrintStream;

public class PisaInterface {

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
	
}

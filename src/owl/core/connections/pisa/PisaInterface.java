package owl.core.connections.pisa;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

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
	 * Returns a map containing 2 {@link PisaRimCore} objects (see getRimAndCore in {@link PisaMolecule})
	 * for each of the 2 members of the interface.
	 * The sum of the residues of the 2 cores is required to be at least minNumResidues. 
	 * If the minimum is not reached with the bsaToAsaSoftCutoff, then the cutoff is 
	 * relaxed in relaxationStep steps until reaching the bsaToAsaHardCutoff.
	 * If either of the 2 molecules of this interface is not a protein, its rimCore 
	 * object in the output map object will be null. If both are not proteins then the map 
	 * will contain null object references for both.
	 * @param bsaToAsaSoftCutoff
	 * @param bsaToAsaHardCutoff
	 * @param relaxationStep
	 * @param minNumResidues
	 * @return
	 */
	public Map<Integer,PisaRimCore> getRimAndCore(double bsaToAsaSoftCutoff, double bsaToAsaHardCutoff, double relaxationStep, int minNumResidues) {
		
		Map<Integer,PisaRimCore> rimcores = new HashMap<Integer, PisaRimCore>();
		if (!firstMolecule.isProtein() && !secondMolecule.isProtein()) {
			rimcores.put(1, null);
			rimcores.put(2, null);
			return rimcores;
		}
		
		// we introduce a margin of relaxationSte*0.10 to be sure we do go all the way down to bsaToAsaHardCutoff (necessary because of rounding)
		for (double cutoff=bsaToAsaSoftCutoff;cutoff>=bsaToAsaHardCutoff-relaxationStep*0.10;cutoff-=relaxationStep) {
			PisaRimCore rimCore1 = null;
			PisaRimCore rimCore2 = null;
			if (this.firstMolecule.isProtein()) rimCore1 = this.firstMolecule.getRimAndCore(cutoff);
			if (this.secondMolecule.isProtein()) rimCore2 = this.secondMolecule.getRimAndCore(cutoff);
			rimcores.put(1,rimCore1);
			rimcores.put(2,rimCore2);
			
			int totalCoreResidues = 0;
			if (firstMolecule.isProtein()) totalCoreResidues+=rimCore1.getCoreSize();
			if (secondMolecule.isProtein()) totalCoreResidues+=rimCore2.getCoreSize();
			if (totalCoreResidues>=minNumResidues) {
				break;
			}
		}
		
		return rimcores;
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
	
	
}

package owl.core.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class InterfaceRimCore implements Serializable {

	private static final long serialVersionUID = -1931015962640626829L;

	private List<Residue> rimResidues;
	private List<Residue> coreResidues;
	private double bsaToAsaCutoff;
	private double minAsaForSurface;
	
	public InterfaceRimCore() {
		this.coreResidues = new ArrayList<Residue>();
		this.rimResidues = new ArrayList<Residue>();
	}
	
	public InterfaceRimCore(List<Residue> rimResidues, List<Residue> coreResidues, double bsaToAsaCutoff, double minAsaForSurface) {
		this.rimResidues = rimResidues;
		this.coreResidues = coreResidues;
		this.bsaToAsaCutoff = bsaToAsaCutoff;
		this.minAsaForSurface = minAsaForSurface;
	}
	
	public void addCoreResidue(Residue residue) {
		coreResidues.add(residue);
	}
	
	public void addRimResidue(Residue residue) {
		rimResidues.add(residue);
	}

	/**
	 * @return the rimResidues
	 */
	public List<Residue> getRimResidues() {
		return rimResidues;
	}
	
	/**
	 * @param rimResidues the rimResidues to set
	 */
	public void setRimResidues(List<Residue> rimResidues) {
		this.rimResidues = rimResidues;
	}

	/**
	 * @return the coreResidues
	 */
	public List<Residue> getCoreResidues() {
		return coreResidues;
	}

	/**
	 * @param coreResidues the coreResidues to set
	 */
	public void setCoreResidues(List<Residue> coreResidues) {
		this.coreResidues = coreResidues;
	}
	
	/**
	 * @return the bsaToAsaCutoff
	 */
	public double getBsaToAsaCutoff() {
		return bsaToAsaCutoff;
	}

	/**
	 * @param bsaToAsaCutoff the bsaToAsaCutoff to set
	 */
	public void setBsaToAsaCutoff(double bsaToAsaCutoff) {
		this.bsaToAsaCutoff = bsaToAsaCutoff;
	}
	
	public double getMinAsaForSurface() {
		return minAsaForSurface;
	}

	public int getCoreSize() {
		return this.coreResidues.size();
	}
	
	public int getRimSize() {
		return this.rimResidues.size();
	}

	public double getAsaRim() {
		double asa = 0;
		for (Residue res:rimResidues) {
			asa+=res.getAsa();
		}
		return asa;
	}
	
	public double getAsaCore() {
		double asa = 0;
		for (Residue res:coreResidues) {
			asa+=res.getAsa();
		}
		return asa;		
	}

	public double getBsaRim() {
		double asa = 0;
		for (Residue res:rimResidues) {
			asa+=res.getBsa();
		}
		return asa;
	}
	
	public double getBsaCore() {
		double asa = 0;
		for (Residue res:coreResidues) {
			asa+=res.getBsa();
		}
		return asa;		
	}

}

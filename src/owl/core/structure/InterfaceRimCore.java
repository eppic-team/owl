package owl.core.structure;

import java.io.Serializable;
import java.util.List;

public class InterfaceRimCore implements Serializable {

	private static final long serialVersionUID = -1931015962640626829L;

	private List<Residue> rimResidues;
	private List<Residue> coreResidues;
	private double bsaToAsaCutoff;
	
	public InterfaceRimCore(List<Residue> rimResidues, List<Residue> coreResidues, double bsaToAsaCutoff) {
		this.rimResidues = rimResidues;
		this.coreResidues = coreResidues;
		this.bsaToAsaCutoff = bsaToAsaCutoff;
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

	public int getCoreSize() {
		return this.coreResidues.size();
	}
	
	public int getRimSize() {
		return this.rimResidues.size();
	}

}

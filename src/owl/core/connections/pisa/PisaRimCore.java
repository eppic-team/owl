package owl.core.connections.pisa;

import java.io.Serializable;
import java.util.List;

public class PisaRimCore implements Serializable {

	private static final long serialVersionUID = -1931015962640626829L;

	private List<PisaResidue> rimResidues;
	private List<PisaResidue> coreResidues;
	private double bsaToAsaCutoff;
	
	public PisaRimCore(List<PisaResidue> rimResidues, List<PisaResidue> coreResidues, double bsaToAsaCutoff) {
		this.rimResidues = rimResidues;
		this.coreResidues = coreResidues;
		this.bsaToAsaCutoff = bsaToAsaCutoff;
	}

	/**
	 * @return the rimResidues
	 */
	public List<PisaResidue> getRimResidues() {
		return rimResidues;
	}

	/**
	 * @param rimResidues the rimResidues to set
	 */
	public void setRimResidues(List<PisaResidue> rimResidues) {
		this.rimResidues = rimResidues;
	}

	/**
	 * @return the coreResidues
	 */
	public List<PisaResidue> getCoreResidues() {
		return coreResidues;
	}

	/**
	 * @param coreResidues the coreResidues to set
	 */
	public void setCoreResidues(List<PisaResidue> coreResidues) {
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

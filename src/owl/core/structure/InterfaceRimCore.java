package owl.core.structure;

import java.io.Serializable;
import java.util.List;

public class InterfaceRimCore implements Serializable {

	private static final long serialVersionUID = -1931015962640626829L;

	private List<AaResidue> rimResidues;
	private List<AaResidue> coreResidues;
	private double bsaToAsaCutoff;
	
	public InterfaceRimCore(List<AaResidue> rimResidues, List<AaResidue> coreResidues, double bsaToAsaCutoff) {
		this.rimResidues = rimResidues;
		this.coreResidues = coreResidues;
		this.bsaToAsaCutoff = bsaToAsaCutoff;
	}

	/**
	 * @return the rimResidues
	 */
	public List<AaResidue> getRimResidues() {
		return rimResidues;
	}

	/**
	 * @param rimResidues the rimResidues to set
	 */
	public void setRimResidues(List<AaResidue> rimResidues) {
		this.rimResidues = rimResidues;
	}

	/**
	 * @return the coreResidues
	 */
	public List<AaResidue> getCoreResidues() {
		return coreResidues;
	}

	/**
	 * @param coreResidues the coreResidues to set
	 */
	public void setCoreResidues(List<AaResidue> coreResidues) {
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

	public String getRimResString() {
		return getResString(rimResidues);
	}
	
	public String getCoreResString() {
		return getResString(coreResidues);
	}
	
	private static String getResString(List<AaResidue> residues) {
		String str = "";
		for (int i=0;i<residues.size();i++) {
			if (i!=residues.size()-1)
				str+=residues.get(i).getSerial()+",";
			else
				str+=residues.get(i).getSerial();
		}
		return str;
	}
	
}

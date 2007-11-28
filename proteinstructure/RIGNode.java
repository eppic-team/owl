package proteinstructure;

/**
 * Class representing a Residue Interaction Graph node, i.e. a residue
 * 
 */
public class RIGNode {

	private int residueSerial;
	private String residueType;
	private SecStrucElement sselem;
	private boolean observed;
	
	
	public RIGNode(int residueSerial, String residueType, SecStrucElement sselem) {
		this.residueSerial = residueSerial;
		this.residueType = residueType;
		this.sselem = sselem;
		this.observed = true;
	}
	
	public RIGNode(int residueSerial, String residueType) {
		this.residueSerial = residueSerial;
		this.residueType = residueType;
		this.sselem = null;
		this.observed = true;
	}
	
	public RIGNode(int residueSerial) {
		this.residueSerial = residueSerial;
		this.residueType = null;
		this.sselem = null;
		this.observed = true;
	}

	public int getResidueSerial() {
		return residueSerial;
	}

	public String getResidueType() {
		return residueType;
	}
	
	public SecStrucElement getSecStrucElement() {
		return sselem;
	}
	
	public boolean isObserved()	{
		return observed;
	}
	
	public void setObserved(boolean observed) {
		this.observed = observed;
	}
	
	public String toString() {
		return "V"+residueSerial;
	}
	
	/**
	 * Deep copies this RIGNode
	 * @return
	 */
	public RIGNode copy() {
		SecStrucElement newsselem = this.sselem==null?null:this.sselem.copy();
		RIGNode copy = new RIGNode(this.residueSerial, this.residueType, newsselem);
		return copy;
	}
	
}

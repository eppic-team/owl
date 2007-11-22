package proteinstructure;

/**
 * Class representing a Residue Interaction Graph node, i.e. a residue
 * 
 */
public class RIGNode {

	private int residueSerial;
	private String residueType;
	private SecStrucElement sselem;
	
	
	public RIGNode(int residueSerial, String residueType, SecStrucElement sselem) {
		this.residueSerial = residueSerial;
		this.residueType = residueType;
		this.sselem = sselem;
	}
	
	public RIGNode(int residueSerial, String residueType) {
		this.residueSerial = residueSerial;
		this.residueType = residueType;
		this.sselem = null;
	}
	
	public RIGNode(int residueSerial) {
		this.residueSerial = residueSerial;
		this.residueType = null;
		this.sselem = null;
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
	
	public String toString() {
		return "V"+residueSerial;
	}
	
	/**
	 * Deep copies this RIGNode
	 * @return
	 */
	public RIGNode copy() {
		RIGNode copy = new RIGNode(this.residueSerial, this.residueType, this.sselem.copy());
		return copy;
	}
	
}

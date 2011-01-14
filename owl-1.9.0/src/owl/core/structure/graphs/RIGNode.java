package owl.core.structure.graphs;

import owl.core.structure.AAinfo;
import owl.core.structure.features.SecStrucElement;

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
		return residueSerial+AAinfo.threeletter2oneletter(residueType);
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
	
	public boolean equals(Object obj) {
		if (obj instanceof RIGNode) {
			RIGNode v = (RIGNode) obj;
			return ((this.residueSerial == v.getResidueSerial()) && (this.residueType.equals(v.getResidueType())));
		} else {
			return false;
		}
	}
	
}

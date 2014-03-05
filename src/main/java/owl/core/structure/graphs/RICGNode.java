package owl.core.structure.graphs;

import java.io.Serializable;

import owl.core.structure.AminoAcid;

public class RICGNode implements Serializable {

	private static final long serialVersionUID = 1L;
	
	private String pdbChainCode;
	private int residueSerial;
	private String residueType;
	
	public RICGNode(int residueSerial, String residueType, String pdbChainCode) {
		this.residueSerial = residueSerial;
		this.residueType = residueType;
		this.pdbChainCode = pdbChainCode;
	}

	public int getResidueSerial() {
		return residueSerial;
	}

	public void setResidueSerial(int residueSerial) {
		this.residueSerial = residueSerial;
	}

	public String getResidueType() {
		return residueType;
	}

	public void setResidueType(String residueType) {
		this.residueType = residueType;
	}

	public String getPdbChainCode() {
		return pdbChainCode;
	}

	public void setPdbChainCode(String pdbChainCode) {
		this.pdbChainCode = pdbChainCode;
	}
	
	public String toString() {
		return pdbChainCode+"-"+ residueSerial+""+AminoAcid.three2one(residueType);
	}
	
	public boolean equals(Object o) {
		if (o==null) return false;
		if (!(o instanceof RICGNode)) return false;
		
		RICGNode other = (RICGNode) o;
		if (other.residueSerial!=this.residueSerial) return false;
		if (!other.pdbChainCode.equals(this.pdbChainCode)) return false;
		if (!other.residueType.equals(this.residueType)) return false;
		
		return true;
	}
	
	public int hashCode() {
		int hash = residueSerial;
	    hash = hash * 31 + pdbChainCode.hashCode();
	    hash = hash * 31 + residueType.hashCode();
	    return hash;
	}
	
}

package owl.core.structure;

import javax.vecmath.Point3i;

public class SubunitId implements Comparable<SubunitId>{

	private char pdbChainCode;
	private CrystalTransform transform;
	
	public SubunitId(char pdbChainCode, CrystalTransform transform) {
		this.pdbChainCode = pdbChainCode;
		this.transform = transform;
	}
	
	public char getPdbChainCode() {
		return pdbChainCode;
	}
	
	public int getTransformId() {
		return transform.getTransformId();
	}
	
	public Point3i getTransl() {
		return transform.getCrystalTranslation();
	}
	
	public boolean isSymRelatedEquivalent(SubunitId o) {
		if (this.pdbChainCode==o.pdbChainCode && this.transform.getTransformId()==o.transform.getTransformId()) return true;
		return false;
	}
	
	@Override
	public boolean equals(Object o) {
		if (!(o instanceof SubunitId)) return false;
		SubunitId other = (SubunitId) o;
		if (this.pdbChainCode==other.pdbChainCode && 
				this.transform.getTransformId()==other.transform.getTransformId() && 
				this.transform.getCrystalTranslation().equals(other.transform.getCrystalTranslation()))
			return true;
		
		return false;
	}
	
	@Override
	public int hashCode() {
		int hash = pdbChainCode;
		hash = hash * 31 + transform.getTransformId();
		hash = hash * 31 + transform.getCrystalTranslation().hashCode();
	    return hash;
	}
	
	@Override
	public String toString() {
		return this.pdbChainCode+""+getTransformId()+getTransl().toString();
	}

	@Override
	public int compareTo(SubunitId o) {
		if (this.getTransformId()<o.getTransformId()) {
			return -1;
		} 
		if (this.getTransformId()>o.getTransformId()) {
			return 1;
		}
		
		if (this.pdbChainCode<o.pdbChainCode) {
			return -1;
		}
		if (this.pdbChainCode>o.pdbChainCode) {
			return 1;
		}
		if (this.getTransl().x<o.getTransl().x) {
			return -1;
		}
		if (this.getTransl().x>o.getTransl().x) {
			return 1;
		}
		if (this.getTransl().y<o.getTransl().y) {
			return -1;
		}
		if (this.getTransl().y>o.getTransl().y) {
			return 1;
		}
		if (this.getTransl().z<o.getTransl().z) {
			return -1;
		}
		if (this.getTransl().z>o.getTransl().z) {
			return 1;
		}
		
		return 0;
	}
}

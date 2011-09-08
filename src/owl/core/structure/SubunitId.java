package owl.core.structure;

import javax.vecmath.Point3i;

public class SubunitId implements Comparable<SubunitId>{

	private char pdbChainCode;
	private int transformId;
	private Point3i transl;
	
	public SubunitId(char pdbChainCode, int transformId, Point3i transl) {
		this.pdbChainCode = pdbChainCode;
		this.transformId = transformId;
		this.transl = transl;
	}
	
	public char getPdbChainCode() {
		return pdbChainCode;
	}
	
	public int getTransformId() {
		return transformId;
	}
	
	public Point3i getTransl() {
		return transl;
	}
	
	public boolean isSymRelatedEquivalent(SubunitId o) {
		if (this.pdbChainCode==o.pdbChainCode && this.transformId==o.transformId) return true;
		return false;
	}
	
	@Override
	public boolean equals(Object o) {
		if (!(o instanceof SubunitId)) return false;
		SubunitId other = (SubunitId) o;
		if (this.pdbChainCode==other.pdbChainCode && this.transformId==other.transformId && this.transl.equals(other.transl)) return true;
		return false;
	}
	
	@Override
	public int hashCode() {
		int hash = pdbChainCode;
		hash = hash * 31 + transformId;
		hash = hash * 31 + transl.hashCode();
	    return hash;
	}
	
	@Override
	public String toString() {
		return this.pdbChainCode+""+transformId+transl.toString();
	}

	@Override
	public int compareTo(SubunitId o) {
		if (this.transformId<o.transformId) {
			return -1;
		} 
		if (this.transformId>o.transformId) {
			return 1;
		}
		
		if (this.pdbChainCode<o.pdbChainCode) {
			return -1;
		}
		if (this.pdbChainCode>o.pdbChainCode) {
			return 1;
		}
		if (this.transl.x<o.transl.x) {
			return -1;
		}
		if (this.transl.x>o.transl.x) {
			return 1;
		}
		if (this.transl.y<o.transl.y) {
			return -1;
		}
		if (this.transl.y>o.transl.y) {
			return 1;
		}
		if (this.transl.z<o.transl.z) {
			return -1;
		}
		if (this.transl.z>o.transl.z) {
			return 1;
		}
		
		return 0;
	}
}

package owl.core.structure;

public class SubunitId implements Comparable<SubunitId>{

	private char pdbChainCode;
	private int transformId;
	
	public SubunitId(char pdbChainCode, int transformId) {
		this.pdbChainCode = pdbChainCode;
		this.transformId = transformId;
	}
	
	@Override
	public boolean equals(Object o) {
		if (!(o instanceof SubunitId)) return false;
		SubunitId other = (SubunitId) o;
		if (this.pdbChainCode==other.pdbChainCode && this.transformId==other.transformId) return true;
		return false;
	}
	
	@Override
	public int hashCode() {
	    return pdbChainCode * 31 + transformId;
	}
	
	@Override
	public String toString() {
		return this.pdbChainCode+""+transformId;
	}

	@Override
	public int compareTo(SubunitId o) {
		if (this.transformId<o.transformId) {
			return -1;
		} 
		if (this.transformId>o.transformId) {
			return 1;
		}
		if (this.transformId==o.transformId) {
			if (this.pdbChainCode<o.pdbChainCode) {
				return -1;
			}
			if (this.pdbChainCode>o.pdbChainCode) {
				return 1;
			}
		}
		return 0;
	}
}

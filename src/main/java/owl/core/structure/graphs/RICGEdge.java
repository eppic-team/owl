package owl.core.structure.graphs;

import java.io.Serializable;

public class RICGEdge implements Serializable {

	private static final long serialVersionUID = 1L;

	private int nAtoms;
	private int nHBonds;
	private boolean isDisulfide;
	private boolean isClash;
	
	
	
	public int getnAtoms() {
		return nAtoms;
	}
	public void setnAtoms(int nAtoms) {
		this.nAtoms = nAtoms;
	}
	public int getnHBonds() {
		return nHBonds;
	}
	public void setnHBonds(int nHBonds) {
		this.nHBonds = nHBonds;
	}
	public boolean isDisulfide() {
		return isDisulfide;
	}
	public void setDisulfide(boolean isDisulfide) {
		this.isDisulfide = isDisulfide;
	}
	public boolean isClash() {
		return isClash;
	}
	public void setClash(boolean isClash) {
		this.isClash = isClash;
	}
	
	
}

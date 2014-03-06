package owl.core.structure.graphs;

import java.io.Serializable;

public class RICGEdge implements Serializable {

	private static final long serialVersionUID = 1L;

	private int nAtoms;
	private int nHBonds;
	private boolean isDisulfide;
	private boolean isClash;
	private double minDistance;
	
	public RICGEdge() {
		minDistance = Double.MAX_VALUE;
	}
	
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
	
	public double getMinDistance() {
		return minDistance;
	}
	
	public void setMinDistance(double minDistance) {
		this.minDistance = minDistance;
	}
	
	public void addDistance(double distance) {
		if (distance<minDistance) {
			minDistance = distance;
		}
	}
	
}

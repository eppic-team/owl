package proteinstructure;

import javax.vecmath.Point3d;

/**
 * A very basic class representing a residue as part of a Polymer 3-dimensional 
 * structure. At the moment only admits one atom per residue.
 * @author duarte
 *
 */
public class Residue {
	
	private String type;
	private int serial;
	private String chainCode;
	private String pdbChainCode;
	
	//TODO eventually there should be an atom class
	//TODO here we restrict a residue to one atom only (CA), it should allow multiple atoms
	private Point3d coords;
	
	
	public Residue(String type, int serial, String chainCode, String pdbChainCode, Point3d coords) {
		this.type = type;
		this.serial = serial;
		this.chainCode = chainCode;
		this.pdbChainCode = pdbChainCode;
		this.coords = coords;
	}

	public String getType() {
		return type;
	}



	public void setType(String type) {
		this.type = type;
	}



	public int getSerial() {
		return serial;
	}



	public void setSerial(int serial) {
		this.serial = serial;
	}



	public String getChainCode() {
		return chainCode;
	}



	public void setChainCode(String chainCode) {
		this.chainCode = chainCode;
	}



	public String getPdbChainCode() {
		return pdbChainCode;
	}



	public void setPdbChainCode(String pdbChainCode) {
		this.pdbChainCode = pdbChainCode;
	}



	public Point3d getCoords() {
		return coords;
	}



	public void setCoords(Point3d coords) {
		this.coords = coords;
	}

	public String toString() {
		return chainCode+serial;
	}

	
}

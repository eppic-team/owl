package owl.core.structure;

import java.util.HashMap;

import javax.vecmath.Point3d;

/**
 * A very basic class representing a residue as part of a Polymer 3-dimensional 
 * structure. At the moment only admits one atom per residue.
 * @author duarte
 *
 */
public class PolResidue {

	// masses taken from http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
	// masses measured in unified atomic mass units (u)
	private static double C_MASS = 12.0;
	private static double P_MASS = 30.974;	
	private static HashMap<String, Double> masses = initialiseMassesMap();
	
	
	private String type;
	private int serial;
	private String chainCode;
	private String pdbChainCode;
	
	//TODO eventually there should be an atom class
	//TODO here we restrict a residue to one atom only (CA), it should allow multiple atoms
	private String atomType;
	private Point3d coords;
	
	
	public PolResidue(String type, int serial, String chainCode, String pdbChainCode, String atomType, Point3d coords) {
		this.type = type;
		this.serial = serial;
		this.chainCode = chainCode;
		this.pdbChainCode = pdbChainCode;
		this.atomType = atomType;
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



	public String getAtomType() {
		return atomType;
	}

	public void setAtomType(String atomType) {
		this.atomType = atomType;
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

	public double getMass() {
		return getMassForAtom(this.atomType.substring(0, 1));
	}

	/**
	 * Deep copies this Residue
	 * @return
	 */
	public PolResidue copy() {
		Point3d coords = this.getCoords();
		PolResidue residue = new PolResidue(this.getType(),this.getSerial(),this.getChainCode(),this.getPdbChainCode(),this.getAtomType(),new Point3d(coords.x,coords.y,coords.z));
		return residue;
	}
	
	private static HashMap<String,Double> initialiseMassesMap() {
		HashMap<String,Double> masses = new HashMap<String, Double>();
		masses.put("C", C_MASS);
		masses.put("P", P_MASS);
		return masses;
	}
	
	public static double getMassForAtom(String atom) {
		return masses.get(atom);
	}
	
}

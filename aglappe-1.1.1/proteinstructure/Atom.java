package proteinstructure;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Point3d;

/**
 * A class representing an atom within a residue of a PDB protein structure
 * @author duarte
 *
 */
public class Atom {

	private static final double DEFAULT_B_FACTOR  = 0.00;		// default value if no b-factor is given
	private static final double DEFAULT_OCCUPANCY = 1.00;		// default value if no occupancy is given
	
	
	private String type;
	private String code;
	private int serial;
	private Point3d coords;
	
	private Residue parentResidue; 
	
	private double bfactor;
	private double occupancy;

	public Atom(int serial, String code, Point3d coords, Residue parentResidue) {
		this.serial = serial;
		this.code = code;
		this.coords = coords;
		this.parentResidue = parentResidue;
		this.bfactor = DEFAULT_B_FACTOR;
		this.occupancy = DEFAULT_OCCUPANCY;
		// getting atom type from code. This works for C, H, S, N and O atoms, not necessarily for others
		Pattern p = Pattern.compile("^\\d?(\\w)\\w*$");
		Matcher m = p.matcher(code);
		if (m.matches()) {
			this.type = m.group(1);
		} else {
			System.err.println("Warning! The atom code "+code+" with serial "+serial+" is not of the standard form! Can't assign an atom type for it!");
		}
	}
	
	public String getType() {
		return type;
	}
	
	public void setType(String type) {
		this.type = type;
	}
	
	public String getCode() {
		return code;
	}

	public void setCode(String code) {
		this.code = code;
	}

	public int getSerial() {
		return serial;
	}
	
	public void setSerial(int serial) {
		this.serial = serial;
	}
	
	public Residue getParentResidue() {
		return parentResidue;
	}
	
	public int getParentResSerial() {
		return this.parentResidue.getSerial();
	}

	public void setParentResidue(Residue parentResidue) {
		this.parentResidue = parentResidue;
	}

	public double getBfactor() {
		return bfactor;
	}

	public void setBfactor(double bfactor) {
		this.bfactor = bfactor;
	}

	public double getOccupancy() {
		return occupancy;
	}

	public void setOccupancy(double occupancy) {
		this.occupancy = occupancy;
	}

	public Point3d getCoords() {
		return coords;
	}
	
	public void setCoords(Point3d coords) {
		this.coords = coords;
	}
	
	public String toString() {
		return serial+code;
	}
	
}

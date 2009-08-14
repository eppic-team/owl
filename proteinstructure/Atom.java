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

	public static final double DEFAULT_B_FACTOR  = 0.00;		// default value if no b-factor is given
	public static final double DEFAULT_OCCUPANCY = 1.00;		// default value if no occupancy is given
	
	
	private AtomType type;
	private String code;
	private int serial;
	private Point3d coords;
	
	private Residue parentResidue; 
	
	private double bfactor;
	private double occupancy;

	/**
	 * Constructs a new Atom given a serial, an atom code (standard PDB name, e.g. CA), 
	 * its coordinates and the parentResidue.
	 * @param serial
	 * @param code
	 * @param coords
	 * @param parentResidue
	 */
	public Atom(int serial, String code, Point3d coords, Residue parentResidue) {
		this.serial = serial;
		this.code = code;
		this.coords = coords;
		this.parentResidue = parentResidue;
		this.bfactor = DEFAULT_B_FACTOR;
		this.occupancy = DEFAULT_OCCUPANCY;
		// getting atom type from code. This works for C, H, S, N and O atoms (all that is needed for standard aas), not necessarily for others
		Pattern p = Pattern.compile("^(\\w)\\w*$");
		Matcher m = p.matcher(code);
		if (m.matches()) {
			this.type = AtomType.getBySymbol(m.group(1));
			if (type==null) {
				System.err.println("Warning! Can not recognise atom type "+m.group(1)+" parsed from atom code "+code+" serial "+serial);
			}
		} else {
			System.err.println("Warning! The atom code "+code+" with serial "+serial+" is not of the standard form! Can't assign an atom type for it!");
		}
	}
	
	/**
	 * Returns the chemical type of this Atom 
	 * @return
	 */
	public AtomType getType() {
		return type;
	}
	
	/**
	 * Sets this Atom's chemical type
	 * @param type the AtomType to set, e.g. AtomType.C
	 */
	public void setType(AtomType type) {
		this.type = type;
	}
	
	/**
	 * Returns the atom code, i.e. the PDB standard atom name (e.g. CA, N, C, O, CB, ...)
	 * @return
	 */
	public String getCode() {
		return code;
	}

	/**
	 * Sets the atom code, i.e. the PDB standard atom name (e.g. CA, N, C, O, CB, ...)
	 * @param code
	 */
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

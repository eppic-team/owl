package owl.core.structure;

import java.io.Serializable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Point3d;

/**
 * A class representing an atom within a residue of a PDB protein structure
 * @author duarte
 *
 */
public class Atom implements Serializable {

	private static final long serialVersionUID = 1L;

	public static final double DEFAULT_B_FACTOR  = 0.00;		// default value if no b-factor is given
	public static final double DEFAULT_OCCUPANCY = 1.00;		// default value if no occupancy is given
		
	private static final Pattern ATOM_TYPE_PATTERN = Pattern.compile("^(\\w)\\w*$");
	
	
	private AtomType type;
	private String code;
	private int serial;
	private Point3d coords;
	
	private Residue parentResidue; 
	
	private double bfactor;
	private double occupancy;
	
	private double asa;
	private double bsa;
	private double vdwradius;

	/**
	 * Constructs a new Atom given a serial, an atom code (standard PDB name, e.g. CA), 
	 * its coordinates and the parentResidue.
	 * @param serial
	 * @param code
	 * @param element the element as it appears in last column of PDB file (1 or 2 letters both capitals)
	 * @param coords
	 * @param parentResidue
	 * @param occupancy
	 * @param bfactor
	 */
	public Atom(int serial, String code, String element, Point3d coords, Residue parentResidue, double occupancy, double bfactor) {
		this.serial = serial;
		this.code = code;
		this.coords = coords;
		this.parentResidue = parentResidue;
		this.bfactor = bfactor;
		this.occupancy = occupancy;
		if (element!=null) {
			this.type = AtomType.getBySymbol(element);
			if (type==null) {
				System.err.println("Warning! Can not recognise atom element "+element+" with atom code "+code+" and serial "+serial);
				this.type = AtomType.X;
			}
		} else { 
			// If we couldn't parse an atom type (e.g. not present in PDB file) we still try to get the element from the atom code
			// Getting atom type from code: it works for C, H, S, N and O atoms (all that is needed for standard aas), not necessarily for others
			Matcher m = ATOM_TYPE_PATTERN.matcher(code);
			if (m.matches()) {
				String elementGuess = null;
				if (Character.isDigit(m.group(1).charAt(0))) {
					// if first is a digit we go for second as the element guess
					if (code.length()>1) {
						elementGuess = Character.toString(code.charAt(1));
						this.type = AtomType.getBySymbol(elementGuess);
					}
				} else {
					// if first is not a digit, that's our element guess
					elementGuess = m.group(1);
					this.type = AtomType.getBySymbol(elementGuess);
				}
				// if not of the above is a good guess then we have no type and we warn
				if (type==null) { 
					System.err.println("Warning! Can't guess chemical element for atom code "+code+". Element guess was "+elementGuess+". Does the code follow the PDB standard?");
					this.type = AtomType.X;
				}									
			} else {
				System.err.println("Warning! Can't guess chemical element for atom code "+code+". Does the code follow the PDB standard?");
				this.type = AtomType.X;
			}
		}
	}
	
	/**
	 * Constructs an empty atom. 
	 */
	private Atom() {
		
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

	public void setParentResidue(AaResidue parentResidue) {
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

	public double getAsa() {
		return asa;
	}
	
	public void setAsa(double asa) {
		this.asa = asa;
	}
	
	public double getBsa() {
		return this.bsa;
	}
	
	public void setBsa(double bsa) {
		this.bsa = bsa;
	}
	
	public double getRadius() {
		return vdwradius;
	}
	
	public void setRadius(double vdwradius) {
		this.vdwradius = vdwradius;
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
	
	/**
	 * Deep copies this Atom object
	 * @return
	 */
	public Atom copy(Residue parentResidue) {
		Atom newAtom = new Atom();
		newAtom.serial = this.serial;
		newAtom.code = this.code;
		newAtom.type = this.type;
		newAtom.coords = new Point3d(this.coords);
		newAtom.parentResidue = parentResidue;
		newAtom.bfactor = this.bfactor;
		newAtom.occupancy = this.occupancy;
		return newAtom;
	}
	
	/**
	 * Equality based on same atom serial and same atom code, same type 
	 * of parent residue and same serial of parent residue.
	 */
	public boolean equals(Object other) {
		if (!(other instanceof Atom)) return false;
		Atom o = (Atom) other;
		if (this.serial!=o.serial) { 
			return false;
		}
		if (!this.code.equals(o.code)) {
			return false;
		}
		if (!this.parentResidue.getLongCode().equals(o.parentResidue.getLongCode())) {
			return false;
		}
		if (this.parentResidue.getSerial()!=o.parentResidue.getSerial()) {
			return false;
		}
		
		return true;
	}
	
	/**
	 * Hash code based on serial, atom code, parent residue type and parent residue serial.
	 */
	public int hashCode() {
		int hash = serial;
	    hash = hash * 31 + code.hashCode();
	    hash = hash * 31 + parentResidue.getLongCode().hashCode();
	    hash = hash * 31 + parentResidue.getSerial();
	    return hash;
	}
}

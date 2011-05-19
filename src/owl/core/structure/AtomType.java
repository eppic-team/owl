package owl.core.structure;

import java.util.HashMap;

/**
 * An atom 
 * 
 * Atomic masses from:
 * http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
 * 
 * vdW radii of standard aa elements from:
 * http://en.wikipedia.org/wiki/Van_der_Waals_radius
 * and for other elements: 
 * http://www.ccdc.cam.ac.uk/products/csd/radii/table.php4
 * 
 * See {@link AtomRadii} class for more accurate vdw radii values
 * 
 * @author duarte
 *
 */

public enum AtomType {

	// in standard aas
	H ( 1,  1.008, 1.20, "H", "Hydrogen",   true),
	C ( 6, 12.011, 1.70, "C", "Carbon",     true),
	O ( 8, 15.999, 1.52, "O", "Oxygen",     true),
	N ( 7, 14.007, 1.55, "N", "Nitrogen",   true),
	S (16, 32.065, 1.80, "S", "Sulfur",     true),
	P (15, 30.974, 1.80, "P", "Phosphorus", true),
	// in non-standard aas and hets
	Li( 3,  6.941, 1.82, "LI", "Lithium",   false),
	B ( 5, 10.811, 2.00,  "B", "Boron",     false),
	F ( 9, 18.998, 1.47,  "F", "Fluorine",  false),
	Na(11, 22.990, 2.27, "NA", "Sodium",    false),
	Mg(12, 24.305, 1.73, "MG", "Magnesium", false),
	Cl(17, 35.453, 1.75, "CL", "Chlorine",  false),
	K (19, 39.098, 2.75,  "K", "Potassium", false),
	Ca(20, 40.078, 2.00, "CA", "Calcium",   false),
	V (23, 50.942, 2.00,  "V", "Vanadium",  false),
	Mn(25, 54.938, 2.00, "MN", "Manganese", false),
	Fe(26, 55.845, 2.00, "FE", "Iron",      false),
	Co(27, 58.933, 2.00, "CO", "Cobalt",    false),
	Ni(28, 58.693, 1.63, "NI", "Nickel",    false),
	Cu(29, 63.546, 1.40, "CU", "Copper",    false),
	Zn(30, 65.382, 1.39, "ZN", "Zinc",      false),
	As(33, 74.922, 1.85, "AS", "Arsenic",   false),
	Se(34, 78.963, 1.90, "SE", "Selenium",  false),
	Br(35, 79.904, 1.85, "BR", "Bromine",   false),
	Sr(38, 87.621, 2.00, "SR", "Strontium", false),
	Y (39, 88.906, 2.00,  "Y", "Yttrium",   false),
	Mo(42, 95.962, 2.00, "MO", "Molybdenum",false),
	Ru(44,101.072, 2.00, "RU", "Ruthenium", false),
	Cd(48,112.412, 1.58, "CD", "Cadmium",   false),
	I (53,126.904, 1.98,  "I", "Iodine",    false),
	Xe(54,131.294, 2.16, "XE", "Xenon",     false),
	Cs(55,132.905, 2.00, "CS", "Caesium",   false),
	Sm(62,150.362, 2.00, "SM", "Samarium",  false),
	Ho(67,164.930, 2.00, "HO", "Holmium",   false),
	W (74,183.841, 2.00,  "W", "Tungsten",  false),
	Pt(78,195.085, 1.72, "PT", "Platinum",  false),
	Au(79,196.967, 1.66, "AU", "Gold",      false),
	Hg(80,200.592, 1.55, "HG", "Mercury",   false),
	U (92,238.029, 1.86,  "U", "Uranium",   false),
	// unknown atom (we treat it as a nitrogen in terms of mass and radius)
	X ( 0, 14.007, 1.55,  "X", "Unknown",   false);
	
	private int atomicNumber;
	private double atomicMass;
	private double radius; // see AtomRadii class for more accurate vdw radii values
	private String symbol;
	private String name;
	private boolean isInStandardAA;
	
	private static final HashMap<String,AtomType> symbol2AtomType = initSymbol2AtomType();

	private AtomType(int atomicNumber, double atomicMass, double radius, String symbol, String name, boolean isInStandardAA) {
		this.atomicMass = atomicMass;
		this.atomicNumber = atomicNumber;
		this.radius = radius;
		this.name = name;
		this.symbol = symbol;
		this.isInStandardAA = isInStandardAA;
	}

	public int getAtomicNumber() {
		return atomicNumber;
	}

	public double getAtomicMass() {
		return atomicMass;
	}

	public double getRadius() {
		return radius;
	}
	
	public String getSymbol() {
		return symbol;
	}

	public String getName() {
		return name;
	}

	public boolean isInStandardAA() {
		return isInStandardAA;
	}
	
	public static AtomType getBySymbol(String symbol) {
		if (symbol2AtomType.containsKey(symbol)) {
			return symbol2AtomType.get(symbol);
		} else {
			return null;
		}
		
	}
	
	private static HashMap<String,AtomType> initSymbol2AtomType() {
		HashMap<String,AtomType> map = new HashMap<String, AtomType>();
		for (AtomType type:AtomType.values()) {
			map.put(type.getSymbol(), type);
		}
		return map;
	}
}

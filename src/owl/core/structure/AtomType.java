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
	H ( 1,  1.008, 1.20, "H", "Hydrogen",   false, true),
	C ( 6, 12.011, 1.70, "C", "Carbon",     false, true),
	N ( 7, 14.007, 1.55, "N", "Nitrogen",   true,  true),
	O ( 8, 15.999, 1.52, "O", "Oxygen",     true,  true),
	P (15, 30.974, 1.80, "P", "Phosphorus", false, true),
	S (16, 32.065, 1.80, "S", "Sulfur",     false, true),
	// in non-standard aas and hets
	D ( 1,  2.014, 1.20,  "D", "Deuterium", false, false), // couldn't find a value for vdw radius, using Hydrogen's
	Li( 3,  6.941, 1.82, "LI", "Lithium",   false, false),
	Be( 4,  9.012, 2.00, "BE", "Beryllium", false, false),
	B ( 5, 10.811, 2.00,  "B", "Boron",     false, false),
	F ( 9, 18.998, 1.47,  "F", "Fluorine",  true,  false),
	Ne(10, 20.178, 1.54, "NE", "Neon",      false, false),
	Na(11, 22.990, 2.27, "NA", "Sodium",    false, false),
	Mg(12, 24.305, 1.73, "MG", "Magnesium", false, false),
	Al(13, 26.982, 2.00, "AL", "Aluminium", false, false),
	Si(14, 28.086, 2.10, "SI", "Silicon",   false, false),
	Cl(17, 35.453, 1.75, "CL", "Chlorine",  false, false),
	K (19, 39.098, 2.75,  "K", "Potassium", false, false),
	Ca(20, 40.078, 2.00, "CA", "Calcium",   false, false),
	V (23, 50.942, 2.00,  "V", "Vanadium",  false, false),
	Mn(25, 54.938, 2.00, "MN", "Manganese", false, false),
	Fe(26, 55.845, 2.00, "FE", "Iron",      false, false),
	Co(27, 58.933, 2.00, "CO", "Cobalt",    false, false),
	Ni(28, 58.693, 1.63, "NI", "Nickel",    false, false),
	Cu(29, 63.546, 1.40, "CU", "Copper",    false, false),
	Zn(30, 65.382, 1.39, "ZN", "Zinc",      false, false),
	Ga(31, 69.723, 1.87, "Ga", "Gallium",   false, false),
	As(33, 74.922, 1.85, "AS", "Arsenic",   false, false),
	Se(34, 78.963, 1.90, "SE", "Selenium",  false, false),
	Br(35, 79.904, 1.85, "BR", "Bromine",   false, false),
	Rb(37, 85.468, 2.00, "RB", "Rubidium",  false, false),
	Sr(38, 87.621, 2.00, "SR", "Strontium", false, false),
	Y (39, 88.906, 2.00,  "Y", "Yttrium",   false, false),
	Mo(42, 95.962, 2.00, "MO", "Molybdenum",false, false),
	Ru(44,101.072, 2.00, "RU", "Ruthenium", false, false),
	Ag(47,107.868, 1.72, "AG", "Silver",    false, false),
	Cd(48,112.412, 1.58, "CD", "Cadmium",   false, false),
	I (53,126.904, 1.98,  "I", "Iodine",    false, false),
	Xe(54,131.294, 2.16, "XE", "Xenon",     false, false),
	Cs(55,132.905, 2.00, "CS", "Caesium",   false, false),
	Ba(56,137.327, 2.00, "BA", "Barium",    false, false),
	La(57,138.905, 2.00, "LA", "Lanthanum", false, false),
	Pr(59,140.908, 2.00, "PR", "Praseodymium",false,false),
	Sm(62,150.362, 2.00, "SM", "Samarium",  false, false),
	Eu(63,151.964, 2.00, "EU", "Europium",  false, false),
	Gd(64,157.253, 2.00, "GD", "Gadolinium",false, false),
	Tb(65,158.925, 2.00, "TB", "Terbium",   false, false),
	Ho(67,164.930, 2.00, "HO", "Holmium",   false, false),
	Er(68,167.259, 2.00, "ER", "Erbium",    false, false),
	Yb(70,173.054, 2.00, "YB", "Ytterbium", false, false),
	Ta(73,180.948, 2.00, "TA", "Tantalum",  false, false),
	W (74,183.841, 2.00,  "W", "Tungsten",  false, false),
	Ir(77,192.217, 2.00, "IR", "Iridium",   false, false),
	Pt(78,195.085, 1.72, "PT", "Platinum",  false, false),
	Au(79,196.967, 1.66, "AU", "Gold",      false, false),
	Hg(80,200.592, 1.55, "HG", "Mercury",   false, false),
	Tl(81,204.383, 1.96, "TL", "Thallium",  false, false),
	Pb(82,207.200, 2.02, "PB", "Lead",      false, false),
	U (92,238.029, 1.86,  "U", "Uranium",   false, false),
	// unknown atom (we treat it as a nitrogen in terms of mass and radius)
	X ( 0, 14.007, 1.55,  "X", "Unknown",   false, false);
	
	private int atomicNumber;
	private double atomicMass;
	private double radius; // see AtomRadii class for more accurate vdw radii values
	private String symbol;
	private String name;
	private boolean isInStandardAA;
	private boolean isHbondAcceptor;
	
	private static final HashMap<String,AtomType> symbol2AtomType = initSymbol2AtomType();

	private AtomType(int atomicNumber, double atomicMass, double radius, String symbol, String name, boolean isHbondAcceptor, boolean isInStandardAA) {
		this.atomicMass = atomicMass;
		this.atomicNumber = atomicNumber;
		this.radius = radius;
		this.name = name;
		this.symbol = symbol;
		this.isHbondAcceptor = isHbondAcceptor;
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

	public boolean isHbondAcceptor() {
		return isHbondAcceptor;
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

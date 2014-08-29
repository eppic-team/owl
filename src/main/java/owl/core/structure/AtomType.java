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
 * The Cambridge Crystallography Data Center page above seems to be gone now (Aug 2013).
 * Anyway a better source should be this:
 * http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
 * which contains radii from Bondi et al 1964 and Mantina et al 2009 (same sources as webelements)
 * 
 * TODO radii in here are not yet corrected to be those in the wikipedia article above (based on Bondi and Mantina)
 * 
 * See {@link AtomRadii} class for more accurate vdw radii values
 * 
 * @author duarte
 *
 */

public enum AtomType {

	// in standard aas
	H ( 1,  1.008, 1.20, "H", "Hydrogen",	true),
	C ( 6, 12.011, 1.70, "C", "Carbon",     true),
	N ( 7, 14.007, 1.55, "N", "Nitrogen",   true),
	O ( 8, 15.999, 1.52, "O", "Oxygen",     true),
	P (15, 30.974, 1.80, "P", "Phosphorus", true),
	S (16, 32.065, 1.80, "S", "Sulfur",     true),
	// in non-standard aas and hets
	D ( 1,  2.014, 1.20,  "D", "Deuterium", false), // couldn't find a value for vdw radius, using Hydrogen's
	He( 2,  4.002, 1.40, "HE", "Helium",    false),
	Li( 3,  6.941, 1.82, "LI", "Lithium",   false),
	Be( 4,  9.012, 2.00, "BE", "Beryllium", false),
	B ( 5, 10.811, 2.00,  "B", "Boron",     false),
	F ( 9, 18.998, 1.47,  "F", "Fluorine",  false),
	Ne(10, 20.178, 1.54, "NE", "Neon",      false),
	Na(11, 22.990, 2.27, "NA", "Sodium",    false),
	Mg(12, 24.305, 1.73, "MG", "Magnesium", false),
	Al(13, 26.982, 2.00, "AL", "Aluminium", false),
	Si(14, 28.086, 2.10, "SI", "Silicon",   false),
	Cl(17, 35.453, 1.75, "CL", "Chlorine",  false),
	Ar(18, 39.948, 1.88, "AR", "Argon",     false),
	K (19, 39.098, 2.75,  "K", "Potassium", false),
	Ca(20, 40.078, 2.00, "CA", "Calcium",   false),
	Sc(21, 44.955, 2.00, "SC", "Scandium",  false),
	Ti(22, 47.867, 2.00, "TI", "Titanium",  false),
	V (23, 50.942, 2.00,  "V", "Vanadium",  false),
	Cr(24, 51.996, 2.00, "CR", "Chromium",  false),
	Mn(25, 54.938, 2.00, "MN", "Manganese", false),
	Fe(26, 55.845, 2.00, "FE", "Iron",      false),
	Co(27, 58.933, 2.00, "CO", "Cobalt",    false),
	Ni(28, 58.693, 1.63, "NI", "Nickel",    false),
	Cu(29, 63.546, 1.40, "CU", "Copper",    false),
	Zn(30, 65.382, 1.39, "ZN", "Zinc",      false),
	Ga(31, 69.723, 1.87, "GA", "Gallium",   false),
	Ge(32, 72.641, 2.00, "GE", "Germanium", false),
	As(33, 74.922, 1.85, "AS", "Arsenic",   false),
	Se(34, 78.963, 1.90, "SE", "Selenium",  false),
	Br(35, 79.904, 1.85, "BR", "Bromine",   false),
	Kr(36, 83.798, 2.02, "KR", "Krypton",   false),
	Rb(37, 85.468, 2.00, "RB", "Rubidium",  false),
	Sr(38, 87.621, 2.00, "SR", "Strontium", false),
	Y (39, 88.906, 2.00,  "Y", "Yttrium",   false),
	Zr(40, 91.224, 2.00, "ZR", "Zirconium", false),
	Nb(41, 92.906, 2.00, "NB", "Niobium",   false),
	Mo(42, 95.962, 2.00, "MO", "Molybdenum",false),
	Tc(43, 98.000, 2.00, "TC", "Technetium",false),
	Ru(44,101.072, 2.00, "RU", "Ruthenium", false),
	Rh(45,102.905, 2.00, "RH", "Rhodium",   false),
	Pd(46,106.421, 1.63, "PD", "Palladium", false),
	Ag(47,107.868, 1.72, "AG", "Silver",    false),
	Cd(48,112.412, 1.58, "CD", "Cadmium",   false),
	In(49,114.818, 1.93, "IN", "Indium",    false),
	Sn(50,118.711, 2.17, "SN", "Tin",       false),
	Sb(51,121.760, 2.00, "SB", "Antimony",  false),
	Te(52,127.603, 2.06, "TE", "Tellurium", false),
	I (53,126.904, 1.98,  "I", "Iodine",    false),
	Xe(54,131.294, 2.16, "XE", "Xenon",     false),
	Cs(55,132.905, 2.00, "CS", "Caesium",   false),
	Ba(56,137.327, 2.00, "BA", "Barium",    false),
	La(57,138.905, 2.00, "LA", "Lanthanum", false),
	Ce(58,140.116, 2.00, "CE", "Cerium",    false),
	Pr(59,140.908, 2.00, "PR", "Praseodymium",false),
	Nd(60,144.242, 2.00, "ND", "Neodymium", false),
	Pm(61,145.000, 2.00, "PM", "Promethium",false),
	Sm(62,150.362, 2.00, "SM", "Samarium",  false),
	Eu(63,151.964, 2.00, "EU", "Europium",  false),
	Gd(64,157.253, 2.00, "GD", "Gadolinium",false),
	Tb(65,158.925, 2.00, "TB", "Terbium",   false),
	Dy(66,162.500, 2.00, "DY", "Dysprosium",false),
	Ho(67,164.930, 2.00, "HO", "Holmium",   false),
	Er(68,167.259, 2.00, "ER", "Erbium",    false),
	Tm(69,168.934, 2.00, "TM", "Thulium",   false),
	Yb(70,173.054, 2.00, "YB", "Ytterbium", false),
	Lu(71,174.967, 2.00, "LU", "Lutetium",  false),
	Hf(72,178.492, 2.00, "HF", "Hafnium",   false),
	Ta(73,180.948, 2.00, "TA", "Tantalum",  false),
	W (74,183.841, 2.00,  "W", "Tungsten",  false),
	Re(75,186.207, 2.00, "RE", "Rhenium",   false),
	Os(76,190.230, 2.00, "OS", "Osmium",    false),
	Ir(77,192.217, 2.00, "IR", "Iridium",   false),
	Pt(78,195.085, 1.72, "PT", "Platinum",  false),
	Au(79,196.967, 1.66, "AU", "Gold",      false),
	Hg(80,200.592, 1.55, "HG", "Mercury",   false),
	Tl(81,204.383, 1.96, "TL", "Thallium",  false),
	Pb(82,207.200, 2.02, "PB", "Lead",      false),
	Bi(83,208.980, 2.00, "BI", "Bismuth",   false),
	Po(84,209.000, 2.00, "PO", "Polonium",  false),
	At(85,210.000, 2.00, "AT", "Astatine",  false),
	Rn(86,222.000, 2.00, "RN", "Radon",     false),
	Fr(87,223.000, 2.00, "FR", "Francium",  false),
	Ra(88,226.000, 2.00, "RA", "Radium",    false),
	Ac(89,227.000, 2.00, "AC", "Actinium",  false),
	Th(90,232.038, 2.00, "TH", "Thorium",   false),
	Pa(91,231.036, 2.00, "PA", "Protactinium",false),
	U (92,238.029, 1.86,  "U", "Uranium",   false),
	Np(93,237.000, 2.00, "NP", "Neptunium",   false),
	Pu(94,244.000, 2.00, "PU", "Plutonium",   false),
	Am(95,243.000, 2.00, "AM", "Americium",   false),
	Cm(96,247.000, 2.00, "CM", "Curium",      false),
	Bk(97,247.000, 2.00, "BK", "Berkelium",   false),
	Cf(98,251.000, 2.00, "CF", "Californium", false),
	Es(99,252.000, 2.00, "ES", "Einsteinium", false),
	Fm(100,257.000, 2.00, "FM", "Fermium",    false),
	Md(101,258.000, 2.00, "MD", "Mendelevium",false),
	No(102,259.000, 2.00, "NO", "Nobelium",   false),
	Lr(103,262.000, 2.00, "LR", "Lawrencium", false),
	Rf(104,265.000, 2.00, "RF", "Rutherfordium",false),
	Db(105,268.000, 2.00, "DB", "Dubnium",    false),
	Sg(106,271.000, 2.00, "SG", "Seaborgium", false),
	Bh(107,271.000, 2.00, "BH", "Bohrium",    false),
	Hs(108,270.000, 2.00, "HS", "Hassium",    false),
	Mt(109,276.000, 2.00, "MT", "Meitnerium", false),
	Ds(110,281.000, 2.00, "DS", "Darmstadtium", false),
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

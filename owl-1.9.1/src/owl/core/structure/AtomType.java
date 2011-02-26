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
	H( 1,  1.008, 1.20, "H", "Hydrogen",   true),
	C( 6, 12.000, 1.70, "C", "Carbon",     true),
	O( 8, 15.995, 1.52, "O", "Oxygen",     true),
	N( 7, 14.003, 1.55, "N", "Nitrogen",   true),
	S(16, 31.972, 1.80, "S", "Sulfur",     true),
	P(15, 30.974, 1.80, "P", "Phosphorus", true),
	// in non-standard aas and hets
	Na(11, 22.990, 2.27, "NA", "Sodium", false),
	Mg(12, 23.985, 1.73, "MG", "Magnesium",false),
	Cl(17, 34.969, 1.75, "CL", "Chlorine",false),
	Ca(20, 39.962, 2.00, "CA", "Calcium",false),
	Mn(25, 54.938, 2.00, "MN", "Manganese",false),
	Fe(26, 53.940, 2.00, "FE", "Iron",false),
	Co(27, 58.933, 2.00, "CO", "Cobalt",false),
	Cu(29, 62.930, 1.40, "CU", "Copper",false),
	Zn(30, 63.929, 1.39, "ZN", "Zinc",false),
	Se(34, 73.922, 1.90, "SE", "Selenium",false),
	Br(35, 78.918, 1.85, "BR", "Bromine",false),
	Cd(48,105.906, 1.58, "CD", "Cadmium",false),
	I (53,126.904, 1.98,  "I", "Iodine",false);
	
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

	public boolean iInStandardAA() {
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

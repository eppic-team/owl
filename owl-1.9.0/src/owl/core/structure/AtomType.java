package owl.core.structure;

import java.util.HashMap;

/**
 * An atom 
 * 
 * Atomic masses from:
 * http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=html&isotype=some
 * @author duarte
 *
 */

public enum AtomType {

	H( 1,  1.008, 1.20, "H", "Hydrogen"),
	C( 6, 12.000, 1.70, "C", "Carbon"),
	O( 8, 15.995, 1.52, "O", "Oxygen"),
	N( 7, 14.003, 1.55, "N", "Nitrogen"),
	S(16, 31.972, 1.80, "S", "Sulfur"),
	P(15, 30.974, 1.80, "P", "Phosphorus");
	
	
	private int atomicNumber;
	private double atomicMass;
	private double radius;
	private String symbol;
	private String name;
	
	private static final HashMap<String,AtomType> symbol2AtomType = initSymbol2AtomType();

	private AtomType(int atomicNumber, double atomicMass, double radius, String symbol, String name) {
		this.atomicMass = atomicMass;
		this.atomicNumber = atomicNumber;
		this.radius = radius;
		this.name = name;
		this.symbol = symbol;
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

	public static AtomType getBySymbol(String symbol) {
		return symbol2AtomType.get(symbol);
		
	}
	
	private static HashMap<String,AtomType> initSymbol2AtomType() {
		HashMap<String,AtomType> map = new HashMap<String, AtomType>();
		for (AtomType type:AtomType.values()) {
			map.put(type.getSymbol(), type);
		}
		return map;
	}
}

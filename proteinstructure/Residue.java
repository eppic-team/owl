package proteinstructure;

import java.util.TreeMap;

/**
 * A class representing a residue as part of a protein structure 
 * @author duarte
 *
 */
public class Residue {

	private String type;
	private int serial;
	private String chainCode;

	private double consurfScore;
	private int consurfColor;
	private double rsa;
	private double scrsa;
	
	private TreeMap<String, Atom> atoms;
	
	public Residue(String type, int serial, String chainCode) {
		atoms = new TreeMap<String, Atom>();
		this.type = type;
		this.serial = serial;
		this.chainCode = chainCode;
	}

	public void addAtom(Atom atom) {
		atoms.put(atom.getType(), atom);
	}
	
	public Atom getAtom(String atomType) {
		return atoms.get(atomType);
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

	public double getConsurfScore() {
		return consurfScore;
	}

	public void setConsurfScore(double consurfScore) {
		this.consurfScore = consurfScore;
	}

	public int getConsurfColor() {
		return consurfColor;
	}

	public void setConsurfColor(int consurfColor) {
		this.consurfColor = consurfColor;
	}

	public double getRsa() {
		return rsa;
	}

	public void setRsa(double rsa) {
		this.rsa = rsa;
	}

	public double getScrsa() {
		return scrsa;
	}

	public void setScrsa(double scrsa) {
		this.scrsa = scrsa;
	}

	public String toString() {
		return chainCode+serial;
	}
	
}

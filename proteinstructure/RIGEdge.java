package proteinstructure;

/**
 * Class representing a Residue Interaction Graph edge
 *
 */
public class RIGEdge {

	private static double DEFAULT_WEIGHT = 1.0;
	private static int DEFAULT_ATOM_WEIGHT = 1;
	
	private double weight;
	private int atomWeight;  
	private double distance;
	
	public RIGEdge(double weight, double distance, int atomWeight) {
		this.weight = weight;
		this.distance = distance;
		this.atomWeight = atomWeight;
	}
	
	public RIGEdge(int atomWeight){
		this.atomWeight = atomWeight;
		this.distance = 0;
		this.weight = DEFAULT_WEIGHT;
	}
	
	public RIGEdge(double weight) {
		this.weight=weight;
		this.distance = 0;
		this.atomWeight = DEFAULT_ATOM_WEIGHT;
	}
	
	public RIGEdge() {
		this.weight=DEFAULT_WEIGHT;
		this.distance = 0;
		this.atomWeight = DEFAULT_ATOM_WEIGHT;
	}
	
	public int getAtomWeight() {
		return atomWeight;
	}
	
	public double getWeight() {
		return weight;
	}
	
	public double getDistance() {
		return distance;
	}
	
	public String toString() {
		return "E";
	}
	
}

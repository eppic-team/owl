package owl.core.structure;

/**
 * Class representing a Residue Interaction Graph edge.
 * The RIGEdge stores properties of the edge (like weights) but not the 
 * end points. Those can be taken from the RIGraph with getEndPoints(edge)
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
	
	public void setAtomWeight(int atomWeight) {
		this.atomWeight = atomWeight;
	}

	public double getWeight() {
		return weight;
	}
	
	public void setWeight(double weight) {
		this.weight = weight;
	}
	
	public double getDistance() {
		return distance;
	}
	
	public void setDistance(double distance) {
		this.distance = distance;
	}

	public String toString() {
		return "E";
	}
	
	public RIGEdge copy() {
		return new RIGEdge(this.weight,this.distance,this.atomWeight);
	}
	
}

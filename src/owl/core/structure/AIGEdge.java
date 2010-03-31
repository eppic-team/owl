package owl.core.structure;

public class AIGEdge {

	private static final double DEFAULT_WEIGHT = 1; 
	private double weight;
	private double distance;
	
	public AIGEdge(double distance) {
		this.weight= DEFAULT_WEIGHT;
		this.distance = distance;
	}
	
	public AIGEdge(double weight, double distance) {
		this.weight = weight;
		this.distance = distance;
	}
	
	public double getWeight() {
		return this.weight;
	}
	
	public double getDistance() {
		return this.distance;
	}

	public AIGEdge copy() {
		return (new AIGEdge(weight, distance));
	}
}

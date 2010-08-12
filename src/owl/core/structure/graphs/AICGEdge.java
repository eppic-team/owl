package owl.core.structure.graphs;

public class AICGEdge {

	private double distance;
	
	public AICGEdge(double distance) {
		this.distance = distance;
	}

	public double getDistance() {
		return this.distance;
	}

	public AICGEdge copy() {
		return (new AICGEdge(distance));
	}

}

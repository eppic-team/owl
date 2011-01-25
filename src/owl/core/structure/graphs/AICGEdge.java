package owl.core.structure.graphs;

import java.io.Serializable;

public class AICGEdge implements Serializable {

	private static final long serialVersionUID = 1L;

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

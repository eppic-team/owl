package embed;

/**
 * Class encapsulating a distance range with lower and upper bounds
 * @author duarte
 *
 */
public class Bound {

	public double lower;
	public double upper;
	
	public Bound(double lower, double upper) {
		this.lower = lower;
		this.upper = upper;
	}
	
	public String toString() {
		return String.format("[%4.1f %4.1f]", lower, upper);
	}
	
}




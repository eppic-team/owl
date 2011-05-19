package owl.core.util;

import javax.vecmath.Vector3d;


/**
 * A class to represent a moment of inertia together with its associated axis.
 * Implements Comparable, being the comparison based on the values of the moment of 
 * inertia.
 * @author duarte
 *
 */
public class InertiaMomentAndAxis implements Comparable<InertiaMomentAndAxis> {
	public double moment;
	public Vector3d axis;
	public InertiaMomentAndAxis(double moment, Vector3d axis) {
		this.moment = moment;
		this.axis = axis;
	}
	public int compareTo(InertiaMomentAndAxis o) {
		if (this.moment<o.moment) return -1;
		if (this.moment>o.moment) return 1;
		return 0;
	}
}

package vecmath;


/**
 * Our replacement for javax.vecmath.Vector3d
 */
public class Vector3d extends Tuple3d {
	
	public Vector3d(double x, double y, double z) {
		super(x,y,z);
	}
	
	public Vector3d() {
		super();
	}
	
	public double lengthSquared() {
		return (x*x + y*y + z*z);
	}
	
}

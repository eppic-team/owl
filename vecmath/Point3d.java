package vecmath;


/**
 * Our replacement for javax.vecmath.Point3d.
 * @author stehr
 *
 */
public class Point3d extends Tuple3d {

	public Point3d(double x, double y, double z) {
		super(x,y,z);
	}
	
	public Point3d() {
		super();
	}
	
	public double distance(Point3d p) {
		return Math.sqrt((x-p.x)*(x-p.x)+(y-p.y)*(y-p.y)+(z-p.z)*(z-p.z));
	}
}

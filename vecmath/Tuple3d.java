package vecmath;

/**
 * Our replacement for javax.vecmath.Tuple3d
 */
public class Tuple3d {
	
	public double x,y,z;
	
	public Tuple3d(double x, double y, double z) {
		this.x=x;
		this.y=y;
		this.z=z;
	}
	
	public Tuple3d() {
		this.x = 0.0;
		this.y = 0.0;
		this.z = 0.0;
	}

	public boolean equals(Object o) {
		try {
			Tuple3d p = (Tuple3d) o;
			return (x==p.x && y==p.y && z==p.z);
		} catch (NullPointerException e) {
			return false;
		} catch (ClassCastException e) {
			return false;
		}
	}
	
	public void add(Tuple3d p) {
		this.x += p.x;
		this.y += p.y;
		this.z += p.z;
	}
	
	public void sub(Tuple3d p) {
		this.x -= p.x;
		this.y -= p.y;
		this.z -= p.z;
	}
	
	public void scale(double s) {
		this.x *= s;
		this.y *= s;
		this.z *= s;
	}
	
	/**
	 * Copies the tuple to the given array.
	 */
	public void get(double[] a) {
		a[0] = x;
		a[1] = y;
		a[2] = z;
	}

}

package vecmath;

public class Point3i {
	
	public int x,y,z;
	
	public Point3i(int x, int y, int z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public boolean equals(Object o) {
		try {
			Point3i p = (Point3i) o;
			return (x==p.x && y==p.y && z==p.z);
		}
		catch (NullPointerException e) {
			return false;
		}
		catch (ClassCastException e) {
			return false;
		}
	}

	public int hashCode() {
		long bits = 1L;
		bits = 31L * bits + (long)x;
		bits = 31L * bits + (long)y;
		bits = 31L * bits + (long)z;
		return (int) (bits ^ (bits >> 32));
	}
	
}

package owl.core.structure;

//import javax.vecmath.Vector3d;

public class InterfaceGraphEdge implements Comparable<InterfaceGraphEdge> {

	private ChainInterface interf;
	private int id;
	
	public InterfaceGraphEdge(ChainInterface interf) {
		this.interf = interf;
		this.id = interf.getId();
	}
	
	@Override
	public int compareTo(InterfaceGraphEdge arg0) {
		return new Integer(id).compareTo(arg0.id);
	}
	
	/**
	 * Compares the two InterfaceGraphEdge based in their ids only.
	 * @param o
	 */
	@Override
	public boolean equals(Object o) {
		if (!(o instanceof InterfaceGraphEdge)) return false;
		InterfaceGraphEdge other = (InterfaceGraphEdge) o;
		return (other.id==this.id);
	}

	@Override
	public int hashCode() {
		return id;
	}
	
	@Override
	public String toString() {
		String str = id+"";
		return str;
	}
	
	public int getId() {
		return id;
	}
	
	public ChainInterface getInterface() {
		return interf;
	}
	
//	public Vector3d getConnectionVector() {
//		return interf.getConnectionVector();
//	}
	
}

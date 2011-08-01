package owl.core.structure;

public class InterfaceGraphEdge implements Comparable<InterfaceGraphEdge> {

	private int id;
	private boolean withinUnitCell;
	
	public InterfaceGraphEdge(int id, boolean withinUnitCell) {
		this.id = id;
		this.withinUnitCell = withinUnitCell;
	}
	
	@Override
	public int compareTo(InterfaceGraphEdge arg0) {
		return new Integer(id).compareTo(arg0.id);
	}
	
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
		if (!withinUnitCell) str+="'";
		return str;
	}
	
	public int getId() {
		return id;
	}
	
	public boolean isWithinUnitCell() {
		return withinUnitCell;
	}
}

package owl.core.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class InterfaceCluster implements Serializable, Comparable<InterfaceCluster> {

	private static final long serialVersionUID = 1L;

	private int id;
	
	private ChainInterface representative;
	
	private List<ChainInterface> members;
	
	public InterfaceCluster(int clusterId) {
		this.id = clusterId;
		this.members = new ArrayList<ChainInterface>();
	}
	
	public boolean addMember(ChainInterface interf) {
		if (members.contains(interf)) {
			return false;
		}
		
		this.members.add(interf);
		
		return true;
	}
	
	public int getId() {
		return id;
	}
	
	public void setId(int clusterId) {
		this.id = clusterId;
	}
	
	public ChainInterface getRepresentative() {
		return representative;
	}
	
	public List<ChainInterface> getMembers() {
		return members;
	}
	
	public int size() {
		return members.size();
	}
	
	public double getMeanArea() {
		double area=0.0;
		for (ChainInterface member:members){
			area += member.getInterfaceArea();
		}
		return area/((double)size());
	}
	
	/**
	 * Sorts the member interfaces by area descending and
	 * assigns the largest area member as the representative 
	 * @throws NullPointerException if members are not initialised or empty
	 */
	public void initialiseRepresentative() {
		// first we sort by area descending (see ChainInterface.comparteTo)
		Collections.sort(members);
		// and assign the largest area member as representative
		this.representative = members.get(0);
	}
	
	public String toString() {
		String str = ""+getId();

		str += " (";
		for (int i=0;i<members.size();i++) {
			str+= members.get(i).getId();
			if (i!=members.size()-1)
				str+=",";
			else str+=")";
		}
		
		return str;
	}

	@Override
	public int compareTo(InterfaceCluster o) {
		// this will sort descending on mean interface areas
		return Double.compare(o.getMeanArea(), this.getMeanArea()); 		
	}
}

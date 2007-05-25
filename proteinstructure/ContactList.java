package proteinstructure;

import java.util.ArrayList;

public class ContactList extends ArrayList<Contact> {

	private static final long serialVersionUID = 1L;

	public ContactList() {
		super();
	}
	
	/**
	 * Gets the maximum range of this ContactList
	 * i.e. the sequence separation for the pair with maximum sequence separation
	 * @return
	 */
	public int getMaxRange() {
		int max=0;
		for (Contact cont:this){
			Math.max(max, cont.getRange());
		}
		return max;
	}

	/**
	 * Gets the maximum node serial in this ContactList
	 * @return
	 */
	public int getMaxNode(){
		int max=0;
		for (Contact cont:this){
			int contactMax=Math.max(cont.i, cont.j);
			max = Math.max(max,contactMax);
		}
		return max;
	}
}

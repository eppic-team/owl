package owl.core.structure;

import java.util.Set;
import java.util.TreeMap;


/**
 * Class representing a contact type, i.e. the atoms that belong to each residue for a given contact type
 * Data hold in a TreeMap structure with keys the three letter residue types and values sets of atom names
 *
 */
public class ContactType extends TreeMap<String,Set<String>> {

	private static final long serialVersionUID = 1L;
	
	private String name;
	private boolean multiAtom; // true when multiAtom type, false when singleAtom type
	
	public ContactType(String name, boolean multiAtom) {
		super();
		this.name = name;
		this.multiAtom = multiAtom;
	}
	
	public boolean isMultiAtom() {
		return multiAtom;
	}
	
	public String getName() {
		return name;
	}
	
}

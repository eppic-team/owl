package owl.core.structure;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * A list of all the interfaces of a crystal structure (a PdbAsymUnit)
 * 
 * 
 * @author duarte_j
 *
 */
public class ChainInterfaceList implements Iterable<ChainInterface>{

	private List<ChainInterface> list;
	
	private PdbAsymUnit pdb; //the PDB entry that this list of interfaces refers to
	
	public ChainInterfaceList() {
		list = new ArrayList<ChainInterface>();
	}
	
	public void addInterface(ChainInterface interf) {
		list.add(interf);
	}
	
	public int getNumProtProtInterfaces() {
		int count = 0;
		for (ChainInterface interf:this) {
			if (interf.isProtein()) {
				count++;
			}
		}
		return count;
	}
	
	public void sort() {
		Collections.sort(list);
		int i=1;
		for (ChainInterface interf:list) {
			interf.setId(i);
			i++;
		}
	}
	
	public int size() {
		return list.size();
	}
	
	public ChainInterface get(int i) {
		return list.get(i);
	}
	
	public PdbAsymUnit getPdb() {
		return pdb;
	}
	
	public void setPdb(PdbAsymUnit pdb) {
		this.pdb = pdb;
	}
	
	@Override
	public Iterator<ChainInterface> iterator() {
		return list.iterator();
	}
}

package owl.core.structure;

import java.util.ArrayList;
import java.util.Iterator;

public class PdbxPolySeqLineList implements Iterable<PdbxPolySeqLine> {

	private ArrayList<PdbxPolySeqLine> list;

	public PdbxPolySeqLineList(){
		list = new ArrayList<PdbxPolySeqLine>();
	}
	
	public void add(PdbxPolySeqLine line) {
		list.add(line);
	}
	
	public boolean isEmpty() {
		return list.isEmpty();
	}
	
	@Override
	public Iterator<PdbxPolySeqLine> iterator() {
		return list.iterator();
	}
}

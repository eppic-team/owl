package owl.core.connections.pisa;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class PisaAsmSetList implements Iterable<PisaAsmSet> {

	private String pdbCode;
	private List<PisaAsmSet> list;
	private String status;
	
	public PisaAsmSetList() {
		list = new ArrayList<PisaAsmSet>();
	}
	
	public PisaAsmSet get(int i) {
		return list.get(i);
	}
	
	public boolean add(PisaAsmSet pisaAsmSet) {
		return list.add(pisaAsmSet);
	}
	
	public String getPdbCode() {
		return pdbCode;
	}
	
	public void setPdbCode(String pdbCode) {
		this.pdbCode = pdbCode;
	}

	public String getStatus() {
		return status;
	}
	
	public void setStatus(String status) {
		this.status = status;
	}

	@Override
	public Iterator<PisaAsmSet> iterator() {
		return list.iterator();
	}
}

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
	
	public int size() {
		return list.size();
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
	
	public OligomericPrediction getOligomericPred() {
		
		// no AsmSet at all, PISA didn't find any possible assemblies (stable or not) and so we call monomer 
		if (size()==0) {
			return new OligomericPrediction(1);
		} 
		PisaAsmSet pas = get(0);
		
		// NOTE pisa has a grey zone (anything with deltaGdiss between -2 and 2 I believe) that we don't use here.
		// we simply call either a valid assembly when deltaGdiss>0, not valid deltaGdiss<0
		
		// we check now for cases when the AsmSet is composed entirely of assemblies with deltaGdiss<0, e.g. 3hzl
		// those will be called monomers directly
		boolean validAsmExist = false;
		for (PisaAssembly pa:pas) {
			if (pa.getDissEnergy()>0 && pa.isMacromolecular()) validAsmExist = true;
		}
		if (!validAsmExist) 
			return new OligomericPrediction(1);

		// all other cases: one or more stable assemblies in the AsmSet
		OligomericPrediction op = new OligomericPrediction(1);
		boolean sizeSet = false;
		for (PisaAssembly pa:pas) {
			if (!pa.isMacromolecular()) continue;

			if (!sizeSet && pa.getDissEnergy()>0) {
				op.setMmSize(pa.getMmsize());
				sizeSet = true;
			}
			op.addAssembly(pa);

			//System.out.printf("\t%2d\t%2d\t%5.1f\t%20s\t%s\n",
			//		mmsizePred,pa.getMmsize(),
			//		pa.getDissEnergy(),
			//		pa.getFormula(),
			//		pa.getInterfaceIdsString());
		}
		return op;
	}
	
}

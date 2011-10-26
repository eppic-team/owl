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
		
		// NOTE that in principle PISA uses deltaGdiss to classify between stable, unstable and gray:
		//  stable >0
		//  unstable ~<-2
		//  gray ~>-2 && <0
		// But the rule is not always so clear, see for instance 1ibr where first assembly is deltaGdiss>0 but gray
		// or 1eer where deltaGdiss=-0.9 and it's considered stable, while for 1bam deltaGdiss=-1.6 is gray
		
		// check whether all assemblies of the set are unstable or all gray
		boolean allGray = true;
		boolean allUnstable = true;
		
		for (PisaAssembly pa:pas) {
			if (pa.isMacromolecular()) {
				if (pa.getPredictionType()!=PisaAssembly.PredictionType.GRAY) allGray = false;
				if (pa.getPredictionType()!=PisaAssembly.PredictionType.UNSTABLE) allUnstable = false;
			}
		}

		if (allGray) {
			return new OligomericPrediction(-1);
		}
		if (allUnstable) {
			// then it's a monomer
			return new OligomericPrediction(1);				
		}

			

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

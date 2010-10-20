package owl.core.structure;

import java.io.PrintStream;
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
	
	public enum AsaCalcMethod {
		INTERNAL("internal"),NACCESS("naccess"),PISA("pisa");
		
		private String name;
		private AsaCalcMethod(String name) {
			this.name = name;
		}
		
		public String getName() {
			return name;
		}
	}

	private List<ChainInterface> list;
	
	private PdbAsymUnit pdb; //the PDB entry that this list of interfaces refers to
	
	private AsaCalcMethod asaCalcMethod;
	private int asaCalcAccuracyParam;
	
	public ChainInterfaceList(AsaCalcMethod asaCalcMethod) {
		this.list = new ArrayList<ChainInterface>();
		this.asaCalcMethod = asaCalcMethod;
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
	
	public AsaCalcMethod getAsaCalcType() {
		return asaCalcMethod;
	}
	
	public void setAsaCalcAccuracyParam(int asaCalcAccuracyParam) {
		this.asaCalcAccuracyParam = asaCalcAccuracyParam;
	}
	
	public int getAsaCalcAccuracyParam() {
		return asaCalcAccuracyParam;
	}
	
	/**
	 * Calculates the rims and cores of all member interfaces for each of the
	 * bsaToAsa cutoffs given.
	 * @param bsaToAsaCutoffs
	 */
	public void calcRimAndCores(double[] bsaToAsaCutoffs) {
		for (ChainInterface interf:list){
			interf.calcRimAndCore(bsaToAsaCutoffs);
		}
	}
	
	@Override
	public Iterator<ChainInterface> iterator() {
		return list.iterator();
	}
	
	public void printTabular(PrintStream ps, String pdbName) {
		ps.println("Interfaces for "+pdbName);
		ps.print("ASAs values from "+this.asaCalcMethod.getName());
		if (asaCalcMethod==AsaCalcMethod.INTERNAL) {
			ps.println(" (accuracy="+this.asaCalcAccuracyParam+")");
		} else {
			ps.println();
		}
		for (ChainInterface interf:this) {
			interf.printTabular(ps);
		}
	}
}

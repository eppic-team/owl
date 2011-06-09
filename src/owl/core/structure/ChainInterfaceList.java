package owl.core.structure;

import java.io.PrintStream;
import java.io.Serializable;
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
public class ChainInterfaceList implements Iterable<ChainInterface>, Serializable {
	
	private static final long serialVersionUID = 1L;

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
	
	private AsaCalcMethod asaCalcMethod;
	private int asaCalcAccuracyParam;
	
	public ChainInterfaceList(AsaCalcMethod asaCalcMethod) {
		this.list = new ArrayList<ChainInterface>();
		this.asaCalcMethod = asaCalcMethod;
	}
	
	public void addInterface(ChainInterface interf) {
		list.add(interf);
	}
	
	public int getNumInterfaces() {
		return list.size();
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
	
	/**
	 * Sorts the interface list descending on total interface areas assigning identifiers 
	 * from 1 to n in that order.
	 */
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
	
	/**
	 * Gets the interface corresponding to given index i.
	 * Indices go from 0 to n-1.
	 * If {@link #sort()} was called then the order is descendent by area.
	 * Note that the interface id is actually i+1
	 * @param i
	 * @return
	 */
	public ChainInterface get(int i) {
		return list.get(i);
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
	
	public int getNumInterfacesAboveArea(double area) {
		int count = 0;
		for (ChainInterface interf:this) {
			if (interf.getInterfaceArea()>area) {
				count++;
			}
		}
		return count;
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
	
	/**
	 * Calculates the rims and cores of all member interfaces with a zooming approach.
	 * @param bsaToAsaSoftCutoff
	 * @param bsaToAsaHardCutoff
	 * @param relaxationStep
	 * @param minNumResidues
	 */
	public void calcRimAndCores(double bsaToAsaSoftCutoff, double bsaToAsaHardCutoff, double relaxationStep, int minNumResidues) {
		for (ChainInterface interf:list){
			interf.calcRimAndCore(bsaToAsaSoftCutoff, bsaToAsaHardCutoff, relaxationStep, minNumResidues);
		}		
	}
	
	public boolean hasInterfacesWithClashes(double clashDistance) {
		for (ChainInterface interf:list){
			if (interf.hasClashes(clashDistance)) return true;
		}
		return false;
	}
	
	public List<ChainInterface> getInterfacesWithClashes(double clashDistance) {
		List<ChainInterface> clashyInterfs = new ArrayList<ChainInterface>();
		for (ChainInterface interf:list){
			if (interf.hasClashes(clashDistance)) {
				clashyInterfs.add(interf);
			}
		}
		return clashyInterfs;
	}
	
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

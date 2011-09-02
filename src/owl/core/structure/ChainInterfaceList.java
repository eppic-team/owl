package owl.core.structure;

import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

//import owl.core.util.CombinationsGenerator;


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
	
	private InterfaceGraph graph; 
	
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
	 * Gets the interface corresponding to given id.
	 * The ids go from 1 to n
	 * If {@link #sort()} was called then the order is descendent by area.
	 * @param id
	 * @return
	 */
	public ChainInterface get(int id) {
		return list.get(id-1);
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
	 * Calculates the rims and cores of all interfaces in list for the
	 * bsaToAsaCutoff given.
	 * @param bsaToAsaCutoff
	 */
	public void calcRimAndCores(double bsaToAsaCutoff) {
		for (ChainInterface interf:list){
			interf.calcRimAndCore(bsaToAsaCutoff);
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
	
	public boolean hasInterfacesWithClashes() {
		for (ChainInterface interf:list){
			if (interf.hasClashes()) return true;
		}
		return false;
	}
	
	public List<ChainInterface> getInterfacesWithClashes() {
		List<ChainInterface> clashyInterfs = new ArrayList<ChainInterface>();
		for (ChainInterface interf:list){
			if (interf.hasClashes()) {
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
	
	/**
	 * Given a PDB chain code returns the List of residues that are in the surface but belong to NO
	 * interface (above given minInterfArea) 
	 * @param pdbChainCode
	 * @param minInterfArea
	 * @return
	 */
	public List<Residue> getResiduesNotInInterfaces(String pdbChainCode, double minInterfArea) {

		List<Residue> surfResidues = new ArrayList<Residue>();
		PdbChain chain = null;
		for (ChainInterface interf:this) {
			if (interf.getFirstMolecule().getPdbChainCode().equals(pdbChainCode)) {
				chain = interf.getFirstMolecule();
				break;
			}
		}
		
		for (Residue res:chain) {
			if (res.getAsa()>0) surfResidues.add(res);
		}
		//System.out.println(" in surface: "+surfResidues.size());
		
		Set<Integer> interfResSerials = new HashSet<Integer>();
		
		for (ChainInterface interf:this) {
			if (interf.getInterfaceArea()>minInterfArea) {
				if (interf.getFirstMolecule().getPdbChainCode().equals(pdbChainCode)) {
					for (Residue res:interf.getFirstRimCore().getCoreResidues()) {
						interfResSerials.add(res.getSerial());
					}
					for (Residue res:interf.getFirstRimCore().getRimResidues()) {
						interfResSerials.add(res.getSerial());
					}
				}
				if (interf.getSecondMolecule().getPdbChainCode().equals(pdbChainCode)) {
					for (Residue res:interf.getSecondRimCore().getCoreResidues()) {
						interfResSerials.add(res.getSerial());
					}
					for (Residue res:interf.getSecondRimCore().getRimResidues()) {
						interfResSerials.add(res.getSerial());
					}
				}
				
			}
		}

		Iterator<Residue> it = surfResidues.iterator();
		while (it.hasNext()) {
			Residue res = it.next();
			if (interfResSerials.contains(res.getSerial())) {
				it.remove();
			}
		}
		//System.out.println(" in no interface: "+surfResidues.size());
		
		return surfResidues;
	}
	
	/**
	 * Returns true if this interface list contains parallel interfaces
	 * A parallel interface is an interface between two equivalent subunits in different unit cells 
	 * related by a translation. Creating an assembly with it would result in an infinite assembly. 
	 * @return
	 */
	public boolean hasParallelInterfaces() {
		for (ChainInterface interf:this) {
			if (interf.isParallel()) {
				return true;
			}
		}
		return false;
	}
	
	/**
     * Get a list of ids of all interfaces that are parallel.
	 * A parallel interface is an interface between two equivalent subunits in different unit cells 
	 * related by a translation. Creating an assembly with it would result in an infinite assembly.
	 * @return
	 */
	public ArrayList<Integer> getParallelInterfacesIds() {
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (ChainInterface interf:this) {
			if (interf.isParallel()) {
				list.add(interf.getId());
			}
		}		
		return list;
	}
	
	/**
	 * Get a list of ids of all interfaces that are not parallel.
	 * A parallel interface is an interface between two equivalent subunits in different unit cells 
	 * related by a translation. Creating an assembly with it would result in an infinite assembly.
	 * @return
	 */
	public ArrayList<Integer> getNonParallelInterfacesIds() {
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (ChainInterface interf:this) {
			if (!interf.isParallel()) {
				list.add(interf.getId());
			}
		}		
		return list;
	}
	
	public InterfaceGraph getInterfacesGraph() {
		if (graph==null) {
			graph = new InterfaceGraph(this);
		}
		return graph;
	}
	
//	public List<Assembly> getAllAssemblies() {
//		List<Assembly> assemblies = new ArrayList<Assembly>();
//		
//		getInterfacesGraph();
//		
//		ArrayList<Integer> candidates = getNonParallelInterfacesIds();
//		System.out.println("Total "+candidates.size()+" interfaces to use ("+(size()-candidates.size())+" parallel interfaces discarded)");
//		
//		// we then enumerate all assemblies with 1 interface, 2 interfaces, 3 interfaces .... up to n
//		System.out.println("Theoretical total assemblies: "+(((int)Math.pow(2, candidates.size()))-1));
//		for (int n=1;n<=candidates.size();n++) {
//			List<Assembly> sizenassemblies = enumerateAssembliesSizeN(n,candidates);
//			
//			assemblies.addAll(sizenassemblies);
//		}
//		
//		return assemblies;
//	}
//	
//	private List<Assembly> enumerateAssembliesSizeN(int n, ArrayList<Integer> candidates) {
//		List<Assembly> assemblies = new ArrayList<Assembly>();
//		CombinationsGenerator x = new CombinationsGenerator (candidates.size(), n);
//		while (x.hasMore()) {
//			Assembly ass = new Assembly(this);
//			int[] indices = x.getNext();
//			for (int i = 0; i < indices.length; i++) {
//				ass.add(candidates.get(indices[i]));
//			}
//			assemblies.add(ass);
//		}
//		System.out.println("size "+n+": "+assemblies.size());
//		return assemblies;
//	}
}

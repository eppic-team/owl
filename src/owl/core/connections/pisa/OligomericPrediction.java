package owl.core.connections.pisa;

import java.util.ArrayList;
import java.util.List;

/**
 * Our interpretation of PISA assembly predictions: represents a single oligomeric prediction
 * extracted from the PISA assemblies output.
 * 
 * Note that each OligomericPrediction can be composed of one or more 
 * equivalent PisaAssemblies (because of equivalent NCS related components) 
 *  
 * @author duarte_j
 *
 */
public class OligomericPrediction {

	private int mmSize;
	private List<PisaAssembly> assemblies;
	
	public OligomericPrediction(int mmSize) {
		this.mmSize = mmSize;
		assemblies = new ArrayList<PisaAssembly>();
	}
	
	public int getMmSize() {
		return mmSize;
	}
	
	public void setMmSize(int mmSize) {
		this.mmSize = mmSize;
	}
	
	public void addAssembly(PisaAssembly assembly) {
		if (assembly.getMmsize()!=mmSize) {
			System.err.println("WARNING! Grouping a PISA assembly of size "+assembly.getMmsize()+" with one of size "+mmSize+". Pisa oligomeric prediction will be ambiguous!");
		}
		this.assemblies.add(assembly);
	}
	
	public List<PisaAssembly> getAssemblies() {
		return assemblies;
	}
	
	/**
	 * Returns the PISA interface ids of all engaged protein-protein interfaces of this 
	 * PISA oligomeric prediction.
	 * @return
	 */
	public List<Integer> getProtInterfacesIds(PisaInterfaceList pil) {
		List<Integer> ids = new ArrayList<Integer>();
		for (PisaAssembly ass:assemblies) {
			for (int id:ass.getInterfaceIds()) {
				if (pil.getById(id).isProtein()) ids.add(id);
			}
		}
		return ids;
	}
	
	/**
	 * Returns true if given PISA interface id is one of the interfaces engaged in the 
	 * assemblies represented by this prediction
	 * @param pil
	 * @param id
	 * @return
	 */
	public boolean containProtInterface(PisaInterfaceList pil, int id) {
		for (PisaAssembly ass:assemblies) {
			for (int assid:ass.getInterfaceIds()) {
				if (assid==id) return true;
			}
		}
		return false;
		
	}
	
}

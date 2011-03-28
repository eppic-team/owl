package owl.core.connections.pisa;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import owl.core.structure.ChainInterfaceList;
import owl.core.structure.ChainInterfaceList.AsaCalcMethod;
import owl.core.structure.PdbAsymUnit;

public class PisaInterfaceList implements Iterable<PisaInterface> {

	private String pdbCode;
	private List<PisaInterface> list;
	
	public PisaInterfaceList() {
		list = new ArrayList<PisaInterface>();
	}
	
	public boolean add(PisaInterface interf) {
		return list.add(interf);
	}

	@Override
	public Iterator<PisaInterface> iterator() {
		return list.iterator();
	}
	
	public String getPdbCode() {
		return pdbCode; 
	}
	
	public void setPdbCode(String pdbCode) {
		this.pdbCode = pdbCode;
	}
	
	/**
	 * Converts this list of Pisa interfaces to one of our ChainInterfaceList containing interfaces
	 * including coordinates appropriately transformed as indicated by transformations read from Pisa.
	 * @param pdb the PDB entry data corresponding to this PisaInterfaceList
	 * @return
	 * @throws IllegalArgumentException if the pdb code of given pdb does not match this PisaInterface's code
	 */
	public ChainInterfaceList convertToChainInterfaceList(PdbAsymUnit pdb) {
		if (!pdb.getPdbCode().equals(this.pdbCode)) {
			throw new IllegalArgumentException("Pdb code of given PdbAsymUnit does not match pdb code from PISA");
		}
		ChainInterfaceList interfList = new ChainInterfaceList(AsaCalcMethod.PISA);
		for (PisaInterface interf:this){
			interfList.addInterface(interf.convertToChainInterface(pdb));
		}
		return interfList;
	}
}

package owl.core.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * A cluster of protein chains: a group of identical in sequence (except for 
 * unobserved residues) NCS-related protein chains (also called entities in PDB)
 * 
 * @author duarte_j
 *
 */
public class ChainCluster implements Serializable {

	private static final long serialVersionUID = 1L;

	/**
	 * The chain representing the cluster: by convention the one with the
	 * lowest (alphabetical order) PDB chain code
	 */
	private PdbChain representative;
	
	/**
	 * The members of the chain cluster, including the representative
	 */
	private List<PdbChain> members;
	
	public ChainCluster(PdbAsymUnit pdb, List<String> pdbChainCodes) {
		this.members = new ArrayList<PdbChain>();
		
		for (String pdbChainCode:pdbChainCodes) {
			this.members.add(pdb.getChain(pdbChainCode));
		}
		
		Collections.sort(pdbChainCodes);
		this.representative = pdb.getChain(pdbChainCodes.get(0));
	}
	
	/**
	 * Returns the chain representing the cluster: by convention the one with the
	 * lowest (alphabetical order) PDB chain code
	 */
	public PdbChain getRepresentative() {
		return representative;
	}
	
	/**
	 * Returns a List with all members of the chain cluster, including the representative
	 */
	public List<PdbChain> getMembers() {
		return members;
	}
	
	/**
	 * Returns a string containing the representative PDB chain code followed by the
	 * PDB chain codes of all member chains of this cluster in brackets. 
	 * The string is formatted like: A (B,C,D,E)  
	 * @param repPdbChainCode
	 * @return
	 */
	public String getClusterString() {
		
		String str = representative.getPdbChainCode();
		
		if (members.size()>1) str+=" (";
		
		for (int i=0;i<members.size();i++) {
			
			if (members.get(i)==representative) {
				continue;
			}
			
			str+= members.get(i).getPdbChainCode()+",";
 
		}
		
		if (members.size()>1)
			str = str.substring(0, str.length()-1)+")";
		
		return str;
	}
	
	/**
	 * Finds the common list of residue serials of all members for which the
	 * given atomName (of a standard aminoacid) is observed (has coordinates) in every member of this cluster.
	 * 
	 * @param atomName the atom name of a standard aminoacid residue, e.g. CA
	 * @return
	 */
	public List<Integer> getCommonObservedSet(String atomName) {
		List<Integer> common = new ArrayList<Integer>();

		// first we get the list of residue serials from first member
		for (Residue res:members.get(0)) {
			if (res instanceof AaResidue && res.containsAtom(atomName)) {
				common.add(res.getSerial()); 
			}
		}

		// then we check all others and take intersection
		if (members.size()>1) {
			for (int i=1;i<members.size();i++) {
				PdbChain member = members.get(i);
				Iterator<Integer> it = common.iterator();
				while (it.hasNext()) {
					int resser = it.next();
					if (!member.containsResidue(resser) || 
							!member.getResidue(resser).containsAtom(atomName)) {
						 
						it.remove();
					} 
				}
			}
		}
		
		return common;
	}
}

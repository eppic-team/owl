package owl.core.connections.pisa;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class PisaMolecule implements Iterable<PisaResidue> {

	private int id;
	private String chainId;
	private String molClass;
	
	private List<PisaResidue> residues;
	
	public PisaMolecule() {
		residues = new ArrayList<PisaResidue>();
	}
	
	public void addResidue(PisaResidue residue) {
		residues.add(residue);
	}

	public void printTabular(PrintStream ps) {
		ps.println(this.getId()+"\t"+this.getChainId()+"\t"+this.getMolClass());
		for (PisaResidue residue:this) {
			residue.printTabular(ps);
		}
	}
	
	/**
	 * @return the id
	 */
	public int getId() {
		return id;
	}

	/**
	 * @param id the id to set
	 */
	public void setId(int id) {
		this.id = id;
	}

	/**
	 * @return the chainId
	 */
	public String getChainId() {
		return chainId;
	}

	/**
	 * @param chainId the chainId to set
	 */
	public void setChainId(String chainId) {
		this.chainId = chainId;
	}

	/**
	 * @return the molClass
	 */
	public String getMolClass() {
		return molClass;
	}

	/**
	 * @param molClass the molClass to set
	 */
	public void setMolClass(String molClass) {
		this.molClass = molClass;
	}

	/**
	 * @return the residues
	 */
	public List<PisaResidue> getResidues() {
		return residues;
	}

	/**
	 * @param residues the residues to set
	 */
	public void setResidues(List<PisaResidue> residues) {
		this.residues = residues;
	}

	@Override
	public Iterator<PisaResidue> iterator() {
		return residues.iterator();
	}
}

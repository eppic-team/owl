package owl.core.connections.pisa;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.vecmath.Matrix4d;

import owl.core.structure.Pdb;
import owl.core.structure.Residue;

public class PisaMolecule implements Iterable<PisaResidue> {
	
	public static final String CLASS_PROTEIN = "Protein";

	private int id;
	private String chainId;
	private String molClass;
	
	private Matrix4d transf;
	private Matrix4d transfOrth;
		
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
	 * Returns true if this is a protein molecule, false otherwise
	 * @return
	 */
	public boolean isProtein() {
		return getMolClass().equals(CLASS_PROTEIN);	
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
	 * Gets the symmetry operator used to generate this molecule (in crystal coordinates).
	 * @return
	 */
	public Matrix4d getTransf() {
		return transf;
	}
	
	/**
	 * Sets the symmetry operator used to generate this molecule (in crystal coordinates).
	 * @param transf
	 */
	public void setTransf(Matrix4d transf) {
		this.transf = transf;
	}

	/**
	 * Gets the symmetry operator used to generate this molecule (in orthonormal coordinates).
	 * @return
	 */
	public Matrix4d getTransfOrth() {
		return transfOrth;
	}
	
	/**
	 * Sets the symmetry operator used to generate this molecule (in orthonormal coordinates).
	 * @param transfOrth
	 */
	public void setTransfOrth(Matrix4d transfOrth) {
		this.transfOrth = transfOrth;
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
	
	/**
	 * Sets the asa/bsa values of the passed pdb object to the ones present in this PisaMolecule's 
	 * PisaResidues
	 * @param pdb
	 */
	public void setAsaAndBsas(Pdb pdb) {
		for (PisaResidue pisaRes:this) {
			int resSerial = pdb.getResSerFromPdbResSer(pisaRes.getPdbResSer());
			if (pdb.containsResidue(resSerial)) {
				Residue res = pdb.getResidue(resSerial);
				res.setAsa(pisaRes.getAsa());
				res.setBsa(pisaRes.getBsa());
			}
		}
	}
	
}

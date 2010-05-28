package owl.core.connections.pisa;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class PisaMolecule implements Iterable<PisaResidue> {
	
	public static final String CLASS_PROTEIN = "Protein";

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
	 * Returns 2 lists of residues as a {@link PisaRimCore} object: core residues are those for 
	 * which the bsa/asa ratio is above the given cut-off, rim those with bsa>0 and with 
	 * bsa/asa below the cut-off.
	 * @param bsaToAsaCutoff
	 * @return
	 */
	public PisaRimCore getRimAndCore(double bsaToAsaCutoff) {
		List<PisaResidue> core = new ArrayList<PisaResidue>();
		List<PisaResidue> rim = new ArrayList<PisaResidue>();
		for (PisaResidue residue:this) {
			if (residue.getBsa()>0) {
				if (residue.getBsaToAsaRatio()<bsaToAsaCutoff) {
					rim.add(residue);
				} else {
					core.add(residue);
				}
			}
		}
		return new PisaRimCore(rim,core,bsaToAsaCutoff);
	}
	
	/**
	 * Returns 2 list of residues as a {@link PisaRimCore} object (see {@link #getRimAndCore(double)})
	 * The core is required to have a minimum of minNumResidues. If the minimum is not 
	 * reached with the bsaToAsaSoftCutoff, then the cutoff is relaxed in relaxationStep steps 
	 * until reaching the bsaToAsaHardCutoff.
	 * @param bsaToAsaSoftCutoff
	 * @param bsaToAsaHardCutoff
	 * @param relaxationStep
	 * @param minNumResidues
	 * @return
	 */
	public PisaRimCore getRimAndCore(double bsaToAsaSoftCutoff, double bsaToAsaHardCutoff, double relaxationStep, int minNumResidues) {
		PisaRimCore rimCore = null;
		for (double cutoff=bsaToAsaSoftCutoff;cutoff>=bsaToAsaHardCutoff;cutoff-=relaxationStep) {
			rimCore = getRimAndCore(cutoff);
			//System.out.printf("cutoff %4.2f, core size: %d\n",cutoff,rimCore.getCoreSize());
			if (rimCore.getCoreSize()>=minNumResidues) {
				break;
			}
		}
		return rimCore;
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

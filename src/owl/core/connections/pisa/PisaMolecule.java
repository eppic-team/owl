package owl.core.connections.pisa;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

public class PisaMolecule implements Iterable<PisaResidue> {
	
	public static final String CLASS_PROTEIN = "Protein";

	private int id;
	private String chainId;
	private String molClass;
	
	private double rxx;
	private double rxy;
	private double rxz;
	private double ryx;
	private double ryy;
	private double ryz;
	private double rzx;
	private double rzy;
	private double rzz;
	private double tx;
	private double ty;
	private double tz;
	
	private String transf; // the symop transformation in algebraic notation (we upper case it)
	
	
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
		// we introduce a margin of relaxationSte*0.10 to be sure we do go all the way down to bsaToAsaHardCutoff (necessary because of rounding)
		for (double cutoff=bsaToAsaSoftCutoff;cutoff>=bsaToAsaHardCutoff-relaxationStep*0.10;cutoff-=relaxationStep) {
			rimCore = getRimAndCore(cutoff);
			//System.out.printf("cutoff %4.2f, core size: %d\n",cutoff,rimCore.getCoreSize());
			if (rimCore.getCoreSize()>=minNumResidues) {
				break;
			}
		}
		return rimCore;
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
	 * Gets the symmetry operator used to generate this molecule. 
	 * @return
	 */
	public Matrix4d getSymOp() {
		return new Matrix4d(new Matrix3d(rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz),new Vector3d(tx,ty,tz),1.0);
	}
	
	/**
	 * Gets the symmetry operator (algebraic notation, upper-case) used to generate this molecule.
	 * @return
	 */
	public String getTransf() {
		return transf;
	}
	
	/**
	 * Sets the symmetry operator (algebraic notation, upper-case) used to generate this molecule.
	 * @param transf
	 */
	public void setTransf(String transf) {
		this.transf = transf;
	}
	
	/**
	 * @param rxx the rxx to set
	 */
	public void setRxx(double rxx) {
		this.rxx = rxx;
	}

	/**
	 * @param rxy the rxy to set
	 */
	public void setRxy(double rxy) {
		this.rxy = rxy;
	}

	/**
	 * @param rxz the rxz to set
	 */
	public void setRxz(double rxz) {
		this.rxz = rxz;
	}

	/**
	 * @param ryx the ryx to set
	 */
	public void setRyx(double ryx) {
		this.ryx = ryx;
	}

	/**
	 * @param ryy the ryy to set
	 */
	public void setRyy(double ryy) {
		this.ryy = ryy;
	}

	/**
	 * @param ryz the ryz to set
	 */
	public void setRyz(double ryz) {
		this.ryz = ryz;
	}

	/**
	 * @param rzx the rzx to set
	 */
	public void setRzx(double rzx) {
		this.rzx = rzx;
	}

	/**
	 * @param rzy the rzy to set
	 */
	public void setRzy(double rzy) {
		this.rzy = rzy;
	}

	/**
	 * @param rzz the rzz to set
	 */
	public void setRzz(double rzz) {
		this.rzz = rzz;
	}

	/**
	 * @param tx the tx to set
	 */
	public void setTx(double tx) {
		this.tx = tx;
	}

	/**
	 * @param ty the ty to set
	 */
	public void setTy(double ty) {
		this.ty = ty;
	}

	/**
	 * @param tz the tz to set
	 */
	public void setTz(double tz) {
		this.tz = tz;
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

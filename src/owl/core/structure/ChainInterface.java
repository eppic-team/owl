package owl.core.structure;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

import javax.vecmath.Matrix4d;

import owl.core.runners.NaccessRunner;
import owl.core.structure.graphs.AICGraph;

public class ChainInterface implements Comparable<ChainInterface> {
	
	public static final String TYPE_PROTEIN = "Protein";

	private int id;
	private double interfaceArea;
	private AICGraph graph;
	
	private Pdb firstMolecule;
	private Pdb secondMolecule;
	
	private Matrix4d firstTransf; 		// the transformation applied to first molecule expressed in crystal axes coordinates
	private Matrix4d firstTransfOrth;	// the transformation applied to first molecule expressed in orthonormal axes coordinates
	private Matrix4d secondTransf; 		// the transformation applied to second molecule expressed in crystal axes coordinates
	private Matrix4d secondTransfOrth;	// the transformation applied to second molecule expressed in orthonormal axes coordinates
	
	private String firstMolType;	// interface descriptions taken from PISA can have ligands as members, 
									// our own interfaces can only be protein (TYPE_PROTEIN constant)
	private String secondMolType;
	
	private double score; 			// a score value assigned to the interface (if from PISA this is the solvation energy) 
	
	/**
	 * Constructs an empty ChainInterface. Use the setters to set the values.
	 */
	public ChainInterface() {
		
	}
	
	public ChainInterface(Pdb firstMolecule, Pdb secondMolecule, AICGraph graph, Matrix4d firstTransf, Matrix4d secondTransf) {
		this.firstMolecule = firstMolecule;
		this.secondMolecule = secondMolecule;
		this.graph = graph;
		this.firstTransf = firstTransf;
		this.secondTransf = secondTransf;
		this.firstMolType = TYPE_PROTEIN;
		this.secondMolType = TYPE_PROTEIN;
	}
	
	public int getId() {
		return id;
	}
	public void setId(int id) {
		this.id = id;
	}
	public double getInterfaceArea() {
		return interfaceArea;
	}
	public void setInterfaceArea(double interfaceArea) {
		this.interfaceArea = interfaceArea;
	}
	public Pdb getFirstMolecule() {
		return firstMolecule;
	}
	public void setFirstMolecule(Pdb firstMolecule) {
		this.firstMolecule = firstMolecule;
	}
	public Pdb getSecondMolecule() {
		return secondMolecule;
	}
	public void setSecondMolecule(Pdb secondMolecule) {
		this.secondMolecule = secondMolecule;
	}
	
	public double getCutoff() {
		return graph.getDistCutoff();
	}
	
	public AICGraph getAICGraph() {
		return graph;
	}
	
	public Matrix4d getFirstTransf() {
		return firstTransf;
	}
	
	public void setFirstTransf(Matrix4d firstTransf) {
		this.firstTransf = firstTransf;
	}
	
	/**
	 * Returns the transformation used to generate the first molecule in orthonormal axes coordinates.
	 * @return
	 */
	public Matrix4d getFirstTransfOrth(){
		if (firstTransfOrth==null) {
			firstTransfOrth = firstMolecule.getCrystalCell().transfToOrthonormal(firstTransf);
		}
		return firstTransfOrth;
	}
	
	public void setFirstTransfOrth(Matrix4d firstTransfOrth) {
		this.firstTransfOrth = firstTransfOrth;
	}
	
	public Matrix4d getSecondTransf() {
		return secondTransf;
	}
	
	public void setSecondTransf(Matrix4d secondTransf) {
		this.secondTransf = secondTransf;
	} 
	
	public Matrix4d getSecondTransfOrth() {
		if (secondTransfOrth==null) {
			secondTransfOrth = secondMolecule.getCrystalCell().transfToOrthonormal(secondTransf);
		}
		return secondTransfOrth;
	}
	
	public void setSecondTransfOrth(Matrix4d secondTransfOrth) {
		this.secondTransfOrth = secondTransfOrth;
	}
	
	public double getScore() {
		return score;
	}
	
	public void setScore(double score) {
		this.score = score;
	}
	
	public String getFirstMolType() {
		return firstMolType;
	}
	
	public void setFirstMolType(String firstMolType) {
		this.firstMolType = firstMolType;
	}
	
	public String getSecondMolType() {
		return secondMolType;
	}
	
	public void setSecondMolType(String secondMolType) {
		this.secondMolType = secondMolType;
	}
	
	/**
	 * Sets the absolute surface accessibility values of this interface's two members from the given map 
	 * of pdb chain codes to maps of residue serials to ASA values.
	 * @param asas
	 */
	private void setAbsSurfaceAccessibilities(HashMap<String, HashMap<Integer,Double>> asas) {
		this.getFirstMolecule().setAbsSurfaceAccessibilities(asas.get(getFirstMolecule().getPdbChainCode()));
		this.getSecondMolecule().setAbsSurfaceAccessibilities(asas.get(getSecondMolecule().getPdbChainCode()));
	}
	
	/**
	 * Runs the NACCESS program to calculate accessible surface areas of the complex of 
	 * the two molecules making up this interface. 
	 * A HashMap with the ASA values of the uncomplexed molecules must be given (pdb chain codes to 
	 * HashMaps of residue serials to ASA values).
	 * Then the BSAs are set from the uncomplexed and complexed ASA values and finally the total 
	 * interface area is calculated (use {@link #getInterfaceArea()} to get it)
	 * @param naccessExecutable
	 * @throws IOException
	 */
	public void calcBSAnaccess(File naccessExecutable, HashMap<String, HashMap<Integer,Double>> asas) throws IOException {
		this.setAbsSurfaceAccessibilities(asas);
		NaccessRunner nar = new NaccessRunner(naccessExecutable, "");
		
		PdbAsymUnit complex = new PdbAsymUnit(firstMolecule.getPdbCode(), 1, null, null, null);
		complex.setChain("A", firstMolecule);
		complex.setChain("B", secondMolecule);
 
		nar.runNaccess(complex); // this will set the bsa members of the Residues of firstMolecule and secondMolecule
	
		double totBuried = 0.0;
		for (Residue residue:firstMolecule.getResidues().values()) {
			totBuried+=residue.getBsa();
		}
		for (Residue residue:secondMolecule.getResidues().values()) {
			totBuried+=residue.getBsa();
		}

		this.interfaceArea = totBuried;
	}
	
	public int getNumAtomsInContact() {
		return this.graph.getVertexCount();
	}
	
	public int getNumContacts() {
		return this.graph.getEdgeCount();
	}
	
	public boolean isProtein() {
		return (isFirstProtein() && isSecondProtein());
	}
	
	public boolean isFirstProtein() {
		return this.firstMolType.equals(TYPE_PROTEIN);
	}
	
	public boolean isSecondProtein() {
		return this.secondMolType.equals(TYPE_PROTEIN);
	}
	
	@Override
	public int compareTo(ChainInterface o) {
		// this will sort descending on interface areas
		return (Double.compare(o.interfaceArea,this.interfaceArea));
	}
	
	public String toString() {
		return firstMolecule.getPdbChainCode()+"-"+secondMolecule.getPdbChainCode()+" "+
		getAICGraph().getEdgeCount()+" "+SpaceGroup.getAlgebraicFromMatrix(secondTransf);
	}
	
	public boolean equals(Object o) {
		if (!(o instanceof ChainInterface)) return false;
		ChainInterface other = (ChainInterface) o;
		String tcc1 = firstMolecule.getPdbChainCode();
		String tcc2 = secondMolecule.getPdbChainCode();
		String occ1 = other.firstMolecule.getPdbChainCode();
		String occ2 = other.secondMolecule.getPdbChainCode();
		
		if ( !(occ1.equals(tcc1) && occ2.equals(tcc2)) &&
			 !(occ2.equals(tcc1) && occ1.equals(tcc2)) ) {
			return false;
		}
		return this.graph.equals(other.graph);
	}
	
	public int hashCode() {
		int hash = graph.getEdgeCount();
	    hash = hash * 31 + firstMolecule.getPdbChainCode().hashCode()+secondMolecule.getPdbChainCode().hashCode();
	    return hash; 
	}
	
	public void printTabular(PrintStream ps) {
		ps.print("# ");
		ps.printf("%d\t%9.2f\t%5.2f\n",this.getId(),this.getInterfaceArea(),this.getScore());
		ps.print("## ");
		this.printFirstMolInfoTabular(ps);
		ps.print("## ");
		this.printSecondMolInfoTabular(ps);
	}
	
	private void printFirstMolInfoTabular(PrintStream ps) {
		ps.println("1\t"+firstMolecule.getPdbChainCode()+"\t"+this.getFirstMolType());
		for (Residue residue:firstMolecule.getResidues().values()) {
			residue.printTabular(ps);
		}
	}

	private void printSecondMolInfoTabular(PrintStream ps) {
		ps.println("2\t"+secondMolecule.getPdbChainCode()+"\t"+this.getFirstMolType());
		for (Residue residue:secondMolecule.getResidues().values()) {
			residue.printTabular(ps);
		}
	}

	/**
	 * Returns a map containing 2 {@link InterfaceRimCore} objects (see getRimAndCore in {@link PisaMolecule})
	 * for each of the 2 members of the interface.
	 * The sum of the residues of the 2 cores is required to be at least minNumResidues. 
	 * If the minimum is not reached with the bsaToAsaSoftCutoff, then the cutoff is 
	 * relaxed in relaxationStep steps until reaching the bsaToAsaHardCutoff.
	 * If either of the 2 molecules of this interface is not a protein, its rimCore 
	 * object in the output map object will be null. If both are not proteins then the map 
	 * will contain null object references for both.
	 * @param bsaToAsaSoftCutoff
	 * @param bsaToAsaHardCutoff
	 * @param relaxationStep
	 * @param minNumResidues
	 * @return
	 */
	public Map<Integer,InterfaceRimCore> getRimAndCore(double bsaToAsaSoftCutoff, double bsaToAsaHardCutoff, double relaxationStep, int minNumResidues) {
		
		Map<Integer,InterfaceRimCore> rimcores = new HashMap<Integer, InterfaceRimCore>();
		if (!isFirstProtein() && !isSecondProtein()) {
			rimcores.put(1, null);
			rimcores.put(2, null);
			return rimcores;
		}
		
		// we introduce a margin of relaxationSte*0.10 to be sure we do go all the way down to bsaToAsaHardCutoff (necessary because of rounding)
		for (double cutoff=bsaToAsaSoftCutoff;cutoff>=bsaToAsaHardCutoff-relaxationStep*0.10;cutoff-=relaxationStep) {
			InterfaceRimCore rimCore1 = null;
			InterfaceRimCore rimCore2 = null;
			if (isFirstProtein()) rimCore1 = this.firstMolecule.getRimAndCore(cutoff);
			if (isSecondProtein()) rimCore2 = this.secondMolecule.getRimAndCore(cutoff);
			rimcores.put(1,rimCore1);
			rimcores.put(2,rimCore2);
			
			int totalCoreResidues = 0;
			if (isFirstProtein()) totalCoreResidues+=rimCore1.getCoreSize();
			if (isSecondProtein()) totalCoreResidues+=rimCore2.getCoreSize();
			if (totalCoreResidues>=minNumResidues) {
				break;
			}
		}
		
		return rimcores;
	}

}

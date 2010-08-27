package owl.core.structure;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import owl.core.runners.NaccessRunner;
import owl.core.structure.graphs.AICGraph;

public class ChainInterface implements Comparable<ChainInterface> {

	private int id;
	private double interfaceArea;
	private AICGraph graph;
	
	private Pdb firstMolecule;
	private Pdb secondMolecule;
	
	private String firstTransf; // the transformation applied to first molecule expressed in algebraic notation 
	private String secondTransf; // the transformation applied to second molecule expressed in algebraic notation
	
	public ChainInterface(Pdb firstMolecule, Pdb secondMolecule, AICGraph graph, String firstTransf, String secondTransf) {
		this.firstMolecule = firstMolecule;
		this.secondMolecule = secondMolecule;
		this.graph = graph;
		this.firstTransf = firstTransf;
		this.secondTransf = secondTransf;
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
	
	public String getFirstTransf() {
		return firstTransf;
	}
	
	public String getSecondTransf() {
		return secondTransf;
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
	
		double firstTotSurface = 0.0;
		double secondTotSurface = 0.0;
		double complexTotSurface = 0.0;
		for (Residue residue:firstMolecule.getResidues().values()) {
			firstTotSurface+=residue.getAsa();
			complexTotSurface+=residue.getBsa();
		}
		for (Residue residue:secondMolecule.getResidues().values()) {
			secondTotSurface+=residue.getAsa();
			complexTotSurface+=residue.getBsa();
		}

		this.interfaceArea = firstTotSurface+secondTotSurface-complexTotSurface;
	}
	
	@Override
	public int compareTo(ChainInterface o) {
		return (Double.compare(this.interfaceArea,o.interfaceArea));
	}
}

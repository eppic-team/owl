package owl.core.structure;

import java.io.File;
import java.io.IOException;

import owl.core.runners.NaccessRunner;
import owl.core.structure.graphs.AICGraph;

public class ChainInterface implements Comparable<ChainInterface> {

	private int id;
	private double interfaceArea;
	private AICGraph graph;
	
	private Pdb firstMolecule;
	private Pdb secondMolecule;
	
	private String transf; // the transformation applied expressed in algebraic notation 
	
	public ChainInterface(Pdb firstMolecule, Pdb secondMolecule, AICGraph graph, String transf) {
		this.firstMolecule = firstMolecule;
		this.secondMolecule = secondMolecule;
		this.graph = graph;
		this.transf = transf;
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
	
	public String getTransf() {
		return transf;
	}
	
	/**
	 * Runs the NACCESS program to calculate accessible surface areas of both interface 
	 * partners and of the complex (the Buried Surface Area). Calculates also the total
	 * interface area, use {@link #getInterfaceArea()} to get it.
	 * @param naccessExecutable
	 * @throws IOException
	 */
	public void calcBSAnaccess(File naccessExecutable) throws IOException {
		NaccessRunner nar = new NaccessRunner(naccessExecutable, "");
		nar.runNaccess(firstMolecule);
		nar.runNaccess(secondMolecule);
		
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
	
	/**
	 * Equality based on same transformation string and same chain identifiers (independent of order)
	 */
	public boolean equals(Object other) {
		if (!(other instanceof ChainInterface)) {
			return false;
		}
		ChainInterface o = (ChainInterface) other;
		if (!this.getTransf().equals(o.getTransf())) {
			return false;
		}
		if (!(firstMolecule.getPdbChainCode().equals(o.firstMolecule.getPdbChainCode()) && 
			  secondMolecule.getPdbChainCode().equals(o.secondMolecule.getPdbChainCode())) ||
			!(firstMolecule.getPdbChainCode().equals(o.secondMolecule.getPdbChainCode()) &&
			  secondMolecule.getPdbChainCode().equals(o.firstMolecule.getPdbChainCode()))
		) {
			return false;
		}
		return true;
	}
	
	/**
	 * Hash code based on transformation string and chain identifiers (independent of order)
	 */
	public int hashCode() {
		int hash = 1;
		hash = hash*31 + transf.hashCode();
		hash = hash*31 + firstMolecule.getPdbChainCode().hashCode()+secondMolecule.getPdbChainCode().hashCode();
		return hash;
	}
	
	@Override
	public int compareTo(ChainInterface o) {
		return (Double.compare(this.interfaceArea,o.interfaceArea));
	}
}

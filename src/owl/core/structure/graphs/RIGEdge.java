package owl.core.structure.graphs;

/**
 * Class representing a Residue Interaction Graph edge.
 * The RIGEdge stores properties of the edge (like weights) but not the 
 * end points. Those can be taken from the RIGraph with getEndPoints(edge)
 */
public class RIGEdge {

	private static double DEFAULT_WEIGHT = 1.0;
	private static int DEFAULT_ATOM_WEIGHT = 1;
	private static double DEFAULT_PHI_PSI = 0.0;
	private static double DEFAULT_DISTANCE = 0.0;
	
	private double weight;
	private int atomWeight;  
	private double distance;
	private double phifrom;
	private double phito;
	private double psifrom;
	private double psito;

	
	public RIGEdge(double weight, double distance, int atomWeight, double phifrom, double phito, double psifrom, double psito) {
		this.weight = weight;
		this.distance = distance;
		this.atomWeight = atomWeight;
		this.phifrom = phifrom;
		this.phito = phito;
		this.psifrom = psifrom;
		this.psito = psito;
	}
	
	public RIGEdge(double weight, double distance, int atomWeight) {
		this.weight = weight;
		this.distance = distance;
		this.atomWeight = atomWeight;
		this.phifrom = DEFAULT_PHI_PSI;
		this.psifrom = DEFAULT_PHI_PSI;
		this.phito = DEFAULT_PHI_PSI;
		this.psito = DEFAULT_PHI_PSI;
	}
	
	public RIGEdge(int atomWeight){
		this.atomWeight = atomWeight;
		this.distance = DEFAULT_DISTANCE;
		this.weight = DEFAULT_WEIGHT;
		this.phifrom = DEFAULT_PHI_PSI;
		this.psifrom = DEFAULT_PHI_PSI;
		this.phito = DEFAULT_PHI_PSI;
		this.psito = DEFAULT_PHI_PSI;
	}
	
	public RIGEdge(double weight) {
		this.weight=weight;
		this.distance = DEFAULT_DISTANCE;
		this.atomWeight = DEFAULT_ATOM_WEIGHT;
		this.phifrom = DEFAULT_PHI_PSI;
		this.psifrom = DEFAULT_PHI_PSI;
		this.phito = DEFAULT_PHI_PSI;
		this.psito = DEFAULT_PHI_PSI;
	}
	
	public RIGEdge() {
		this.weight=DEFAULT_WEIGHT;
		this.distance = DEFAULT_DISTANCE;
		this.atomWeight = DEFAULT_ATOM_WEIGHT;
		this.phifrom = DEFAULT_PHI_PSI;
		this.psifrom = DEFAULT_PHI_PSI;
		this.phito = DEFAULT_PHI_PSI;
		this.psito = DEFAULT_PHI_PSI;
	}
	
	public int getAtomWeight() {
		return atomWeight;
	}
	
	public void setPhiPsi(double phifrom, double phito, double psifrom, double psito) {
		this.phifrom = phifrom;
		this.phito = phito;
		this.psifrom = psifrom;
		this.psito = psito;
	}
	
	public void setAtomWeight(int atomWeight) {
		this.atomWeight = atomWeight;
	}

	public double getWeight() {
		return weight;
	}
	
	public void setWeight(double weight) {
		this.weight = weight;
	}
	
	public double getDistance() {
		return distance;
	}
	
	public void setDistance(double distance) {
		this.distance = distance;
	}
	
	public double getPhiFrom() {
		return phifrom;
	}
	public double getPhiTo() {
		return phito;
	}
	public double getPsiFrom() {
		return psifrom;
	}
	public double getPsito() {
		return psito;
	}
	
	public String toString() {
		return "E";
	}
	
	public RIGEdge copy() {
		return new RIGEdge(this.weight,this.distance,this.atomWeight,this.phifrom, this.phito, this.psifrom, this.psito);
	}
	
}

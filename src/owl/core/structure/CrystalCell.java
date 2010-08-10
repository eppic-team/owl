package owl.core.structure;


/**
 * A crystal cell's parameters.
 * 
 * @author duarte_j
 *
 */
public class CrystalCell {

	private double a;
	private double b;
	private double c;
	
	private double alpha;
	private double beta;
	private double gamma;
	
	public CrystalCell(double a, double b, double c, double alpha, double beta, double gamma){
		this.a = a;
		this.b = b;
		this.c = c;
		this.alpha = alpha;
		this.beta = beta;
		this.gamma = gamma;
	}

	public double getA() {
		return a;
	}

	public void setA(double a) {
		this.a = a;
	}

	public double getB() {
		return b;
	}

	public void setB(double b) {
		this.b = b;
	}

	public double getC() {
		return c;
	}

	public void setC(double c) {
		this.c = c;
	}

	public double getAlpha() {
		return alpha;
	}

	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	public double getBeta() {
		return beta;
	}

	public void setBeta(double beta) {
		this.beta = beta;
	}

	public double getGamma() {
		return gamma;
	}

	public void setGamma(double gamma) {
		this.gamma = gamma;
	}
	
	/**
	 * Returns the volume of this unit cell.
	 * See http://en.wikipedia.org/wiki/Parallelepiped
	 * @return
	 */
	public double getVolume() {
		double alphaRad = Math.toRadians(alpha);
		double betaRad  = Math.toRadians(beta);
		double gammaRad = Math.toRadians(gamma);
		
		return a*b*c*
		Math.sqrt(1-Math.cos(alphaRad)*Math.cos(alphaRad)-Math.cos(betaRad)*Math.cos(betaRad)-Math.cos(gammaRad)*Math.cos(gammaRad)
				+2.0*Math.cos(alphaRad)*Math.cos(betaRad)*Math.cos(gammaRad));
	}
}

package owl.core.structure;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;


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
	
	private double alphaRad;
	private double betaRad;
	private double gammaRad;
	
	private Matrix3d transfToOrthonormal; // cached transformation crystal to orthonormal matrix
	private Matrix3d inverseTrasnfToOrthonormal; // cached inverse transformation crystal to orthonormal
	
	public CrystalCell(double a, double b, double c, double alpha, double beta, double gamma){
		this.a = a;
		this.b = b;
		this.c = c;
		this.alpha = alpha;
		this.beta = beta;
		this.gamma = gamma;
		
		this.alphaRad = Math.toRadians(alpha);
		this.betaRad  = Math.toRadians(beta);
		this.gammaRad = Math.toRadians(gamma);

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
		
		return a*b*c*
		Math.sqrt(1-Math.cos(alphaRad)*Math.cos(alphaRad)-Math.cos(betaRad)*Math.cos(betaRad)-Math.cos(gammaRad)*Math.cos(gammaRad)
				+2.0*Math.cos(alphaRad)*Math.cos(betaRad)*Math.cos(gammaRad));
	}
	
	/**
	 * Returns the unit cell transformation, given 3 integers for each of the directions of the unit cell.
	 * See "Fundamentals of Crystallography" C. Giacovazzo, section 2.5, eq 2.30
	 * @param direction
	 * @return
	 */
	public Matrix4d getTransform(Point3i direction) {
		Matrix4d mat = new Matrix4d();
		mat.setIdentity();
		
		mat.setTranslation(new Vector3d(getXtranslation(direction),getYtranslation(direction),getZtranslation(direction)));
		return mat;
	}

	public double getXtranslation(Point3i direction) {
		return  (double)direction.x*a+
				(double)direction.y*b*Math.cos(gammaRad)+
				(double)direction.z*c*Math.cos(betaRad);
	}
	
	public double getYtranslation(Point3i direction) {
		// see Table 2.1 of chapter 2 of Giacovazzo
		double cosAlphaStar = (Math.cos(betaRad)*Math.cos(gammaRad)-Math.cos(alphaRad))/(Math.sin(betaRad)*Math.sin(gammaRad));
		return  (double)direction.y*b*Math.sin(gammaRad)+
				(double)direction.z*(-1.0)*c*Math.sin(betaRad)*cosAlphaStar;
	}
	
	public double getZtranslation(Point3i direction) {
		// see Table 2.1 of chapter 2 of Giacovazzo
		double cStar = (this.a*this.b*Math.sin(gammaRad))/getVolume();
		return (double)direction.z*(1.0/cStar);
	}
	
	public Matrix3d getTransfToOrthonormal() {
		if (transfToOrthonormal!=null) {
			return transfToOrthonormal;
		}
		// see Table 2.1 of chapter 2 of Giacovazzo
		double cosAlphaStar = (Math.cos(betaRad)*Math.cos(gammaRad)-Math.cos(alphaRad))/(Math.sin(betaRad)*Math.sin(gammaRad));
		double cStar = (this.a*this.b*Math.sin(gammaRad))/getVolume();
		// see eq. 2.30 Giacovazzo
		double m32 = -this.c*Math.sin(betaRad)*cosAlphaStar;
		double m33 = 1.0/cStar;
		transfToOrthonormal =  new Matrix3d(                   this.a,                         0,    0,
							  				this.b*Math.cos(gammaRad), this.b*Math.sin(gammaRad),    0,
							  				this.c*Math.cos(betaRad) ,                       m32,  m33);
		return transfToOrthonormal;
	}
	
	public Matrix3d getInverseTransfToOrthonormal() {
		if (inverseTrasnfToOrthonormal!=null){
			return inverseTrasnfToOrthonormal;
		}
		inverseTrasnfToOrthonormal = new Matrix3d();
		inverseTrasnfToOrthonormal.invert(getTransfToOrthonormal());
		return inverseTrasnfToOrthonormal;
	}
}

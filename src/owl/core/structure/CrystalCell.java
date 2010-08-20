package owl.core.structure;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
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

	private Matrix3d M; 	// cached basis change transformation matrix
	private Matrix3d Minv;  // cached basis change transformation matrix
	private Matrix3d Mtransp;
	private Matrix3d MtranspInv;
	
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
	 * Returns the unit cell translation transformation matrix, given 3 integers for each 
	 * of the directions of the unit cell.
	 * See "Fundamentals of Crystallography" C. Giacovazzo, section 2.5, eq 2.30
	 * @param direction
	 * @return
	 */
	public Matrix4d getTransform(Vector3d direction) {
		Matrix4d mat = new Matrix4d();
		mat.setIdentity();
		mat.setTranslation(getTranslationVector(direction));
		return mat;
	}
	
	private Vector3d getTranslationVector(Vector3d direction) {
		Vector3d translationVec = new Vector3d();
		// see Giacovazzo section 2.E, eq. 2.E.1
		getMTransposeInv().transform(direction, translationVec);
		return translationVec;
	}
	
	/**
	 * Transform given Matrix4d in crystal basis to the orthonormal basis using
	 * the PDB convention (NCODE=1)
	 * @param m
	 * @return
	 */
	public Matrix4d transfToOrthonormal(Matrix4d m) {
		Vector3d trans = this.getTranslationVector(new Vector3d(m.m03,m.m13,m.m23));
		
		Matrix3d rot = new Matrix3d();
		m.getRotationScale(rot);
		// see Giacovazzo section 2.E, eq. 2.E.1
		// Rprime = MT-1 * R * MT
		rot.mul(this.getMTranspose());
		rot.mul(this.getMTransposeInv(),rot);

		return new Matrix4d(rot,trans,1.0);
	}

	/**
	 * Returns the change of basis (crystal to orthonormal) transform matrix, that is 
	 * M inverse in the notation of Giacovazzo using the PDB convention (NCODE=1). 
	 * The matrix is only calculated upon first call of this method, thereafter it is cached.
	 * See "Fundamentals of Crystallography" C. Giacovazzo, section 2.5 
	 * @return
	 */
	private Matrix3d getMInv() {
		if (Minv!=null) {
			return Minv;
		}
		// see Table 2.1 of chapter 2 of Giacovazzo
		double cosAlphaStar = (Math.cos(betaRad)*Math.cos(gammaRad)-Math.cos(alphaRad))/(Math.sin(betaRad)*Math.sin(gammaRad));
		double cStar = (this.a*this.b*Math.sin(gammaRad))/getVolume();
		// see eq. 2.30 Giacovazzo
		double m21 = -this.c*Math.sin(betaRad)*cosAlphaStar;
		double m22 = 1.0/cStar;
		Minv =  new Matrix3d(                    this.a,                         0,    0,
							  this.b*Math.cos(gammaRad), this.b*Math.sin(gammaRad),    0,
							  this.c*Math.cos(betaRad) ,                       m21,  m22);
		return Minv;
	}
	
	/**
	 * Returns the change of basis (orthonormal to crystal) transform matrix, that is
	 * M in the notation of Giacovazzo using the PDB convention (NCODE=1).
	 * The matrix is only calculated upon first call of this method, thereafter it is cached. 
	 * See "Fundamentals of Crystallography" C. Giacovazzo, section 2.5 
	 * @return
	 */
	private Matrix3d getM() {
		if (M!=null){
			return M;
		}
		M = new Matrix3d();
		M.invert(getMInv());
		return M;
	}
	
	private Matrix3d getMTranspose() {
		if (Mtransp!=null){
			return Mtransp;
		}
		Matrix3d M = getM();
		Mtransp = new Matrix3d();
		Mtransp.transpose(M);
		return Mtransp;
	}
	
	private Matrix3d getMTransposeInv() {
		if (MtranspInv!=null){
			return MtranspInv;
		}
		MtranspInv = new Matrix3d();
		MtranspInv.invert(getMTranspose());
		return MtranspInv;
	}
}

package owl.core.structure;

import java.io.Serializable;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3i;

/**
 * Representation of a transformation in a crystal: 
 * - a transformation id (each of the transformations in a space group, 0 to m)
 * - a crystal translation
 * The transformation matrix in crystal basis is stored, representing the basic 
 * transformation together with the crystal translation. 
 * Contains methods to check for equivalent transformations.
 * 
 * 
 * @author duarte_j
 *
 */
public class CrystalTransform implements Serializable {

	private static final long serialVersionUID = 1L;


	public static final Matrix4d IDENTITY = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);


	/**
	 * The transform id corresponding to the SpaceGroup's transform indices. 
	 * From 0 (identity) to m (m=number of symmetry operations of the space group)
	 * It is unique within the unit cell but equivalent units of different crystal unit cells 
	 * will have same id 
	 */
	private int transformId;
	
	/**
	 * The 4-dimensional matrix transformation in crystal basis.
	 * Note that the translational component of this matrix is not necessarily
	 * identical to crystalTranslation since some operators have fractional 
	 * translations within the cell
	 */
	private Matrix4d matTransform;
	
	/**
	 * The crystal translation (always integer)
	 */
	private Point3i crystalTranslation;
	
	
	/**
	 * Creates a new CrystalTransform representing the identity transform 
	 * in cell (0,0,0)
	 */
	public CrystalTransform() {
		this.transformId = 0;
		this.matTransform = (Matrix4d)IDENTITY.clone();
		this.crystalTranslation = new Point3i(0,0,0);
	}
	
	public CrystalTransform(SpaceGroup sg, int transformId) {
		this.transformId = transformId;
		this.matTransform = (Matrix4d)sg.getTransformation(transformId).clone();
		this.crystalTranslation = new Point3i(0,0,0);

	}
	
	public CrystalTransform(CrystalTransform transform) {
		this.transformId = transform.transformId;
		this.matTransform = new Matrix4d(transform.matTransform);
		this.crystalTranslation = new Point3i(transform.crystalTranslation);
		
	}
	
	public Matrix4d getMatTransform() {
		return matTransform;
	}
	
	public void setMatTransform(Matrix4d matTransform) {
		this.matTransform = matTransform;
	}
	
	public Point3i getCrystalTranslation() {
		return crystalTranslation;
	}
	
	public void translate(Point3i translation) {
		matTransform.m03 = matTransform.m03+(double)translation.x;
		matTransform.m13 = matTransform.m13+(double)translation.y;
		matTransform.m23 = matTransform.m23+(double)translation.z;
		
		crystalTranslation.add(translation); 

	}
	
	/**
	 * Returns true if the given CrystalTransform is equivalent to this one.
	 * Two crystal transforms are equivalent if one is the inverse of the other, i.e.
	 * their transformation matrices multiplication is equal to the identity.
	 * @param other
	 * @return
	 */
	public boolean isEquivalent(CrystalTransform other) {
		Matrix4d mul = new Matrix4d();
		mul.mul(this.matTransform,other.matTransform);

		if (mul.epsilonEquals(IDENTITY, 0.0001)) {
			return true;
		}
		return false;
	}
	
	public boolean isPureCrystalTranslation() {
		return transformId==0;
	}
	
	public int getTransformId() {
		return transformId;
	}
	
	public void setTransformId(int transformId) {
		this.transformId = transformId;
	}
	
	public String toString() {
		return String.format("[%2d-(%2d,%2d,%2d)]",transformId,crystalTranslation.x,crystalTranslation.y,crystalTranslation.z);
	}

}

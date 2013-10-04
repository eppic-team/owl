package owl.core.util;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

/**
 * A class to contain the results of an rmsd calculation using the Kabsch algorithm:
 * a 3x3 optimal superposition matrix, the translation vector of the centroids and the rmsd value. 
 * 
 * 
 * @author duarte_j
 *
 */
public class OptSuperposition {
	
	private double rmsd;
	private Matrix3d supMatrix;
	private Vector3d transl;
	
	public OptSuperposition(double rmsd, Matrix3d supMatrix, Vector3d transl) {
		this.rmsd = rmsd;
		this.supMatrix = supMatrix;
		this.transl = transl;
	}
	
	public double getRmsd() {
		return rmsd;
	}
	
	public Matrix3d getSupMatrix() {
		return supMatrix;
	}
	
	public Vector3d getCentroidsTranslation() {
		return transl;
	}
	
	public Matrix4d getTransformMatrix() {
		return new Matrix4d(supMatrix, transl, 1.0);
	}

}

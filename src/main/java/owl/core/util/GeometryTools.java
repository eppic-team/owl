package owl.core.util;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeMap;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Tuple3d;
import javax.vecmath.Vector3d;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import Jama.SingularValueDecomposition;
import owl.core.structure.Atom;
import owl.core.structure.BackboneCoords;

/**
 * Some tools for molecular geometry calculations.
 * @author stehr
 *
 */
public class GeometryTools {
	
	private static final double DELTA=0.0000001;
	
	/**
	 * Returns the angle between two plains defined by six points.
	 * Plain 1 is defined by points a1,a2,a3 and plain 2 is defined by points b1,b2,b3.
	 * @return the torsion angle in radians
	 */
	public static double torsionAngle(Point3d a1, Point3d a2, Point3d a3, Point3d b1, Point3d b2, Point3d b3) {
		Vector3d vA1 = new Vector3d(a2);
		vA1.sub(a1);
		
		Vector3d vA2 = new Vector3d(a3);
		vA2.sub(a2);
		
		Vector3d vB1 = new Vector3d(b2);
		vB1.sub(b1);
		
		Vector3d vB2 = new Vector3d(b3);
		vB2.sub(b2);
		
		Vector3d normA = new Vector3d();	// norm vector for plane A
		normA.cross(vA1, vA2);
		
		Vector3d normB = new Vector3d();	// norm vector for plane B
		normB.cross(vB1, vB2);
		
		double cos = normA.dot(normB) / (normA.length()*normB.length());	// angle between planes
		
		double sign = Math.signum(normA.dot(vB2)); 				// sign of the angle

		return sign * Math.acos(cos);
	}
	
	/**
	 * Returns the dihedral angle in radians given the coordinates of four consecutive atoms in a polymer.
	 * e.g. for amino acid i
	 * phi: C(i-1), N, C-alpha, C
	 * psi: N, C-alpha, C, N(i+1)
	 * omega: C-alpha, C, N(i+1), C-alpha(i+1)
	 * @return the dihedral angle in radians
	 */
	public static double dihedralAngle(Point3d p1, Point3d p2, Point3d p3, Point3d p4) {
		return torsionAngle(p1,p2,p3,p2,p3,p4);
	}
	
	/**
	 * Read backbone coordinates in pdb format (see http://www.wwpdb.org/documentation/format32/v3.2.html)
	 * @param pdbFile the input file
	 * @param chain the chain for which coordinates are read
	 * @return a map from residue number to Backbone objects holding the atom coordinates
	 * @throws IOException
	 */
	public static TreeMap<Integer, BackboneCoords> readBackbone(File pdbFile, char chain) throws IOException {
		TreeMap<Integer, BackboneCoords> backbone = new TreeMap<Integer, BackboneCoords>();
		BufferedReader in = new BufferedReader(new FileReader(pdbFile));
		String line;
		while((line=in.readLine()) != null) {
			if(!line.startsWith("ATOM")) continue;		// process only atom lines
			if(line.charAt(21) != chain) continue;		// process only given chain
			Integer resNum = new Integer(line.substring(22,26).trim());
			String atomType = line.substring(13, 15).trim();
			double x = Double.parseDouble(line.substring(31, 38).trim());
			double y = Double.parseDouble(line.substring(38, 46).trim());
			double z = Double.parseDouble(line.substring(46, 54).trim());			
			BackboneCoords bb;
			if(backbone.containsKey(resNum)) {
				bb = backbone.get(resNum);
			} else {
				bb = new BackboneCoords();
			}
			if(atomType.equals("CA")) bb.Ca = new Point3d(x,y,z);
                        if(atomType.equals("C")) bb.C = new Point3d(x,y,z);
                        if(atomType.equals("N")) bb.N = new Point3d(x,y,z);
                        if(atomType.equals("O")) bb.O = new Point3d(x,y,z);
			backbone.put(resNum, bb);
		}
		in.close();
		return backbone;
	}
	
	/**
	 * Returns a predicted secondary structure state based on dihedral angles
	 * @param phi the phi angle
	 * @param psi the psi angle
	 * @return one of the states 'H' (helix), 'E' (exteded strand), 'O' (other)
	 */
	public static char classifyByAngle(double phi, double psi) {
		return 'O';
	}
	
	/**
	 * Returns a predicted secondary structure state based on Ca_i-Ca_iplus3 distance
	 * @param distCa the distance between the C-alpha atoms at position i and i+3
	 * @return one of the states 'H' (helix), 'E' (exteded strand), 'O' (other)
	 */
	public static char classifyByDistance(double distCa) {
		return 'O';
	}
	
	/**
	 * Gets the momemt of inertia tensor for the given center point and set of atoms.
	 * @param center
	 * @param atoms
	 * @return a 3x3 array with containing the moments of inertia tensor
	 */
	public static double[][] getMomentInertiaTensor(Point3d center, Atom[] atoms) {
		double[][] tensor = new double[3][3];
		for (Atom atom:atoms) {
			Point3d coords = new Point3d(atom.getCoords());
			coords.sub(center);
			double m = atom.getType().getAtomicMass();
			tensor[0][0] += m*(coords.y*coords.y+coords.z*coords.z); 
			tensor[1][1] += m*(coords.x*coords.x+coords.z*coords.z);
			tensor[2][2] += m*(coords.x*coords.x+coords.y*coords.y);
			tensor[0][1] += -m*coords.x*coords.y;
			tensor[0][2] += -m*coords.x*coords.z;
			tensor[1][2] += -m*coords.y*coords.z;
		}
		tensor[1][0] = tensor[0][1];
		tensor[2][0] = tensor[0][2];
		tensor[2][1] = tensor[1][2];
		return tensor;
	}
	
	/**
	 * Gets the principal moments of inertia, i.e. the result of diagonalizing the
	 * moments of inertia tensor.
	 * @param center
	 * @param atoms
	 * @return
	 */
	public static ArrayList<InertiaMomentAndAxis> getPrincipalMomentsInertia(Point3d center, Atom[] atoms) {
		Matrix tensor = new Matrix(getMomentInertiaTensor(center,atoms));
		EigenvalueDecomposition eig = tensor.eig();
		Matrix eigVals = eig.getD();
		Matrix eigVecs = eig.getV();
		//eigVals.print(12, 3);
		//eigVecs.print(12, 3);
		ArrayList<InertiaMomentAndAxis> momentsAndAxes = new ArrayList<InertiaMomentAndAxis>();
		for (int i=0;i<3;i++) {
			momentsAndAxes.add(new InertiaMomentAndAxis(eigVals.get(i, i), new Vector3d(eigVecs.get(0, i), eigVecs.get(1, i), eigVecs.get(2, i))));
		}
		return momentsAndAxes;
	}
	
	/**
	 * Gets the inertia axis corresponding to the maximum principal moment of 
	 * inertia.
	 * @param center
	 * @param atoms
	 * @return
	 */
	public static Vector3d getBiggestInertiaAxis(Point3d center, Atom[] atoms) {
		ArrayList<InertiaMomentAndAxis> momentsAndAxes = getPrincipalMomentsInertia(center,atoms);
		Collections.sort(momentsAndAxes); // sorts in ascending order: biggest value is last
		return momentsAndAxes.get(2).axis;
	}
	
	/**
	 * Given a transformation matrix containing a rotation and translation returns the
	 * screw component of the rotation. 
	 * See http://www.crystallography.fr/mathcryst/pdf/Gargnano/Aroyo_Gargnano_1.pdf
	 * @param m
	 * @return
	 */
	public static Vector3d getTranslScrewComponent(Matrix4d m) {
		
		int foldType = getRotAxisType(m);
		// For reference see:
		// http://www.crystallography.fr/mathcryst/pdf/Gargnano/Aroyo_Gargnano_1.pdf

		Vector3d transl = null;

		Matrix3d W = 
				new Matrix3d(m.m00,m.m01,m.m02,
						m.m10,m.m11,m.m12,
						m.m20,m.m21,m.m22);

		if (foldType>=0) {

			// the Y matrix: Y = W^k-1 + W^k-2 ... + W + I  ; with k the fold type
			Matrix3d Y = new Matrix3d(1,0,0, 0,1,0, 0,0,1);					
			Matrix3d Wk = new Matrix3d(1,0,0, 0,1,0, 0,0,1);

			for (int k=0;k<foldType;k++) {						
				Wk.mul(W); // k=0 Wk=W, k=1 Wk=W^2, k=2 Wk=W^3, ... k=foldType-1, Wk=W^foldType
				if (k!=foldType-1) Y.add(Wk);
			}

			transl = new Vector3d(m.m03, m.m13, m.m23);
			Y.transform(transl);

			transl.scale(1.0/foldType);

		} else {

			if (foldType==-2) { // there are glide planes only in -2
				Matrix3d Y = new Matrix3d(1,0,0, 0,1,0, 0,0,1);
				Y.add(W);

				transl = new Vector3d(m.m03, m.m13, m.m23);
				Y.transform(transl);

				transl.scale(1.0/2.0);
			} else { // for -1, -3, -4 and -6 there's nothing to do: fill with 0s 
				transl = new Vector3d(0,0,0);
			}
		}

		return transl;
	}
	
	/**
	 * Given a transformation matrix containing a rotation returns the type of rotation:
	 * 1 for identity, 2 for 2-fold rotation, 3 for 3-fold rotation, 4 for 4-fold rotation, 
	 * 6 for 6-fold rotation,
	 * -1 for inversions, -2 for mirror planes, -3 for 3-fold improper rotation, 
	 * -4 for 4-fold improper rotation and -6 for 6-fold improper rotation
	 * @param m
	 * @return
	 */
	public static int getRotAxisType(Matrix4d m) {
		int axisType = 0;
		
		Matrix3d rot = new Matrix3d(m.m00,m.m01,m.m02,
				m.m10,m.m11,m.m12,
				m.m20,m.m21,m.m22);
		
		double determinant = rot.determinant();
		
		if (!deltaComp(determinant,1.0) && !deltaComp(determinant, -1.0)) {
			throw new IllegalArgumentException("Given matrix does not seem to be a rotation matrix.");
		}
		
		int trace = (int)(rot.m00+rot.m11+rot.m22);
		if (determinant>0) {
			switch (trace) {
			case 3:
				axisType=1;
				break;
			case -1:
				axisType=2;
				break;
			case 0:
				axisType=3;
				break;						
			case 1:
				axisType=4;
				break;
			case 2:
				axisType=6;
				break;
			default:
				throw new NullPointerException("Trace of transform does not correspond to one of the expected types. This is most likely a bug");
			}
		} else {
			switch (trace) {
			case -3:
				axisType=-1;
				break;
			case 1:
				axisType=-2;
				break;
			case 0:
				axisType=-3;
				break;						
			case -1:
				axisType=-4;
				break;
			case -2:
				axisType=-6;
				break;
			default:
				throw new NullPointerException("Trace of transform does not correspond to one of the expected types. This is most likely a bug");
			}
		}
		return axisType;
	}
	
	/**
	 * Given a rotation matrix calculates the rotation axis and angle for it.
	 * The angle is calculated from the trace, the axis from the eigenvalue
	 * decomposition.
	 * If given matrix is improper rotation or identity matrix then 
	 * axis (0,0,0) and angle 0 are returned.
	 * @param m
	 * @return
	 * @throws IllegalArgumentException if given matrix is not a rotation matrix (determinant not 1 or -1)
	 */
	public static AxisAngle4d getRotAxisAndAngle(Matrix3d m) {	
		double determinant = m.determinant();
		
		if (!(Math.abs(determinant)-1.0<DELTA)) throw new IllegalArgumentException("Given matrix is not a rotation matrix"); 
		
		AxisAngle4d axisAndAngle = new AxisAngle4d(new Vector3d(0,0,0),0);
		
		double[] d = {m.m00,m.m10,m.m20,
				m.m01,m.m11,m.m21,
				m.m02,m.m12,m.m22};

		Matrix r = new Matrix(d,3);
		
		if (!deltaComp(r.det(), 1.0)) {
			// improper rotation: we return axis 0,0,0 and angle 0
			return axisAndAngle;
		}
		
		EigenvalueDecomposition evd = new EigenvalueDecomposition(r);
		
		Matrix eval = evd.getD();
		if (deltaComp(eval.get(0, 0),1.0) && deltaComp(eval.get(1, 1),1.0) && deltaComp(eval.get(2, 2),1.0)) {
			// the rotation is an identity: we return axis 0,0,0 and angle 0
			return axisAndAngle;
		}
		int indexOfEv1;
		for (indexOfEv1=0;indexOfEv1<3;indexOfEv1++) {
			if (deltaComp(eval.get(indexOfEv1, indexOfEv1),1)) break;
		}
		Matrix evec = evd.getV(); 
		axisAndAngle.set(new Vector3d(evec.get(0,indexOfEv1), evec.get(1, indexOfEv1), evec.get(2, indexOfEv1)), 
						Math.acos((eval.trace()-1.0)/2.0));
		
		return axisAndAngle;
	}
	
	/**
	 * Calculates the optimal superposition (minimal RMSD) between two conformations.      
	 * conformation1: Vector3d array (matrix of dimensions [N,3])       
	 * conformation2: Vector3d array (matrix of dimensions [N,3]) 
	 * 
	 * Beware, if option transform is set to true the coordinates that are transformed are ONLY those passed: 
	 * usually they will be a subset of the structure (e.g. only CA atoms) and thus the rest of the structure
	 * will keep its old coordinates. 
	 * 
	 * The output optimal superposition contains the rotation matrix needed to transform
	 * conformation1 to be superimposed onto conformation2.
	 * 
	 * Implementation taken (python) from: 
	 *  http://boscoh.com/protein/rmsd-root-mean-square-deviation
	 * then ported to java using Jama matrix package  
	 * 
	 * See also:
	 *  http://en.wikipedia.org/wiki/Kabsch_algorithm
	 *  http://cnx.org/content/m11608/latest/
	 * 
	 * @param conformation1
	 * @param conformation2
	 * @param transform if true conformation1 and conformation2 will be transformed to be 
	 * optimally superposed: first brought to a common center (subtracting centroids) and then 
	 * conformation1 rotated to be optimally superimposed to conformation2. If false the given
	 * conformations are not altered at all.
	 * @return
	 * @throws IllegalArgumentException if the 2 given arrays are not of the same size or if they have 0 size
	 */
	public static OptSuperposition calcOptimalSuperposition(Tuple3d[] conformation1, Tuple3d[] conformation2, boolean transform) {
		if (conformation1.length!=conformation2.length) {
			throw new IllegalArgumentException(
					"Given conformations have different size: conformation1: "+conformation1.length+", conformation2: "+conformation2.length);
		}
		
		int n = conformation1.length;
		
		if (n==0) throw new IllegalArgumentException("The given conformations are of 0 size");


		Tuple3d[] conf1 = null;
		Tuple3d[] conf2 = null;
		
		if (transform) {
			conf1 = conformation1;
			conf2 = conformation2;
		} else {
			conf1 = copyConformation(conformation1);
			conf2 = copyConformation(conformation2);
		}
		
		
		// 1st we bring both conformations to the same centre by subtracting their respective centroids
		Point3d center1 = getCentroid(conf1);
		Point3d center2 = getCentroid(conf2);

		// translating our conformations to the same coordinate system by subtracting centers
		for (Tuple3d vec:conf1){
			vec.sub(center1);
		}
		for (Tuple3d vec:conf2){
			vec.sub(center2);
		}

		//E0: initial sum of squared lengths of both conformations
		double sum1 = 0.0;
		double sum2 = 0.0;
		for (int i=0;i<n;i++){			
			sum1 += (new Vector3d(conf1[i])).lengthSquared();
			sum2 += (new Vector3d(conf2[i])).lengthSquared();
		}
		double E0 = sum1 + sum2;


		Matrix covarianceMatrix = new Matrix(3,3);
		// this would be the same operation as we do below but in jama Matrix notation: 
		//Matrix covarianceMatrix = vecs2.transpose().times(vecs1); //gives a 3x3 matrix
		double[][] covMatrixArray = covarianceMatrix.getArray(); // this is just a reference to covarianceMatrix
		for (int k=0;k<n;k++) {
			covMatrixArray[0][0] += conf2[k].x * conf1[k].x;
			covMatrixArray[0][1] += conf2[k].x * conf1[k].y;
			covMatrixArray[0][2] += conf2[k].x * conf1[k].z;

			covMatrixArray[1][0] += conf2[k].y * conf1[k].x;
			covMatrixArray[1][1] += conf2[k].y * conf1[k].y;
			covMatrixArray[1][2] += conf2[k].y * conf1[k].z;
			
			covMatrixArray[2][0] += conf2[k].z * conf1[k].x;
			covMatrixArray[2][1] += conf2[k].z * conf1[k].y;
			covMatrixArray[2][2] += conf2[k].z * conf1[k].z;

		}

		// singular value decomposition
		SingularValueDecomposition svd = covarianceMatrix.svd();
		Matrix U = svd.getU();
		Matrix V_trans = svd.getV().transpose(); 
		double[] singularValues = svd.getSingularValues();

		boolean isReflection = false;
		if (U.det()*V_trans.det()<0.0){ 
			isReflection = true;
		}
		if (isReflection){
			// reflect along smallest principal axis:
			// we change sign of last coordinate (smallest singular value)
			singularValues[singularValues.length-1]=(-1)*singularValues[singularValues.length-1];  			
		}

		// getting sum of singular values
		double sumSV = 0.0;
		for (int i=0;i<singularValues.length;i++){
			sumSV += singularValues[i];
		}

		// rmsd square: Kabsch formula
		double rmsd_sq = (E0 - 2.0*sumSV)/((double) n);
		rmsd_sq = Math.max(rmsd_sq, 0.0);

		if (isReflection) { // first we check if we are in is_reflection case: we need to change sign to last row of U
			for (int j=0;j<U.getColumnDimension();j++){
				// we change sign to last row of U
				int lastRow = U.getRowDimension()-1;
				U.set(lastRow, j, (-1)*U.get(lastRow,j));
			}
		}
		Matrix optimalRotation = U.times(V_trans); 

		Matrix3d supMat = new Matrix3d(optimalRotation.getRowPackedCopy());
		double rmsd = Math.sqrt(rmsd_sq);
		Vector3d transl = new Vector3d(center1);
		transl.sub(center2); 

		if (transform) {
			for (int i=0;i<n;i++) {
				supMat.transform(conf1[i]);
			}
		}
		
		return new OptSuperposition(rmsd, supMat, transl);
	}
	
	/**
	 * Calculates the centroid of given coordinates
	 * @param conformation
	 * @return
	 */
	public static Point3d getCentroid(Tuple3d[] conformation) {
		Point3d center = new Point3d();

		for (int i=0;i<conformation.length;i++){ 
			center.add(conformation[i]);

		}
		
		center.scale((double)1/conformation.length);

		return center;
	}
	
	/**
	 * Deep copies the given conformation
	 * @param conformation
	 * @return
	 */
	private static Tuple3d[] copyConformation(Tuple3d[] conformation) {
		Tuple3d[] newConformation = new Tuple3d[conformation.length];
		for (int i=0;i<conformation.length;i++) {
			newConformation[i] = new Vector3d(conformation[i]);
		}
		return newConformation;
	}
	
	/**
	 * Calculates the root mean square distance of the given coordinates as given
	 * without optimally superposing them.
	 * @param conformation1
	 * @param conformation2
	 * @return
	 */
	public static double getCoordinatesRmsd(Tuple3d[] conformation1, Tuple3d[] conformation2) {
		if (conformation1.length!=conformation2.length) {
			throw new IllegalArgumentException(
					"Given conformations have different size: conformation1: "+conformation1.length+", conformation2: "+conformation2.length);
		}
		
		int n = conformation1.length;
		
		if (n==0) throw new IllegalArgumentException("The given conformations are of 0 size");
		
		double sumSqDistances = 0.0;
		for (int i=0;i<n;i++) {
			sumSqDistances += (new Point3d(conformation1[i])).distanceSquared(new Point3d(conformation2[i]));
		}
		return Math.sqrt(sumSqDistances/(double)n);
	}
	
	private static boolean deltaComp(double d1, double d2) {
		return Math.abs(d1-d2)<DELTA;
	}
	
	/**
	 * Takes the name of a pdb file, reads in the coordinates and prints out the dihedral angles.
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		if(args.length < 1) {
			System.out.println("Usage: " + GeometryTools.class.getName() + " <pdb_file_name>");
			System.exit(1);
		}

		File pdbFile = new File(args[0]);
		char chain = 'A';	// define the chain to be processed
		
		// read coordinates
		TreeMap<Integer,BackboneCoords> myBackbone = readBackbone(pdbFile, chain);

		// calculate angles and predict secondary structure state
		System.out.println("resNum\tphi\tpsi\tomega\tpredSS\tCaDist\tpredSS");
		int maxResNum = myBackbone.lastKey();
		for (int i = 2; i <= maxResNum - 1; i++) {
			BackboneCoords bb_i = myBackbone.get(i);
			BackboneCoords bb_iplus1 = myBackbone.get(i+1);
			BackboneCoords bb_iminus1 = myBackbone.get(i-1);
			double phi = Math.toDegrees(dihedralAngle(bb_iminus1.C, bb_i.N, bb_i.Ca, bb_i.C)); // phi: C(i-1), N, C-alpha, C
			double psi = Math.toDegrees(dihedralAngle(bb_i.N, bb_i.Ca, bb_i.C, bb_iplus1.N));  // psi: N, C-alpha, C, N(i+1)
			double omega = Math.toDegrees(dihedralAngle(bb_i.Ca, bb_i.C, bb_iplus1.N, bb_iplus1.Ca)); // omega: C-alpha, C, N(i+1), C-alpha(i+1)
			char ssPred = classifyByAngle(phi, psi);
			System.out.printf("%3d\t%4.0f\t%4.0f\t%4.0f\t%s", i, phi, psi, omega, ssPred);
			if(myBackbone.containsKey(i+3)) {
				BackboneCoords bb_iplus3 = myBackbone.get(i+3);
				double distCa = bb_i.Ca.distance(bb_iplus3.Ca);
				char ssPred2 = classifyByDistance(distCa);
				System.out.printf("\t%4.1f\t%s", distCa, ssPred2);
			}
			System.out.println();
		}
	}
}

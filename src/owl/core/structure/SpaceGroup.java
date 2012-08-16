package owl.core.structure;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;


/**
 * A crystallographic space group. We store the standard numeric identifier,
 * the international short symbol and the transformations corresponding to
 * each space group (as Matrix4ds and in algebraic notation).
 * The information for all (protein crystallography) space groups can be 
 * parsed from CCP4 package's symop.lib file.
 * 
 * See: http://en.wikipedia.org/wiki/Space_group
 * 
 * @author duarte_j
 * @see SymoplibParser
 */
public final class SpaceGroup implements Serializable {

	private static final long serialVersionUID = 1L;

	public enum BravaisLattice {
		
		TRICLINIC    (1, "TRICLINIC",    new CrystalCell(1.00,1.25,1.50, 60,70,80)), // alpha,beta,gamma!=90
		MONOCLINIC   (2, "MONOCLINIC",   new CrystalCell(1.00,1.25,1.50, 90,60,90)), // beta!=90, alpha=gamma=90
		ORTHORHOMBIC (3, "ORTHORHOMBIC", new CrystalCell(1.00,1.25,1.50, 90,90,90)), // alpha=beta=gamma=90
		TETRAGONAL   (4, "TETRAGONAL",   new CrystalCell(1.00,1.00,1.25, 90,90,90)), // alpha=beta=gamma=90, a=b
		TRIGONAL     (5, "TRIGONAL",     new CrystalCell(1.00,1.00,1.25, 90,90,120)),// a=b!=c, alpha=beta=90, gamma=120
		HEXAGONAL    (6, "HEXAGONAL",    new CrystalCell(1.00,1.00,1.25, 90,90,120)),// a=b!=c, alpha=beta=90, gamma=120
		CUBIC        (7, "CUBIC",        new CrystalCell(1.00,1.00,1.00, 90,90,90)); // a=b=c, alpha=beta=gamma=90
		
		private static HashMap<String, BravaisLattice> name2bl = initname2bl();
		private String name;
		private int id;
		private CrystalCell exampleUnitCell;
		private BravaisLattice(int id, String name, CrystalCell exampleUnitCell) {
			this.name = name;
			this.id = id;
			this.exampleUnitCell = exampleUnitCell;
		}
		public String getName() {return name;}
		public int getId() {return id;}
		public CrystalCell getExampleUnitCell() {return exampleUnitCell;}
		
		private static HashMap<String,BravaisLattice> initname2bl(){
			HashMap<String,BravaisLattice> name2bl = new HashMap<String, SpaceGroup.BravaisLattice>();
			for (BravaisLattice bl:BravaisLattice.values()) {
				name2bl.put(bl.getName(), bl);
			}
			return name2bl;
		}
		public static BravaisLattice getByName(String blName) {
			return name2bl.get(blName);
		}
	}
	
	private static final Pattern splitPat1 = Pattern.compile("((?:[+-]?[XYZ])+)([+-][0-9/.]+)");
	private static final Pattern splitPat2 = Pattern.compile("([+-]?[0-9/.]+)((?:[+-][XYZ])+)");
	private static final Pattern coordPat = Pattern.compile("(?:([+-])?([XYZ]))+?"); // the last +? is for ungreedy matching
	private static final Pattern transCoefPat = Pattern.compile("([-+]?[0-9.]+)(?:/([0-9.]+))?");
	
	private static final Pattern nonEnantPat = Pattern.compile("[-abcmnd]");
	
	private static final double delta=0.0000001;
	
	private final int id;
	private final String shortSymbol;
	private final List<Matrix4d> transformations;
	private final List<String> transfAlgebraic;
	
	private List<Double> angles; // all rotation angles except for the identity transformation (thus one element fewer than transformations)
	private List<Vector3d> axes; // all rotation axes except for the identity transformation (thus one element fewer than transformations)
	
	// we store here the transformIds (the rotational ones) for each of the possible rotation axes (the key of the map)
	private HashMap<Vector3d,ArrayList<Integer>> operatorsPerAxis; // axes of rotation to lists of operator ids
	
	private HashMap<Integer,Integer> transformId2AxisType;
	
	private BravaisLattice bravLattice;
	
	public SpaceGroup(int id, String shortSymbol, BravaisLattice bravLattice) {
		this.id = id;
		this.shortSymbol = shortSymbol;
		transformations = new ArrayList<Matrix4d>();
		transfAlgebraic = new ArrayList<String>();
		this.bravLattice = bravLattice;
	}
	
	public void addTransformation(String transfAlgebraic) {
		this.transfAlgebraic.add(transfAlgebraic);
		this.transformations.add(getMatrixFromAlgebraic(transfAlgebraic));
	}
	
	public static Matrix4d getMatrixFromAlgebraic(String transfAlgebraic) {
		String[] parts = transfAlgebraic.toUpperCase().split(",");
		double[] xCoef = convertAlgebraicStrToCoefficients(parts[0].trim());
		double[] yCoef = convertAlgebraicStrToCoefficients(parts[1].trim());
		double[] zCoef = convertAlgebraicStrToCoefficients(parts[2].trim());
		
		Matrix4d mat = new Matrix4d();
		mat.setIdentity();
		mat.setRotation(new Matrix3d(xCoef[0],xCoef[1],xCoef[2],yCoef[0],yCoef[1],yCoef[2],zCoef[0],zCoef[1],zCoef[2]));
		mat.setTranslation(new Vector3d(xCoef[3],yCoef[3],zCoef[3]));
		return mat;
		//return new Matrix4d(xCoef[0],xCoef[1],xCoef[2],xCoef[3],
		//					yCoef[0],yCoef[1],yCoef[2],yCoef[3],
		//					zCoef[0],zCoef[1],zCoef[2],zCoef[3],
		//					0,0,0,1);
	}
	
	private static double[] convertAlgebraicStrToCoefficients(String algString) {
		String letters = null;
		String noLetters = null;
		Matcher m = splitPat1.matcher(algString);
		if (m.matches()) {
			letters = m.group(1);
			noLetters = m.group(2);
		} else {
			m = splitPat2.matcher(algString);
			if (m.matches()) {
				letters = m.group(2);
				noLetters = m.group(1);
			} else {
				letters = algString;
			}
		}
		double[] coefficients = new double[4];
		m = coordPat.matcher(letters);
		while(m.find()){
			String sign = "";
			if (m.group(1)!=null) {
				sign = m.group(1);
			}
			double s = 1.0;
			if (sign.equals("-")){
				s = -1.0;
			}
			String coord = m.group(2);
			if (coord.equals("X")) {
				coefficients[0] = s;
			} else if (coord.equals("Y")) {
				coefficients[1] = s;
			} else if (coord.equals("Z")) {
				coefficients[2] = s;
			}
		}		
		if (noLetters!=null) {
			m = transCoefPat.matcher(noLetters);
			if (m.matches()) {
				double num = Double.parseDouble(m.group(1));
				double den = 1;
				if (m.group(2)!=null) {
					den = Double.parseDouble(m.group(2));
				}
				coefficients[3] = num/den;
			}
		} else {
			coefficients[3]=0;
		}
		return coefficients;
	}
	
	public int getId() {
		return id;
	}
	
	public String getShortSymbol() {
		return shortSymbol;
	}
	
	/**
	 * Gets all transformations except for the identity in crystal axes basis.
	 * @return
	 */
	public List<Matrix4d> getTransformations() {
		List<Matrix4d> transfs = new ArrayList<Matrix4d>();
		for (int i=1;i<this.transformations.size();i++){
			transfs.add(transformations.get(i));
		}
		return transfs;
	}
	
	private void calcRotAxesAndAngles() {
		// beware this works only for non enantiomorphic groups (those without inversion points/mirror planes)
		angles = new ArrayList<Double>();
		axes = new ArrayList<Vector3d>();
		operatorsPerAxis = new HashMap<Vector3d, ArrayList<Integer>>();
		transformId2AxisType = new HashMap<Integer, Integer>();
		transformId2AxisType.put(0, 1); // for operator transformId=0, axis type is "1-fold"
		
		for (int i=1;i<this.transformations.size();i++){
			double[] d = {transformations.get(i).m00,transformations.get(i).m01,transformations.get(i).m02,
					transformations.get(i).m10,transformations.get(i).m11,transformations.get(i).m12,
					transformations.get(i).m20,transformations.get(i).m21,transformations.get(i).m22};
			Matrix r = new Matrix(d,3);
			EigenvalueDecomposition evd = new EigenvalueDecomposition(r);
			
			Matrix eval = evd.getD();
			if (eval.get(0, 0)==1 && eval.get(1, 1)==1 && eval.get(2, 2)==1) {
				// the rotation is an identity: no axis, we add a (0,0,0)
				axes.add(new Vector3d(0,0,0)); 
				angles.add(0.0);
				transformId2AxisType.put(i, 1);
				continue;
			}
			int indexOfEv1;
			for (indexOfEv1=0;indexOfEv1<3;indexOfEv1++) {
				if (deltaComp(eval.get(indexOfEv1, indexOfEv1),1,0.00001)) break;
			}
			Matrix evec = evd.getV();
			double angle = Math.acos((eval.trace()-1)/2);
			angles.add(angle);
			Vector3d axis = new Vector3d(evec.get(0,indexOfEv1),evec.get(1, indexOfEv1),evec.get(2, indexOfEv1));
			axes.add(axis);
			if (operatorsPerAxis.containsKey(axis)) {
				operatorsPerAxis.get(axis).add(i);
			} else {
				ArrayList<Integer> ops = new ArrayList<Integer>();
				ops.add(i);
				operatorsPerAxis.put(axis, ops);
			}
			if (deltaComp(angle, Math.PI, delta)) { 
				transformId2AxisType.put(i, 2);
			} else if (deltaComp(angle, 2.0*Math.PI/3.0, delta)) {
				transformId2AxisType.put(i, 3);
			} else if (deltaComp(angle, Math.PI/2.0, delta)) {
				transformId2AxisType.put(i, 4);
			} else if (deltaComp(angle, Math.PI/3.0, delta)) {
				transformId2AxisType.put(i, 6);
			}
		}	
	}
	
	/**
	 * Gets all axes (crystal basis) corresponding to the rotations of each of the transformations of 
	 * this space group except for identity.
	 * @return
	 */
	public List<Vector3d> getRotAxes() {
		if (axes!=null) return axes;
		calcRotAxesAndAngles();
		return axes;
	}
	
	/**
	 * Gets all angles of rotation corresponding to each of the transformations of this space group
	 * except for identity.
	 * @return
	 */
	public List<Double> getRotAngles() {
		if (angles!=null) return angles;
		calcRotAxesAndAngles();
		return angles;
	}
	
	/**
	 * Returns true if both given transform ids belong to the same crystallographic axis (a, b or c)
	 * For two non-rotation transformations (i.e. identity operators) it returns true
	 * @param tId1
	 * @param tId2
	 * @return
	 */
	public boolean areInSameAxis(int tId1, int tId2) {
		if (tId1==tId2) return true;
		
		if (operatorsPerAxis== null) calcRotAxesAndAngles();
		
		if (isRotationIdentity(tId1) && isRotationIdentity(tId2)) return true;
		
		for (Vector3d axis:operatorsPerAxis.keySet()) {
			ArrayList<Integer> ops = operatorsPerAxis.get(axis);
			if (ops.contains(tId1) && ops.contains(tId2)) {
				return true;
			}
		}
		
		return false;
	}
	
	/**
	 * Given a transformId returns the type of axis of rotation: 2, 3, 4 or 6 -fold
	 * (or 1-fold for identity transformations)
	 * @param transformId
	 * @return
	 */
	public int getAxisFoldType(int transformId) {
		if (transformId2AxisType== null) calcRotAxesAndAngles();
		return transformId2AxisType.get(transformId);
	}
	
	/**
	 * Gets a transformation by index expressed in crystal axes basis.
	 * Index 0 corresponds always to the identity transformation.
	 * Beware the returned Matrix4d is not a copy but it stays linked 
	 * to the one stored in this SpaceGroup object
	 * @param i
	 * @return
	 */
	public Matrix4d getTransformation(int i) {
		return transformations.get(i);
	}
	
	/**
	 * Gets a transformation algebraic string given its index.
	 * Index 0 corresponds always to the identity transformation. 
	 * @param i
	 * @return
	 */
	public String getTransfAlgebraic(int i) {
		return transfAlgebraic.get(i);
	}
	
	public boolean equals(Object o) {
		if (! (o instanceof SpaceGroup)) {
			return false;
		}
		SpaceGroup other = (SpaceGroup) o;
		if (other.getId()==this.getId()) {
			return true;
		}
		return false;
	}
	
	/**
	 * Gets the number of symmetry operators corresponding to this SpaceGroup (counting 
	 * the identity operator)
	 * @return
	 */
	public int getNumOperators() {
		return this.transformations.size();
	}
	
	public static String getAlgebraicFromMatrix(Matrix4d m) {
		String x = formatAlg(m.m00,m.m01,m.m02,m.m03);
		String y = formatAlg(m.m10,m.m11,m.m12,m.m13);
		String z = formatAlg(m.m20,m.m21,m.m22,m.m23);
		String alg = x+","+y+","+z;
		return alg;
	}
	
	private static String formatAlg(double xcoef, double ycoef, double zcoef, double trans) {
		boolean[] leading = {false,false,false};
		if (xcoef!=0) {
			leading[0] = true;
		} else if (ycoef!=0) {
			leading[1] = true;
		} else if (zcoef!=0) {
			leading[2] = true;
		}
		String x = deltaComp(xcoef,0,delta)?"":formatCoef(xcoef,leading[0])+"X";
		String y = deltaComp(ycoef,0,delta)?"":formatCoef(ycoef,leading[1])+"Y";
		String z = deltaComp(zcoef,0,delta)?"":formatCoef(zcoef, leading[2])+"Z";
		String t = deltaComp(trans,0,delta)?"":formatTransCoef(trans);
		return x+y+z+t;
				
	}
	
	private static String formatCoef(double c, boolean leading) {
		if (leading) {
			return (deltaComp(Math.abs(c),1,delta)?(c>0?"":"-"):String.format("%4.2f",c));
		} else {
			return (deltaComp(Math.abs(c),1,delta)?(c>0?"+":"-"):String.format("%+4.2f",c));
		}
	}
	
	private static String formatTransCoef(double c) {
		if (Math.abs((Math.rint(c)-c))<delta) { // this is an integer
			return String.format("%+d",(int)Math.rint(c));
		} else { // it is a fraction
			int num,den;
			int floor = (int)Math.floor(c);
			double decPart = c - floor;
			if (deltaComp(decPart,0.3333333,delta)) {
				num=1;den=3;
			} else if (deltaComp(decPart,0.6666667,delta)) {
				num=2;den=3;
			} else if (deltaComp(decPart,0.2500000,delta)) {
				num=1;den=4;
			} else if (deltaComp(decPart,0.5000000,delta)) {
				num=1;den=2;
			} else if (deltaComp(decPart,0.7500000,delta)) {
				num=3;den=4;
			} else if (deltaComp(decPart,0.1666667,delta)) {
				num=1;den=6;
			} else if (deltaComp(decPart,0.8333333,delta)) {
				num=5;den=6;
			} else {
				num=0;den=0; // this in an error
			}
			num = floor*den+num;
			return String.format("%+d/%d", num,den);
			//return String.format("%+4.2f",c);
		}
	}
	
	private static boolean deltaComp(double d1, double d2, double delta) {
		return Math.abs(d1-d2)<delta;
	}
	
	public BravaisLattice getBravLattice() {
		return bravLattice;
	}
	
	public boolean isEnantiomorphic() {
		Matcher m = nonEnantPat.matcher(shortSymbol);
		if (m.find()) {
			return false;
		}
		return true;
	}
	
	/**
	 * Returns true if the transform identified by given id is an identity rotation
	 * It happens always for transformId==0 or for face centered space groups (e.g. C 1 2 1)
	 * that have operators that are simply translations in cell unit without rotation 
	 * @param transformId
	 * @return
	 */
	public boolean isRotationIdentity(int transformId) {
		if (transformId==0) return true;
		Vector3d axis = getRotAxes().get(transformId-1);
		if (axis.epsilonEquals(new Vector3d(0,0,0), delta)) {
			return true;
		}
		return false;
	}
}

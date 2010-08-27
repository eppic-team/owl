package owl.core.structure;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;


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
public final class SpaceGroup {

	
	private int id;
	private String shortSymbol;
	private List<Matrix4d> transformations;
	private List<String> transfAlgebraic;
	
	public SpaceGroup(int id, String shortSymbol) {
		this.id = id;
		this.shortSymbol = shortSymbol;
		transformations = new ArrayList<Matrix4d>();
		transfAlgebraic = new ArrayList<String>();
	}
	
	public void addTransformation(Matrix4d transfMatrix) {
		transformations.add(transfMatrix);
		
	}
	
	public void addTransformation(String transfAlgebraic) {
		this.transfAlgebraic.add(transfAlgebraic);
		this.transformations.add(convertAlgebraicToMatrix(transfAlgebraic));
	}
	
	private Matrix4d convertAlgebraicToMatrix(String transfAlgebraic) {
		String[] parts = transfAlgebraic.split(",");
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
	
	private double[] convertAlgebraicStrToCoefficients(String algString) {
		double[] coefficients = new double[4];
		Pattern coordPat = Pattern.compile("([+-])?([XYZ])");
		Matcher m = coordPat.matcher(algString);
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
		Pattern transCoefPat = Pattern.compile(".*([-+]?\\d+)/(\\d+).*");
		m = transCoefPat.matcher(algString);
		if (m.matches()) {
			double num = Double.parseDouble(m.group(1));
			double den = 1;
			if (m.group(2)!=null) {
				den = Double.parseDouble(m.group(2));
			}
			coefficients[3] = num/den;
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
	
	/**
	 * Gets a transformation by index expressed in crystal axes basis.
	 * Index 0 corresponds always to the identity transformation.
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
		double delta=0.0000001;
		return ((deltaComp(xcoef,0,delta)?"":formatCoef(xcoef)+"X")+
				(deltaComp(ycoef,0,delta)?"":formatCoef(ycoef)+"Y")+
				(deltaComp(zcoef,0,delta)?"":formatCoef(zcoef)+"Z")+
				(deltaComp(trans,0,delta)?"":String.format("%+3.1f", trans)));
	}
	
	private static String formatCoef(double c) {
		double delta=0.0000001;
		return (deltaComp(Math.abs(c),1,delta)?(c>0?"+":"-"):String.format("%3.1f",c));
	}
	
	private static boolean deltaComp(double d1, double d2, double delta) {
		return Math.abs(d1-d2)<delta;
	}
}

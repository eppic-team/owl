package owl.core.structure;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.vecmath.Matrix4d;

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
public class SpaceGroup {

	
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
		
		return new Matrix4d(xCoef[0],xCoef[1],xCoef[2],xCoef[3],yCoef[0],yCoef[1],yCoef[2],yCoef[3],zCoef[0],zCoef[1],zCoef[2],zCoef[3],0,0,0,1);
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
	
	public List<Matrix4d> getTransformations() {
		return transformations;
	}
	
	public Matrix4d getTransformation(int i) {
		return transformations.get(i);
	}
	
	public String getTransfAlgebraic(int i) {
		return transfAlgebraic.get(i);
	}
	
}

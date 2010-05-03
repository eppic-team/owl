package owl.core.runners.tinker;

public class TinkerConstraint {
	public static double DEFAULT_FORCE_CONSTANT = 1.0;
	public static enum CONSTRAINT{DISTANCE, GMBPTHETA,GMBPPHI}
	
	private int i;
	private int j;
	private double min;
	private double max;
	private CONSTRAINT type;
	private double forceConstant;
	
	public TinkerConstraint(int i, int j, double min, double max, double forceConstant, CONSTRAINT type) {
		this.type = type;
		this.i = i;
		this.j = j;
		this.min = min;
		this.max = max;
		this.forceConstant = forceConstant;
	}
	
	public double getForceConstant() {
		return this.forceConstant;
	}
	
	public CONSTRAINT getType() {
		return type;
	}
	
	public int getI() {
		return i;
	}
	
	public int getJ() {
		return j;
	}
	
	public double getMin() {
		return min;
	}
	
	public double getMax() {
		return max;
	}
	

}

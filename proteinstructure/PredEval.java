package proteinstructure;


public class PredEval {
	
	public int given;
	public int predicted;
	public int original;
	public int cmtotal;

	public int TruePos;
	public int FalsePos;
	public int TrueNeg;
	public int FalseNeg;
	
	public double sensitivity;
	public double specificity;
	public double precision;

	
	public PredEval(int TruePos, int FalsePos, int TrueNeg, int FalseNeg, int given, int predicted, int original, int cmtotal){
		this.TruePos=TruePos;
		this.FalsePos=FalsePos;
		this.TrueNeg=TrueNeg;
		this.FalseNeg=FalseNeg;
		this.given=given;
		this.predicted=predicted;
		this.original=original;
		this.cmtotal=cmtotal;
		this.sensitivity = (double) TruePos / (TruePos + FalseNeg);
		this.specificity = (double) TrueNeg / (TrueNeg + FalsePos);
		this.precision = (double) TruePos / (TruePos + FalsePos);
	}
	
	public void print(){
		System.out.println("Size of full contact map: "+cmtotal);
		System.out.println("Number original contacts: "+original);
		System.out.println("Number of given contacts: "+given);
		System.out.println("Number of predicted contacts: "+predicted);
		System.out.println();
		System.out.println("True Positives: "+TruePos);
		System.out.println("True Negatives: "+TrueNeg);
		System.out.println("False Positives: "+FalsePos);
		System.out.println("False Negatives: "+FalseNeg);
		System.out.println();
		System.out.println(String.format("Sensitivity: %4.3f",sensitivity));
		System.out.println(String.format("Specificity: %4.3f", specificity));
		System.out.println(String.format("Precision(=accuracy=PPV=TP/TP+FP): %4.3f", precision));
	}
	
}

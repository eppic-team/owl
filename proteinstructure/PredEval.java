package proteinstructure;


public class PredEval {
	
	// user customizable fields
	public String title;
	public int rank;
	
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
	public double accuracy;
	public double coverage;

	
	public PredEval(int TruePos, int FalsePos, int TrueNeg, int FalseNeg, int given, int predicted, int original, int cmtotal){
		this.title="";
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
		this.accuracy = (double) TruePos / (TruePos + FalsePos);
		this.coverage = (double) TruePos / original;
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
		System.out.println(String.format("Accuracy(=precision=PPV=TP/TP+FP): %4.3f", accuracy));
		System.out.println(String.format("Coverage(=TP/native): %4.3f", coverage));		
	}
	
	/**
	 * Print the headers for the rows written by printRow().
	 */
	public static void printHeaders() {
		System.out.print("Title\t");
		System.out.print("orig\t");
		System.out.print("pred\t");
		System.out.print("TP\t");
		System.out.print("TN\t");
		System.out.print("FP\t");
		System.out.print("FN\t");
		System.out.print("Sens\t");
		System.out.print("Spec\t");
		System.out.print("Acc\t");
		System.out.print("Cov\t");
		System.out.println();
	}
	
	/**
	 * Print selected fields of this object in one line with tab separated columns. 
	 */
	public void printRow() {
		System.out.printf("%s\t", title);
		System.out.printf("%d\t", original);
		System.out.printf("%d\t", predicted);
		System.out.printf("%d\t", TruePos);
		System.out.printf("%d\t", TrueNeg);
		System.out.printf("%d\t", FalsePos);
		System.out.printf("%d\t", FalseNeg);
		System.out.printf("%4.2f\t", sensitivity);
		System.out.printf("%4.2f\t", specificity);
		System.out.printf("%4.2f\t", accuracy);
		System.out.printf("%4.2f\n", coverage);
	}
	
}

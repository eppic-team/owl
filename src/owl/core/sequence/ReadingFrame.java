package owl.core.sequence;

public enum ReadingFrame {

	ONE ( 1,false),TWO ( 2,false),THREE ( 3,false),		// forward strand 
	RONE(-1, true),RTWO(-2, true),RTHREE(-3, true);	    // reverse strand
	
	private int number;
	private boolean reverse;
	
	private ReadingFrame(int number, boolean reverse) {
		this.number = number;
		this.reverse = reverse;
	}
	
	public int getNumber(){
		return this.number;
	}
	
	public boolean isReverse(){
		return this.reverse;
	}
}

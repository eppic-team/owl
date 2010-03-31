package owl.core.structure;

/**
 * Objects of this class are returned by the {@link PrositeScanner} to represent matching Prosite motifs.
 * TODO: Make this class implement the Feature interace (and merge with PrositeFeature)
 * @author stehr
 *
 */
public class PrositeHit {
	
	/*------------------------------ constants ------------------------------*/

	public static final int 	UNDEF_START = -1;
	public static final int 	UNDEF_STOP = -1;
	public static final String 	UNDEF_SIGNATURE = null;
	public static final double 	UNDEF_SCORE = -1.0;
	public static final int 	UNDEF_LEVEL = -1;
	
	/*--------------------------- member variables --------------------------*/
	
	public int start;
	public int stop;
	public String signatureAc;
	public double score;
	public int level;
	
	/*----------------------------- constructors ----------------------------*/
	
	public PrositeHit() {
		this.start = UNDEF_START;
		this.stop = UNDEF_STOP;
		this.signatureAc = UNDEF_SIGNATURE;
		this.score = UNDEF_SCORE;
		this.level = UNDEF_LEVEL;	
	}
	
	public PrositeHit(int start, int stop, String signatureAc, double score, int level) {
		this.start = start;
		this.stop = stop;
		this.signatureAc = signatureAc;
		this.score = score;
		this.level = level;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	public String toString() {
		return String.format("%s-%s %s %s %s", start==UNDEF_START?"?":start, 
											   stop==UNDEF_STOP?"?":stop,
											   signatureAc==UNDEF_SIGNATURE?"?":signatureAc,
											   score==UNDEF_SCORE?"?":String.format("%6.3f",score),
											   level==UNDEF_LEVEL?"?":level);
	}
}

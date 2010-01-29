package features;

import proteinstructure.Interval;
import proteinstructure.IntervalSet;
import proteinstructure.PrositeHit;
import proteinstructure.PrositeScanner;

/**
 * A feature based on a matching ProSite motif.
 * @author stehr
 */
public class PrositeFeature implements Feature {

	//;
	
	/*--------------------------- member variables --------------------------*/
	IntervalSet position;
	String description;
	FeatureType type;
	
	// type specific data
	String signatureAc;
	double score;
	int level;
	
	/*----------------------------- constructors ----------------------------*/
	
	public PrositeFeature(PrositeHit hit) {
		this.type = FeatureType.PROSITE;
		this.position = new IntervalSet();
		this.position.add(new Interval(hit.start, hit.stop));
		
		this.signatureAc = hit.signatureAc;
		this.score = hit.score;
		this.level = hit.level;
		
		this.description = String.format("%s %s", signatureAc, PrositeScanner.getPatternDescription(signatureAc));		
	}
	
	/**
	 * Createas a PrositeFeature using individual parameters.
	 */
	public PrositeFeature(int start, int stop, String signatureAc, double score, int level) {
		this.type = FeatureType.PROSITE;
		this.position = new IntervalSet();
		this.position.add(new Interval(start, stop));
		
		this.signatureAc = signatureAc;
		this.score = score;
		this.level = level;
		
		this.description = String.format("%s %s");
	}
	
	/*-------------------------- implemented methods ------------------------*/
	
	public String getDescription() {
		return this.description;
	}

	public IntervalSet getIntervalSet() {
		return this.position;
	}

	public FeatureType getType() {
		return this.type;
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	public String toString() {
		return this.description + " : " + Interval.createSelectionString(this.getIntervalSet().getIntegerSet());
	}

}

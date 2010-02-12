package features;

import java.util.TreeSet;

import proteinstructure.CatalyticSite;
import tools.Interval;
import tools.IntervalSet;

/**
 * A catalytic site from the Catalytic Site Atlas.
 * @author stehr
 */
public class CsaFeature implements Feature {

	/*--------------------------- member variables --------------------------*/
	IntervalSet position;
	String description;
	FeatureType type;
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Createas a CsaFeature based on the legacy class CatalyticSite.
	 */
	public CsaFeature(CatalyticSite cs) {
		this.type = FeatureType.CSA;
		
		TreeSet<Integer> posSet = new TreeSet<Integer>();
		for(int p:cs.getRes()) {
			posSet.add(p);
		}
		this.position = Interval.getIntervals(posSet);
		
		this.description = String.format("Id:%d Ev:%8s Pdb:%s", cs.getSerial(), cs.getEvidence(), cs.getLittEntryPdbCode()+cs.getLittEntryPdbChainCode());
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
		return  this.description  + " : " + Interval.createSelectionString(this.getIntervalSet().getIntegerSet());
	}

}

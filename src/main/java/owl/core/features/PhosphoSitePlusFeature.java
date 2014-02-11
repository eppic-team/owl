package owl.core.features;

import owl.core.structure.AminoAcid;
import owl.core.util.Interval;
import owl.core.util.IntervalSet;

/**
 * A feature based on a modification site annotation from the PhosphoSitePlus database.
 * @author stehr
 */
public class PhosphoSitePlusFeature implements Feature{

	/*--------------------------- member variables --------------------------*/
	IntervalSet position;
	String description;
	FeatureType type;
	
	// type specific data
	int modId;							// the accession number in PhoshoSitePlus
	ProteinModificationType modType;	// the type of modication (e.g. phosphorylation)
	AminoAcid aa;						// the base type of the modified amino acid
	
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Createas a new PhosphoSitePlus feature.
	 */
	public PhosphoSitePlusFeature(int modId, ProteinModificationType modType, int pos, AminoAcid aa) {
		this.type = FeatureType.PHOSPHOSITE;
		this.position = new IntervalSet();
		this.position.add(new Interval(pos, pos));
		this.modId = modId;
		this.modType = modType;
		this.aa = aa;
		this.description = String.format("%s%d-%s",aa.getOneLetterCode(), pos, modType.getSymbol());
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
		int pos = this.getPosition();
		return String.format("%s%d-%s",this.aa.getOneLetterCode(), pos, this.modType.getSymbol());
	}
	
	public int getModId() {
		return this.modId;
	}
	
	public int getPosition() {
		return this.position.first().beg;
	}
	
	public AminoAcid getAminoAcid() {
		return this.aa;
	}
	
	public ProteinModificationType getModType() {
		return this.modType;
	}
}

package owl.core.features;

import owl.core.util.Interval;
import owl.core.util.IntervalSet;

public class SiftsFeature implements Feature {


	private String pdbCode;
	private String pdbChainCode;	
	private IntervalSet cifPosition;
	private IntervalSet uniPosition;
	private String uniprotId;
	private FeatureType type;

	public SiftsFeature(String pdbCode, String pdbChainCode, String uniprotId, int cifBeg, int cifEnd, int uniBeg, int uniEnd) {
		this.type = FeatureType.SIFTS;
		this.cifPosition = new IntervalSet();
		this.cifPosition.add(new Interval(cifBeg, cifEnd));
		this.uniPosition = new IntervalSet();
		this.uniPosition.add(new Interval(uniBeg, uniEnd));
		this.pdbCode = pdbCode;
		this.pdbChainCode = pdbChainCode;
		this.uniprotId = uniprotId;
	}

	@Override
	public String getDescription() {
		return getUniprotId();
	}

	@Override
	public IntervalSet getIntervalSet() {
		return cifPosition;
	}

	@Override
	public FeatureType getType() {
		return type;
	}
	
	public String getUniprotId() {
		return uniprotId;
	}

	public String getPdbCode() {
		return pdbCode;
	}
	
	public String getPdbChainCode() {
		return pdbChainCode;
	}
	
	public IntervalSet getUniprotIntervalSet() {
		return uniPosition;
	}
}

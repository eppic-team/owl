package owl.core.structure.scoring;

import owl.core.sequence.Sequence;
import owl.core.structure.Pdb;
import owl.core.structure.features.SecondaryStructure;
import owl.core.structure.graphs.RIGraph;

public class DRToThree implements ResidueContactScoringFunction {

	@Override
	public String getMethodName() {
		return "Gripps^3";
	}

	@Override
	public double getOverallScore() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getScore(int i, int j) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getScoreForSelection(RIGraph subSet) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void init(Sequence sequence, RIGraph contacts,
			SecondaryStructure ss, Pdb coordinates) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean requiresCoordinates() {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public void updateData(Sequence sequence, RIGraph contacts,
			SecondaryStructure ss, Pdb coordinates) {
		// TODO Auto-generated method stub
		
	}
	

}

package owl.core.structure.scoring;
import owl.core.structure.graphs.RIGMatrix;

import java.sql.SQLException;

import owl.core.sequence.Sequence;
import owl.core.structure.PdbChain;
import owl.core.structure.features.SecondaryStructure;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.MySQLConnection;

public class DRToThree implements ResidueContactScoringFunction {

	private RIGMatrix contactMatrix, scoreMatrix;
	private double score;
	private MySQLConnection conn;
	
	
	public String getMethodName() {
		return "Gripps^3";
	}

	
	public double getOverallScore() {
		return score;
	}

	
	public double getScore(int i, int j) {
		double v  = scoreMatrix.getElement(i, j);
		if (v== -1) {
			return v;
		}
		return 1-v;
	}

	
	public double getScoreForSelection(RIGraph subSet) {
		// TODO Auto-generated method stub
		return 0;
	}

	
	public void init(Sequence sequence, RIGraph contacts, SecondaryStructure ss, PdbChain coordinates, MySQLConnection conn) {
		contactMatrix = new RIGMatrix(contacts);
		this.conn = conn;
		try {
			this.scoreSet(3);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
	}

	
	public boolean requiresCoordinates() {
		return false;
	}

	
	public void updateData(Sequence sequence, RIGraph contacts,
			SecondaryStructure ss, PdbChain coordinates) {
		contactMatrix = new RIGMatrix(contacts);
		try {
			this.scoreSet(3);
		} catch (SQLException e) {
			e.printStackTrace();
		}
		
	}
	
	private void scoreSet(int power) throws SQLException {
		// faster scoring : only need the contacts in the subset^power
		// generate rim raised to power
		RIGMatrix rip = new RIGMatrix(contactMatrix);
		for (int i=1; i<power; i++) { // and multiply subset rim into it
			rip.mul(contactMatrix);  // by doing that <power> times we raise rim^power 
		} // next power 
		
		scoreMatrix = new RIGMatrix(contactMatrix);
		scoreMatrix=contactMatrix.scoreDeltaMul(conn,rip); 
		score = scoreMatrix.getSum(); // score the entire matrix 
		scoreMatrix.reScale(0, 1);
	}
}

package owl.core.structure.scoring;

import owl.core.sequence.Sequence;
import owl.core.structure.AminoAcid;
import owl.core.structure.Pdb;
import owl.core.structure.features.SecondaryStructure;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.MySQLConnection;

/**
 * Example implementation of the ResidueContactScoringFunction interface.
 * This potential tries to capture the hydrophobic effect in proteins by
 * assuming that the hydrophobic interaction propensity of two residues
 * is proportional to the product of the individual hydrophobicity values.
 * Note that hydrophobic interactions are not actually real interactions 
 * but are due to forces implied by the solvent environment.
 * This class is not synchronized, i.e. may lead to unpredictable results
 * when used cuncurrently by multiple threads.
 * TODO: Could make this a template for any residue-type based potential
 * by making initializeScoringMatrix abstract, i.e. any derived class
 * would have to implement only this function.
 * @author stehr
 */
public class HydrophobicPotential implements ResidueContactScoringFunction {

	// constants
	public static final double INVALID_SCORE = Double.NaN;
	
	// member variables
	Sequence sequence;			// sequence of the target protein
	RIGraph contacts;			// contacts of the target protein
	double[][] scoringMatrix;	// the pairwise residue-type scores
	double minScore;			// minimum score for normalization
	double maxScore;			// maximum score for normalization
	
	// private methods
	/**
	 * Normalizes raw scores to the [0;1] range using 
	 * minscore and maxscore calculated in init.
	 */
	private double normalize(double score) {
		return (score - minScore) / (maxScore - minScore);
	}
	
	/**
	 * Initializes the scoring matrix. The indices in this matrix
	 * correspond to the 20 residue types (counted from 0).
	 */
	private void initializeScoringMatrix() {
		int size = AminoAcid.values().length;
		maxScore = Double.MAX_VALUE * (-1.0);
		minScore = Double.MAX_VALUE;
		scoringMatrix = new double[size][size];
		for(AminoAcid aa:AminoAcid.values()) {
			for(AminoAcid bb:AminoAcid.values()) {
				double score = aa.getHydrophobicity() * bb.getHydrophobicity();
				scoringMatrix[aa.ordinal()][bb.ordinal()] = score;
				if(score > maxScore) maxScore = score;
				if(score < minScore) minScore = score;
			}
		}
	}
	
	/**
	 * Returns the raw (i.e. unnormalized) score for residues i and j
	 * @param i index of first residue in contact (counted from 1)
	 * @param j index of second residue in contact (counted from 1)
	 * @return the raw score of a contact between i and j or INVALID_SCORE
	 * if either i or j is an invalid index or maps to a non-standard
	 * amino acid.
	 */
	public double getRawScore(int i, int j) {
		if(i < 1 || i > sequence.getLength() ||
		   j < 1 || j > sequence.getLength()) {
			return INVALID_SCORE;
		}
		char charI = sequence.getSeq().charAt(i-1);
		char charJ = sequence.getSeq().charAt(j-1);
		if(AminoAcid.isStandardAA(charI) && AminoAcid.isStandardAA(charJ)) {
			AminoAcid resI = AminoAcid.getByOneLetterCode(charI);
			AminoAcid resJ = AminoAcid.getByOneLetterCode(charJ);
			return scoringMatrix[resI.ordinal()][resJ.ordinal()];
		}
		return INVALID_SCORE;
	}
	
	// implemented methods
	
	public String getMethodName() {
		return "Hydrophobic potential";
	}

	
	public boolean requiresCoordinates() {
		return false;	// score is purely sequence based
						// so does not require coordinates
	}
	
	
	public void init(Sequence sequence, RIGraph contacts,
			SecondaryStructure ss, Pdb coordinates, MySQLConnection conn) {
		this.sequence = sequence;
		this.contacts = contacts;
		// ignoring secondary structure and coordinates
		initializeScoringMatrix();
	}
	
	
	public void updateData(Sequence sequence, RIGraph contacts,
			SecondaryStructure ss, Pdb coordinates) {
		// if data has changed, simply update the sequence and contacts variables
		// no need to recalculate scores
		this.sequence = sequence;
		this.contacts = contacts;
	}

	
	public double getScore(int i, int j) {
		return normalize(getRawScore(i,j));
	}

	
	public double getScoreForSelection(RIGraph subSet) {
		double sumScore = 0;
		for(RIGEdge e: subSet.getEdges()) {
			// do some Jung magic to get the residue numbers
			int i = subSet.getEndpoints(e).getFirst().getResidueSerial();
			int j = subSet.getEndpoints(e).getSecond().getResidueSerial();
			sumScore += getRawScore(i,j);
		}
		// here, normalization does not make sense, or does it?
		return sumScore;
	}

	
	public double getOverallScore() {
		return getScoreForSelection(contacts);	// whole graph
	}
	
	/**
	 * Simple example for how to use this scoring function.
	 */
	public static void main(String[] args) {
		
		// Building a test sequence of all possible amino acids
		String s = "";
		for(int i = 1; i <= 20; i++) {
			s += AminoAcid.getByNumber(i).getOneLetterCode();
		}
		Sequence testSeq = new Sequence("Test", s);
		
		// initialize scoring function
		HydrophobicPotential scoringFunction = new HydrophobicPotential();
		scoringFunction.init(testSeq, new RIGraph(testSeq.getSeq()), null, null,null);
		
		// printing all pairwise scores
		for(int i = 1; i <= testSeq.getLength(); i++) {
			for(int j = 2; j <= testSeq.getLength(); j++) {
				if(j <= i) {
					System.out.print("     ");
				} else {
					double score = scoringFunction.getScore(i, j);
					System.out.printf("%4.2f ", score);
				}
			}
			System.out.println();
		}
		// Overall score should be zero because we have not defined any contacts
		System.out.println("Overall score: " + scoringFunction.getOverallScore());

	}

}

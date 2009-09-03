package proteinstructure;

import java.io.File;
import java.io.IOException;

/**
 * A class representing a combined type + count scoring.
 * @author duarte
 *
 */
public class CombinedScorer extends Scorer {
	
	private static final double TYPE_WEIGHT = 1;
	private static final double COUNT_WEIGHT = 1;
	
	private Scorer countScorer;
	private Scorer typeScorer;

	public CombinedScorer(File scMatType, File scMatCount) throws IOException, FileFormatError {
		

		this.typeScorer = Scorer.readScoreMatFromFile(scMatType);
		this.countScorer = Scorer.readScoreMatFromFile(scMatCount);
		
		
		if ((typeScorer.getScoringMethod()!=ScoringMethod.ATOMTYPE && typeScorer.getScoringMethod()!=ScoringMethod.RESTYPE) || 
				(countScorer.getScoringMethod()!=ScoringMethod.ATOMCOUNT && countScorer.getScoringMethod()!=ScoringMethod.RESCOUNT)) {
			throw new IllegalArgumentException("Given scoring matrices are not of the right type");
		}
		if ((typeScorer.getScoringMethod()==ScoringMethod.ATOMTYPE && countScorer.getScoringMethod()!=ScoringMethod.ATOMCOUNT) 
				|| (typeScorer.getScoringMethod()==ScoringMethod.RESTYPE && countScorer.getScoringMethod()!=ScoringMethod.RESCOUNT)){
			throw new IllegalArgumentException("Given scoring matrices are not of the right type");
		}
		if (typeScorer.getScoringMethod()==ScoringMethod.RESTYPE) {
			this.scoringMethod = ScoringMethod.RESCOMBINED;  
		} else if (typeScorer.getScoringMethod()==ScoringMethod.ATOMTYPE) {
			this.scoringMethod = ScoringMethod.ATOMCOMBINED;
		}
		
		if (!typeScorer.getContactType().equals(countScorer.getContactType())) {
			throw new IllegalArgumentException("Count scorer and type scorer are based on different contact types");
		}
		this.ct = typeScorer.getContactType();
		if (typeScorer.getCutoff()!=countScorer.getCutoff()) {
			throw new IllegalArgumentException("Count scorer and type scorer are based on different cutoffs");
		}
		this.cutoff = typeScorer.getCutoff();
		
		if (!typeScorer.getListFile().equals(countScorer.getListFile())) {
			throw new IllegalArgumentException("Count scorer and type scorer are based on different training set files");
		}
		this.listFile = typeScorer.getListFile();
		if (typeScorer.sizeOfTrainingSet()!=countScorer.sizeOfTrainingSet()) {
			System.err.println("Warning: count scorer and type scorer are based on training sets of different sizes " +
					"(type: "+typeScorer.sizeOfTrainingSet()+", count: "+countScorer.sizeOfTrainingSet()+"). Using max of both for CombinedScorer");
		}
		this.totalStructures = Math.max(typeScorer.sizeOfTrainingSet(),countScorer.sizeOfTrainingSet());
		
		this.minSeqSep = -1; // this shouldn't be used at all in CombinedScorer (the minSeqSep of the type or count scorers should be taken)
	}
	
	@Override
	public double scoreIt(Pdb pdb) {

		double countScore = this.countScorer.scoreIt(pdb);
		double typeScore = this.typeScorer.scoreIt(pdb);
		
		return TYPE_WEIGHT*typeScore+COUNT_WEIGHT*countScore;
	}

}

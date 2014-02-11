package owl.core.structure.scoring;

import owl.core.sequence.Sequence;
import owl.core.structure.PdbChain;
import owl.core.structure.features.SecondaryStructure;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.MySQLConnection;

/**
 * This interface defines a scoring function for residue-residue contacts in proteins.
 * The scoring can be based on sequence, known contacts, secondary structure
 * and 3D coordinates. The main functionality is to return a normalized score for a
 * given pair of residues or a given set of residue pairs. The latter can also be
 * used to get a score for the structure as a whole for applications like Decoy scoring.
 * 
 * The init function takes a sequence, a contact graph, secondary structure annotation 
 * and a PdbChain object with 3d coordinates. The sequence object may not be null and has to be
 * consistent with the sequence the RIGraph, secondary structure and PdbChain object are based on.
 * The graph object may not be null but may contain no contacts. A scoring function may be
 * explicitly based on atomic coordinates or distances. In this case, the method
 * requiresCoordinates has to return true so that the calling application knows that the PdbChain
 * object may not be null. The secondary structure object may be null, so implementations
 * may use secondary structure information if given but also need to work in the absence of it.
 * 
 * The contacts should be assumed to be C-beta atoms being not more than 8 Angstrom apart unless
 * the implementation specifically takes into account the contact type defined in the RIGraph.
 * 
 * It is encouraged to return scores normalized to values between 0 and 1. We do not make this
 * a strict requirement yet because we need to gain more experience in how for this is a too severe
 * restriction. (TODO).
 * 
 * For sets of contacts, the score can be a simple average of the pairwise scores or something
 * more complex. This makes it possible to implement multi-body scoring functions which are
 * not strictly additive, as long as they can also return scores for single contacts.
 * @author stehr
 */
public interface ResidueContactScoringFunction {

	/**
	 * Returns the name of this scoring function. The idea is that scoring functions can be
	 * used like plug-ins in applications and the name will appear in menus or drop-down lists
	 * to the user can select a scoring function.
	 * @return
	 */
	public String getMethodName();
	
	/**
	 * Returns true if this scoring function requires 3d coordinates to calculate scores.
	 * This allows the calling application to decide what scoring function(s) to include
	 * depending on whether 3d information is available or not. 
	 */
	public boolean requiresCoordinates();
	
	/**
	* Initializes the ScoringFunction object with the data that the scoring can potentially be based on.
	* It is assumed that init may take some time to initialize data and to precaculate and cache results
	* such that scores can be quickly retrieved with getScore. Before this method is called, the results
	* of all other method calls are undefined.
	*/
	public void init(Sequence sequence, RIGraph contacts, SecondaryStructure ss, PdbChain coordinates, MySQLConnection conn);

	/**
	* Notifies the ScoringFunction object that the underlying data has changed and scores have to be recalculated.
	*/
	public void updateData(Sequence sequence, RIGraph contacts, SecondaryStructure ss, PdbChain coordinates);

	/**
	* Returns the normalized score for contact (i,j). If (i,j) is not really a contact, should
	* returns the score assuming that it was a contact.
	*/
	public double getScore(int i, int j);

	/**
	* Returns a score for the whole contact map, as it would be used for decoy scoring.
	* This can be simply the sum of individual scores or something more complex.
	*/
	public double getOverallScore();

	/**
	* Returns a summed score for the given subset of contacts.
	* This allows the implementation of multi-body potentials which are not strictly additive.
	*/
	public double getScoreForSelection(RIGraph subSet);

}
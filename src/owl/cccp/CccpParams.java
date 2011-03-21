package owl.cccp;

/**
 * Stores some parameters used for consensus contact prediction. The intention is to adapt the values
 * in this class as part of an automatical optimization procedure. The contents of this class will
 * change over time as new parameters are introduces and different strategies are attempted.
 * Since this class is mainly to be used internally by ConsensusPredictor (and potentially in the
 * future by something which might be called OptimizeConsensusPrediction) members are directly accessible.
 */
public class CccpParams {
	
	// TODO: Implement the parameters below, note that the default values should be specified here.
	// For benchmarking they should be overwritten and then passed to ConsensusPredictor constructor.
	
	// parameter for contact prediction
	protected boolean	useOnlyFirstModels = false;		// filter out all models > 1
	protected boolean	restrictNumTemplates = false;	// limit the number of templates to be used (for testing)
	protected int		maxNumTemplates = 20;			// if above: use at most this number of templates (ignore rest)
	protected boolean	useConsensusThreshold = false;	// if true use consensus threshold, otherwise pick top n edges
	protected double	consensusThreshold = 0.4;		// keep edges above this consensus score (if useConsThresh=true)
	protected boolean	useConsensusNumContacts = true; // in case of top edges, choose n by consensus (t) or global parameter (f)
	protected double	contactsPerNode = 5;			// if above is false use this number of contacts per residue
	protected boolean	useQuantileEdges = false;		// if useConsensusNumContacts=true, use median or quantile
	protected double	numContactsQuantile = 0.5;		// if using quantiles, choose n such that % models have less contacts 
	protected boolean	filterByMinConsensus = false;	// ignore all models with less than the given ensemble consensus score
	protected double	minConsensus = 2.0;				// keep only models with at least this consensus score (see above)
	protected boolean	filterByBestConsensus = false;	// ignore % models with worst consensus scores
	protected double	keepBestConsensusPercentage = 0.5; // keep this % models with the highest consensus scores
	
	// parameters for 3D prediction
	protected int	 	numTinkerModels = 40;			// number of models to generate out of which one is chosen
	protected boolean	tinkerFastMode = false;			// fast tinker reconstruction without annealing refinement
	protected boolean	forceTransOmega = true;			// force omega angles into trans conformation
	protected int		transOmegaInterval = 5;			// the radius around 178 defining the trans conformation 
	protected boolean	useDsspConsensusSecondaryStructure = false;  // if true, assign SS using DSSP and use consensus for 3D models
	protected double	consensusSecondaryStructureThreshold = 0.5;	// consensus threshold for consensus ss prediction
	protected boolean	useParallelTinker = false;		// whether to run on cluster, has (hopefully) no impact on results
	
	// parameters for quality prediction
	protected double	cons2qualOffset = -31;			// offset of linear fit to map model consensus score to GDT 
	protected double	cons2qualSlope = 27;			// slope of linear fit to map model consensus score to GDT
	
	// parameters for difficulty estimation
	protected double	cons2diffOffset = 0;			// offset of linear fit to map ensemble consensus score to avg GDT	
	protected double	cons2diffSlope = 0;				// slope of linear fit to map ensemble consensus score to avg GDT
	
	// other strategies/parameters to test:
	
	// choose different parameters for easy/medium/hard targets
	// choose different parameters for SR/LR contacts
	// use sum of pairs overlap instead of ensemble consensus score (what's the relation between the two?)
	
	/*---------------------------- public methods ---------------------------*/
	/**
	 * Print the current parameters to stdout
	 */
	public void print() {
		System.out.printf("#first=%b,quantile=%4.2f,filter=%4.2f", 
				useOnlyFirstModels,numContactsQuantile,keepBestConsensusPercentage);
	}
}

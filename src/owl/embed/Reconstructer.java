package owl.embed;

import java.io.File;

import javax.vecmath.Vector3d;

import owl.core.sequence.Sequence;
import owl.core.structure.ContactType;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import Jama.Matrix;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * Reconstructs a given protein contact map using the EMBED algorithm 
 * from Crippen and Havel
 * 
 * 
 * @author duarte
 *
 */
public class Reconstructer {
	
	protected static final double BB_CA_DIST = 3.8;
	private static final String EMBEDDING_ATOM_TYPE = "CA"; // this is the atom for which the embedded coordinates
															// will be set in the PdbChain models produced by reconstruct()
															// It shouldn't be a constant but at the moment we are only supporting 
															// reconstruction of Ca contact maps  
	
	private RIGraph rig;
	
	BoundsSmoother bs;
	
	private Bound[][] initialBounds;			// bounds as they are given from the input graph (i.e. sparse)
	
	/**
	 * Constructs a new Reconstructer object from a RIGraph
	 * The RIGraph is converted into a set of distance ranges (bounds matrix) using the cutoff
	 * as upper limit and hard-spheres as lower limit.  
	 * @param graph
	 */
	public Reconstructer(RIGraph graph) {
		this.rig = graph;
		this.initialBounds = convertRIGraphToBoundsMatrix(this.rig);
	}	
	
	/**
	 * Gets a reference to the initial all pairs bounds matrix 
	 * @return
	 */
	private Bound[][] getInitialBoundsAllPairs() {
		return bs.getInitialBoundsAllPairs();
	}

	/**
	 * Maps from residue serials to indices of the bounds matrices 
	 * @param idx the matrix index
	 * @return
	 */
	public int getResserFromIdx(int idx) {
		return idx+1;
	}
	
	/**
	 * Reconstructs the contact map given in constructor returning the desired number of PdbChain
	 * models, representing a sample of the conformational space of the contact map
	 * @param numModels the desired number of models
	 * @param metrize whether metrization is to be used or not. If not used a simple sampling 
	 * of the bounds matrix is performed.
	 * @param scalingMethod either ScalingMethod.RADGYRATION or ScalingMethod.AVRG_INTER_CA_DIST
	 * @param debug if true the bounds matrices and sampled matrices are printed to stdout
	 * @return
	 */
	public PdbChain[] reconstruct(int numModels, boolean metrize, Embedder.ScalingMethod scalingMethod, boolean debug) {
		
		PdbChain[] models = new PdbChain[numModels];
		
		// BoundsSmoother makes its own copy of initialBounds so we are sure that the copy of initialBounds in this class is not modified
		// and can potentially be reused (by calling again reconstruct)
		bs = new BoundsSmoother(initialBounds);
		
		if (debug) {
			// all pairs bounds after triangle inequality
			System.out.println("Bounds for all pairs after triangle inequality:");
			BoundsSmoother.printBounds(bs.getInitialBoundsAllPairs());
		}
		
		for (int model=0;model<numModels;model++) {
			
			Matrix matrix;
			if (!metrize) {
				matrix = bs.sampleBounds();
			} else {
				matrix = bs.metrize();
			}
			
			if (debug && metrize) {
				// bounds after metrization (no need to print them if sampling because they don't change)
				System.out.println("Bounds after metrization");
				bs.printBounds();
			}
			
			if (debug) {
				// sampled matrix after sampling/metrization
				System.out.println("Sampled matrix:");
				printMatrix(matrix);
			}
			
			Embedder emb = new Embedder(matrix);
			Vector3d[] embedding = emb.embed(scalingMethod);
			models[model] = new PdbChain(new Sequence(rig.getContactType(),this.rig.getSequence()), embedding, EMBEDDING_ATOM_TYPE);
			
		}
		
		return models;
	}
	
	/*------------------------ statics  ------------------------------*/
	
	/**
	 * Convert the given RIGraph to a bounds matrix. The indices of the matrix can be mapped back to residue 
	 * serials through {@link #getResserFromIdx(int)}
	 * At the moment supports only Ca contact type 
	 * @return
	 * @throws IllegalArgumentException if contact type of given RIGraph is not a single atom contact type
	 */
	public static Bound[][] convertRIGraphToBoundsMatrix(RIGraph rig) {
		// sanity check
		if (rig.getSequence().length()!=rig.getFullLength()) {
			throw new IllegalArgumentException("Full length of RIG and length of sequence differ!");
		}
		
		int conformationSize = rig.getFullLength();
		// code cloned from ConstraintsMaker.createDistanceConstraints with some modifications
		Bound[][] bounds = new Bound[conformationSize][conformationSize];
		double cutoff = rig.getCutoff();
		String ct = rig.getContactType();
		String i_ct = ct;
		String j_ct = ct;
		if (ct.contains("/")){
			i_ct = ct.split("/")[0];
			j_ct = ct.split("/")[1];
		}
		
		if (!ContactType.isValidSingleAtomContactType(i_ct) || !ContactType.isValidSingleAtomContactType(j_ct)){
			throw new IllegalArgumentException("Contact type "+i_ct+" or "+j_ct+" is not valid for reconstruction");
		}
		
		for (RIGEdge cont:rig.getEdges()){
			Pair<RIGNode> pair = rig.getEndpoints(cont);
			String i_res = pair.getFirst().getResidueType();
			String j_res = pair.getSecond().getResidueType();

			// as dist_min we take the average of the two dist mins, if i_ct and j_ct are the same then this will be the same as dist_min for ct
			double dist_min = (ContactType.getLowerBoundDistance(i_ct,i_res,j_res)+ContactType.getLowerBoundDistance(j_ct,i_res,j_res))/2;
			// for single atom contact types getUpperBoundDistance and getLowerBoundDistance will return 0 thus for those cases dist_max = cutoff
			double dist_max = ContactType.getUpperBoundDistance(i_ct, i_res, j_res)/2+ContactType.getUpperBoundDistance(i_ct, i_res, j_res)/2+cutoff;
			
			if (pair.getSecond().getResidueSerial()>pair.getFirst().getResidueSerial()+1) { //we don't add the first diagonal, we add it later as contiguous CA constraints 
				bounds[pair.getFirst().getResidueSerial()-1][pair.getSecond().getResidueSerial()-1] = new Bound(dist_min, dist_max);
			}
		}
		
		// adding contiguous CA distance backbone restraints
		// TODO this assumes we are using Ca contact type, thus this method is not general for any contact type!
		// of course for other contact types would be good to have the CA constraints too
		addBackboneRestraints(bounds);
		return bounds;
	}	
	
	/**
	 * Adds backbone restrains to the given bounds matrix. 
	 * At the moment the restraints are the contiguous CA distances only
	 * @param bounds
	 */
	protected static void addBackboneRestraints(Bound[][] bounds) {
		for (int i=0;i<bounds.length-1;i++) {
			bounds[i][i+1]=new Bound(BB_CA_DIST,BB_CA_DIST);
		}
	}
	
	private static void printViolations (Matrix matrixEmbedded, Bound[][] bounds) {
		int count = 0;
		for (int i=0;i<matrixEmbedded.getRowDimension();i++) {
			for (int j=i+1;j<matrixEmbedded.getColumnDimension();j++) {
				if ((matrixEmbedded.get(i,j)<bounds[i][j].lower) || (matrixEmbedded.get(i,j)>bounds[i][j].upper)) {
					System.out.printf("%3d %3d %4.1f %s\n",i,j,matrixEmbedded.get(i,j),bounds[i][j].toString());
					count++;
				}
			}
		}
		int cells = (matrixEmbedded.getRowDimension()*(matrixEmbedded.getRowDimension()-1))/2;
		System.out.println("Number of violations: "+count+" out of "+cells+" cells in half matrix");
	}
	
	private static int[] getViolations(Matrix matrixEmbedded, Bound[][] bounds) {
		int upperViol = 0;
		int lowerViol = 0;
		for (int i=0;i<matrixEmbedded.getRowDimension();i++) {
			for (int j=i+1;j<matrixEmbedded.getColumnDimension();j++) {
				if (bounds[i][j]!=null) { // for a sparse bounds matrix there will be missing cells, we want to only count violations to cells with values
					if ((matrixEmbedded.get(i,j)<bounds[i][j].lower)) {
						lowerViol++;
					}				
					if ((matrixEmbedded.get(i,j)>bounds[i][j].upper)) {
						upperViol++;
					}
				}
			}
		} 
		int[] viols = {lowerViol, upperViol};
		return viols;
	}

	private static void printMatrix(Matrix matrix) {
		for (int i=0;i<matrix.getRowDimension();i++) {
			for (int j=0;j<matrix.getColumnDimension();j++) {
				System.out.printf("%4.1f ",matrix.get(i, j));
			}
			System.out.println();
		}
	}

	public static void benchmarkReconstruction(String pdbCode, String pdbChainCode,
			String ct, double cutoff,
			boolean debug, boolean writeFiles, 
			String outDir, 
			int numModels,Embedder.ScalingMethod scalingMethod, boolean metrize) 
	throws Exception{
		
		File cifFile = new File(System.getProperty("java.io.tmpdir"),pdbCode+".cif");
		cifFile.deleteOnExit();
		PdbAsymUnit.grabCifFile("/path/to/local/mmCIF/gz/all/repo", null, pdbCode, cifFile, false);
		
		PdbAsymUnit fullpdb = new PdbAsymUnit(cifFile);
		PdbChain pdb = fullpdb.getChain(pdbChainCode);
		PdbChain pdbMirror = pdb.copy(pdb.getParent());
		pdbMirror.mirror();

		RIGraph graph = pdb.getRIGraph(ct, cutoff);
		int numberContacts = graph.getEdgeCount();
		int sizeHalfMatrix = (graph.getFullLength()*(graph.getFullLength()-1))/2;
		Reconstructer rec = new Reconstructer(graph);
		PdbChain[] pdbs = rec.reconstruct(numModels, metrize, scalingMethod, debug);

		System.out.println(pdbCode+pdbChainCode);
		System.out.println("Total restraints: "+numberContacts+", total cells half matrix: "+sizeHalfMatrix);
		
		System.out.printf("%6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s\t%6s", "rmsd","rmsdm","low_restr_viols", "upp_restr_viols", "low_bounds_viols", "upp_bounds_viols", "restr_viols", "bounds_viols");
		System.out.println();
		
		int modelnum=1;
		for (PdbChain model:pdbs) {
			
			double rmsd = pdb.rmsd(model, "Ca");
			double rmsdm = pdbMirror.rmsd(model, "Ca");

			if (rmsd>rmsdm ) {
				model.mirror();
			}
			
			Matrix matrixEmbedded = model.calcDistMatrixJamaFormat("Ca");
			int[] restViols = getViolations(matrixEmbedded, rec.initialBounds);
			int[] boundsViols = getViolations(matrixEmbedded, rec.getInitialBoundsAllPairs());
			System.out.printf("%6.3f\t%6.3f\t%6d\t%6d\t%6d\t%6d\t%6d\t%6d",rmsd,rmsdm,restViols[0],restViols[1],boundsViols[0],boundsViols[1],restViols[0]+restViols[1],boundsViols[0]+boundsViols[1]);
			System.out.println();
			
			if (debug) {
				printViolations(matrixEmbedded, rec.getInitialBoundsAllPairs());
			}

			
			if (writeFiles)
				model.writeToPDBFile(new File(outDir,"embed_"+pdbCode+pdbChainCode+"_"+modelnum+".pdb"));
			
			modelnum++;
		}
		

	}

	
	/*-------------------------- main  -------------------------------*/
	
	/**
	 * To test the class
	 */
	public static void main (String[] args) throws Exception {
		boolean debug = false;
		boolean writeFiles = false;
		String outDir = "/project/StruPPi/jose/embed";
		int numModels = 20;
		Embedder.ScalingMethod scalingMethod = Embedder.ScalingMethod.AVRG_INTER_CA_DIST;
		//boolean metrize = true; // if true metrization performed, if false random sampling
		
		//String pdbCode = "1bxy";
		//String pdbChainCode = "A";
		String ct = "Ca";
		double cutoff = 8.0;
		
		//String[] pdbCodes = {"1i1b", "1agd", "1mjc", "2acy", "1sha", "1rbp", "3eca", "1ho4", "1fap"};
		//String[] pdbChainCodes = {"A", "B", "A", "A", "A", "A", "A", "A", "B"};
		String[] pdbCodes = {"1bxy"};
		String[] pdbChainCodes = {"A"};
		
		System.out.println("#### SAMPLING ");
		for (int i=0; i<pdbCodes.length; i++) {
			benchmarkReconstruction(pdbCodes[i], pdbChainCodes[i], ct, cutoff, debug, writeFiles, outDir, numModels, scalingMethod, false);
		}
		
		System.out.println("#### METRIZATION ");
		for (int i=0; i<pdbCodes.length; i++) {
			benchmarkReconstruction(pdbCodes[i], pdbChainCodes[i], ct, cutoff, debug, writeFiles, outDir, numModels, scalingMethod, true);
		}
	}
	
}

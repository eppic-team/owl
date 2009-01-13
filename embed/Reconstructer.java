package embed;

import java.util.TreeMap;

import javax.vecmath.Vector3d;

import Jama.Matrix;

import proteinstructure.AAinfo;
import proteinstructure.ModelPdb;
import proteinstructure.Pdb;
import proteinstructure.PdbasePdb;
import proteinstructure.RIGEdge;
import proteinstructure.RIGNode;
import proteinstructure.RIGraph;
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
	
	private RIGraph rig;
	private TreeMap<Integer,Integer> idx2resser;
	
	private Bound[][] initialBounds;			// bounds as they are given from the input graph (i.e. sparse)
	private Bound[][] initialBoundsAllPairs;	// bounds for all pairs of atoms after triangle inequality
	
	public Reconstructer(RIGraph graph) {
		this.rig = graph;
		this.initialBounds = convertRIGraphToBoundsMatrix();
	}
	
	
	
	/**
	 * Convert the given RIGraph to a bounds matrix. The indices of the matrix can be mapped back to residue 
	 * serials through {@link #getResserFromIdx(int)}
	 * Will only admit single atom contact type RIGraphs
	 * @return
	 * @throws IllegalArgumentException if contact type of given RIGraph is not a single atom contact type
	 */
	private Bound[][] convertRIGraphToBoundsMatrix() {
		//TODO this class works right now only for proteins with no non-observed residues, FIX IT!
		int conformationSize = rig.getObsLength();
		// code cloned from ConstraintsMaker.createDistanceConstraints with some modifications
		Bound[][] bounds = new Bound[conformationSize][conformationSize];
		TreeMap<Integer,Integer> resser2idx = new TreeMap<Integer, Integer>();
		this.idx2resser = new TreeMap<Integer, Integer>();
		int idx = 0;
		for (int resser:rig.getSerials()) {
			resser2idx.put(resser,idx);
			idx2resser.put(idx,resser);
			idx++;
		}
		double cutoff = rig.getCutoff();
		String ct = rig.getContactType();
		String i_ct = ct;
		String j_ct = ct;
		if (ct.contains("/")){
			i_ct = ct.split("/")[0];
			j_ct = ct.split("/")[1];
		}
		
		if (!AAinfo.isValidSingleAtomContactType(i_ct) || !AAinfo.isValidSingleAtomContactType(j_ct)){
			throw new IllegalArgumentException("Contact type "+i_ct+" or "+j_ct+" is not valid for reconstruction");
		}
		
		for (RIGEdge cont:rig.getEdges()){
			Pair<RIGNode> pair = rig.getEndpoints(cont);
			String i_res = pair.getFirst().getResidueType();
			String j_res = pair.getSecond().getResidueType();

			// as dist_min we take the average of the two dist mins, if i_ct and j_ct are the same then this will be the same as dist_min for ct
			double dist_min = (AAinfo.getLowerBoundDistance(i_ct,i_res,j_res)+AAinfo.getLowerBoundDistance(j_ct,i_res,j_res))/2;
			// for single atom contact types getUpperBoundDistance and getLowerBoundDistance will return 0 thus for those cases dist_max = cutoff
			double dist_max = AAinfo.getUpperBoundDistance(i_ct, i_res, j_res)/2+AAinfo.getUpperBoundDistance(i_ct, i_res, j_res)/2+cutoff;
			
			if (pair.getSecond().getResidueSerial()>pair.getFirst().getResidueSerial()+1) { //we don't add the first diagonal, we add it later as contiguous CA constraints 
				bounds[resser2idx.get(pair.getFirst().getResidueSerial())][resser2idx.get(pair.getSecond().getResidueSerial())] = new Bound(dist_min, dist_max);
			}
		}
		// adding contiguous CA distance backbone constraints
		// TODO this assumes we are using Ca contact type, thus this method is not general for any contact type!
		// of course for other contact types would be good to have the CA constraints too
		for (int i=0;i<conformationSize-1;i++) {
			bounds[i][i+1]=new Bound(BB_CA_DIST,BB_CA_DIST);
		}
		return bounds;
	}
	
	/**
	 * Maps from residue serials to indices of the matrices returned by {@link #getBoundsAllPairs()} and
	 * {@link #sampleBounds(Bound[][])}
	 * @param idx
	 * @return
	 */
	public int getResserFromIdx(int idx) {
		return idx2resser.get(idx);
	}
	
	/**
	 * Reconstructs the contact map given in constructor returning the desired number of Pdb
	 * models, representing a sample of the conformational space of the contact map
	 * @param numModels the desired number of models
	 * @param metrize whether metrization is to be used or not. If not used a simple sampling 
	 * of the bounds matrix is performed.
	 * @param scalingMethod either ScalingMethod.RADGYRATION or ScalingMethod.AVRG_INTER_CA_DIST
	 * @param debug if true the bounds matrices and sampled matrices are printed to stdout
	 * @return
	 */
	public Pdb[] reconstruct(int numModels, boolean metrize, Embedder.ScalingMethod scalingMethod, boolean debug) {
		
		Pdb[] models = new ModelPdb[numModels];
		
		// we deep copy before passing to BoundsSmoother, to make sure we keep here a copy of initialBounds that is not modified
		// and can potentially be reused (by calling again reconstruct)
		BoundsSmoother bs = new BoundsSmoother(BoundsSmoother.copyBounds(initialBounds));
		
		initialBoundsAllPairs = bs.getBoundsAllPairs();
		
		if (debug) {
			// all pairs bounds after triangle inequality
			System.out.println("Bounds for all pairs after triangle inequality:");
			BoundsSmoother.printBounds(initialBoundsAllPairs);
		}
		
		for (int model=0;model<numModels;model++) {
			
			Matrix matrix;
			if (!metrize) {
				matrix = bs.sampleBounds();
			} else {
				matrix = bs.metrize();
			}
			
			if (debug) {
				// bounds after sampling/metrization
				System.out.println("Bounds after "+((metrize)?"metrization":"sampling")+":");
				bs.printBounds();
			}
			
			if (debug) {
				// sampled matrix after sampling/metrization
				System.out.println("Sampled matrix:");
				printMatrix(matrix);
			}
			
			Embedder emb = new Embedder(matrix);
			Vector3d[] embedding = emb.embed(scalingMethod);
			models[model] = new ModelPdb(this.rig.getSequence(), embedding, "CA");
			
		}
		
		return models;
	}
	
	/*------------------------ statics  ------------------------------*/
	
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
	
	private static int getNumberViolations(Matrix matrixEmbedded, Bound[][] bounds) {
		int count = 0;
		for (int i=0;i<matrixEmbedded.getRowDimension();i++) {
			for (int j=i+1;j<matrixEmbedded.getColumnDimension();j++) {
				if ((matrixEmbedded.get(i,j)<bounds[i][j].lower) || (matrixEmbedded.get(i,j)>bounds[i][j].upper)) {
					count++;
				}
			}
		}
		return count;
	}
	
	private static void printMatrix(Matrix matrix) {
		for (int i=0;i<matrix.getRowDimension();i++) {
			for (int j=0;j<matrix.getColumnDimension();j++) {
				System.out.printf("%4.1f ",matrix.get(i, j));
			}
			System.out.println();
		}
	}

	
	/*-------------------------- main  -------------------------------*/
	
	/**
	 * To test the class
	 */
	public static void main (String[] args) throws Exception {
		boolean debug = false;
		boolean writeFiles = false;
		int numModels = 10;
		Embedder.ScalingMethod scalingMethod = Embedder.ScalingMethod.AVRG_INTER_CA_DIST;
		boolean metrize = true; // if true metrization performed, if false random sampling
		
		String pdbCode = "1bxy";
		String pdbChainCode = "A";
		String ct = "Ca";
		double cutoff = 8.0;
		
		Pdb pdb = new PdbasePdb(pdbCode);
		pdb.load(pdbChainCode);
		Pdb pdbMirror = new PdbasePdb(pdbCode);
		pdbMirror.load(pdbChainCode);
		pdbMirror.mirror();

		RIGraph graph = pdb.get_graph(ct, cutoff);
		Reconstructer rec = new Reconstructer(graph);
		Pdb[] pdbs = rec.reconstruct(numModels, metrize, scalingMethod, debug);

		
		System.out.printf("%6s\t%6s\t%6s", "rmsd","rmsdm","viols");
		System.out.println();
		
		int modelnum=1;
		for (Pdb model:pdbs) {
			
			double rmsd = pdb.rmsd(model, "Ca");
			double rmsdm = pdbMirror.rmsd(model, "Ca");

			if (rmsd>rmsdm ) {
				model.mirror();
			}
			
			Matrix matrixEmbedded = model.calculateDistMatrix("Ca");
			
			System.out.printf("%6.3f\t%6.3f\t%6d",rmsd,rmsdm,getNumberViolations(matrixEmbedded, rec.initialBoundsAllPairs));
			System.out.println();
			
			if (debug) {
				printViolations(matrixEmbedded, rec.initialBoundsAllPairs);			
			}

			
			if (writeFiles)
				model.dump2pdbfile("/project/StruPPi/jose/embed/embed_"+pdbCode+pdbChainCode+"_"+modelnum+".pdb");
			
			modelnum++;
		}
	}
}

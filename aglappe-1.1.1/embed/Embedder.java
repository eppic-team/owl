package embed;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.TreeMap;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import proteinstructure.ModelPdb;
import proteinstructure.Pdb;
import proteinstructure.PdbasePdb;

import tools.Goodies;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 * Implementation of the embedding part of the EMBED algorithm of Crippen and Havel.
 * Given a complete distance matrix with exact distances it returns the 3dimensional embedding for it.
 * 
 * Taken from "Distance Geometry: Theory, Algorithms, and Chemical Applications" (section 3.2.2) by T.F. Havel, 
 * in Encyclopedia of Computational Chemistry (Wiley, New York, 1998).
 * 
 * @author duarte
 *
 */
public class Embedder {
	
	// scaling method for scaling step: radius of gyration, or averaged inter-Calpha distances
	public enum ScalingMethod {RADGYRATION, AVRG_INTER_CA_DIST};
	
	private Matrix sqDists;
	private double[] masses;
	private double[] weights;
	private int n;
	private double avrgMass;
	

	/**
	 * Constructs a new Embedder object. Get the embedding by calling {@link #embed()}
	 * The weights and masses will all be 1
	 * @param squareDistances
	 */
	public Embedder(Matrix squareDistances) {
		this(squareDistances,createTrivialVector(1.0, squareDistances.getRowDimension()),createTrivialVector(1.0, squareDistances.getRowDimension()));
	}
	
	/**
	 * Constructs a new Embedder object. Get the embedding by calling {@link #embed()}
	 * @param squareDistances complete matrix of squared distances
	 * @param atomMasses the atom masses with indices corresponding to those of the matrix
	 * @param weights the atom weights with indices corresponding to those of the matrix
	 */
	public Embedder(Matrix squareDistances, double[] atomMasses, double[] weights) {
		this.sqDists = squareDistances;
		if (sqDists.getRowDimension()!=sqDists.getColumnDimension()) {
			throw new IllegalArgumentException("Distance matrix must be square!");
		}
		this.n = this.sqDists.getRowDimension();
		this.masses = atomMasses;
		if (n!=masses.length) {
			throw new IllegalArgumentException("Array of atom masses doesn't have same dimensions as matrix of square distances");
		}
		this.weights = weights;
		if (n!=weights.length) {
			throw new IllegalArgumentException("Array of weights doesn't have same dimensions as matrix of square distances");
		}
		this.avrgMass = getAvrgMass();
	}
	
	/**
	 * Returns the 3-dimensional embedding for the data given in constructor. 
	 * The vector will be indexed as the original squareDistances Matrix
	 * @param scalingMethod 
	 * @return
	 */
	public Vector3d[] embed(ScalingMethod scalingMethod) {
		Matrix A = new Matrix(n,n);
		double[] Do = calculateDoFromDij();
		for (int i=0;i<n;i++){
			for (int j=0;j<n;j++){
				A.set(i, j, 0.5*(Do[i]+Do[j]-sqDists.get(i, j)));
			}
		}
		Matrix W = new Matrix(n,n);
		for (int i=0;i<n;i++) {
			W.set(i, i, weights[i]);
		}
		Matrix B = W.times(A.times(W));
		EigenvalueDecomposition eig = B.eig(); 
		Matrix eigValMatrix = eig.getD();
		int[] biggestEigValIndices = getIndexOf3BiggestEigVals(eigValMatrix);
		Matrix eigVecMatrix = eig.getV();
		double[] Ycol0 = scale(getColumn(eigVecMatrix,biggestEigValIndices[0]), 
				Math.sqrt(eigValMatrix.get(biggestEigValIndices[0], biggestEigValIndices[0]))); 
		double[] Ycol1 = scale(getColumn(eigVecMatrix,biggestEigValIndices[1]), 
				Math.sqrt(eigValMatrix.get(biggestEigValIndices[1], biggestEigValIndices[1])));
		double[] Ycol2 = scale(getColumn(eigVecMatrix,biggestEigValIndices[2]), 
				Math.sqrt(eigValMatrix.get(biggestEigValIndices[2], biggestEigValIndices[2])));
		Matrix Y = new Matrix(n,3);
		for (int i=0;i<n;i++){
			Y.set(i, 0, Ycol0[i]);
			Y.set(i, 1, Ycol1[i]);
			Y.set(i, 2, Ycol2[i]);
		}
		Matrix X = W.inverse().times(Y);
		Vector3d[] embedding = new Vector3d[n];
		for (int i=0;i<n;i++){
			embedding[i] = new Vector3d(X.get(i, 0), X.get(i,1), X.get(i,2));
		}
		
		scaleEmbedding(embedding, scalingMethod);
		
		return embedding;
	}
	
	/**
	 * Returns a new array result of scaling input vector by factor.
	 * @param vector
	 * @param factor
	 * @return
	 */
	private double[] scale(double[] vector, double factor){
		double[] scaledVector = new double[vector.length];
		for (int i=0;i<vector.length;i++){
			scaledVector[i]=vector[i]*factor;
		}
		return scaledVector;
	}
	
	/**
	 * Gets column j of Matrix mat in an array.
	 * @param mat
	 * @param j
	 * @return
	 */
	private double[] getColumn(Matrix mat, int j) {
		double[] col = new double[mat.getRowDimension()];
		for (int i=0;i<mat.getRowDimension();i++){
			col[i]=mat.get(i, j);
		}
		return col;
	}
	
	/**
	 * Given an eigenvalue matrix, returns the indices of the 3 biggest eigenvalues.
	 * If any of the 3 biggest eigenvalues are negative it prints a warning.
	 * @param eigValMatrix
	 * @return
	 */
	private int[] getIndexOf3BiggestEigVals(Matrix eigValMatrix) {
		HashMap<Integer, Double> eigVals = new HashMap<Integer, Double>();
		for (int i=0;i<eigValMatrix.getColumnDimension();i++) {
			eigVals.put(i, eigValMatrix.get(i, i));
		}
		LinkedHashMap<Integer, Double> eigValsSorted = Goodies.sortMapByValue(eigVals, Goodies.DESCENDING);
		ArrayList<Integer> indices = new ArrayList<Integer>(eigValsSorted.keySet());
		ArrayList<Double> values = new ArrayList<Double>(eigValsSorted.values());
		
		//System.out.println("3 biggest eigenvalues: "+values.get(0)+", "+values.get(1)+", "+values.get(2));
		// we want the three biggest non-negative eigenvalues, what to do if they are negative? no idea! 
		// anyway we print a warning if one of them is negative
		if (values.get(0)<0 || values.get(1)<0 || values.get(2)<0) {
			System.err.println("Warning one of the 3 biggest eigenvalues is negative!");
		}
		int[] biggestEigIndices = {indices.get(0), indices.get(1), indices.get(2)};
		return biggestEigIndices;
	}
	
	private double getAvrgMass() {
		double avrgMass = 0;
		for (int i=0;i<n;i++) {
			avrgMass+=masses[i];
		}
		return avrgMass/(double)n;
	}
	
	/**
	 * Calculates the vector of square distances to the centre of mass
	 * for all atoms from the square inter-distances of all atoms Dij
	 * @return
	 */
	private double[] calculateDoFromDij() {
		double[] Do = new double[n];
		for(int i=0;i<n;i++){
			Do[i] = calculateDoiFromDij(i);
		}
		return Do;
	}
	
	/**
	 * Calculates the square distance to the centre of mass for atom i
	 * from the square inter-distances of all atoms Dij 
	 * @param i
	 * @return
	 */
	private double calculateDoiFromDij (int i) {
		double summjDij = 0;
		double summjmkDjk = 0;
		for (int j=0;j<n;j++) {
			summjDij+=masses[j]*sqDists.get(i, j);
		}
		for (int j=0;j<n;j++) {
			for (int k=j+1;k<n;k++) {
				summjmkDjk+=masses[j]*masses[k]*sqDists.get(j, k);
			}
		}
		return (1/avrgMass)*summjDij-(1/(avrgMass*avrgMass))*summjmkDjk;
	}
	
	/**
	 * Modifies the given embedding conformation to scale it to match the input
	 * distance matrix.
	 * Two ways of doing this right now: scale through radius of gyration, scale 
	 * through average of the inter-Calpha distances. 
	 * @param embedding
	 * @param method 
	 */
	private void scaleEmbedding(Vector3d[] embedding, ScalingMethod method) {
		
		boolean debug = false;
		
		double scale = 0;
		// getting the scale factor
		if (method==ScalingMethod.RADGYRATION) {
			// radius of gyration 
			double radGyrInput = getRadiusGyration(sqDists);
			double radGyrOutput = getRadiusGyration(embedding);
			if (debug) {
				System.out.println("radius gyration input: "+radGyrInput);
				System.out.println("radius gyration output: "+radGyrOutput);
			}
			scale = radGyrInput/radGyrOutput;
		} 
		else if (method==ScalingMethod.AVRG_INTER_CA_DIST) {
			// average of the inter-Calpha distances
			double sumCAdists = 0;
			for (int i=0;i<n-1;i++){
				sumCAdists+=new Point3d(embedding[i]).distance(new Point3d(embedding[i+1]));
			}
			scale = Reconstructer.BB_CA_DIST/(sumCAdists/(n-1));			
		}

		// scaling the coordinates
		for (int i=0;i<n;i++) {
			embedding[i].scale(scale);
		}
		
		if (debug) {
			System.out.println("radius gyration after scaling: "+getRadiusGyration(embedding));
		}

	}
	
	/*----------------- statics -----------------*/
	
	private static double getRadiusGyration(Matrix sqDistMatrix) {
		double squareSum = 0;
		for (int i=0;i<sqDistMatrix.getRowDimension();i++){
			for (int j=i+1; j<sqDistMatrix.getColumnDimension();j++){
				squareSum+=sqDistMatrix.get(i,j);
			}
		}
		return Math.sqrt(squareSum/(2.0*sqDistMatrix.getRowDimension()));
	}
	
	private static double getRadiusGyration(Vector3d[] conformation) {
		Point3d centerMass = new Point3d(0,0,0);
		for (int i=0;i<conformation.length;i++) {
			centerMass.add(conformation[i]);
		}
		centerMass.scale(1.0/conformation.length);
		double squareSum = 0;
		for (int i=0;i<conformation.length;i++) {
			squareSum+=new Point3d(conformation[i]).distanceSquared(centerMass);
		}
		return Math.sqrt(squareSum/conformation.length);
	}	
	
	/**
	 * Creates a vector (double[]) of given size with all values set to given value.
	 * @param value
	 * @param size
	 * @return
	 */
	public static double[] createTrivialVector(double value, int size) {
		double[] vector = new double[size];
		for (int i=0;i<size;i++) {
			vector[i] = value;
		}
		return vector;
	}
	
	/*------------------- main ------------------*/
	
	/**
	 * To test the class.
	 */
	public static void main(String[] args) throws Exception {
		
		String pdbCode = "1i1b";
		String pdbChainCode = "A";
		File outPdbFile = new File("/project/StruPPi/jose/embed_"+pdbCode+pdbChainCode+".pdb");
		
		Pdb pdb = new PdbasePdb(pdbCode);
		pdb.load(pdbChainCode);
		int n = pdb.get_length();
		int ind = 0;
		TreeMap<Integer, Integer> resser2ind = new TreeMap<Integer,Integer>();
		TreeMap<Integer, Integer> ind2resser = new TreeMap<Integer, Integer>();
		for (int resser:pdb.getAllSortedResSerials()) {
			resser2ind.put(resser, ind);
			ind2resser.put(ind, resser);
			ind++;
		}
		
		// the returned matrix from calculateDistMatrix is not yet squared
		Matrix sqDistMatrix = pdb.calculateDistMatrix("Ca");
		for (int i =0;i<sqDistMatrix.getRowDimension();i++) {
			for (int j=0;j<sqDistMatrix.getColumnDimension();j++){
				sqDistMatrix.set(i, j, sqDistMatrix.get(i, j)*sqDistMatrix.get(i,j));
			}
		}

		double[] masses = createTrivialVector(1.0, n);
		double[] weights = createTrivialVector(1.0, n);
		
		System.out.println("Embedding...");
		Embedder embedder = new Embedder(sqDistMatrix, masses, weights);
		Vector3d[] embedding = embedder.embed(ScalingMethod.RADGYRATION);
		Pdb pdbEmbedded = new ModelPdb(pdb.getSequence(),embedding, "CA");
		
		double rmsd = pdb.rmsd(pdbEmbedded, "Ca");
		pdb.mirror();
		double rmsdm = pdb.rmsd(pdbEmbedded, "Ca");
		System.out.println("rmsd of embedded to original conformation: "+rmsd);
		System.out.println("rmsd of embedded to mirrored original conformation: "+rmsdm);

		// of the 2 enantiomers we take the one with lowest rmsd
		if (rmsdm<rmsd) {
			pdbEmbedded.mirror();
		}
		
		System.out.println("Writing out embedding as CA trace pdb file "+outPdbFile);
		pdbEmbedded.dump2pdbfile(outPdbFile.getAbsolutePath());

	}
}

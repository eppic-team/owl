package owl.embed;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import javax.vecmath.GMatrix;

import Jama.Matrix;

import edu.uci.ics.jung.graph.util.Pair;
/**
 * This class represents adjacency matrices of undirected graphs. To some extend, this class
 * is a copy of {@link SparseMatrix}, except no subclasses can be implemented.
 * Additionally to the functionality provided by the superclass, this class contains methods
 * for normalization. Note, that all objects of this class are square matrices.
 * @author gmueller
 *
 */
public final class AdjacencyMatrix
{
	/**the adjacency matrix*/
	private AdMat mat;
	/**
	 * Default constructor
	 */
	public AdjacencyMatrix (){
		mat = new AdMat();
	}
	/**
	 * Constructs an zero <tt>dim x dim</tt> adjacency matrix.
	 * @param dim the dimension (number of nodes)
	 */
	public AdjacencyMatrix (int dim){
		mat = new AdMat(dim,dim);
	}
	/**
	 * Constructs an adjacency matrix by converting the argument <tt>mat</tt>.
	 * @param mat the matrix (square and symmetric)
	 * @throws SparseMatrixException if the argument is neither square nor symmetric
	 */
	public AdjacencyMatrix (GMatrix mat) throws SparseMatrixException {
		GMatrix tr = new GMatrix (mat);
		tr.transpose();
		boolean test1 = mat.getNumCol()!=mat.getNumRow();
		boolean test2 = !mat.equals(tr);
		if(test1||test2){
			if(test1&&!test2)
				throw new SparseMatrixException ("Only square matrices accepted!");
			else{
				if(test2)
					throw new SparseMatrixException ("Only symmetric matrices accepted!");
				else
					throw new SparseMatrixException ("Argument neither symmmetric nor square!");
			}
		}
		this.mat = new AdMat(mat);
	}
	/**
	 * Constructs an adjacency matrix by converting the argument <tt>mat</tt>.
	 * @param mat the matrix (square and symmetric)
	 * @throws SparseMatrixException if the argument is neither square nor symmetric
	 */
	public AdjacencyMatrix (Matrix mat){
		Matrix tr = mat.transpose();
		boolean test1 = mat.getColumnDimension()!=mat.getRowDimension();
		boolean test2 = !mat.equals(tr);
		if(test1||test2){
			if(test1&&!test2)
				throw new SparseMatrixException ("Only square matrices accepted!");
			else{
				if(test2)
					throw new SparseMatrixException ("Only symmetric matrices accepted!");
				else
					throw new SparseMatrixException ("Argument neither symmmetric nor square!");
			}
		}
		this.mat = new AdMat(mat);
	}
	/**
	 * Constructs an adjacency matrix by deep copying the argument
	 * <tt>mat</tt>.
	 * @param mat the argument to copy
	 */
	public AdjacencyMatrix (AdjacencyMatrix mat){
		this.mat = new AdMat(mat.mat);
	}
	/**
	 * Constructs an adjacency matrix by convert the argument matrix.
	 * @param mat a <code>SparseMatrix</code> object
	 * @throws SparseMatrixException if the argument is neither symmetric nor square
	 */
	public AdjacencyMatrix (SparseMatrix mat) throws SparseMatrixException {
		boolean test1 = mat.getColumnDimension()!=mat.getRowDimension();
		boolean test2 = !mat.isSymmetric();
		if(test1||test2){
			if(test1&&!test2)
				throw new SparseMatrixException ("Only square matrices accepted!");
			else{
				if(test2)
					throw new SparseMatrixException ("Only symmetric matrices accepted!");
				else
					throw new SparseMatrixException ("Argument neither symmmetric nor square!");
			}
		}
		this.mat = new AdMat(mat);
	}
	/**
	 * Adds the <tt>value</tt> to this adjacency matrix. Note, in order to keep
	 * the matrix symmetric, both entries <tt>a_ij</tt> and <tt>a_ji</tt> contain
	 * the same entry.
	 * @param i the first index
	 * @param j the second index
	 * @param value the value
	 */
	public void addEntry(int i, int j, double value){
		mat.addEntry(i, j, value);
		mat.addEntry(j, j, value);
	}
	/**
	 * Adds the <tt>value</tt> to this adjacency matrix. Note, in order to keep
	 * the matrix symmetric, both entries <tt>a_ij</tt> and <tt>a_ji</tt> contain
	 * the same entry.
	 * @param index the index pair
	 * @param value the value
	 */
	public void addEntry(Pair<Integer> index, double value){
		mat.addEntry(index, value);
		mat.addEntry(new Pair<Integer> (index.getSecond(), index.getFirst()),value);
	}
	/**
	 * Adds a new entry to the matrix, if the row and column indices do not
	 * exceed the matrix dimensions. Note, that refreshing is suppressed, so this
	 * method should be used, when multiple adding is needed, followed by explicit
	 * refreshing.
	 * @param rowIndex the row index
	 * @param colIndex the column index
	 * @param value the entry
	 * @throws SparseMatrixException entry indices exceed the matrix dimensions
	 * @see addEntrySuppressRefresh(Pair,double) 
	 */
	public void addEntrySuppressRefresh(int rowIndex, int colIndex, double value) throws SparseMatrixException {
		mat.addEntrySuppressRefresh(rowIndex, colIndex, value);
		mat.addEntrySuppressRefresh(colIndex,rowIndex,value);
	}
	/**
	 * Returns true, if the entry <tt>a_ij</tt> is not zero and does not exceed the
	 * matrix dimensions. 
	 * @param i the first index
	 * @param j the second index
	 * @return true if the entry is not zero
	 */
	public boolean containsIndexPair (int i, int j){
		return mat.containsIndexPair(i,j)?true:false;
	}
	/**
	 * Returns true, if the entry <tt>a_ij</tt> is not zero and does not exceed the
	 * matrix dimensions. 
	 * @param pair the index pair <tt>(i,j)</tt>
	 * @return true if the entry is not zero
	 */
	public boolean containsIndexPair(Pair<Integer> pair){
		return mat.containsIndexPair(pair);
	}
	/**
	 * Constructs the <tt>dim x dim</tt> identity matrix. 
	 * @param dim the dimension
	 * @return the <tt>dim x dim</tt> identity matrix
	 */
	public AdjacencyMatrix constructIdentity(int dim){
		return new AdjacencyMatrix(mat.constructSquareIdentity(dim));
	}
	/**
	 * Returns the dimension of this adjacency matrix.
	 * @return the dimension
	 */
	public int getDimension (){
		return mat.getColumnDimension();
	}
	/**
	 * Returns a <code>HashSet</code> of all index pairs <tt>(i,j)</tt>,
	 * mapping on non-zero entries.
	 * @return index pairs mapping on non-zero entries
	 */
	public HashSet<Pair<Integer>> getIndexPairs (){
		return mat.getIndexPairs();
	}
	/**
	 * Returns a <code>HashSet</code> of all index pairs <tt>(i',j')</tt>,
	 * where <tt>i = i'</tt>.
	 * @param i the first index
	 * @return a set of all index pairs
	 */
	public HashSet<Pair<Integer>> getPairSetFromIndex (int i){
		return mat.getPairSetFromFirstIndex(i);
	}
	/**
	 * Returns this adjacency matrix as a <code>HashMap</code>
	 * @return a <code>HashMap</code> representing this matrix
	 */
	public HashMap<Pair<Integer>,Double> getMatrix (){
		return mat.getMatrix();
	}
	/**
	 * Returns the entry <tt>a_ij</tt>
	 * @param i the first index
	 * @param j the second index
	 * @return the entry
	 */
	public double getMatrixEntry(int i, int j){
		return mat.getMatrixEntry(i, j);
	}
	/**
	 * Returns the entry <tt>a_pair</tt>
	 * @param pair the index pair
	 * @return the entry
	 */
	public double getMatrixEntry(Pair<Integer> pair){
		return mat.getMatrixEntry(pair);
	}
	/**
	 * Returns the maximal entry of this matrix
	 * @return the maximal entry
	 */
	public double getMax (){
		return mat.getMax();
	}
	/**
	 * Returns the minimal non-zero entry of this matrix
	 * @return the minimal non-zero entry
	 */
	public double getMinNonZero (){
		return mat.getMinNonZero();
	}
	/**
	 * Returns the number of non-zero entries
	 * @return the number of non-zero entries
	 */
	public int getNumberOfEntries (){
		return mat.getNumberOfEntries();
	}
	/**
	 * Sets this adjacency matrix to identity
	 * @param dim the dimension
	 */
	public void identity(int dim){
		mat.identity(dim);
	}
	/**
	 * Returns the inner product of <tt>this</tt> an <tt>another</tt>
	 * adjacency matrix
	 * @param another another adjacency matrix 
	 * @return the inner product
	 */
	public double innerProduct (AdjacencyMatrix another){
		return mat.innerProduct(another.mat);
	}
	/**
	 * Returns the product of <tt>this</tt> and <tt>mat</tt>:
	 * <p><tt>X = this * mat</tt>
	 * @param mat the right hand factor
	 * @return the product <tt>X</tt>
	 */
	public AdjacencyMatrix multiply(AdjacencyMatrix mat){
		return new AdjacencyMatrix (mat.multiply(mat));
	}
	/**
	 * Returns the scalar multiple of this adjacency matrix
	 * @param scalar the scalar
	 * @return the scalar multiple
	 */
	public AdjacencyMatrix multiply(double scalar){
		return new AdjacencyMatrix (mat.multiply(scalar));
	}
	/**
	 * Returns the product of <tt>this</tt> and <tt>mat</tt>:
	 * <p><tt>X = this * another</tt>
	 * @param another the right hand factor
	 * @return the product <tt>X</tt>
	 */
	public AdjacencyMatrix multiply(SparseMatrix another){
		return new AdjacencyMatrix (mat.multiply(another));
	}
	/**
	 * Normalizes this adjacency matrix by the absolute value of
	 * its greatest entry.
	 * @throws ArithmeticException if this matrix is equivalent to zero matrix
	 */
	public void normalizeByGreatestValue () throws ArithmeticException {
		normalizeByValue(getMax());
	}
	/**
	 * Normalizes this adjacency matrix by the <tt>value</tt>
	 * @param value the normalizing value, must not be zero
	 * @throws ArithmeticException if <tt>value == 0</tt> returns true
	 */
	public void normalizeByValue (double value) throws ArithmeticException {
		if(value!=0d)
			normalizeByValue(value,1E-14);
		else throw new ArithmeticException ("Argument value cannot be zero!");
	}
	/**
	 * Normalizes this adjacency matrix by the <tt>value</tt>, where <tt>thresh</tt>
	 * is a marginal value. 
	 * @param value the normalizing value
	 * @param thresh the marginal value 
	 * @throws ArithmeticException if<ol>
	 * <li>value is equal to zero</li>
	 * <li>thresh is less than zero</li></ol>
	 */
	public void normalizeByValue (double value, double thresh) throws ArithmeticException {
		if(value==0||thresh<0)
			throw new ArithmeticException ("Argument value cannot be zero!");		
		HashMap<Pair<Integer>,Double> cp = mat.getMatrix();
		HashMap<Pair<Integer>,Double> map = new HashMap<Pair<Integer>,Double>();
		for (Pair<Integer> indices:cp.keySet()){
			double entry = cp.get(indices);
			double ratio = entry/value; 
			if(Math.abs(ratio)>=thresh) map.put(indices, entry/value);
		}
		mat = new AdMat();
		mat.setMatrix(map);
	}
	/**
	 * Returns <tt>this^exp</tt>, if exp is positive.
	 * @param exp the exponent
	 * @return the power <tt>this^exp</tt>
	 * @throws IllegalArgumentException if the exponent is not positive
	 */
	public AdjacencyMatrix pow (int exp) throws IllegalArgumentException {
		
		return new AdjacencyMatrix (mat.pow(exp));
	}
	/**
	 * Refreshes all index mapping
	 * @see owl.embed.SparseMatrix.refreshMatrix()
	 */
	public void refreshMatrix (){
		mat.refreshMatrix();
	}
	/**
	 * Removes the entries <tt>a_ij</tt> and <tt>a_ji</tt>
	 * @param i the first index
	 * @param j the second index
	 */
	public void removeEntry(int i, int j){
		mat.removeEntry(i, j);
		mat.removeEntry(j, i);
	}
	/**
	 * Removes the entries <tt>a_pair</tt> and <tt>a_riap</tt>,
	 * where <tt>riap</tt> is inverted index pair
	 * @param pair the index pair <tt>(i,j)</tt>
	 */
	public void removeEntry(Pair<Integer> pair){
		mat.removeEntry(pair);
		mat.removeEntry(new Pair<Integer>(pair.getSecond(),pair.getFirst()));
	}
	/**
	 * Removes the entries <tt>a_ij</tt> and <tt>a_ji</tt>, but suppresses
	 * the refreshing of the index mapping 
	 * @param i the first index
	 * @param j the second index
	 */
	public void removeEntrySuppressRefresh(int i, int j){
		mat.removeEntrySuppressRefresh(i, j);
		mat.removeEntrySuppressRefresh(j, i);
	}
	/**
	 * Removes all entries <tt>a_ij</tt> from this matrix, if <tt>threshold > a_ij</tt>
	 * returns true. Note, if <tt>threshold</tt> is less than the return value of
	 * {@link getMinNonZero()}, this matrix is equal the zero matrix
	 * @param threshold the threshold
	 */
	public void removeEntriesBelowThreshold(double threshold){
		mat.removeEntriesBelowThreshold(threshold);
	}
	/**
	 * Sets the dimension of this matrix object.
	 * @param dim the dimension
	 */
	public void setDimension (int dim){
		mat.setDimension(dim,dim);
	}
	/**
	 * Inner auxiliary class: subclass of {@link SparseMatrix}
	 * @author gmueller
	 *
	 */
	private class AdMat extends SparseMatrix 
	{
		private AdMat (){
			super();
		}
		private AdMat (int rowDim, int colDim){
			super(rowDim,colDim);
		}
		private AdMat (GMatrix mat){
			super(mat);
		}
		private AdMat (Matrix mat){
			super(mat);
		}
		private AdMat (AdMat mat){
			super(mat);
		}
		private AdMat (SparseMatrix mat){
			super(mat);
		}
		/**
		 * Removes all entries the are below the given threshold <tt>thresh</tt>.
		 * @param thresh a threshold value
		 */
		private void removeEntriesBelowThreshold (double thresh){
			HashMap<Pair<Integer>,Double> map = getMatrix();
			Iterator<Pair<Integer>> it = map.keySet().iterator();
			while(it.hasNext()){
				Pair<Integer> index = it.next();
				double entry = map.get(index);
				if(entry<thresh)
					it.remove();
			}
			setMatrix(map);
		}
		
	}
}

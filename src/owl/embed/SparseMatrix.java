package owl.embed;

import java.util.*;

import javax.vecmath.GMatrix;

import Jama.Matrix;

import edu.uci.ics.jung.graph.util.*; 
//import embed.contactmaps.SparseMatrix;

/**
 * Class to deal with sparse double precision matrices, i.e. having only a few non-zero entries. This class provides
 * essential methods as addition and subtraction, multiplication (matrix or scalar), transposition
 * and fundamental setter and getter methods for this class.
 * @author gmueller
 *TODO inner product method might be helpful
 */
public class SparseMatrix {
		
	/** the actual representation of the matrix is a HashMap, where the key is a Pair of integers
	 * and the value is a Double
	 */
	private HashMap<Pair<Integer>,Double> matrix;
	/** a HashMap mapping the first (row) index <tt>i</tt> to all index pairs <tt>(i,j)</tt>, where it is present*/ 
	private HashMap<Integer,HashSet<Pair<Integer>>> first;
	/** a HashMap mapping the second (column) index <tt>j</tt> to all index pairs <tt>(i,j)</tt>, where it is present*/
	private HashMap<Integer,HashSet<Pair<Integer>>> second;
	/** the row dimension*/
	private int rowDim;
	/** the column dimension*/ 
	private int colDim;
	/**flag, indicating whether this matrix is considered sparse, i.e. number of entries does not exceed
	 * <tt>0.2 * row_dim * col_dim</tt>*/ 
	private boolean isSparse;
	/**flag, indicating whether this matrix is symmetric*/
	private boolean isSymmetric;
	/**flag, indicating whether this matrix object was changed without refreshing index mapping*/
	private boolean wasModified;
	/**
	 * Constructs an empty <tt>SparseMatrix</tt>. Note, that invocation of instance methods
	 * causes <tt>NullPointerException</tt>, if no initial setters are invoked on the instance,
	 * constructed by this constructor.
	 */
	public SparseMatrix (){
		matrix = new HashMap<Pair<Integer>,Double>();
		first  = new HashMap<Integer,HashSet<Pair<Integer>>>();
		second = new HashMap<Integer,HashSet<Pair<Integer>>>();
	}
	public SparseMatrix (int row, int col){
		this();
		setDimension(row,col);
	}
	/**
	 * Constructs a <tt>SparseMatrix</tt> from the <tt>double</tt> matrix parameter <tt>mat</tt>.
	 * If at least one row dimension does not match the other dimensions, an <tt>SparseMatrxException</tt>
	 * is thrown.
	 * @param mat a <tt>double</tt> matrix
	 * @throws SparseMatrixException if at least one dimension does not match the others
	 */
	public SparseMatrix (double[][] mat) throws SparseMatrixException {
		setMatrix(mat);
	}
	/**
	 * Constructs a <tt>SparseMatrix</tt> from a <tt>HashMap</tt>.
	 * @param map the <tt>HashMap</tt>
	 * @param row the row dimension
	 * @param col the column dimension
	 * @throws SparseMatrixException if at least one index in the <tt>HashMap</tt>
	 * exceeds the dimensions
	 */
	public SparseMatrix (HashMap<Pair<Integer>,Double> map, int row, int col) throws SparseMatrixException {
		setDimension(row, col);
		setMatrix(map);
	}
	/**
	 * Constructs a <tt>SparseMatrix</tt> as a shallow copy of <tt>mat</tt>.
	 * @param mat an other <tt>SparseMatrix</tt>
	 */
	public SparseMatrix (SparseMatrix mat){
		setMatrix(mat);
	}
	/**
	 * Creates a <code>SparseMatrix</code> instance, by converting
	 * a <code>javax.vecmath.GMatrix</code>
	 * @param mat a GMatrix
	 */
	public SparseMatrix (GMatrix mat){
		this();
		setDimension(mat.getNumRow(),mat.getNumCol());
		for(int i = 0; i < rowDim; i++){
			for(int j = 0; j < colDim; j++){
				double entry = mat.getElement(i, j);
				if(entry!=0){
					Pair<Integer> indPair = new Pair<Integer>(new Integer(i),new Integer(j));
					matrix.put(indPair, new Double (entry));
				}
			}
		}
		refreshMatrix();
	}
	/**
	 * Creates a <code>SparseMatrix</code> instance, by converting
	 * a <code>Jama.Matrix</code>.
	 * @param mat a Matrix
	 */
	public SparseMatrix (Matrix mat){
		setMatrix(mat.getArray());
	}
	/**
	 * Initial setter: converts a <tt>double</tt> matrix to an instance of this class
	 * by iterating over all entries and storing only non-zero entries <tt>a_ij</tt> and their corresponding
	 * index pair <tt>(i,j)</tt>.
	 * @param mat a <tt>double</tt> matrix to be converted
	 * @throws SparseMatrixException if the inner dimensions do not match
	 */
	private void setMatrix (double[][] mat) throws SparseMatrixException {
		int length1  = mat.length, length2 = mat[0].length;
		matrix  = new HashMap<Pair<Integer>,Double> ();
		first   = new HashMap<Integer,HashSet<Pair<Integer>>> ();
		second  = new HashMap<Integer,HashSet<Pair<Integer>>> ();
		for(int i = 0; i < length1; i++){
			for(int j = 0; j < mat[i].length; j++){
				if(mat[i].length == length2){
					if(mat[i][j] != 0.0){
						Integer f_ind = new Integer(i), s_ind = new Integer (j);
						Pair<Integer> pair = new Pair<Integer> (f_ind,s_ind);
						matrix.put(pair, new Double (mat[i][j]));						
					}
				}
				else throw new SparseMatrixException ("\nAt least one dimension did not match. Please, " +
						"make sure, only matrices with matching dimensions are converted.");
			}
		}
		rowDim = length1;
		colDim = length2;
		refreshMatrix();
	}
	/**
	 * Initial setter: converts a <tt>HashMap</tt> instance to a <tt>SparseMatrix</tt> instance.
	 * Note, that this method refreshes all index mappings.
	 * @param map a <tt>HashMap</tt>
	 * @throws SparseMatrixException if at least one index exceeds the specified dimensions
	 */
	protected void setMatrix (HashMap<Pair<Integer>,Double> map) throws SparseMatrixException {
		matrix  = new HashMap<Pair<Integer>,Double> (map);
		if(first==null)
			first   = new HashMap<Integer,HashSet<Pair<Integer>>> ();
		if(second==null)
			second  = new HashMap<Integer,HashSet<Pair<Integer>>> ();
		refreshMatrix();
	}
	/**
	 * Initial setter: copies <tt>mat</tt> to a new <tt>SparseMatrix</tt>. Note, that this
	 * method refreshes all index mappings.
	 * @param mat a <tt>SparseMatrix</tt>
	 */
	private void setMatrix (SparseMatrix mat){
		matrix  = mat.getMatrix();
		first   = mat.getFirstIndexSet();
		second  = mat.getSecondIndexSet();
		rowDim = mat.rowDim;
		colDim = mat.colDim;
		refreshMatrix();
	}
	/**
	 * Initial setter: sets the row and column dimension fields.
	 * @param row the row dimension
	 * @param col the column dimension
	 * @throws SparseMatrixException if at least one dimension parameter is less than or equal to zero
	 */
	protected void setDimension (int row, int col) throws SparseMatrixException {
		if(row > 0 && col > 0){
			colDim = col;
			rowDim = row;
		}
		else throw new SparseMatrixException ("Any dimension field must be greater than zero!");
	}
	/**
     * Adds the <tt>mat</tt> matrix to this matrix.
     * @param subMat the sub matrix
     * @param rowBeg the row onset
     * @param colBeg the column onset
     * @throws util.SparseMatrixException if the onset exceeds the matrix dimension
     */
    public void setSubMatrix (SparseMatrix subMat, int rowBeg, int colBeg)
            throws SparseMatrixException {
        if(rowBeg>rowDim-subMat.rowDim) throw new SparseMatrixException ("Onset greater than/equal to offset!");
        if(colBeg>colDim-subMat.colDim) throw new SparseMatrixException ("Onset greater than/equal to offset!");
        for (Pair<Integer> matIndex : subMat.matrix.keySet()){
            int currRow = matIndex.getFirst(), currCol = matIndex.getSecond();
            Pair<Integer> nPair = new Pair<Integer> (new Integer (currRow+rowBeg),new Integer(currCol+colBeg));
            matrix.put(nPair,subMat.matrix.get(matIndex));
        }
    }
	/**
	 * Adds a column vector (an <tt>n x 1</tt> <code>SparseMatrix</code>) to the 
	 * column, specified by <tt>index</tt>. 
	 * @param vec the vector as an <tt>n x 1</tt> <code>SparseMatrix</code>
	 * @param index the column index
	 * @throws SparseMatrixException if either <li><tt>index</tt> exceeds the column dimension</li>
	 * <li>the row dimensions of both instances are unequal, or</li>
	 * <li>the column dimension of <tt>vec</tt> is not 1</li>
	 */
	public void addColumnVector (SparseMatrix vec, int index) throws SparseMatrixException {
		if(vec.rowDim!=rowDim)
			throw new SparseMatrixException ("Row dimension mismatch!");
		if(vec.colDim!=1)
			throw new SparseMatrixException ("A vector must have column dimension = 1!");
		if(index>=colDim)
			throw new SparseMatrixException ("Index i = "+index+" exceeds matrix column dimension!");
		for (Pair<Integer> pairs : vec.matrix.keySet()){
			Pair<Integer> nPair = new Pair<Integer>(pairs.getFirst(),new Integer(index));
			addEntrySuppressRefresh(nPair,vec.matrix.get(pairs));
		}
		refreshMatrix();
	}
	/**
	 * Adds a row vector (an <tt>1 x n</tt> <code>SparseMatrix</code>) to the 
	 * row, specified by <tt>index</tt>. 
	 * @param vec the vector as an <tt>1 x n</tt> <code>SparseMatrix</code>
	 * @param index the row index
	 * @throws SparseMatrixException if either <li><tt>index</tt> exceeds the row dimension</li>
	 * <li>the column dimensions of both instances are unequal, or</li>
	 * <li>the row dimension of <tt>vec</tt> is not 1</li>
	 */
	public void addRowVector (SparseMatrix vec, int index) throws SparseMatrixException {
		if(vec.colDim!=colDim)
			throw new SparseMatrixException ("Row dimension mismatch!");
		if(vec.rowDim!=1)
			throw new SparseMatrixException ("A vector must have column dimension = 1!");
		if(index>=rowDim)
			throw new SparseMatrixException ("Index i = "+index+" exceeds matrix column dimension!");
		for (Pair<Integer> pairs : vec.matrix.keySet()){
			Pair<Integer> nPair = new Pair<Integer>(new Integer(index),pairs.getSecond());
			addEntrySuppressRefresh(nPair,vec.matrix.get(pairs));
		}
		refreshMatrix();
	}
	/**
	 * Adds a new entry to the matrix, if the row and column indices do not
	 * exceed the matrix dimensions. Note, that refreshing is suppressed, so this
	 * method should be used, when multiple adding is needed, followed by explicit
	 * refreshing. Example:
	 * <p>
	 * <p>SparseMatrix matrix = new SparseMatrix(...);
	 * <p><tt>while(...){
	 * <p>matrix.addEntrySuppressRefresh(//some entry index and value);
	 * <p>}
	 * <p>matrix.refresh();</tt>
	 * <p>Calling some method, that depends on the row or column index maps without refreshing
	 * will cause a SparseMatrixException  
	 * @param index the index
	 * @param value the entry
	 * @throws SparseMatrixException entry indices exceed the matrix dimensions
	 */
	public void addEntrySuppressRefresh(Pair<Integer> index, double value) throws SparseMatrixException {
		int rowIndex = index.getFirst(), colIndex = index.getSecond();
		if(rowDim>rowIndex&&rowIndex>=0&&colDim>colIndex&&colIndex>=0){
			if(value!=0d){
				matrix.put(index,new Double(value));
				wasModified = true;
			}
		}
		else throw new SparseMatrixException ("Index values must be within the range of " +
				rowDim+" and "+colDim+" and equal to or greater than zero!");		
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
		if(rowDim>rowIndex&&rowIndex>=0&&colDim>colIndex&&colIndex>=0){
			Pair<Integer> pair = new Pair<Integer>(new Integer(rowIndex),new Integer(colIndex));
			//if (matrix.containsKey(pair)) System.err.println("Matrix contains entry for " +
				//	firstIndex+", "+secondIndex+"!");
			if(value!=0d){
				matrix.put(pair,new Double(value));
				wasModified = true;
			}
		}
		else throw new SparseMatrixException ("Index values must be within the range of " +
				rowDim+" and "+colDim+" and equal to or greater than zero!");		
	}
	/**
	 * Method adding a new entry to the matrix. If the matrix already contained an entry in the specified
	 * position the old value is replaced.
	 * @param firstIndex the first index (row index)
	 * @param secondIndex the second index (column index)
	 * @param value the entry value
	 * @throws SparseMatrixException if the index values exceed the matrix dimensions
	 */
	public void addEntry(int firstIndex, int secondIndex, double value) throws SparseMatrixException
	{
		if(rowDim>firstIndex&&firstIndex>=0&&colDim>secondIndex&&secondIndex>=0){
			Pair<Integer> pair = new Pair<Integer>(new Integer(firstIndex),new Integer(secondIndex));
			if(value!=0d){
				matrix.put(pair,new Double(value));
				refreshMatrix();
			}
		}
		else throw new SparseMatrixException ("Index values must be within the range of " +
				rowDim+" and "+colDim+" and equal to or greater than zero!");		
	}
	/**
	 * Method adding a new entry to the matrix. If the matrix already contained an entry in the specified
	 * position the old value is replaced.
	 * @param pair the index pair
	 * @param value the entry value
	 * @throws SparseMatrixException if the index values exceed the matrix dimensions
	 */
	public void addEntry(Pair<Integer> pair, double value) throws SparseMatrixException
	{
		int firstIndex=pair.getFirst().intValue(), secondIndex=pair.getSecond().intValue();
		if(rowDim>firstIndex&&firstIndex>=0&&colDim>secondIndex&&secondIndex>=0){
			if(value!=0d){
				matrix.put(pair,new Double(value));
				refreshMatrix();
			}
		}
		else throw new SparseMatrixException ("Index values must be within the range of " +
				rowDim+" and "+colDim+" and equal to or greater than zero!");		
	}
	/**
	 * Removes the entry associated to the <tt>index</tt> from this matrix, if the row and column indices do not
	 * exceed the matrix dimensions. Note, that refreshing is suppressed, so this
	 * method should be used, when multiple removing is needed, followed by explicit
	 * refreshing.
	 * @param index the index
	 * @throws SparseMatrixException entry indices exceed the matrix dimensions
	 */
	public void removeEntrySuppressRefresh(Pair<Integer> index){
		if(containsIndexPair(index)) {matrix.remove(index); wasModified = true;}
	}
	/**
	 * Removes the entry associated to the indices from this matrix, if they do not
	 * exceed the matrix dimensions. Note, that refreshing is suppressed, so this
	 * method should be used, when multiple removing is needed, followed by explicit
	 * refreshing.
	 * @param rowIndex the row index
	 * @param colIndex the column index
	 * @throws SparseMatrixException entry indices exceed the matrix dimensions
	 */
	public void removeEntrySuppressRefresh(int rowIndex, int colIndex){
		Pair<Integer> index = new Pair<Integer> (new Integer(rowIndex),new Integer(colIndex));
		if(containsIndexPair(index)) {matrix.remove(index); wasModified = true;}
	}
	/**
	 * Removes the entry of the matrix specified by <tt>firstIndex</tt> and <tt>secondIndex</tt>.
	 * If no such entry is present, this object is left unchanged. 
	 * @param firstIndex the row index
	 * @param secondIndex the column index
	 */
	public void removeEntry(int firstIndex, int secondIndex){
		Pair<Integer> pair = new Pair<Integer>(new Integer(firstIndex),new Integer(secondIndex));
		if(containsIndexPair(pair)) {matrix.remove(pair); refreshMatrix();}
	}
	/**
	 * Removes the entry of the matrix specified by <tt>pair</tt>.
	 * If no such entry is present, this object is left unchanged. 
	 * @param pair the index pair
	 */
	public void removeEntry(Pair<Integer> pair){
		if(containsIndexPair(pair)) {matrix.remove(pair); refreshMatrix();}
	}
	/**
	 * Refreshes all index mapping, that is, each instance of this class contains a
	 * mapping of all row and column indices onto a set of index pairs, respectively. So after any modification
	 * of a SparseMatrix object, this method is called before returning from most of the
	 * modifying methods. Though, some methods explicitly suppress refreshing, so
	 * the user has to call this method, before some index map may be used.
	 * <p>Additionally, this method sets a flag, whether or not
	 * this matrix is considered 'sparse', that is if and only if the number
	 * of entries is less than 20%.
	 * <p><b>Note</b>, that this method must be called before calling any method dealing
	 * with indexing operations, otherwise a <code>SparseMatrixException</code> will be thrown
	 */
	public final void refreshMatrix (){
		if(this==null||matrix==null) throw new NullPointerException ("Null object!");

		first.clear();
		second.clear();
		for (Pair<Integer> pairs : matrix.keySet()){
			Integer firstInd  = pairs.getFirst();
			Integer secondInd = pairs.getSecond();
			if(first.containsKey(firstInd)){
				HashSet<Pair<Integer>> helper1 = first.get(firstInd);
				helper1.add(pairs);
				first.put(firstInd,helper1);
			}
			else{
				HashSet<Pair<Integer>> helper1 = new HashSet<Pair<Integer>> (); 
				helper1.add(pairs);
				first.put(firstInd,helper1);
			}
			if(second.containsKey(secondInd)){
				HashSet<Pair<Integer>> helper1 = second.get(secondInd);
				helper1.add(pairs);
				second.put(secondInd,helper1);
			}
			else{
				HashSet<Pair<Integer>> helper1 = new HashSet<Pair<Integer>> (); 
				helper1.add(pairs);
				second.put(secondInd,helper1);
			}
		}
		double frac = ((double)matrix.size())/((double) colDim*rowDim);
		if(frac<.2) isSparse = true;
		wasModified = false;
		symmetric();			
	}
	/**
	 * Returns the eigenvector of the absolute greatest eigenvalue of this <code>SparseMatrix</code>.
	 * @return an <tt>n x 1</tt> <code>SparseMatrix</code> representing the eigenvector associated to the
	 * absolute greatest eigenvalue
	 * @throws SparseMatrixException if this matrix is not a square matrix
	 */
	public SparseMatrix computeCentralityVector () throws SparseMatrixException {
		if(colDim!=rowDim)
			throw new SparseMatrixException ("Intended only for square matrices!");
		SparseMatrix vec = new SparseMatrix (rowDim,1);
		for(int i = 0; i < rowDim; i++){
			vec.addEntrySuppressRefresh(i, 0, 1/Math.sqrt(rowDim));
		}
		vec.refreshMatrix();
		SparseMatrix old = new SparseMatrix (rowDim,1);
		double diff = (vec.subtract(old)).norm2();
		double mean = .5*(vec.add(old)).norm2();
		while(diff/mean>1e-12){
			SparseMatrix nVec = multiply(vec);
			old = new SparseMatrix (vec);
			vec = nVec.multiply(1/nVec.norm2());
			diff = (vec.subtract(old)).norm2();
			mean = .5*(vec.add(old)).norm2();
			System.out.println("||x_n - x_n-1||2 = "+diff);
		}
		return vec;
	}
	
	/**
	 * Constructs the <tt>n x n</tt> dimensional identity matrix, specified by the dimension
	 * parameter <tt>dim</tt>.
	 * @param dim the dimension
	 * @return the identity matrix
	 */
	public SparseMatrix constructSquareIdentity (int dim){
		SparseMatrix id = new SparseMatrix (dim,dim);
		for (int i = 0; i < dim; i++){
			id.matrix.put(new Pair<Integer>(new Integer(i),new Integer(i)),new Double(1d));
		}
		id.refreshMatrix();
		return id;
	}
	/**
	 * Sets this matrix object to identity.
	 */
	public void identity (){
		if(this==null) throw new NullPointerException ("Null object!");
		else{
			if(matrix==null) matrix = new HashMap<Pair<Integer>,Double> ();
			else matrix.clear();
			int dim = Math.min(rowDim, colDim);
			for(int i = 0; i < dim; i++){
				matrix.put(new Pair<Integer>(new Integer(i),new Integer(i)),new Double(1d));
			}
			refreshMatrix();
		}
	}
	/**
	 * Sets this SparseMatrix to the square (<tt>dim x dim</tt>) identity matrix
	 * @param dim the dimension
	 */
	public void identity (int dim){
		if(this==null) throw new NullPointerException ("Null object!");
		else{
			if(rowDim==0&&colDim==0) {rowDim = dim; colDim = dim;}
			if(matrix==null) matrix = new HashMap<Pair<Integer>,Double> ();
			else matrix.clear();
			for(int i = 0; i < dim; i++){
				matrix.put(new Pair<Integer>(new Integer(i),new Integer(i)),new Double(1d));
			}
			refreshMatrix();
		}
	}
	/**
	 * Sets this matrix to an <tt>row x col</tt> identity matrix.
	 * @param row the row dimension
	 * @param col the column dimension
	 */
	public void identity (int row, int col){
		setDimension(row,col);
		if(matrix==null)
			matrix = new HashMap<Pair<Integer>,Double>();
		int minDim = Math.min(row, col);
		for(int i = 0; i < minDim; i++){
			Pair<Integer> index = new Pair<Integer>(new Integer(i),new Integer(i));
			matrix.put(index, new Double(1d));
		}
		refreshMatrix();
	}
	/**
	 * Checks whether the <code>Pair</code> instance <tt>pair</tt> is present in the <code>SparseMatrix</code>.
	 * Returns true, iff the entry for the given index pair is non-zero.
	 * @param pair an index pair
	 * @return true, iff the entry for the index pair is non-zero
	 */
	public boolean containsIndexPair (Pair<Integer> pair){
		return matrix.containsKey(pair)?true:false;
	}
	/**
	 * Checks whether the index pair <tt>(i,j)</tt> is present in the <tt>SparseMatrix</tt>.
	 * Returns true, iff the entry for the given index pair is non-zero.
	 * @param i the first (row) index
	 * @param j the second (column) index
	 * @return true, iff the entry for the index pair is non-zero
	 */
	public boolean containsIndexPair (int i, int j){
		Pair<Integer> pair = new Pair<Integer> (i,j);
		return containsIndexPair (pair)?true:false;
	}
	/**
	 * Returns a <tt>HashMap</tt>, where the keys are <tt>Pair</tt>s of <tt>Integer</tt>s and the value
	 * is a <tt>Double</tt> instance. Note, only non-zero entries are stored in this <tt>HashMap</tt>. 
	 * @return a <tt>HashMap</tt> with all non-zero entries
	 */
	public HashMap<Pair<Integer>,Double> getMatrix (){
		if(matrix != null) return new HashMap<Pair<Integer>,Double>(matrix);
		else throw new NullPointerException ("\nThe matrix field was not initialized before calling this method.");
	}
	/**
	 * Returns the entry at position <tt>(i,j)</tt>. If the method <tt>containsIndexPair(Pair)</tt> returns
	 * true, the entry is return, otherwise a zero, as long as both parameters do not exceed the dimensions.
	 * @param i the row index
	 * @param j the column index
	 * @return the entry
	 * @throws SparseMatrixException if at least one index exceeds the dimensions
	 */
	public double getMatrixEntry(int i, int j) throws SparseMatrixException {
		if(i < rowDim && j < colDim&&i>=0&&j>=0){
			Pair<Integer> pair = new Pair<Integer>  (new Integer (i), new Integer (j));
			if(containsIndexPair(pair)) return matrix.get(pair).doubleValue();
			else return 0.0;
		}
		else throw new SparseMatrixException ("\nAt least one index exceeded the matrix dimensions, or was negative.");
	}
	//position<0||position>=seq.length){		if(position>=seq.length)
	/**
	 * Returns the entry specified by the <tt>indexPair</tt> parameter, if they dont exceed
	 * the matrix dimensions or are negative. If the index pair is within the range, but no
	 * entry present, zero is returned.
	 * @param indexPair the index pair <tt>(i,j)</tt>
	 * @return the matrix entry at position <tt>(i,j)</tt>, if present
	 * @throws SparseMatrixException if at least one index is negative or exceeds the matrix dimensions
	 */
	public double getMatrixEntry(Pair<Integer> indexPair) throws SparseMatrixException {
		int firstInd = indexPair.getFirst().intValue(), secondInd = indexPair.getSecond().intValue();
		if((firstInd>=rowDim||firstInd<0)||(secondInd>=colDim||secondInd<0)){
			if(firstInd>=rowDim&&secondInd>=colDim) throw new SparseMatrixException("At least one" +
					" index exceeds the matrix dimensions!");
			else throw new SparseMatrixException("At least one index is negative!");			
		}
		else{
			if(this.containsIndexPair(indexPair)) return matrix.get(indexPair).doubleValue();
			else return 0d;
		}
	}
	/**
	 * Returns the maximal entry in this matrix
	 * @return the maximal entry or zero
	 */
	public double getMax (){
		if(matrix==null) throw new NullPointerException ("Null object!");
		else{
			int counter = 0;
			double max = Double.NEGATIVE_INFINITY;
			for(Pair<Integer> pairs : matrix.keySet()){
				double currMax = getMatrixEntry(pairs);
				if(currMax>max) max = currMax;
				counter++;
			}
			boolean test = counter<colDim*rowDim;
			return max>=0&&test?max:test?0d:max;
		}
	}
	/**
	 * Returns the minimal entry in this matrix
	 * @return the minimal entry or zero
	 */
	public double getMin (){
		if(matrix==null) throw new NullPointerException ("Null object!");
		else{
			double min = Double.POSITIVE_INFINITY;
			int counter = 0;
			for(Pair<Integer> pairs : matrix.keySet()){
				double currMin = getMatrixEntry(pairs);
				if(currMin<min) min = currMin;
				counter++;
			}
			boolean test = counter<colDim*rowDim;
			return min<=0&&test?min:test?0:min;
		}
	}
	/**
	 * Returns the minimal non-zero entry in this matrix
	 * @return the minimal non-zero entry
	 */
	public double getMinNonZero (){
		if(matrix==null) throw new NullPointerException ("Null object!");
		double min = Double.POSITIVE_INFINITY;
		for(Pair<Integer> pairs : matrix.keySet()){
			double currMin = getMatrixEntry(pairs);
			if(currMin<min) min = currMin;
		}
		return min;
	}
	/**
	 * Returns the row dimension
	 * @return the row dimension
	 */
	public int getRowDimension (){
		return rowDim;
	}
	/**
	 * Returns the column dimension
	 * @return the column dimension
	 */
	public int getColumnDimension (){
		return colDim;
	}
	/**
	 * Returns the <tt>i</tt>-th column vector of this matrix as a <tt>n x 1</tt> <code>SparseMatrix</code>
	 * instance.
	 * @param i the column index
	 * @return the <tt>i</tt>-th column vector
	 * @throws SparseMatrixException if <tt>i</tt> exceeds the column dimension or
	 * <tt>i</tt> is negative 
	 */
	public SparseMatrix getColumnVector (int i) throws SparseMatrixException {
		if(i>=colDim||i<0) throw new SparseMatrixException ("Column dimension exceeded or negative index: i = "+i+"!");
		SparseMatrix col = new SparseMatrix (rowDim,1);
		Integer colInd = new Integer(i);
		if(!second.containsKey(colInd))
			return col;
		for (Pair<Integer> pairs : second.get(colInd)){
			int index = pairs.getFirst();
			col.addEntrySuppressRefresh(index,0, matrix.get(pairs));
		}
		col.refreshMatrix();
		return col;
	}
	/**
	 * Returns the <tt>i</tt>-th row vector of this matrix as a <tt>n x 1</tt> <code>SparseMatrix</code>
	 * instance.
	 * @param i the row index
	 * @return the <tt>i</tt>-th row vector
	 * @throws SparseMatrixException if <tt>i</tt> exceeds the row dimension or
	 * <tt>i</tt> is negative 
	 */
	public SparseMatrix getRowVector (int i) throws SparseMatrixException {
		if(i>=rowDim||i<0) throw new SparseMatrixException ("Row dimension exceeded or negative index: i = "+i+"!");
		SparseMatrix row = new SparseMatrix (colDim,1);
		Integer rowInd = new Integer(i);
		if(!second.containsKey(rowInd))
			return row;
		for (Pair<Integer> pairs : first.get(rowInd)){
			int index = pairs.getSecond();
			row.addEntrySuppressRefresh(0,index, matrix.get(pairs));
		}
		row.refreshMatrix();
		return row;
	}
	/**
	 * Converts an object of this class to a <tt>double</tt> matrix, including zero entries.
	 * @return a <tt>double</tt> matrix
	 */
	public double[][] getFullMatrix (){
		if(matrix == null || matrix.size() == 0){
			if(matrix==null)
				throw new NullPointerException ("The matrix field was not initialized before calling this method.");
			else
				return new double [rowDim][colDim];
		}
		else{
			double[][] fullmat = new double[rowDim][colDim];
			Set<Pair<Integer>> keys = getMatrix().keySet();
			Iterator<Pair<Integer>> it = keys.iterator();
			while(it.hasNext()){
				Pair<Integer> pair  = it.next();
				int f_ind = pair.getFirst().intValue(), s_ind = pair.getSecond().intValue();
				fullmat[f_ind][s_ind] = getMatrixEntry(f_ind, s_ind);
			}
			return fullmat;
		}
	}
	/**
	 * Returns a set of all index pairs. Note, the return object is a deep copy of the key set, so
	 * any change of the Set while an Iterator is running over the collection causes no ConcurrentModificationException!  
	 * @return a set of all index pairs
	 */
	public HashSet<Pair<Integer>> getIndexPairs (){
		if(matrix != null) return new HashSet<Pair<Integer>>(matrix.keySet());
		else throw new NullPointerException ("The matrix field was not initialized before calling this method.");
	}
	

    /**
     * Getter: returns a HashMap representation of all first (row) indices. The index
     * is used as key mapping onto a HashSet of Pairs of Integers.
     * @return a HashSet of the first index of the matrix
     * @throws SparseMatrixException matrix was modified, but not refreshed before calling this method
     */
    public HashMap<Integer,HashSet<Pair<Integer>>> getFirstIndexSet () throws SparseMatrixException {
    	if(wasModified) throw new SparseMatrixException ("Matrix was not refreshed after modification!");
    	else return new HashMap<Integer,HashSet<Pair<Integer>>> (first);
    }

    /**
     * Getter: returns a HashMap representation of all second indices. The index
     * is used as key mapping onto a HashSet of Pairs of Integers.
     * @return a HashSet of the second index of the matrix
     * @throws SparseMatrixException matrix was modified, but not refreshed before calling this method 
     */
    public HashMap<Integer,HashSet<Pair<Integer>>> getSecondIndexSet () throws SparseMatrixException {
    	if(wasModified) throw new SparseMatrixException ("Matrix was not refreshed after modification!");
        return new HashMap<Integer,HashSet<Pair<Integer>>>  (second);
    }

    /**
     * Returns a HashSet of Pairs of Integers, that is all index pairs sharing the
     * first index, or null, if no such entry exists.
     * @param i the first index in an index pair <tt>(i,j)</tt>
     * @return a HashSet of all index pairs sharing the first index or null
     * @throws SparseMatrixException matrix was modified, but not refreshed before calling this method 
     */
    public HashSet<Pair<Integer>> getPairSetFromFirstIndex (int i) throws SparseMatrixException {
    	Integer index = new Integer (i);
    	if(wasModified) throw new SparseMatrixException ("Matrix was not refreshed after modification!");    	
    	else return first.get(index)==null?null:new HashSet<Pair<Integer>>(first.get(index));
    }

    /**
     * Returns a HashSet of Pairs of Integers, that is all index pairs sharing a
     * common second index, or null, if no such entry exists.
     * @param j the second index in an index pair <tt>(i,j)</tt>
     * @return a HashSet of all index pairs sharing a common second index
     * @throws SparseMatrixException matrix was modified, but not refreshed before calling this method 
     */
    public HashSet<Pair<Integer>> getPairSetFromSecondIndex (int j) throws SparseMatrixException {
    	Integer index = new Integer (j);
    	if(wasModified) throw new SparseMatrixException ("Matrix was not refreshed after modification!");    	
    	else return second.get(index)==null?null:new HashSet<Pair<Integer>>(second.get(index));
    }
    
    /**
     * Lists all entries and returns the largest entry of this SparseMatrix.
     * @return the largest entry as <tt>double</tt>
     * @deprecated replaced by {@link getMax()}
     */
    @Deprecated
    public double getLargestEntrie (){
    	if(matrix.size() >= 0){
    		double[] array = new double[matrix.size()];
    		Set<Pair<Integer>> key = getIndexPairs();
    		Iterator<Pair<Integer>> it = key.iterator();
    		int counter = 0;
    		while(it.hasNext()){
    			Pair<Integer> pair = it.next();
    			array[counter] = matrix.get(pair).doubleValue();
    			counter++;
    		}
    		Arrays.sort(array);
    		return array[counter - 1];
    	}
    	else return 0.0;
    }
    /**
     * Returns the matrix of this <code>SparseMatrix.
     * @param rowBegin the row onset
     * @param rowEnd the row offset
     * @param colBegin the column onset
     * @param colEnd the column offset
     * @return the matrix
     * @throws util.SparseMatrixException<li>if onset is greater than or equal to offset, or</tt>
     * <li>the offset exceeds the matrix dimensions</li>
     */
    public SparseMatrix getSubMatrix (int rowBegin, int rowEnd, int colBegin,
            int colEnd) throws SparseMatrixException {
        if(rowBegin>rowEnd) throw new SparseMatrixException ("Onset greater than/equal to offset!");
        if(colBegin>colEnd) throw new SparseMatrixException ("Onset greater than/equal to offset!");
        if(rowEnd>rowDim||colEnd>colDim)
             throw new SparseMatrixException ("Offset exceeds matrix dimension!");
        int newRow = rowEnd-rowBegin+1;
        int newCol = colEnd-colBegin+1;
        SparseMatrix sub = new SparseMatrix (newRow,newCol);
        for (int i = 0; i < newRow; i++){
            for (int j = 0; j < newCol; j++){
                Pair<Integer> thisPair = new Pair<Integer> (new Integer (i+rowBegin), new Integer (j+colBegin));
                if(containsIndexPair(thisPair)){
                    Pair<Integer> subIndex = new Pair<Integer> (new Integer (i), new Integer (j));
                    sub.matrix.put(subIndex,matrix.get(thisPair));
                }
            }
        }
        sub.refreshMatrix();
        return sub;
    }
    /*------------------------Operation methods------------------------*/
    /**
     * Addition method: adds <tt>mat</tt> to <tt>this</tt> an returns
     * it as a new instance of this class.
     * @param mat one summand <tt>mat</tt> in <tt>B = this + mat</tt>
     * @return the sum <tt>B</tt> in <tt>B = this + mat</tt>
     * @throws SparseMatrixException if at least one dimension pair does not match
     */
	public SparseMatrix add (SparseMatrix mat) throws SparseMatrixException {
		boolean test1 = getRowDimension() == mat.getRowDimension();
		boolean test2 = getColumnDimension() == mat.getColumnDimension();
		if(test1 && test2){
			HashMap<Pair<Integer>,Double> map1 = getMatrix();
			HashMap<Pair<Integer>,Double> map2 = mat.getMatrix();
			SparseMatrix sum = new SparseMatrix();
			sum.setDimension(rowDim, colDim);
			Set<Pair<Integer>> key1 = map1.keySet();
			Iterator<Pair<Integer>> it1 = key1.iterator();
			while(it1.hasNext()){
				Pair<Integer> pair = it1.next();
				if(map2.containsKey(pair)){
					double sEntry = matrix.get(pair).doubleValue()+mat.matrix.get(pair).doubleValue();
					sum.matrix.put(pair,new Double(sEntry));
					map2.remove(pair);
				}
				else sum.matrix.put(pair, matrix.get(pair));
			}
			sum.matrix.putAll(map2);
			sum.refreshMatrix();
			return sum;
		}
		else throw new SparseMatrixException ("The row and column dimensions must match!");
	}
	/**
	 * Returns true, iff the difference <tt>X = this - mat</tt> equals zero.
	 * @param mat the subtrahend
	 * @return true iff the difference is equal to zero
	 */
	public boolean equals (SparseMatrix mat){
		SparseMatrix diff = subtract(mat);
		if(diff.getMatrix().size() == 0) return true;
		else return false;
	}
	/**
	 * Method to compute the inner product of two SparseMatrices:
	 * <tt>this</tt> and <tt>mat</tt>:
	 * <p><tt>sum(i,j) this(i,j) * mat(i,j)</tt>
	 * <p>Note, calling this method on <tt>this</tt> is equivalent to 
	 * computing the Frobenius-norm ({@link SparseMatrix#norm2()} of a given matrix.
	 * @param mat another matrix
	 * @return the inner product of <tt>this</tt> and <tt>mat</tt> 
	 */
	public double innerProduct (SparseMatrix mat){
		double innerProd = 0d;
		HashSet<Pair<Integer>> key1 = new HashSet<Pair<Integer>>(matrix.keySet());
		HashSet<Pair<Integer>> key2 = new HashSet<Pair<Integer>>(mat.matrix.keySet());
		key1.retainAll(key2);
		for(Pair<Integer> pair : key1){
			if(containsIndexPair(pair)&&mat.containsIndexPair(pair)) innerProd += ((matrix.get(pair)).doubleValue())*((mat.matrix.get(pair)).doubleValue());
		}
		return innerProd;
	}
	/**
	 * Indicates, whether or not this object is considered sparse,
	 * that is if and only if the number of entries is less 20%
	 * of <tt>row_dim * col_dim</tt>, this method returns true.
	 * @return returns true if this object as less than 20% occupancy
	 * rate
	 */
	public boolean isSparse (){
		return isSparse?true:false;
	}
	/**
	 * Setter: sets the symmetry flag if and only if this matrix
	 * is symmetric, that is if <tt>a_ij == a_ji</tt> returns true, for all
	 * <tt>n-1 >= i,j >= 0</tt>, where <tt>n</tt> is the dimension
	 * of this square (!) matrix. A non-square matrix will always be
	 * considered a non-symmetric matrix.
	 * @throws SparseMatrixException matrix was modified, but not refreshed before calling this method 
	 */
	private void symmetric () throws SparseMatrixException {
		if(first==null&&second==null) throw new NullPointerException ("Null object!");
		else{
			if(wasModified) throw new SparseMatrixException ("Matrix was not refreshed after modification!");
			else{
				if(rowDim==colDim){
					for(Pair<Integer> index : matrix.keySet()){
						Double entry = matrix.get(index);
						int fIndex = index.getFirst(), sIndex = index.getSecond();
						Pair<Integer> xedni = new Pair<Integer>(new Integer(sIndex),new Integer(fIndex));
						if(matrix.containsKey(xedni)){
							if(entry.equals(matrix.get(xedni))) isSymmetric = true;
							else {isSymmetric = false;break;}
						}
						else {isSymmetric = false;break;}
					}
				}
			}
		}
	}
	/**
	 * Returns true, if and only if all entries
	 * <tt>a_ij == a_ji</tt> returns true. 
	 * @return true, if this matrix is symmetric
	 */
	public boolean isSymmetric (){
		symmetric();
		return isSymmetric?true:false;
	}
	/**
	 * Returns the kernel vector <tt>v</tt> of this matrix,
	 * which satisfies:
	 * <p><tt>A v = o</tt>, where <tt>o</tt> is the zero vector
	 * @return the kernel vector as a SparseMatrix
	 * TODO should be removed
	 */
	@Deprecated
	public SparseMatrix kernelVector (){
		if(matrix==null) throw new NullPointerException ("Null object!");
		else{
			if(convert2Matrix().rank()==rowDim) return new SparseMatrix (rowDim,1);
			//if this matrix has full rank, the kernel is trivial!
			
			else{
				SparseMatrix kerVec = new SparseMatrix (colDim,1);
				SparseMatrix operator = new SparseMatrix(rowDim,colDim);
				operator.identity(rowDim,colDim);
				for(int i = 0; i < colDim; i++){
					Pair<Integer> index = new Pair<Integer>(new Integer(i),new Integer(0));
					kerVec.matrix.put(index, 1/Math.sqrt(colDim));
				}
				kerVec.refreshMatrix();
				operator = operator.add(this);
				//int pot = 2;
				double normOP = operator.norm2();
				/*while(normOP>=1){
					operator = operator.subtract(this.multiply(2d));
					normOP = operator.norm2();
					pot++;
				}*/
				//operator = operator.add(this);
				operator = operator.multiply(1/(normOP+1));
				SparseMatrix latter = new SparseMatrix (colDim,1);
				double currNorm = kerVec.norm2(), diff = 0;
				while(Math.abs(currNorm-diff)>1e-13){
					latter = new SparseMatrix (kerVec);
					kerVec = operator.multiply(kerVec);
					double norm = kerVec.norm2();
					kerVec = kerVec.multiply(1/norm);
					diff = currNorm;
					currNorm = (kerVec.subtract(latter)).norm2();
					System.out.println(currNorm);
				}
				return kerVec;
			}
		}
	}
	/**
	 * Scalar multiplication method: multiplies this instance with <tt>scalar</tt> and
	 * returns the result as a new instance.
	 * @param scalar the scalar
	 * @return the product of <tt>scalar * this</tt>
	 */
	public SparseMatrix multiply (double scalar){
		HashMap<Pair<Integer>, Double> map  = getMatrix();
		HashMap<Pair<Integer>, Double> nmap = new HashMap<Pair<Integer>, Double>();
		Set<Pair<Integer>> keys = map.keySet();
		Iterator<Pair<Integer>> it = keys.iterator();
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			double newval      = map.get(pair).doubleValue()*scalar;
			if(newval != 0.0) nmap.put(pair, new Double (newval));
		}
		return new SparseMatrix (nmap,getRowDimension(),getColumnDimension());
	}
	
	/**
	 * Multiplication method: multiplies <tt>mat</tt> to <tt>this</tt> and returns the
	 * product as a new instance.
	 * @param mat the right hand factor in <tt>X = this * mat</tt>
	 * @return the product <tt>X</tt> in <tt>X = this * mat</tt>
	 * @throws SparseMatrixException if the outer dimensions do not match
	 */
	public SparseMatrix multiply (SparseMatrix mat) throws SparseMatrixException {
		if(getColumnDimension() == mat.getRowDimension()){
			HashMap<Pair<Integer>,Double> map1 = getMatrix();
			HashMap<Pair<Integer>,Double> map2 = mat.getMatrix();
			HashMap<Pair<Integer>,Double> newmap = new HashMap<Pair<Integer>,Double>(); 
			Set<Pair<Integer>> keys = map1.keySet();
			Iterator<Pair<Integer>> it = keys.iterator();
			while(it.hasNext()){
				Pair<Integer> pair1 = it.next();
				Set<Pair<Integer>> keys2 = map2.keySet();
				Iterator<Pair<Integer>> it2 = keys2.iterator();
				while(it2.hasNext()){
					Pair<Integer> pair2 = it2.next();
					if(pair1.getSecond().intValue() == pair2.getFirst().intValue()){
						Pair<Integer> newpair = new Pair<Integer> (pair1.getFirst(),pair2.getSecond());
						double pro = map1.get(pair1).doubleValue()*map2.get(pair2).doubleValue();
						if(newmap.containsKey(newpair)){
							double sum = newmap.get(newpair).doubleValue();
							if(sum != 0.0) newmap.put(newpair,new Double (sum+pro));	
						}
						else if(pro != 0.0) newmap.put(newpair,new Double (pro));
					}
				}
			}
			return new SparseMatrix (newmap,getRowDimension(),mat.getColumnDimension());
		}
		else throw new SparseMatrixException ("The number of columns of the first matrix has to equal the number of rows of the second matrix.");
	}
	/**
	 * Returns the 1-norm of this matrix.
	 * @return the 1-norm
	 * @throws SparseMatrixException if refreshing was not performed
	 * before calling this method
	 */
	public double norm1 () throws SparseMatrixException {
		if(wasModified) throw new SparseMatrixException ("Matrix was modified without refreshing index mapping!");
		if(matrix.size()==0) return 0d;
		double maxSum = 0d;
		for (Integer colIndex : second.keySet()){
			HashSet<Pair<Integer>> indices = second.get(colIndex);
			double colSum = 0d;
			for (Pair<Integer> index: indices)
				colSum =+ Math.abs(matrix.get(index));
			if(colSum>=maxSum)
				maxSum = colSum;
		}
		return maxSum;
	}
	/**
	 * Returns the 2-norm of this sparse matrix
	 * @return the 2-norm (frobenius-norm)
	 */
	public double norm2 (){
		if(this==null) throw new NullPointerException("Null object!");
		else{
			return Math.sqrt(innerProduct(this));
		}
	}
	
	/**
	 * Returns the infinity norm of this matrix.
	 * @return the infinity norm
	 * @throws SparseMatrixException if refreshing was not performed
	 * before calling this method
	 */
	public double normInf () throws SparseMatrixException {
		if(wasModified) throw new SparseMatrixException ("Matrix was modified without refreshing index mapping!");
		if(matrix.size()==0) return 0d;
		double maxSum = 0d;
		for(Integer rowIndex : first.keySet()){
			HashSet<Pair<Integer>> indices = first.get(rowIndex);
			double rowSum = 0d;
			for (Pair<Integer> index : indices)
				rowSum =+ Math.abs(matrix.get(index));
			if(rowSum>=maxSum)
				maxSum = rowSum;
		}
		return maxSum;
	}
	/**
	 * Takes this SparseMatrix instance and and raises it to the power defined by the
	 * parameter <tt>exp</tt>. So the return object is the result <tt>X</tt> of the equation
	 * <p><tt>X = this^exp</tt>. <b>Note, that in case of a non-square matrix</b>, an intermediate
	 * (power) result will be a SparseMatrix with smaller dimension, causing SparseMatrixException, since
	 * multiplication of matrices with mismatching dimensions causes the same exception. So, square
	 * matrices are the only acceptable arguments. Additionally, the exponent <tt>exp</tt> must be an integer
	 * greater than zero, otherwise an SparseMatrixException is thrown.
	 * @param exp the power to which the SparseMatrix is raised (strictly positive)
	 * @return the result <tt>X</tt> of <tt>X = this^exp</tt>
	 * @throws SparseMatrixException if this matrix is not square
	 * @throws IllegalArgumentException if the exponent is not positive
	 */
	public SparseMatrix pow (int exp) throws SparseMatrixException, IllegalArgumentException  {
		if(getColumnDimension()==getRowDimension() && exp > 0){
			if(exp == 1)
				return new SparseMatrix(this);
			else{
				SparseMatrix cp1 = new SparseMatrix(this);
				SparseMatrix cp2 = new SparseMatrix(this);
				for(int i = 1; i < exp; i++){
					cp1 = cp1.multiply(cp2);
				}
				return cp1;
			}
		}
		else {if(exp <= 0) throw new IllegalArgumentException ("Only positive integers are accepted!");

		else throw new SparseMatrixException("Only square matrices are acceptable arguments!");}
	}
	/**
	 * Subtracts <tt>mat</tt> from <tt>this</tt> and returns the difference as a
	 * new instance. A similar result is obtained by invocation of
	 * <p><tt>SparseMatrix mat1 = new SparseMatrix(...);</tt>//first arbitrary matrix</p>
	 * <p><tt>SparseMatrix mat2 = new SparseMatrix(...);</tt>//second arbitrary matrix</p>
	 * <p><tt>SparseMatrix diff1 = mat1.add(mat2.multiply(-1.0));</tt>//first difference</p>
	 * <p><tt>SparseMatrix diff1 = mat1.subtract(mat2);</tt>//second difference</p>
	 * <p><tt>diff1.equals(diff2)</tt>//returns true</p>
	 * @param mat the subtrahend in <tt>X = this - mat</tt>
	 * @return the difference <tt>X</tt> in <tt> X = this - mat</tt>
	 * @throws SparseMatrixException if at least one dimension does not match
	 */
	public SparseMatrix subtract (SparseMatrix mat) throws SparseMatrixException {
		return (add(mat.multiply(-1.0)));
	}
	
	/**
	 * Returns the number of non-zero entries of this matrix.
	 * @return the number of non-zeros entries
	 */
	public int getNumberOfEntries (){
		if(matrix==null) throw new NullPointerException ("Null object!");
		else return matrix.size();
	}	
	
	/**
     * Returns the transposed matrix of this instance as a new SparseMatrix instance.
     * @return the transposed matrix
     */
    public SparseMatrix transpose (){
        HashMap<Pair<Integer>,Double> cp = getMatrix();
        HashMap<Pair<Integer>,Double> nm = new HashMap<Pair<Integer>,Double> ();
        Set<Pair<Integer>> key = cp.keySet();
        Iterator<Pair<Integer>> it = key.iterator();
        while(it.hasNext()){
            Pair<Integer> pair1 = it.next();
            Pair<Integer> pair2 = new Pair<Integer> (pair1.getSecond(),pair1.getFirst());
            nm.put(pair2,cp.get(pair1));
        }
        SparseMatrix tr = new SparseMatrix();
        tr.setMatrix(nm);
        tr.setDimension(getColumnDimension(), getRowDimension());
        return tr;
    }
	/**
	 * Converts this instance to a String.
	 */
	public String toString (){
		if(matrix != null){
			String str = "";
			double[][] ar = getFullMatrix();
			int length = ar.length;
			for(int i = 0; i < length; i++){
				for(int j = 0; j < ar[i].length; j++){
					if(j < ar[i].length - 1) str += ar[i][j] + "\t";
					else str += ar[i][j] + "\n";
				}
			}
			return str;
		}
		else throw new NullPointerException ("The matrix field was not initialized before calling this method.");
	}
	/*----------------------conversion methods------------------*/
	/**
	 * Converts a set of Pairs to a <code>SparseMatrix</code> 
	 */
	void convertToSparseMatrix (Set<Pair<Integer>> set, int col, int row){
		Iterator<Pair<Integer>> it = set.iterator();
		HashMap<Pair<Integer>,Double> map = new HashMap<Pair<Integer>,Double>();
		while(it.hasNext()){
			Pair<Integer> pair  = it.next();
			Pair<Integer> npair = new Pair<Integer>(pair.getSecond(),pair.getFirst());
			map.put(pair, new Double(1.0));
			map.put(npair, new Double(1.0));
		}
		setDimension(row,col);
		setMatrix(map);
	}
	/*----------------------converting methods-------------------------------*/
	/**
	 * Returns a this matrix as a <code>GMatrix</code> object
	 * @return a GMatrix object
	 */
	public GMatrix convert2GMatrix (){
		GMatrix gMatr = new GMatrix(rowDim,colDim);
		gMatr.setZero();
		if(matrix==null){
			return gMatr;
		}
		else{
			for(Pair<Integer> pair: matrix.keySet()){
				int rowInd = pair.getFirst().intValue(), colInd = pair.getSecond().intValue();
				gMatr.setElement(rowInd, colInd, getMatrixEntry(rowInd,colInd));
			}
			return gMatr;
		}
	}
	/**
	 * Returns a this matrix as a <code>Matrix</code> object.
	 * @return a Matrix object
	 */
	public Matrix convert2Matrix (){
		Matrix gMatr = new Matrix(rowDim,colDim);
		if(matrix==null){
			return gMatr;
		}
		else{
			for(Pair<Integer> pair: matrix.keySet()){
				int rowInd = pair.getFirst().intValue(), colInd = pair.getSecond().intValue();
				gMatr.set(rowInd, colInd, getMatrixEntry(rowInd,colInd));
			}
			return gMatr;
		}
	}
	/**
	 * Returns a randomly picked index pair of this matrix
	 * @return an index pair of this matrix
	 */
	public Pair<Integer> getRandomEntry (){
		Random rand = new Random();
		int randInt = rand.nextInt(matrix.size()), counter = 0;
		Pair<Integer> pair = null;
		for(Pair<Integer> index: matrix.keySet()){
			if(randInt==counter){
				pair = new Pair<Integer>(index);
				break;
			}
			else
				counter++;
		}
		return pair;
	}
	public static void main (String[] args){
		double[][] s1 = {{1d/3d,0d,0d},{0d,0d,1d},{0d,1d,0d}};
		double[][] s2 = {{3d,0d,0d},{0d,0d,1d},{0d,1d,0d}};
		double[][] a  = {{0d,0d,3d},{0d,0d,0d},{0d,0d,0d}};
		SparseMatrix S1 = new SparseMatrix(s1);
		SparseMatrix S2 = new SparseMatrix(s2);
		SparseMatrix m  = new SparseMatrix(a);
		System.out.println((S1.multiply(m.multiply(S2))).toString());
	}

}
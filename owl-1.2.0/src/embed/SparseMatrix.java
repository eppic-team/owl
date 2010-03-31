package embed;

import java.util.*;
import edu.uci.ics.jung.graph.util.*; 
//import embed.contactmaps.SparseMatrix;

/**
 * Class to deal with sparse double precision matrices, i.e. having only a few non-zero entries. This class provides
 * essential methods as addition and subtraction, multiplication (matrix or scalar), transposition
 * and fundamental setter and getter methods for this class.
 * @author gmueller
 *
 */
public class SparseMatrix {
	/**
	 * A class indicating that some exceptional situation occurred during
	 * computation of instances of the <code>{@link SparseMatrix}</code> class.
	 * Typically, addition or multiplication of instances with non-matching dimensions
	 * will cause such an exception, as well as <tt>double</tt> matrix instances, with
	 * various row dimension or index pairs exceeding the dimensions.
	 * @author gmueller
	 *
	 */
	class SparseMatrixException extends RuntimeException {
		private static final long serialVersionUID = 1L;
		/**
		 * Constructs a new sparse matrix exception indicating a serious error
		 * or exceptional situation occurred during computation, with a specified detail
		 * message.
		 * @param message the specified detail message
		 */
		public SparseMatrixException (String message){
			super(message);
		}
	}
	/** the actual representation of the matrix is a HashMap, where the key is a Pair of integers
	 * and the value is a Double
	 */
	private HashMap<Pair<Integer>,Double> matrix;
	/** a HashMap mapping the first index <tt>i</tt> to all index pairs <tt>(i,j)</tt>, where it is present*/ 
	private HashMap<Integer,HashSet<Pair<Integer>>> first;
	/** a HashMap mapping the second index <tt>j</tt> to all index pairs <tt>(i,j)</tt>, where it is present*/
	private HashMap<Integer,HashSet<Pair<Integer>>> second;
	/** the row dimension*/
	private int row_dim;
	/** the column dimension*/ 
	private int col_dim;
	/**
	 * Constructs an empty <tt>SparseMatrix</tt>. Note, that invocation of instance methods
	 * causes <tt>NullPointerException</tt>, if no initial setters are invoked on the instance,
	 * constructed by this constructor.
	 */
	public SparseMatrix (){};
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
	 * Initial setter: converts a <tt>double</tt> matrix to an instance of this class
	 * by iterating over all entries and storing only non-zero entries <tt>a_ij</tt> and their corresponding
	 * index pair <tt>(i,j)</tt>.
	 * @param mat a <tt>double</tt> matrix to be converted
	 * @throws SparseMatrixException if the inner dimensions do not match
	 */
	public void setMatrix (double[][] mat) throws SparseMatrixException {
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
						if(first.containsKey(f_ind)){
							HashSet<Pair<Integer>> sub = first.get(f_ind);
							sub.add(pair);
							first.put(f_ind, sub);
						}
						else{
							HashSet<Pair<Integer>> sub = new HashSet<Pair<Integer>>();
							sub.add(pair);
							first.put(f_ind, sub);
						}
						if(second.containsKey(s_ind)){
							HashSet<Pair<Integer>> sub = second.get(s_ind);
							sub.add(pair);
							second.put(s_ind, sub);
						}
						else{
							HashSet<Pair<Integer>> sub = new HashSet<Pair<Integer>>();
							sub.add(pair);
							second.put(s_ind, sub);
						}
					}
				}
				else throw new SparseMatrixException ("\nAt least one dimension did not match. Please, make sure, only matrices with matching dimensions are converted.");
			}
		}
		setDimension(length1,length2);
	}
	/**
	 * Initial setter: converts a <tt>HashMap</tt> instance to a <tt>SparseMatrix</tt> instance.
	 * @param map a <tt>HashMap</tt>
	 * @throws SparseMatrixException if at least one index exceeds the specified dimensions
	 */
	public void setMatrix (HashMap<Pair<Integer>,Double> map) throws SparseMatrixException {
		matrix  = new HashMap<Pair<Integer>,Double> (map);
		first   = new HashMap<Integer,HashSet<Pair<Integer>>> ();
		second  = new HashMap<Integer,HashSet<Pair<Integer>>> ();
		Set<Pair<Integer>> keys = map.keySet();
		int err1 = 0, err2 = 0;
		Iterator<Pair<Integer>> it = keys.iterator();
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			Integer f_ind = pair.getFirst(), s_ind = pair.getSecond();
			err1 = f_ind.intValue(); err2 = s_ind.intValue();
			if(err1 < row_dim && err2 < col_dim){
				if(first.containsKey(f_ind)){
					HashSet<Pair<Integer>> sub = first.get(f_ind);
					sub.add(pair);
					first.put(f_ind, sub);
				}
				else{
					HashSet<Pair<Integer>> sub = new HashSet<Pair<Integer>>();
					sub.add(pair);
					first.put(f_ind, sub);
				}
				if(second.containsKey(s_ind)){
					HashSet<Pair<Integer>> sub = second.get(s_ind);
					sub.add(pair);
					second.put(s_ind, sub);
				}
				else{
					HashSet<Pair<Integer>> sub = new HashSet<Pair<Integer>>();
					sub.add(pair);
					second.put(s_ind, sub);
				}
			}
			else throw new SparseMatrixException ("At least one index of the pair (i,j) = ("+err1+","+err2+") exceeded the specified dimensions (n,m) = ("+row_dim+","+col_dim+").");
		}
	}
	/**
	 * Initial setter: copies <tt>mat</tt> to a new <tt>SparseMatrix</tt>.
	 * @param mat a <tt>SparseMatrix</tt>
	 */
	public void setMatrix (SparseMatrix mat){
		matrix  = mat.getMatrix();
		first   = mat.getFirstIndexSet();
		second  = mat.getSecondIndexSet();
		setDimension (mat.getRowDimension(),mat.getColumnDimension());
	}
	/**
	 * Initial setter: sets the row and column dimension fields.
	 * @param row the row dimension
	 * @param col the column dimension
	 * @throws SparseMatrixException if at least one dimension parameter is less than or equal to zero
	 */
	public void setDimension (int row, int col) throws SparseMatrixException {
		if(row > 0 && col > 0){
			col_dim = col;
			row_dim = row;
		}
		else throw new SparseMatrixException ("Any dimension field must be greater than zero!");
	}
	/**
	 * Constructs the <tt>n x n</tt> dimensional identity matrix, specified by the dimension
	 * parameter <tt>dim</tt>.
	 * @param dim the dimension
	 * @return the identity matrix
	 */
	public SparseMatrix constructSquareIdentity (int dim){
		HashMap<Pair<Integer>,Double> map = new HashMap<Pair<Integer>,Double>();
		for(int i = 0; i < dim; i++){
			Integer index = new Integer(i);
			Pair<Integer> pair = new Pair<Integer>(index,index);
			map.put(pair, new Double(1.0));			
		}
		SparseMatrix id = new SparseMatrix();
		id.setDimension(dim,dim);
		id.setMatrix(map);
		return id;
	}
	/**
	 * Checks whether the <tt>Pair</tt> instance <tt>pair</tt> is present in the <tt>SparseMatrix</tt>.
	 * Returns true, iff the entry for the given index pair is non-zero.
	 * @param pair an index pair
	 * @return true, iff the entry for the index pair is non-zero
	 */
	public boolean containsIndexPair (Pair<Integer> pair){
		if(matrix.containsKey(pair)) return true;
		else return false;
	}
	/**
	 * Returns a <tt>HashMap</tt>, where the keys are <tt>Pair</tt>s of <tt>Integer</tt>s and the value
	 * is a <tt>Double</tt> instance. Note, only non-zero entries are stored in this <tt>HashMap</tt>. 
	 * @return a <tt>HashMap</tt> with all non-zero entries
	 * @throws NullPointerException if the matrix field is not initialized before calling
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
		if(i < row_dim && j < col_dim){
			Pair<Integer> pair = new Pair<Integer>  (new Integer (i), new Integer (j));
			if(containsIndexPair(pair)) return matrix.get(pair).doubleValue();
			else return 0.0;
		}
		else throw new SparseMatrixException ("\nAt least one index exceeded the matrix dimensions.");
	}
	/**
	 * Returns the row dimension
	 * @return the row dimension
	 */
	public int getRowDimension (){
		return row_dim;
	}
	/**
	 * Returns the column dimension
	 * @return the column dimension
	 */
	public int getColumnDimension (){
		return col_dim;
	}
	/**
	 * Converts an object of this class to a <tt>double</tt> matrix, including zero entries.
	 * @return a <tt>double</tt> matrix
	 * @throws NullPointerException if the matrix field is not initialized before calling
	 */
	public double[][] getFullMatrix (){
		if(matrix != null && matrix.size() > 0){
			int row = getRowDimension(), col = getColumnDimension();
			double[][] fullmat = new double[row][col];
			Set<Pair<Integer>> keys = getMatrix().keySet();
			Iterator<Pair<Integer>> it = keys.iterator();
			while(it.hasNext()){
				Pair<Integer> pair  = it.next();
				int f_ind = pair.getFirst().intValue(), s_ind = pair.getSecond().intValue();
				fullmat[f_ind][s_ind] = getMatrixEntry(f_ind, s_ind);
			}
			return fullmat;
		}
		else{
			if(matrix.size() == 0) return new double [getRowDimension()][getColumnDimension()];
			else throw new NullPointerException ("The matrix field was not initialized before calling this method.");
		}
	}
	/**
	 * Returns a set of all index pairs
	 * @return a set of all index pairs
	 * @throws NullPointerException if the matrix field is not initialized before calling
	 */
	public Set<Pair<Integer>> getIndexPairs () throws NullPointerException {
		if(matrix != null) return getMatrix().keySet();
		else throw new NullPointerException ("The matrix field was not initialized before calling this method.");
	}
	

    /**
     * Getter: returns a HashMap representation of all first indices. The index
     * is used as key mapping onto a HashSet of Pairs of Integers.
     * @return a HashSet of the first index of the matrix
     */
    public HashMap<Integer,HashSet<Pair<Integer>>> getFirstIndexSet (){
        return new HashMap<Integer,HashSet<Pair<Integer>>> (first);
    }

    /**
     * Getter: returns a HashMap representation of all second indices. The index
     * is used as key mapping onto a HashSet of Pairs of Integers.
     * @return a HashSet of the second index of the matrix
     */
    public HashMap<Integer,HashSet<Pair<Integer>>> getSecondIndexSet (){
        return new HashMap<Integer,HashSet<Pair<Integer>>>  (second);
    }

    /**
     * Returns a HashSet of Pairs of Integers, that is all index pairs sharing a
     * common first index.
     * @param i the first index in an index pair <tt>(i,j)</tt>
     * @return a HashSet of all index pairs sharing a common first index
     */
    public HashSet<Pair<Integer>> getPairSetFromFirstIndex (int i){
        return first.get(new Integer (i));
    }

    /**
     * Returns a HashSet of Pairs of Integers, that is all index pairs sharing a
     * common second index.
     * @param j the second index in an index pair <tt>(i,j)</tt>
     * @return a HashSet of all index pairs sharing a common second index
     */
    public HashSet<Pair<Integer>> getPairSetFromSecondIndex (int j){
        return second.get(new Integer (j));
    }
    
    /**
     * Lists all entries and returns the largest entry of this SparseMatrix.
     * @return the largest entry as <tt>double</tt>
     */
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
			HashMap<Pair<Integer>,Double> nmap = new HashMap<Pair<Integer>,Double>();
			Set<Pair<Integer>> key1 = map1.keySet();
			Iterator<Pair<Integer>> it1 = key1.iterator();
			while(it1.hasNext()){
				Pair<Integer> pair = it1.next();
				if(map2.containsKey(pair)){
					double sum = map1.get(pair).doubleValue()+map2.get(pair).doubleValue();
					if(sum != 0.0) nmap.put(pair, new Double(sum));
					map2.remove(pair);
				}
				else nmap.put(pair, map1.get(pair));
			}
			nmap.putAll(map2);
			return new SparseMatrix(nmap,getRowDimension(),getColumnDimension());
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
			return new SparseMatrix (newmap,mat.getRowDimension(),getColumnDimension());
		}
		else throw new SparseMatrixException ("The number of columns of the first matrix has to equal the number of rows of the second matrix.");
	}
	/**
	 * Takes this SparseMatrix instance and and raises it to the power defined by the
	 * parameter <tt>exp</tt>. So the return object is the result <tt>X</tt> of the equation
	 * <p><tt>X = this^exp</tt>. <b>Note, that in case of a dimension mismatch</b>, an intermediate
	 * (power) result will be a SparseMatrix with smaller dimension, causing SparseMatrixException, since
	 * multiplication of matrices with mismatching inner dimensions causes the same exception. So, square
	 * matrices are the only acceptable arguments. Additionally, the exponent <tt>exp</tt> must be an integer
	 * greater than or equal to zero, otherwise an SparseMatrixException is thrown.
	 * @param exp the power to which the SparseMatrix is raised
	 * @return the result <tt>X</tt> of <tt>X = this^exp</tt>
	 * @throws SparseMatrixException if the dimension don't fit
	 */
	public SparseMatrix pow (int exp) throws SparseMatrixException {
		if(getColumnDimension()==getRowDimension() && exp >=0){
			if(exp > 1){
				SparseMatrix cp1 = new SparseMatrix(this);
				SparseMatrix cp2 = new SparseMatrix(this);
				for(int i = 1; i < exp; i++){
					cp1 = cp1.multiply(cp2);
				}
				return cp1;
			}
			else return new SparseMatrix(this);
		}
		else {if(exp < 0) throw new SparseMatrixException("Only integers greater than or equal to zero are accepted!");

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
     * Returns the transposed matrix of this instance as a new SparseMatrix instance.
     * @return the transposed matrix
     * @throws java.lang.NullPointerException if the matrix field was not initialized before
     * calling this method
     */
    public SparseMatrix transpose () throws NullPointerException {
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
	
	public void convertToSparseMatrix (Set<Pair<Integer>> set, int col, int row){
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
	
	public static void main (String[] args){
		double[][] ar1 = {{0.0,1.0,1.0,1.0},{1.0,0.0,1.0,0.0},{1.0,1.0,0.0,0.0},{1.0,0.0,0.0,0.0}};
		//double[][] ar2 = {{2.0,0.0,1.0},{1.0,2.0,0.0},{0.0,0.0,2.0}};
		//double[][] ar3 = {{-1.0,1.0,0.0},{1.0,0.0,-1.0},{0.0,0.0,1.0}};
		SparseMatrix mat = new SparseMatrix(ar1);
		SparseMatrix s   = new SparseMatrix(ar1);
		//SparseMatrix s1  = new SparseMatrix(ar3);
		s = mat.multiply(s);
		System.out.println(s.toString());
	}

}
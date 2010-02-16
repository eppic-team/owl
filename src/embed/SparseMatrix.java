package embed;

import java.util.*;

import Jama.Matrix;
import edu.uci.ics.jung.graph.util.*; 
//import embed.contactmaps.SparseMatrix;

/**
 * Class to deal with sparse matrices, i.e. having only a few non-zero entries. This class provides
 * essential methods 
 * @author gmueller
 *
 */
public class SparseMatrix {
	
	private HashMap<Pair<Integer>,Double> matrix;
	
	private HashMap<Integer,HashSet<Pair<Integer>>> first;
	
	private HashMap<Integer,HashSet<Pair<Integer>>> second;
	
	private int row_dim;
	
	private int col_dim;
	
	public SparseMatrix (){};
	
	public SparseMatrix (double[][] mat){
		setMatrix(mat);
	}
	
	public SparseMatrix (HashMap<Pair<Integer>,Double> map, int row, int col){
		setMatrix(map);
		setDimension(row, col);
	}
	
	public SparseMatrix (SparseMatrix mat){
		setMatrix(mat);
	}
	
	public void setMatrix (double[][] mat){
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
				else throw new IllegalArgumentException ("\nAt least one dimension did not match. Please, make sure, only matrices with matching dimensions are converted.");
			}
		}
		setDimension(length1,length2);
	}
	
	public void setMatrix (HashMap<Pair<Integer>,Double> map){
		matrix  = new HashMap<Pair<Integer>,Double> (map);
		first   = new HashMap<Integer,HashSet<Pair<Integer>>> ();
		second  = new HashMap<Integer,HashSet<Pair<Integer>>> ();
		Set<Pair<Integer>> keys = map.keySet();
		Iterator<Pair<Integer>> it = keys.iterator();
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			Integer f_ind = pair.getFirst(), s_ind = pair.getSecond();
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
	
	public void setMatrix (SparseMatrix mat){
		matrix  = mat.getMatrix();
		first   = new HashMap<Integer,HashSet<Pair<Integer>>> (mat.getFirstIndexSet());
		second  = new HashMap<Integer,HashSet<Pair<Integer>>> (mat.getSecondIndexSet());
		setDimension (mat.getRowDimension(),mat.getColumnDimension());
	}
	
	public void setDimension (int row, int col){
		col_dim = col;
		row_dim = row;
	}
	
	public boolean containsIndexPair (Pair<Integer> pair){
		if(matrix.containsKey(pair)) return true;
		else return false;
	}
	
	public HashMap<Pair<Integer>,Double> getMatrix (){
		if(matrix != null) return new HashMap<Pair<Integer>,Double> (matrix);
		else throw new NullPointerException ("\nThe matrix field was not initialized before calling this method.");
	}
	
	public double getMatrixEntry(int i, int j){
		if(i < col_dim && j < row_dim){
			Pair<Integer> pair = new Pair<Integer>  (new Integer (i), new Integer (j));
			if(matrix.containsKey(pair)) return matrix.get(pair).doubleValue();
			else return 0.0;
		}
		else throw new IndexOutOfBoundsException ("\nAt least one index exceeded the matrix dimensions.");
	}
	
	public int getRowDimension (){
		return row_dim;
	}
	
	public int getColumnDimension (){
		return col_dim;
	}
	
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
	
	public Set<Pair<Integer>> getIndexPairs (){
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
        return new HashMap<Integer,HashSet<Pair<Integer>>> (second);
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

	
	public SparseMatrix add (SparseMatrix mat){
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
					nmap.put(pair, new Double(sum));
					map2.remove(pair);
				}
				else nmap.put(pair, map1.get(pair));
			}
			nmap.putAll(map2);
			return new SparseMatrix(nmap,getRowDimension(),getColumnDimension());
		}
		else throw new IllegalArgumentException ("The row and column dimensions must match!");
	}
	
	public SparseMatrix multiply (double scalar){
		HashMap<Pair<Integer>, Double> map  = getMatrix();
		HashMap<Pair<Integer>, Double> nmap = new HashMap<Pair<Integer>, Double>();
		Set<Pair<Integer>> keys = map.keySet();
		Iterator<Pair<Integer>> it = keys.iterator();
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			double newval      = map.get(pair).doubleValue()*scalar;
			nmap.put(pair, new Double (newval));
		}
		return new SparseMatrix (nmap,getRowDimension(),getColumnDimension());
	}
	
	public SparseMatrix multiply (SparseMatrix mat){
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
							newmap.put(newpair,new Double (sum+pro));	
						}
						else newmap.put(newpair,new Double (pro));
					}
				}
			}
			return new SparseMatrix (newmap,mat.getRowDimension(),getColumnDimension());
		}
		else throw new IllegalArgumentException ("The number of columns of the first matrix has to equal the number of rows of the second matrix.");
	}
	
	public SparseMatrix subtract (SparseMatrix mat){
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
	
	public static void main (String[] args){
		double[][] ar = {{0.0,1.0,0.0},{0.0,0.0,1.0},{0.0,0.0,0.0}};
		double[][] id = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
		SparseMatrix mat = new SparseMatrix (ar);
		SparseMatrix id1 = new SparseMatrix (id);
		Matrix mat1 = new Matrix(ar);
		Matrix id2  = new Matrix(id);
		for(int i = 0; i < 4; i++){
			id1 = id1.multiply(mat);
			System.out.println(id1.toString());
			id2 = id2.times(mat1);
			System.out.println((new SparseMatrix (id2.getArray())).toString());
		}
	}

}

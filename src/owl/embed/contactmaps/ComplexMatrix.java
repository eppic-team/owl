package owl.embed.contactmaps;


import edu.uci.ics.jung.graph.util.*;
import java.util.*;

public class ComplexMatrix {
	
	private HashMap<Pair<Integer>,ComplexDouble> matrix;
	
	private int row;
	
	private int col;
	
	private TreeMap<Integer,HashSet<Pair<Integer>>> first;
	
	private TreeMap<Integer,HashSet<Pair<Integer>>> second;
	
	public ComplexMatrix (){};
	
	public ComplexMatrix (double[][] mat){
		setMatrix(mat);
	}
	
	public ComplexMatrix (ComplexDouble[][] mat){
		setMatrix(mat);
	}
	
	public ComplexMatrix (ComplexMatrix mat){
		setDimension(mat.getRowDimension(),mat.getColumnDimension());
		setMatrix(mat.getMatrix());
		setFirstIndexSet(mat.getFirstIndexSet());
		setSecondIndexSet(mat.getSecondIndexSet());
	}
	
	public void setMatrix (ComplexDouble[][] mat){
		setDimension(mat.length,mat[0].length);
		matrix = new HashMap<Pair<Integer>,ComplexDouble> ();
		first  = new TreeMap<Integer,HashSet<Pair<Integer>>> ();
		second = new TreeMap<Integer,HashSet<Pair<Integer>>> ();
		for(int i = 0; i < row; i++){
			for(int j = 0; j < col; j++){
				if(row == mat[i].length && mat[i][j].getAbs() != 0.0){
					Integer f_val = new Integer (i), s_val = new Integer (j);
					Pair<Integer> pair = new Pair<Integer>(f_val,s_val);
					matrix.put(pair,mat[i][j]);					
					if(first.containsKey(f_val)){
						HashSet<Pair<Integer>> set = first.get(f_val);
						set.add(pair);
						first.put(f_val, set);
					}
					else{
						HashSet<Pair<Integer>> set = new HashSet<Pair<Integer>>();
						set.add(pair);
						first.put(f_val, set);
					}
					if(second.containsKey(s_val)){
						HashSet<Pair<Integer>> set = second.get(s_val);
						set.add(pair);
						second.put(s_val, set);
					}
					else{
						HashSet<Pair<Integer>> set = new HashSet<Pair<Integer>>();
						set.add(pair);
						second.put(s_val, set);
					}
				}
			}
		}
	}
	
	public void setMatrix(double[][] mat){
		setDimension(mat.length,mat[0].length);
		matrix = new HashMap<Pair<Integer>,ComplexDouble> ();
		first  = new TreeMap<Integer,HashSet<Pair<Integer>>> ();
		second = new TreeMap<Integer,HashSet<Pair<Integer>>> ();
		for(int i = 0; i < row; i++){
			for(int j = 0; j < col; j++){
				if(row == mat[i].length && mat[i][j] != 0.0){
					Integer f_val = new Integer (i), s_val = new Integer (j);
					Pair<Integer> pair = new Pair<Integer>(f_val,s_val);
					matrix.put(pair,new ComplexDouble(mat[i][j],0.0));					
					if(first.containsKey(f_val)){
						HashSet<Pair<Integer>> set = first.get(f_val);
						set.add(pair);
						first.put(f_val, set);
					}
					else{
						HashSet<Pair<Integer>> set = new HashSet<Pair<Integer>>();
						set.add(pair);
						first.put(f_val, set);
					}
					if(second.containsKey(s_val)){
						HashSet<Pair<Integer>> set = second.get(s_val);
						set.add(pair);
						second.put(s_val, set);
					}
					else{
						HashSet<Pair<Integer>> set = new HashSet<Pair<Integer>>();
						set.add(pair);
						second.put(s_val, set);
					}
				}
			}
		}
	}
	
	public void setMatrix(HashMap<Pair<Integer>,ComplexDouble> map){
		matrix = new HashMap<Pair<Integer>,ComplexDouble>(map);
	}
	
	public void setFirstIndexSet(TreeMap<Integer,HashSet<Pair<Integer>>> first){
		this.first = new TreeMap<Integer,HashSet<Pair<Integer>>> (first);
	}
	
	public void setSecondIndexSet(TreeMap<Integer,HashSet<Pair<Integer>>> second){
		this.second = new TreeMap<Integer,HashSet<Pair<Integer>>> (second);
	} 
	
	public void setDimension (int row, int col){
		this.row = row;
		this.col = col;
	}
	
	public boolean containsFirstIndex (int i){
		if(first.containsKey(new Integer(i))) return true;
		else return false;
	}
	
	public boolean containsSecondIndex (int i){
		if(second.containsKey(new Integer(i))) return true;
		else return false;
	}
	
	public HashMap<Pair<Integer>,ComplexDouble> getMatrix (){
		return new HashMap<Pair<Integer>,ComplexDouble> (matrix);
	}
	
	public ComplexDouble getMatrix (int i, int j){
		if(i < row && j < col){
			Pair<Integer> pair = new Pair<Integer>(new Integer(i),new Integer(j));
			if(matrix.containsKey(pair)) return matrix.get(pair);
			else return new ComplexDouble (0.0,0.0);
		}
		else throw new ComplexMatrixException ("Indices exceed dimensions!");
	}
	
	public ComplexDouble[][] getFullMatrix (){
		if(matrix != null){
			ComplexDouble[][] mat = new ComplexDouble[row][col];
			Set<Pair<Integer>> key = matrix.keySet();
			Iterator<Pair<Integer>> it = key.iterator();
			while(it.hasNext()){
				Pair<Integer> pair = it.next();
				int f_ind = pair.getFirst().intValue(), s_ind = pair.getSecond().intValue();
				mat[f_ind][s_ind]= matrix.get(pair);
			}
			ComplexDouble.fillMatrix(mat);
			return mat;
		}
		else throw new NullPointerException ("Matrix field not initialized!");
	}
	
	public int getRowDimension (){
		return row;
	}
	
	public int getColumnDimension (){
		return col;
	}
	
	public TreeMap<Integer,HashSet<Pair<Integer>>> getFirstIndexSet (){
		return new TreeMap<Integer,HashSet<Pair<Integer>>> (first);
	}
	
	public TreeMap<Integer,HashSet<Pair<Integer>>> getSecondIndexSet (){
		return new TreeMap<Integer,HashSet<Pair<Integer>>> (second);
	}
	
	public HashSet<Pair<Integer>> getFirstIndexSet (int i){
		Integer f_ind = new Integer(i);
		if(first.containsKey(f_ind)) return first.get(f_ind);
		else throw new ComplexMatrixException ("No such entry!");
	}
	
	public HashSet<Pair<Integer>> getSecondIndexSet (int i){
		Integer f_ind = new Integer(i);
		if(first.containsKey(f_ind)) return first.get(f_ind);
		else throw new ComplexMatrixException ("No such entry!");
	}
	
	public ComplexMatrix add(ComplexMatrix mat){
		if(col == mat.col && row == mat.row){
			HashMap<Pair<Integer>,ComplexDouble> map1 = getMatrix();
			HashMap<Pair<Integer>,ComplexDouble> map2 = mat.getMatrix();
			HashMap<Pair<Integer>,ComplexDouble> sum = new HashMap<Pair<Integer>,ComplexDouble>();
			Set<Pair<Integer>> key1 = map1.keySet();
			Set<Pair<Integer>> key2 = map2.keySet();
			Set<Pair<Integer>> help = map1.keySet();
			help.retainAll(key2);
			key1.removeAll(help);
			key2.removeAll(help);
			Iterator<Pair<Integer>> it = help.iterator();
			while(it.hasNext()){
				Pair<Integer> pair = it.next();
				ComplexDouble val1 = map1.get(pair), val2 = map2.get(pair);
				sum.put(pair, val1.add(val2));
			}
			it = key1.iterator();
			while(it.hasNext()){
				Pair<Integer> pair = it.next();
				sum.put(pair, map1.get(pair));
			}
			it = key2.iterator();
			while(it.hasNext()){
				Pair<Integer> pair = it.next();
				sum.put(pair, map1.get(pair));
			}
			ComplexMatrix sm = new ComplexMatrix ();
			sm.setMatrix(sum);
			sm.setDimension(row, col);
			return sm;
		}
		else throw new ComplexMatrixException ("Dimension mismatch!");
	}
	
	public ComplexMatrix multiply (ComplexDouble z){
		if(matrix != null){
			HashMap<Pair<Integer>,ComplexDouble> map = new HashMap<Pair<Integer>,ComplexDouble>(); 
			Set<Pair<Integer>> keys = matrix.keySet();
			Iterator<Pair<Integer>> it = keys.iterator();
			while(it.hasNext()){
				Pair<Integer> pair = it.next();
				ComplexDouble mult = matrix.get(pair).multiply(z);
				map.put(pair, mult);
			}
			ComplexMatrix mat = new ComplexMatrix();
			mat.setDimension(row, col);
			mat.setMatrix(map);
			mat.setFirstIndexSet(first);
			mat.setSecondIndexSet(second);
			return mat;
		}
		else throw new NullPointerException ("Matrix field not initialized!");
	}
	
	public ComplexMatrix multiply (ComplexMatrix mat){
		if(getColumnDimension() == mat.getRowDimension()){
			HashMap<Pair<Integer>,ComplexDouble> map1 = getMatrix();
			HashMap<Pair<Integer>,ComplexDouble> map2 = mat.getMatrix();
			HashMap<Pair<Integer>,ComplexDouble> prod = new HashMap<Pair<Integer>,ComplexDouble> ();
			TreeMap<Integer,HashSet<Pair<Integer>>> f = new TreeMap<Integer,HashSet<Pair<Integer>>>();
			TreeMap<Integer,HashSet<Pair<Integer>>> s = new TreeMap<Integer,HashSet<Pair<Integer>>>();
			Set<Pair<Integer>> key1 = map1.keySet();
			Iterator<Pair<Integer>> it1 = key1.iterator();
			while(it1.hasNext()){
				Pair<Integer> pair1 = it1.next();
				Set<Pair<Integer>> key2 = map2.keySet();
				Iterator<Pair<Integer>> it2 = key2.iterator();
				while(it2.hasNext()){
					Pair<Integer> pair2 = it2.next();
					if(pair1.getSecond().equals(pair2.getFirst())){
						Integer f_val = pair1.getFirst(), s_val = pair2.getSecond();
						Pair<Integer> npair = new Pair<Integer>(f_val,s_val);
						ComplexDouble pr = map1.get(pair1).multiply(map2.get(pair2));
						if(f.containsKey(f_val)){
							HashSet<Pair<Integer>> set = f.get(f_val);
							set.add(npair);
							f.put(f_val, set);
						}
						else {
							HashSet<Pair<Integer>> set = new HashSet<Pair<Integer>>();
							set.add(npair);
							f.put(f_val, set);
						}
						if(s.containsKey(s_val)){
							HashSet<Pair<Integer>> set = s.get(s_val);
							set.add(npair);
							s.put(s_val, set);
						}
						else {
							HashSet<Pair<Integer>> set = new HashSet<Pair<Integer>>();
							set.add(npair);
							s.put(s_val, set);
						}
						if(prod.containsKey(npair)){
							ComplexDouble su = pr.add(prod.get(npair));
							prod.put(npair, su);
						}
						else prod.put(npair, pr);
					}
				}
			}
			ComplexMatrix prd = new ComplexMatrix();
			prd.setDimension(row, mat.col);
			prd.setMatrix(prod);
			prd.setFirstIndexSet(f);
			prd.setSecondIndexSet(s);
			return prd;
		}
		else throw new ComplexMatrixException("Dimension mismatch!");
	}
	
	public TreeSet<Integer> getUnPresentRow (){
		TreeSet<Integer> tree = new TreeSet<Integer>();
		for(int i = 0; i < row; i++){
			if(!containsFirstIndex(i)) tree.add(new Integer (i));
		}
		return tree;
	}
	
	public TreeSet<Integer> getUnPresentColumn (){
		TreeSet<Integer> tree = new TreeSet<Integer>();
		for(int i = 0; i < col; i++){
			if(!containsSecondIndex(i)) tree.add(new Integer (i));
		}
		return tree;
	}
	
	public ComplexMatrix scalarAddRow (int i, int j, ComplexDouble z){
		HashSet<Pair<Integer>> rows1 = getFirstIndexSet(i);
		HashSet<Pair<Integer>> rows2 = getFirstIndexSet(j);
		HashSet<Pair<Integer>> sur   = new HashSet<Pair<Integer>>();
		HashMap<Pair<Integer>,ComplexDouble> map = new HashMap<Pair<Integer>,ComplexDouble>();
		Iterator<Pair<Integer>> it = rows1.iterator();
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			Iterator<Pair<Integer>> it2 = rows2.iterator();
			while(it2.hasNext()){
				Pair<Integer> pair2 = it2.next();
				if(pair.getSecond().equals(pair2.getSecond())){
					ComplexDouble val1 = matrix.get(pair);
					ComplexDouble val2 = matrix.get(pair2);
					map.put(pair2, val2.add(val1.multiply(z)));
					map.put(pair, val1);
					sur.add(pair);
					sur.add(pair2);
				}
			}	
		}
		Set<Pair<Integer>> key = matrix.keySet();
		key.removeAll(sur);
		it = key.iterator();
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			map.put(pair, matrix.get(pair));
		}
		ComplexMatrix mat = new ComplexMatrix ();
		mat.setMatrix(map);
		mat.setDimension(row, col);
		return mat;
	}
	
	public String toString (){
		if(matrix != null){
			String str = "";
			ComplexDouble[][] ar = getFullMatrix();
			for(int i = 0; i < row; i++){
				for(int j = 0; j < col; j++){
					if(j < col-1) str += ar[i][j].toString() + "\t";
					else str += ar[i][j].toString() + "\n";
				}
			}
			return str;
		}
		else throw new NullPointerException ("Matrix field not initialized!");
	}
	
	public static void main (String[] args){
		double[][] mat1 = {{0.0,1.0,0.0},{1.0,1.0,1.0},{0.0,1.0,0.0}};
		double[][] tras = {{1.0,1.0,1.0},{0.0,-1.0,2.0},{-1.0,1.0,1.0}};
		double[][] tra1 = {{3.0,0.0,-3.0},{2.0,-2.0,2.0},{1.0,2.0,1.0}};
		//Matrix mat = new Matrix(mat1);
		ComplexMatrix tra = new ComplexMatrix(tras);
		ComplexMatrix tr1 = (new ComplexMatrix(tra1)).multiply(new ComplexDouble(1.0/6.0,0.0));
		ComplexMatrix ma  = new ComplexMatrix(mat1);
		//ComplexDouble zer = new ComplexDouble(0.0,0.0);
		ComplexDouble sq2 = new ComplexDouble(Math.sqrt(2.0),0.0);
		ComplexDouble ima = new ComplexDouble(0.0,1.0);
		ComplexDouble[][] ar = new ComplexDouble[3][3];
		ComplexDouble.fillMatrix(ar);
		for(int i = 0; i < 4; i++){
			ar[1][1] = ima.multiply(Math.pow(-1.0,i));
			if(i < 2) ar[2][2] = sq2;
			else ar[2][2] = sq2.multiply(-1.0);
			ComplexMatrix mat = new ComplexMatrix (ar);
			ComplexMatrix pr  = tra.multiply(mat.multiply(tr1));
			ComplexMatrix d   = tr1.multiply(ma.multiply(tra));
			System.out.println(pr.toString());
			System.out.println(pr.multiply(pr).toString());
			System.out.println(tra.multiply(tr1).toString());
			System.out.println(d.toString());
		}
	}
	
	/*public Matrix leibnizTransform (){
		TreeSet<Integer> rows = getUnPresentRow();
		TreeSet<Integer> cols = getUnPresentColumn();
		if(rows.size() > 0 && cols.size() > 0){
			Iterator<Integer> it = null;
			if(rows.size() > cols.size()){
				it = rows.iterator();
			}
			else it = cols.iterator();
			while
		}
	}*/
	
	class ComplexMatrixException extends RuntimeException {
		private static final long serialVersionUID = 1L;
		public ComplexMatrixException (String message){
			super(message);
		}
	}

}


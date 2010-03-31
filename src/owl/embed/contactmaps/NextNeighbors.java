package owl.embed.contactmaps;

import java.util.*;

import owl.embed.SparseMatrix;
import edu.uci.ics.jung.graph.util.*;

public class NextNeighbors {
	
	private SparseMatrix matrix;
	
	private HashMap<Integer,HashSet<Integer>> neighbors1;
	
	public NextNeighbors(SparseMatrix matrix){
		setFields(matrix);
	}
	
	public void setFields(SparseMatrix matrix){
		this.matrix = new SparseMatrix(matrix);
		setNeighbors();
	}
	
	public void setNeighbors (){
		if(matrix != null){
			neighbors1 = new HashMap<Integer,HashSet<Integer>>();
			HashMap<Integer,HashSet<Pair<Integer>>> index1 = matrix.getFirstIndexSet();
			Set<Integer> keys = index1.keySet();
			Iterator<Integer> it = keys.iterator();
			while(it.hasNext()){
				Integer index = it.next();
				HashSet<Pair<Integer>> indexset = index1.get(index);
				Iterator<Pair<Integer>> it2 = indexset.iterator();
				HashSet<Integer> subhash = new HashSet<Integer>();
				while(it2.hasNext()){
					Pair<Integer> pair = it2.next();
					subhash.add(pair.getSecond());
				}
				neighbors1.put(index, subhash);
			}
		}
	}
	
	public static void main (String[] args){
		double[][] ar1 = {{0.0,1.0,1.0,1.0},{1.0,0.0,1.0,0.0},{1.0,1.0,0.0,0.0},{1.0,0.0,0.0,0.0}};
		int length = 7;
		SparseMatrix mat = new SparseMatrix(ar1);
		SparseMatrix s   = new SparseMatrix(ar1);
		NextNeighbors[] nn = new NextNeighbors[length];
		for(int i = 1; i < length; i++){
			nn[i-1] = new NextNeighbors(s);
			s = mat.multiply(s);
			System.out.println(s.toString());
		}
	}

}

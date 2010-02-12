package tools;

import java.util.Comparator;

import edu.uci.ics.jung.graph.util.Pair;

/**
 * A Comparator for Pair<Integer>, it can be used in constructing sorted 
 * collections of Pair<Integer>. Pair doesn't implement comparable, that's why
 * this is needed 
 * 
 * Example use: 
 * 
 * TreeMap<Pair<Integer>> map = new TreeMap<Pair<Integer>>(new IntPairComparator());
 *
 */
public class IntPairComparator implements Comparator<Pair<Integer>> {

	public int compare(Pair<Integer> pair1, Pair<Integer> pair2) {
		if (pair1.getFirst()>pair2.getFirst()){
			return 1;
		} 
		else if (pair1.getFirst()<pair2.getFirst()){
			return -1;
		}
		else {
			if (pair1.getSecond()>pair2.getSecond()){
				return 1;
			}
			else if (pair1.getSecond()<pair2.getSecond()){
				return -1;
			}
		}				
		return 0;
	}

}

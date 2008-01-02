package proteinstructure;

import java.util.TreeSet;

/**
 * Class representing a unique, ordered (according to 
 * ordering define in Interval.compareTo) set of intervals
 *
 */
public class IntervalSet extends TreeSet<Interval>{

	private static final long serialVersionUID = 1L;
	
	/**
	 * Returns an ordered set of integers result of the intersection of all Intervals in this set
	 * @return
	 */
	public TreeSet<Integer> getIntegerSet() {
		TreeSet<Integer> set = new TreeSet<Integer>();
		for (Interval interv:this) {
			for (int i=interv.beg;i<=interv.end;i++) {
				set.add(i);
			}
		}
		return set; 
	}
	
}

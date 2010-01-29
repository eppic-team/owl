package proteinstructure;

import java.util.Iterator;
import java.util.TreeSet;

/**
 * Class representing a unique, ordered (according to 
 * ordering defined in Interval.compareTo) set of intervals
 */
public class IntervalSet extends TreeSet<Interval> implements Comparable<IntervalSet> {

	private static final long serialVersionUID = 1L;
		
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Returns an ordered set of integers resulting from the intersection of all Intervals in this set
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
	
	/*-------------------------- implemented methods ------------------------*/
	
	/**
	 * Compares two IntevalSets according to a lexicographical-like ordering.
	 */
	public int compareTo(IntervalSet other) {
		// whichever set has a smaller first element is smaller
		int compare = this.first().compareTo(other.first());
		if(compare != 0) return compare;
		// otherwise, we have to go through all elements until one is missing
		TreeSet<Integer> thisElements = this.getIntegerSet();
		TreeSet<Integer> otherElements = other.getIntegerSet();
		Iterator<Integer> thisIt = thisElements.iterator();
		Iterator<Integer> otherIt = otherElements.iterator();
		while(thisIt.hasNext() && otherIt.hasNext()) {
			int thisNext = thisIt.next();
			int otherNext = otherIt.next();
			if(thisNext != otherNext) {
				return new Integer(thisNext).compareTo(new Integer(otherNext));
			}
		}
		if(thisIt.hasNext() && !otherIt.hasNext()) return 1;	// this is longer, other ends
		if(!thisIt.hasNext() && otherIt.hasNext()) return -1;	// this ends, other is longer
		// if none of the above checks could break the tie, then we are equal
		return 0;
	}
	
}

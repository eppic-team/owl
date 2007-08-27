package proteinstructure;

import java.util.TreeSet;
import java.util.Vector;

public class NodeSet extends TreeSet<Integer> {

	private static final long serialVersionUID = 1L;

	/**
	 * Returns the set as a vector of intervals of consecutive elements.
	 * @return
	 */
	public Vector<Interval> getIntervals() {
		Vector<Interval> intervals = new Vector<Interval>();
		if(this.size() == 0) return intervals;
		// assuming that this set is sorted, otherwise sort it
		int last = this.first();	// previous element
		int start = last;			// start if current interval
		for(int i:this) {
			if(i > last+1) {
				// output interval and start new one
				intervals.add(new Interval(start, last));
				start = i;
				last = i;
			} else
			if(i == last) {
				// can that even happen?
			} else
			if(i == last+1) {
				last = i;
			}
		}
		// output last interval
		intervals.add(new Interval(start, last));
		return intervals;
	}
	
}

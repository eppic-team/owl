package proteinstructure;

import java.util.TreeSet;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
	
	/**
	 * Returns true if selStr is a valid selection string in 'comma-hyphen' syntax, e.g. 1-3,5,7-8.
	 * @return true if selStr is a syntactically correct selection string, false otherwise
	 */
	public static boolean isValidSelectionString(String selStr) {
		Pattern p = Pattern.compile("\\d+(-\\d+)?(,\\d+(-\\d+)?)*");
		Matcher m = p.matcher(selStr);
		return m.matches();
	}
	
	/**
	 * Create a new NodeSet from a selection string in 'comma-hyphen' nomenclature, e.g. 1-3,5,7-8.
	 * The validity of a selection string can be checked by isValidSelectionString().
	 * @return A node set corresponding to the given selection string or null of string is invalid.
	 */
	public static NodeSet parseSelectionString(String selStr) {
		if(!isValidSelectionString(selStr)) return null;
		NodeSet newSet = new NodeSet();
		String[] tokens = selStr.split(",");
		for(String t:tokens) {
			if(t.contains("-")) {
				String[] range = t.split("-");
				int from = Integer.parseInt(range[0]);
				int to = Integer.parseInt(range[1]);
				for(int i=from; i <= to; i++) {
					newSet.add(i);
				}
				
			} else {
				int num = Integer.parseInt(t);
				newSet.add(num);
			}
		}
		return newSet;
	}
	
}

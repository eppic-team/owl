package proteinstructure;

import java.util.TreeSet;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// TODO: Create a Node class (which can contain node coloring, class, residue type, etc.)
public class NodeSet extends TreeSet<Node> {

	/*------------------------------ constants ------------------------------*/
	
	private static final long serialVersionUID = 1L;
	
	/*----------------------------- constructors ----------------------------*/
	
	public NodeSet() {
		super();
	}	
	
	/*---------------------------- public methods ---------------------------*/
	
	public boolean contains(Object o) {
		if(!(o instanceof Node)) {
			System.err.println("");
			throw new ClassCastException("Trying to call contains() on a node set with a non-Node parameter.");
		}
		return super.contains(o);
	}
	
	/**
	 * Returns a deep copy of this NodeSet
	 */
	public NodeSet copy() {
		NodeSet copy = new NodeSet();
		for(Node n:this) {
			copy.add(n.copy());
		}
		return copy;
	}
	
	/**
	 * Returns true if this NodeSet equals the given other NodeSet. Nodes are considered
	 * equal if they have the same num.
	 * @param other the NodeSet to compare with
	 * @return true if the two sets are equal, false otherwise
	 */
	public boolean equals(NodeSet other) {
		if(this.size() != other.size()) return false;
		for(Node n:this) {
			if(!other.contains(n)) {
				return false;
			}
		}
		// just to be safe:
		for(Node n:other) {
			if(!this.contains(n)) {
				return false;
			}
		}
		return true;
	}
	
	/**
	 * Returns the set as a vector of intervals of consecutive elements.
	 * @return
	 */
	public Vector<Interval> getIntervals() {
		Vector<Interval> intervals = new Vector<Interval>();
		if(this.size() == 0) return intervals;
		// assuming that this set is sorted, otherwise sort it
		int last = this.first().getNum();	// previous element
		int start = last;			// start if current interval
		for(Node n:this) {
			int i = n.getNum();
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
					newSet.add(new Node(i));
				}
				
			} else {
				int num = Integer.parseInt(t);
				newSet.add(new Node(num));
			}
		}
		return newSet;
	}
	
	/**
	 * Returns the first node with the given property set to the given value.
	 * @return the first matching node or null if none found
	 */
	public Node getNodeByProperty(String property, String value) {
		Node result = null;
		for(Node n:this) {
			if(n.hasProperty(property) && n.getProperty(property).equals(value)) {
				result = n;
				break;
			}
		}
		return result;
	}
	
	/*--------------------------------- main --------------------------------*/
	
	/* Test this class */
	public static void main(String[] args) {
		NodeSet nodes = new NodeSet();
		// fill up with nodes from 0 to 9
		for(int i=0; i<10; i++) {
			Node n = new Node(i);
			nodes.add(n);
		}
		System.out.println("Nodes: " + nodes);
		// add some properties to first node
		Node n = new Node(11);
		n.setProperty("Hallo", "Otto");
		Node n2 = new Node(11);
		System.out.println("n == n2:\t" + (n == n2));
		System.out.println("n.equals(n2):\t" + (n.equals(n2)));
		System.out.println("n.compareTo(n2):" + (n.compareTo(n2)));
		System.out.println("nodes.contains(new Node(3)):" + (nodes.contains(new Node(3))));
		nodes.remove(new Node(3));
		System.out.println("Nodes: " + nodes);		
		System.out.println("nodes.contains(new Node(3)):" + (nodes.contains(new Node(3))));
	}
	
}

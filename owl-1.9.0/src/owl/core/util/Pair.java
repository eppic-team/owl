package owl.core.util;

import java.util.HashMap;
import java.util.HashSet;

import javax.swing.JComponent;
import javax.swing.JMenuItem;

/**
 * Simple data structure for the storage of two object of arbitrary types in a 
 * pair. This class is inspired by the STL utility class <code>pair</code>.
 * @author Lars Petzold
 *
 * @param <F> type of the first object
 * @param <S> type of the second object
 */
public class Pair<F,S> {
    
    F first;
    S second;
    
    /**
     * Create a new pair.
     * @param first  the first object
     * @param second  the secon object
     */
    public Pair(F first, S second) {
	this.first = first;
	this.second = second;
    }
    
    /**
     * Gets the F type object.
     * @return the first object
     */
    public F getFirst() {
	return first;
    }
    
    /**
     * Gets the S type object.
     * @return the second object
     */
    public S getSecond() {
	return second;
    }
    
    /**
     * Sets the F type object.
     * @param first  the first object
     */
    public void setFirst(F first) {
	this.first = first;
    }
    
    /**
     * Sets the S type object.
     * @param second  the second object
     */
    public void setSecond(S second) {
	this.second = second;
    }
    
    /**
     * Gets the string representation of the pair which is the comma-separated 
     * string representation of the types F and S encapsulated in squared 
     * brackets.
     * @return the Pair's string representation
     */
    public String toString() {
	return "[" + first.toString() + "," + second.toString() + "]";
    }
    
    public static void main(String args[]) {
	HashSet<Pair<Integer,Integer>> h1 = new HashSet<Pair<Integer,Integer>>();
	h1.add(new Pair<Integer, Integer>(1,2));
	h1.add(new Pair<Integer, Integer>(4,2));
	
	JMenuItem m1 = new JMenuItem();
	JMenuItem m2 = new JMenuItem();
	JMenuItem m3 = new JMenuItem();
	HashSet<Pair<JComponent, Boolean>> h2 = new HashSet<Pair<JComponent, Boolean>>();
	h2.add(new Pair<JComponent, Boolean>(m1,true));
	h2.add(new Pair<JComponent, Boolean>(m2,true));
	h2.add(new Pair<JComponent, Boolean>(m3,false));
	
	HashMap<JComponent,Boolean> h3 = new HashMap<JComponent, Boolean>();
	h3.put(m1,true);
	h3.put(m2,true);
	h3.put(m3,false);
	
    }
}

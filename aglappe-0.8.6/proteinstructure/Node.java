package proteinstructure;
import java.util.*;

/**
 * A generic node in a graph. Special kinds of nodes (like residues in a RIG or proteins in a PPI)
 * can be derived from this class. Alternatively, custom properties can be used to store additional data.
 * See also: NodeSet, Edge, NodesAndEdges
 * @author Henning Stehr
 */
public class Node  implements Comparable {
	
	/*--------------------------- member variables --------------------------*/
	public int num;		// the node id, has to be unique in the graph, public for convenience
	private HashMap<String,String> properties;	// custom properties of this node
	
	/*----------------------------- constructors ----------------------------*/
	
	public Node(int num) {
		this.num = num;
		this.properties = new HashMap<String,String>();
	}
	
	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Returns the number of this node
	 * @return the num
	 */
	public int getNum() {
		return num;
	}
	
	/**
	 * Sets the property with the given name to the given value.
	 * @param name
	 * @param value
	 */
	public void setProperty(String name, String value) {
		this.properties.put(name, value);
	}
	
	/**
	 * Gets the value of the property with the given name. Return value is of type Object
	 * and has to be casted back to the original type.
	 * @return the value of the property or null if the node has no property of that name
	 */
	public String getProperty(String name) {
		return properties.get(name);
	}
	
	/**
	 * Returns true if the node has a property with this name
	 * @param name name of the property
	 * @return true if the node has this property, false otherwise
	 */
	public boolean hasProperty(String name) {
		return properties.containsKey(name);
	}
	
	/**
	 * Returns an (almost deep) copy of this node. The copy has its own identity and its own properties object but
	 * the property values refer to the same objects as the original.
	 * @return A copy of this node
	 */
	public Node copy() {
		Node newNode = new Node(this.num);
		for(String propName:properties.keySet()){
			String val = properties.get(propName);
			newNode.setProperty(propName, val);
		}
		return newNode;
	}
	
	/**
	 * Return true if the nodes are the same. Two nodes are defined to be equal if they have the same num.
	 * @param otherNode
	 * @return true if the two nodes are equal.
	 */
	public boolean equals(Object o) {
		Node other = (Node) o;
		return this.num == other.num;
	}
	
	/**
	 * Standard comparison method.
	 * @param o
	 * @return
	 */
	public int compareTo(Object o) {
		Node other = (Node) o;
		if (this.num>other.num){
			return 1;
		} 
		else if (this.num<other.num){
			return -1;
		}
		return 0; // if none of the conditions before returned, nums are equal 
	}
	
	/**
	 * Returns a string representation of this node.
	 * @return this node as a string 
	 */
	public String toString() {
		return Integer.toString(this.getNum());
	}
	
	/**
	 * Print the properties of this node to stdout.
	 *
	 */
	public void printProperties() {
		System.out.println("Properties of node " + this.getNum() + ":");
		System.out.println(properties);
	}
	
}

package proteinstructure;

import java.util.ArrayList;
import java.util.Collection;
import java.util.TreeMap;



public class RIGNbhood extends TreeMap<Integer,RIGNode> {
	

	private static final long serialVersionUID = 1L;
	
	public static final String centralLetter="x";
	public static final String gapLetter="_";

	// central residue
	private RIGNode centralResidue;
	
	/**
	 * Construct a RIGNbhood by passing the central residue as a RIGNode object
	 * and all neighbors as a Collection<RIGNode>
	 * @param centralResidue
	 * @param neighbors
	 */
	public RIGNbhood(RIGNode centralResidue, Collection<RIGNode> neighbors){
		super();
		this.centralResidue = centralResidue;
		this.put(centralResidue.getResidueSerial(), centralResidue);
		for (RIGNode node:neighbors) {
			this.put(node.getResidueSerial(),node);
		}
	}
	
	/**
	 * Constructs a RIGNbhood with just a central residue and no neighbors
	 * @param centralResidue
	 */
	public RIGNbhood(RIGNode centralResidue) {
		super();
		this.centralResidue = centralResidue;
		this.put(centralResidue.getResidueSerial(), centralResidue);
	}

	public RIGNode getCentralResidue() {
		return centralResidue;
	}
	
	/**
	 * Returns true if one of the neighbors or the central residue is of the given residue type
	 * @param resType
	 * @return
	 */
	public boolean containsResType(String resType) {
		for (RIGNode node:this.values()) {
			if (node.getResidueType().equals(resType)) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Returns whether this RIGNbhood is equal to the given one
	 * Equality defined as: same residue types in the same order, both to the left and right of central residue
	 * @param other
	 * @return
	 */
	public boolean equals(Object other) { 

		
		RIGNbhood otherNbhood = (RIGNbhood) other;
		if (this.size()!=otherNbhood.size()) {
			return false;
		}
		if (!this.centralResidue.getResidueType().equals(otherNbhood.centralResidue.getResidueType())) {
			return false;
		}
		
		// we order the RIGNodes in 2 groups: left and right of the central residues (both for this and other)
		ArrayList<RIGNode> thisLeft = new ArrayList<RIGNode>();
		ArrayList<RIGNode> thisRight = new ArrayList<RIGNode>();
		ArrayList<RIGNode> otherLeft = new ArrayList<RIGNode>();
		ArrayList<RIGNode> otherRight = new ArrayList<RIGNode>();
		for (RIGNode node:this.values()) {
			if (node.getResidueSerial()<this.centralResidue.getResidueSerial()) {
				thisLeft.add(node);
			} else if (node.getResidueSerial()!=this.centralResidue.getResidueSerial()){
				thisRight.add(node);
			}
		}
		for (RIGNode node:otherNbhood.values()) {
			if (node.getResidueSerial()<otherNbhood.centralResidue.getResidueSerial()) {
				otherLeft.add(node);
			} else if (node.getResidueSerial()!=otherNbhood.centralResidue.getResidueSerial()){
				otherRight.add(node);
			}			
		}
		
		// if sizes don't match they are different
		if (thisLeft.size()!=otherLeft.size()) {
			return false;
		}
		if (thisRight.size()!=otherRight.size()) {
			return false;
		}
		
		// as the RIGNodes are ordered as originally (sequence order) we simply compare if types match (left and then right)
		for (int i=0;i<thisLeft.size();i++) {
			if (!thisLeft.get(i).getResidueType().equals(otherLeft.get(i).getResidueType())) {
				return false;
			}
		}
		for (int i=0;i<thisRight.size();i++) {
			if (!thisRight.get(i).getResidueType().equals(otherRight.get(i).getResidueType())) {
				return false;
			}
		}

		// if we are here, we passed all tests: they are equal!
		return true;
	}

	
	/**
	 * Returns true if this neighborhood matches the given one.
	 * A match is defined as: both left and right sides of the this central 
	 * residue are substrings of other's left and right sides.
	 * e.g.: 
	 * 		this: ABxEF matches: 
	 * 		other: --A-----B---xE---F--- (with '-' being other residues in between)
	 * @param other
	 * @return
	 */
	public boolean match(RIGNbhood other) {
		// we order the RIGNodes in 2 groups: left and right of the central residues (both for this and other)
		ArrayList<RIGNode> thisLeft = new ArrayList<RIGNode>();
		ArrayList<RIGNode> thisRight = new ArrayList<RIGNode>();
		ArrayList<RIGNode> otherLeft = new ArrayList<RIGNode>();
		ArrayList<RIGNode> otherRight = new ArrayList<RIGNode>();
		for (RIGNode node:this.values()) {
			if (node.getResidueSerial()<this.centralResidue.getResidueSerial()) {
				thisLeft.add(node);
			} else if (node.getResidueSerial()!=this.centralResidue.getResidueSerial()){
				thisRight.add(node);
			}
		}
		for (RIGNode node:other.values()) {
			if (node.getResidueSerial()<other.centralResidue.getResidueSerial()) {
				otherLeft.add(node);
			} else if (node.getResidueSerial()!=other.centralResidue.getResidueSerial()){
				otherRight.add(node);
			}			
		}
		
		// as the RIGNodes are ordered as originally (sequence order) we simply compare if types match (left and then right)
		int i = 0;
		int j = 0;
		while (i<thisLeft.size()) { 
			if (j==otherLeft.size()) {
				return false;
			}
			if (thisLeft.get(i).getResidueType().equals(otherLeft.get(j).getResidueType())) {
				i++;
				j++;				
			} else {
				j++;
			}
		}
		
		i = 0;
		j = 0;
		while (i<thisRight.size()) {
			if (j==otherRight.size()) {
				return false;
			}
			if (thisRight.get(i).getResidueType().equals(otherRight.get(j).getResidueType())) {
				i++;
				j++;
			} else {
				j++;
			}
		}

		// if we are here, we passed all tests: they match!
		return true;
	}
	
	public String getMotifFullGaps(){
		String motif="";
		for (int i=this.firstKey();i<=this.lastKey();i++) {
			if (this.containsKey(i)){
				if (i!=centralResidue.getResidueSerial()){
					motif+=AAinfo.threeletter2oneletter(this.get(i).getResidueType());
				} else {
					motif+=centralLetter;
				}
			} else {
				motif+=gapLetter;
			}
		}
		return motif;
	}
	
	public String getMotif(){
		String motif="";
		int gapSize = 0;
		String gap = "";
		for (int i=this.firstKey();i<=this.lastKey();i++) {
			if (this.containsKey(i)){
				if (i!=centralResidue.getResidueSerial()){
					motif+=gap;
					motif+=AAinfo.threeletter2oneletter(this.get(i).getResidueType());
					gapSize=0;
					gap="";
				} else {
					motif+=gap;
					motif+=centralLetter;
					gapSize=0;
					gap="";
				}
			} else {
				gapSize++;
				gap=gapLetter+"{"+gapSize+"}";
			}
		}
		return motif;
	}

	public String getMotifNoGaps(){
		String motif="";
		for (int i:this.keySet()) {
			if (i!=centralResidue.getResidueSerial()){
				motif+=AAinfo.threeletter2oneletter(this.get(i).getResidueType());
			} else {
				motif+=centralLetter;
			}
		}
		return motif;
	}

	public String getMotifReducedAlphabet(RIGraph graph) {
		String motif="";
		for (int i:this.keySet()) {
			if (i!=centralResidue.getResidueSerial()){
				// as long as the graph is undirected, JUNG finds the edge correctly independently of the order in which we give the nodes 
				RIGEdge edge = graph.findEdge(centralResidue,this.get(i));
				double weight = edge.getWeight();
				if (weight>0) { // SC dominated
					motif+=AAinfo.threeletter2oneletter(this.get(i).getResidueType());
				} else if (weight<0) { //BB dominated
					// for the sec structure we take the "seen" node (so not the central but the other), following Michael's convention
					char ssType = this.get(i).getSecStrucElement().getType();
					char ssLetter = 'o';
					if (ssType == 0 || ssType == 'O') ssLetter='o';
					if (ssType=='S') ssLetter = 'b';
					else if (ssType=='H') ssLetter= 'z';
					motif+=ssLetter;
				}
			} else {
				motif+=centralLetter;
			}
		}
		return motif;
	}
	
	public String getCommaSeparatedResSerials(){
		String ressers="";
		for (int resser:this.keySet()){
			if (resser!=centralResidue.getResidueSerial()) {
				ressers += resser+",";
			} 
		}
		// we chop off the last comma
		if (ressers.length()!=0) { 
			ressers = ressers.substring(0, ressers.length()-1);
		}
		return ressers;
	}

	public String toString(){
		if (this.isEmpty()) return "";
		else return this.getMotif();
	}
	
	/**
	 * Returns a copy (deep) of this RIGNbhood 
	 * @return
	 */
	public RIGNbhood copy(){
		RIGNbhood copy = new RIGNbhood(centralResidue.copy());
		for (int residueSerial:this.keySet()) {
			if (residueSerial!=centralResidue.getResidueSerial()) {
				copy.put(residueSerial, this.get(residueSerial).copy());
			}
		}
		return copy;
	}
}

package proteinstructure;

import java.util.Set;
import java.util.TreeMap;

import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * Class representing a protein structure graph,   
 * subclasses are RIGraph and AIGraph
 * 
 * NOTE: ProtStructGraph extends SparseGraph from jung-impl and thus 
 * 		it accepts parallel and loop edges.
 * 		It also accepts a mixed of directed/undirected edges.
 * 		In our protein structure graphs we want neither 
 * 		parallel/loop edges nor mixed directed/undirected edges 
 * TODO: override the addEdge function of SparseGraph so that ProtStructGraph 
 * 		doesn't accept any of the above
 * 
 * @param <V>
 * @param <E>
 */
public abstract class ProtStructGraph<V,E> extends SparseGraph<V,E> {

	
	protected static final int DEFAULT_MODEL = 1;
	protected final static String GRAPHFILEFORMATVERSION = "1.0";
	public static final String NO_CONTACT_TYPE = "";
	public static final double NO_CUTOFF = 0.0;
	public static final int	   NO_SEQ_SEP_VAL =	-1; // default seq sep value indicating that no seq sep has been specified

 	
	protected String sequence;		// the full sequence (with unobserved residues and non-standard aas ='X')
	protected String pdbCode;		// the lower-case pdb code
	protected String pdbChainCode;	// The pdb chain code (upper case), i.e. the classic (author's) pdb code (Pdb.NULL_CHAIN_CODE if it is blank in original pdb file)
	protected String chainCode;		// Our internal chain identifier (upper case)
									// - in reading from pdbase or from msdsd it will be set to the internal chain id (asym_id field for pdbase, pchain_id for msdsd)
    								// - in reading from pdb file it coincides with pdbChainCode except for Pdb.NULL_CHAIN_CODE where we use "A"	
	protected int model;			// model serial number (for NMR structures, for all others is 1)
	protected String sid;			// the scop id if this ProtStructGraph comes from a Pdb object restricted to a SCOP domain
	protected String comment;		// a user-defined comment field
	
	// optional fields for graphs based on casp predictions
	protected int targetNum;
	protected int caspModelNum;
	protected int groupNum;
	protected String authorStr;
	protected String methodStr;
	protected String[] caspParents;	// list of parents used for modelling if this graph is derived from a Casp prediction, may be null
	
	protected int fullLength;		// full length of the protein with all aas (same as sequence length)

	protected int minSeqSep;
	protected int maxSeqSep;
	protected boolean interSSE = false;
	
	protected SecondaryStructure secondaryStructure;
	
	protected TreeMap<Integer,V> serials2nodes;
	
	public ProtStructGraph() {
		super();
		this.sequence=null;
		this.pdbCode=null;
		this.chainCode=null;
		this.pdbChainCode=null;
		this.model=DEFAULT_MODEL;
		this.fullLength=0;
		this.minSeqSep = NO_SEQ_SEP_VAL;
		this.maxSeqSep = NO_SEQ_SEP_VAL;
		this.secondaryStructure = null;
		this.comment = null;
	}

	/**
	 * Returns the contact range (sequence separation) of the given edge
	 * @param edge
	 * @return
	 */
	public abstract int getContactRange(E edge);
	
	/**
	 * Returns the residue serial of the given node
	 * @param node
	 * @return
	 */
	public abstract int getResidueSerial(V node);

	protected void setSerials2NodesMap(TreeMap<Integer,V> serials2nodes) {
		this.serials2nodes = serials2nodes;
	}
	
	/**
	 * Returns the RIG/AIG Node object given a residue/atom serial
	 * @param serial
	 * @return
	 */
	public V getNodeFromSerial(int serial) {
		return serials2nodes.get(serial);
	}
	
	/**
	 * Returns the edge connecting the two nodes with the given serials.
	 * @param i
	 * @param j
	 * @return The edge between nodes with serials i and j or null if no such edge exists.
	 */
	public E getEdgeFromSerials(int i, int j) {
		return findEdge(this.getNodeFromSerial(i), this.getNodeFromSerial(j));
	}
	
	/**
	 * Returns true if the graph contains an edge between the nodes with serials i and j.
	 * @param i
	 * @param j
	 * @return True if the graph contains an edge from i to j, false otherwise.
	 */
	public boolean containsEdgeIJ(int i, int j) {
		return (getEdgeFromSerials(i,j) != null);
	}
	
	/**
	 * Returns true if the graph contains a vertex with serial i
	 * @param i
	 * @return
	 */
	public boolean containsVertexI(int i) {
		return this.containsVertex(getNodeFromSerial(i));
	}
	
	/**
	 * Returns all ordered residue serials (if RIGraph) or all atom serials (if AIGraph)
	 * @return
	 */
	public Set<Integer> getSerials() {
		return serials2nodes.keySet();
	}
	
	public String getPdbCode() {
		return pdbCode;
	}
	
	public void setPdbCode(String pdbCode) {
		this.pdbCode = pdbCode;
	}

	public String getChainCode() {
		return chainCode;
	}
	
	public void setChainCode(String chainCode) {
		this.chainCode = chainCode;
	}
	
	public String getPdbChainCode() {
		return pdbChainCode;
	}
	
	public void setPdbChainCode(String pdbChainCode) {
		this.pdbChainCode = pdbChainCode;
	}

	public int getModel() {
		return model;
	}
	
	public void setModel(int model) {
		this.model = model;
	}
	
	/**
	 * @return the caspModelNum
	 */
	public int getCaspModelNum() {
		return caspModelNum;
	}

	/**
	 * @return the groupNum
	 */
	public int getGroupNum() {
		return groupNum;
	}

	/**
	 * @return the targetNum
	 */
	public int getTargetNum() {
		return targetNum;
	}
	
	public String getAuthorStr() {
		return authorStr;
	}
	
	public String getMethodStr() {
		return methodStr;
	}
	
	public void setParents(String[] parents) {
		this.caspParents = parents;
	}
	
	public String[] getParents() {
		return this.caspParents;
	}

	/**
	 * @param caspModelNum the caspModelNum to set
	 */
	public void setCaspModelNum(int caspModelNum) {
		this.caspModelNum = caspModelNum;
	}

	/**
	 * @param groupNum the groupNum to set
	 */
	public void setGroupNum(int groupNum) {
		this.groupNum = groupNum;
	}

	/**
	 * @param targetNum the targetNum to set
	 */
	public void setTargetNum(int targetNum) {
		this.targetNum = targetNum;
	}
	
	/**
	 * @param authorStr the author string to set
	 */
	public void setAuthorStr(String authorStr) {
		this.authorStr = authorStr;
	}
	
	/**
	 * @param methodStr the method string to set
	 */
	public void setMethodStr(String methodStr) {
		this.methodStr = methodStr;
	}
	
	public String getSid() {
		return sid;
	}
	
	public void setSid(String sid) {
		this.sid = sid;
	}
	
	public boolean isRestrictedToScopDomain() {
		return sid!=null;
	}

	public int getFullLength() {
		return fullLength;
	}

	public String getSequence() {
		return sequence;
	}
	
	public void setSequence(String sequence) {
		this.sequence = sequence;
		this.fullLength = sequence.length();
	}
	
	public SecondaryStructure getSecondaryStructure() {
		return this.secondaryStructure;
	}
	
	public void setSecondaryStructure(SecondaryStructure secondaryStructure) {
		this.secondaryStructure = secondaryStructure;
	}
	
	public boolean hasSecondaryStructure() {
		return this.secondaryStructure!=null;
	}

	/**
	 * Sets the user comment for this graph.
	 * @param comment
	 */
	public void setComment(String comment) {
		this.comment = comment;
	}
	
	/**
	 * Returns the user comment.
	 * @return the user comment string
	 */
	public String getComment() {
		return this.comment;
	}
	
	/**
	 * True if this protein structure graph is directed
	 * @return
	 */
	public boolean isDirected() {
		// this is valid here because for our AIG/RIG graphs either all edges are directed or undirected
		if (!this.directed_edges.isEmpty()) return true;
		return false;
	}
	
	/**
	 * True if this ProtStructGraph has the sequence field set to not blank
	 * @return
	 */
	public boolean hasSequence() {
		if (sequence==null) return false;
		return !sequence.equals("");
	}

	/**
	 * Returns the minimum sequence separation filter last applied to this graph or NO_SEQ_SEP_VAL if none was applied.
	 * @return the minimum sequence separation or NO_SEQ_SEP_VAL
	 */
	public int getMinSeqSep() {
		return minSeqSep;
	}
	
	/**
	 * Returns the maximum sequence separation filter last applied to this graph or NO_SEQ_SEP_VAL if none was applied.
	 * @return the maximum sequence separation or NO_SEQ_SEP_VAL
	 */
	public int getMaxSeqSep() {
		return maxSeqSep;
	}
	
	
	/**
	 * Removes edges strictly above the given range (sequence distance).
	 * Does nothing if range is not positive.
	 * @param range
	 */
	public void restrictContactsToMaxRange(int range){
		if(range > 0) {
			for (E edge:this.getEdges()) {
				if (this.getContactRange(edge)>range) {
					this.removeEdge(edge);
				}
			}
			maxSeqSep = range;
		}
	}

	/**
	 * Removes edges strictly below the given range (sequence distance).
	 * Does nothing if range is not positive.
	 * @param range
	 */
	public void restrictContactsToMinRange(int range){
		if(range > 0) {
			for (E edge:this.getEdges()) {
				if (this.getContactRange(edge)<range) {
					this.removeEdge(edge);
				}
			}
			minSeqSep = range;
		}
	}
	
	/**
	 * Removes edges strictly within the same secondary structure element
	 * @param range
	 */
	public void restrictContactsBetweenSs() {
		for(E edge:this.getEdges()) {
			Pair<V> pair = this.getEndpoints(edge);
			SecStrucElement ss1 = secondaryStructure.getSecStrucElement(this.getResidueSerial(pair.getFirst()));
			SecStrucElement ss2 = secondaryStructure.getSecStrucElement(this.getResidueSerial(pair.getSecond()));
			if(ss1 != null && ss2 != null && ss1 == ss2) {
				this.removeEdge(edge);
			}
		}
		interSSE = true;		
	}

	/**
	 * Removes all edges from this ProtStructGraph
	 *
	 */
	public void removeAllEdges() {
		for (E edge:this.getEdges()){
			this.removeEdge(edge);
		}
	}
			
}

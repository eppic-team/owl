package proteinstructure;

import java.util.Set;
import java.util.TreeMap;

import edu.uci.ics.jung.graph.SimpleGraph;
import edu.uci.ics.jung.graph.SparseGraph;

/**
 * Class representing a protein structure graph, known extending 
 * classes are RIGraph and AIGraph
 * 
 * TODO we implement SimpleGraph intending to mark the graph as not 
 * 		accepting parallel/loop edges but that doesn't work:
 * 		SimpleGraph is simply a marker interface, we would have to 
 * 		implement ourselves a SparseGraph that doesn't accept 
 * 		parallel/loop edges
 *
 * @param <V>
 * @param <E>
 */
public abstract class ProtStructGraph<V,E> extends SparseGraph<V,E> implements SimpleGraph<V,E> {

	
	protected static final int DEFAULT_MODEL = 1;
	protected final static String GRAPHFILEFORMATVERSION = "1.0";
	

	protected String sequence;		// the full sequence (with unobserved residues and non-standard aas ='X')
	protected String pdbCode;		// the lower-case pdb code
	protected String pdbChainCode;	// The pdb chain code (upper case), i.e. the classic (author's) pdb code ("NULL" if it is blank in original pdb file)
	protected String chainCode;		// Our internal chain identifier (upper case)
									// - in reading from pdbase or from msdsd it will be set to the internal chain id (asym_id field for pdbase, pchain_id for msdsd)
    								// - in reading from pdb file it coincides with pdbChainCode except for "NULL" where we use "A"	
	protected int model;			// model serial number (for NMR structures, for all others is 1)
	protected int fullLength;		// full length of the protein with all aas (same as sequence length)

	protected int minSeqSep;
	protected int maxSeqSep;
	
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
		this.minSeqSep = -1;
		this.maxSeqSep = -1;
		this.secondaryStructure = null;
	}

	/**
	 * Returns the contact range (sequence separation) of the given edge
	 * @param edge
	 * @return
	 */
	public abstract int getContactRange(E edge);
	
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
	 * Removes edges strictly above the given range (sequence distance)
	 * @param range
	 */
	public void restrictContactsToMaxRange(int range){
		for (E edge:this.getEdges()) {
			if (this.getContactRange(edge)>range) {
				this.removeEdge(edge);
			}
		}
		maxSeqSep = range;
	}

	/**
	 * Removes edges strictly below the given range (sequence distance)
	 * @param range
	 */
	public void restrictContactsToMinRange(int range){
		for (E edge:this.getEdges()) {
			if (this.getContactRange(edge)<range) {
				this.removeEdge(edge);
			}
		}
		minSeqSep = range;
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

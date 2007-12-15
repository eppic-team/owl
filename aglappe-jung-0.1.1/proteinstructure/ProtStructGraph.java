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
	

	protected String sequence;		// the full sequence (with unobserved residues and non-standard aas ='X')
	protected String pdbCode;		// the lower-case pdb code
	protected String pdbChainCode;	// The pdb chain code (upper case), i.e. the classic (author's) pdb code ("NULL" if it is blank in original pdb file)
	protected String chainCode;		// Our internal chain identifier (upper case)
									// - in reading from pdbase or from msdsd it will be set to the internal chain id (asym_id field for pdbase, pchain_id for msdsd)
    								// - in reading from pdb file it coincides with pdbChainCode except for "NULL" where we use "A"	
	protected int model;			// model serial number (for NMR structures, for all others is 1)
	protected String sid;			// the scop id if this ProtStructGraph comes from a Pdb object restricted to a SCOP domain
	
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
	 * Returns the edge connecting the two nodes with the given serials.
	 * @param i
	 * @param j
	 * @return The edge between nodes with serials i and j or null if no such edge exists.
	 */
	public E getEdgeFromSerials(int i, int j) {
		return findEdge(this.getNodeFromSerial(i), this.getNodeFromSerial(j));
	}
	
	/**
	 * Returns true iff the graph contains an edge between the nodes with serials i and j.
	 * @param i
	 * @param j
	 * @return True if the graph contains an edge from i to j, false otherwise.
	 */
	public boolean containsEdgeIJ(int i, int j) {
		return (getEdgeFromSerials(i,j) != null);
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
	
	public String getSid() {
		return sid;
	}
	
	public void setSid(String sid) {
		this.sid = sid;
	}
	
	public boolean isRestrictedToScopDomain() {
		return sid!=null;
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
	
	/**
	 * Removes the given edge from this graph
	 * Overrides removeEdge from SparseGraph just to fix a bug (filed in bug tracker with request ID: 1844767)
	 * TODO: when bug is fixed in the next JUNG2 release get rid of this
	 * @return true if the edge is present and thus can be removed, false if the edge is not present
	 */
	@Override
    public boolean removeEdge(E edge)
    {
        if (!containsEdge(edge)) 
            return false;
        
        Pair<V> endpoints = getEndpoints(edge);
        V v1 = endpoints.getFirst();
        V v2 = endpoints.getSecond();
        
        // remove edge from incident vertices' adjacency maps
        if (getEdgeType(edge) == EdgeType.DIRECTED)
        {
            vertex_maps.get(v1)[OUTGOING].remove(v2);
            vertex_maps.get(v2)[INCOMING].remove(v1);
            directed_edges.remove(edge);
        }
        else
        {
            vertex_maps.get(v1)[INCIDENT].remove(v2);
            vertex_maps.get(v2)[INCIDENT].remove(v1);
            undirected_edges.remove(edge);
        }

        return true;
    }
		
}

package proteinstructure;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.io.IOException;
import java.util.TreeMap;
import java.util.HashMap;

/**
 * A residue interaction graph derived from a single chain pdb protein structure
 * 
 * @author 		Jose Duarte
 * Class:		Graph
 * Package:		proteinstructure
 */
public class Graph {

	public final static String GRAPHFILEFORMATVERSION = "1.0";
	
	public EdgeSet contacts; // we keep it public to be able to re-reference the object directly (getContacts() copies it)
	
	protected TreeMap<Integer,String> nodes; // nodes is a TreeMap of residue serials to residue types (3 letter code)
	protected String sequence; 				// the full sequence (with unobserved residues and non-standard aas ='X')
	protected String pdbCode;
	protected String chainCode;
	protected String pdbChainCode;
	protected double cutoff;
	protected String ct;					// the contact type
	protected boolean directed;
	
	// fullLength is length of full sequence or:
	// -if sequence not provided (when reading from db): length of everything except possible unobserved residues at end of chain
	// -if sequence and nodes not provided (when reading from file and sequence field missing): length except possible unobserved residues at end of chain and possible nodes without contacts at end of chain
	protected int fullLength; 
	protected int obsLength;  // length without unobserved, non standard aas 
	
	protected int numContacts;
	
	protected boolean modified;
	
	public Graph() {
		
	}
		
	/**
	 * Constructs Graph object by passing ArrayList with contacts and TreeMap with nodes (res serials and types)
	 * Must also pass contact type, cutoff, pdbCode and chainCode
	 * @param contacts
	 * @param nodes
	 * @param sequence
	 * @param cutoff
	 * @param ct
	 * @param pdbCode
	 * @param chainCode
	 */
	protected Graph (EdgeSet contacts, TreeMap<Integer,String> nodes, String sequence, double cutoff,String ct, String pdbCode, String chainCode, String pdbChainCode) {
		this.contacts=contacts;
		this.cutoff=cutoff;
		this.nodes=nodes;
		this.sequence=sequence;
		this.pdbCode=pdbCode;
		this.chainCode=chainCode;
		this.pdbChainCode=pdbChainCode;
		this.ct=ct;
		this.fullLength=sequence.length();
		this.obsLength=nodes.size();
		this.numContacts=contacts.size();
		this.modified=false;
		this.directed=false;
		if (ct.contains("/")){
			directed=true;
		}
		
		assert(this.pdbCode.equals(this.pdbCode.toLowerCase()));				// pdb codes should be always lower case 
		assert(this.pdbChainCode.equals(this.pdbChainCode.toUpperCase()));		// pdb chain codes should be always upper case
	}
	

	//TODO implement (from python) write_graph_to_db, do we really need it here??
		
	public void write_graph_to_file (String outfile) throws IOException {
		PrintStream Out = new PrintStream(new FileOutputStream(outfile));
		Out.println("#AGLAPPE GRAPH FILE ver: "+GRAPHFILEFORMATVERSION);
		Out.println("#SEQUENCE: "+sequence);
		Out.println("#PDB: "+pdbCode);
		Out.println("#PDB CHAIN CODE: "+pdbChainCode);
		Out.println("#CHAIN: "+chainCode);
		Out.println("#CT: "+ct);
		Out.println("#CUTOFF: "+cutoff);
		for (Edge pair:contacts){
			int i_resser=pair.i;
			int j_resser=pair.j;
			Out.println(i_resser+"\t"+j_resser);
		}
		Out.close();		
	}
	
	/**
	 * Gets list of contacts as a new EdgeSet (deep copied)
	 * 
	 */
	public EdgeSet getContacts(){
		EdgeSet newContacts = new EdgeSet();
		for (Edge cont:contacts){
			newContacts.add(new Edge(cont.i,cont.j));
		}
		return newContacts;
	}
	
	/**
	 * Gets TreeMap of nodes, deep copying  
	 * 
	 */
	public TreeMap<Integer,String> getNodes(){
		TreeMap<Integer,String> newNodes = new TreeMap<Integer,String>();
		for (int resser:nodes.keySet()){
			newNodes.put(resser, nodes.get(resser));
		}
		return newNodes;
	}
	
	/**
	 * Deep copies this Graph object returning new one
	 * @return
	 */
	public Graph copy(){
		return new Graph(getContacts(),getNodes(),sequence,cutoff,ct,pdbCode,chainCode,pdbChainCode);		
	}
	
	/**
	 * Gets a reference to this Graph deep copying contacts but re-referencing nodes
	 * @return
	 */
	public Graph copyKeepingNodes(){
		return new Graph(getContacts(),nodes,sequence,cutoff,ct,pdbCode,chainCode,pdbChainCode);		
	}	
	
	/**
	 * Returns an int matrix with 1s for contacts and 0s for non contacts, i.e. the contact map
	 * In non-crossed cases this should give us the upper half matrix (contacts are only j>i)
	 * In crossed cases this gives us a full matrix (contacts are both j>i and i>j since they are directed)
	 * @return
	 */
	public int[][] getIntMatrix(){
		// this initialises the matrix to 0 (i.e. no contact)
		int[][] cm = new int[fullLength][fullLength];
		// we put a 1 for all given contacts
		for (Edge cont:contacts){
			int i_resser = cont.i;
			int j_resser = cont.j;
			cm[i_resser-1][j_resser-1]=1;
		}
		return cm;
	}

	/**
	 * Gets a node's residue type given the residue serial
	 * @param resser
	 * @return
	 */
	public String getResType(int resser){
		return nodes.get(resser);
	}
	
	/**
	 * Gets node neighbourhood given a residue serial
	 * @param resser
	 * @return
	 */
	public NodeNbh getNodeNbh(int resser){
		NodeNbh nbh = new NodeNbh(resser, getResType(resser));
		//this could be implemented using the contact map matrix and scanning through 1 column/row
		//it would be just slightly faster, here we do 2*numContacts iterations, using matrix would be only fullLength iterations
		//however we would then have the overhead of creating the matrix
		for (Edge cont:contacts){
			if (cont.i==resser) nbh.put(cont.j, nodes.get(cont.j));
			if (cont.j==resser) nbh.put(cont.i, nodes.get(cont.i));
		}
		return nbh;
	}
	
	/**
	 * Gets edge neighbourhood (common neighbourhood) given a residue serial pair
	 * @param i_resser
	 * @param j_resser
	 * @return
	 */
	public EdgeNbh getEdgeNbh(int i_resser, int j_resser){
		EdgeNbh nbh = new EdgeNbh(i_resser, getResType(i_resser), j_resser, getResType(j_resser));
		NodeNbh i_nbhd = getNodeNbh(i_resser);
		NodeNbh j_nbhd = getNodeNbh(j_resser);
		if (j_nbhd.size()>=i_nbhd.size()) { //with this we will be slightly faster, always iterating through smallest TreeMap
			for (int resser:i_nbhd.keySet()) {
				if (j_nbhd.containsKey(resser)) nbh.put(resser, i_nbhd.get(resser));
			}
		} else {
			for (int resser:j_nbhd.keySet()) {
				if (i_nbhd.containsKey(resser)) nbh.put(resser, j_nbhd.get(resser));			
			}
		}
		return nbh;
	}
	
	/**
	 * Gets 2nd shell node neighbourhood
	 * @param resser
	 */
	public NodeNbh get2ndshellNodeNbh(int resser){
		// first we create a NodeNbh object for the second shell, central residue is given resser
		NodeNbh nbh2ndshell = new NodeNbh(resser,getResType(resser));
		// we get 1st neighbourhood
		NodeNbh nbh = this.getNodeNbh(resser);
		for (int nb:nbh.keySet()){
			NodeNbh nbh2 = this.getNodeNbh(nb); // for each first neighbour we take its neighbourhood
			for (int nb2:nbh2.keySet()){
				if (nb2!=resser && !nbh.containsKey(nb2)){ // if the 2nd neighbour nb2 is not the given resser or is not a 1st neighbour
					nbh2ndshell.put(nb2, getResType(nb2));
				}
			}
		}
		return nbh2ndshell;
	}
	
	public void addEdge(Edge cont){
		if (!directed && cont.i>cont.j){
			// we invert in case of undirected and i>j because in undirected we have only the half of the matrix j>i
			// if we added an edge i>j it could happen that the edge was already there but inverted and wouldn't be detected as a duplicate
			cont = new Edge(cont.j,cont.i); 
		}
		contacts.add(cont); // contacts is a TreeSet and thus takes care of duplicates
		int oldNumContacts = numContacts;
		numContacts=getNumContacts();
		// if number of contacts changed that means we actually added a new contact and thus we modified the graph
		if (numContacts!=oldNumContacts) modified=true;
		
	}
	
	public void delEdge(Edge cont){
		if (!directed && cont.i>cont.j){
			// we invert in case of undirected and i>j because in undirected we have only the half of the matrix j>i
			// if we try to delete an edge i>j it won't be there, we have to invert it and then try to delete
			cont = new Edge(cont.j,cont.i); 
		}
		contacts.remove(cont);
		int oldNumContacts = numContacts;
		numContacts=getNumContacts();
		// if number of contacts changed that means we actually added a new contact and thus we modified the graph
		if (numContacts!=oldNumContacts) modified=true;
	}
	
	public void restrictContactsToMaxRange(int range){
		EdgeSet edgesToDelete = new EdgeSet();
		for (Edge cont:contacts){
			if (cont.getRange()>range) edgesToDelete.add(cont);
		}
		for (Edge cont:edgesToDelete){
			delEdge(cont);
		}
	}
	
	public void restrictContactsToMinRange(int range){
		EdgeSet edgesToDelete = new EdgeSet();
		for (Edge cont:contacts){
			if (cont.getRange()<range) edgesToDelete.add(cont);
		}
		for (Edge cont:edgesToDelete){
			delEdge(cont);
		}
	}

	/**
	 * Returns a HashMap with all edge neighbourhood sizes (if they are >0) for each cell in the contact map
	 * @return
	 */
	public HashMap<Edge,Integer> getAllEdgeNbhSizes() {
		HashMap<Edge,Integer> sizes = new HashMap<Edge, Integer>();
		if (!directed) {
			for (int i=1; i<fullLength;i++){
				for (int j=i+1; j<fullLength;j++){
					int size = getEdgeNbh(i, j).size();
					if (size>0)	sizes.put(new Edge(i,j), size);
				}
			}			
		} else {
			for (int i=1; i<fullLength;i++){
				for (int j=1; j<fullLength;j++){
					if (i!=j){
						int size = getEdgeNbh(i, j).size();
						if (size>0) sizes.put(new Edge(i,j), size);
					}
				}
			}
		}
		return sizes;
	}

	//TODO not sure what kind of return we want, for now is a HashMap with three graph objects 
	public HashMap<String,Graph> compare(Graph other) throws Exception{
		//first check that other has same sequence than this, otherwise throw exception
		if (!this.sequence.equals(other.sequence)){
			//TODO throw specific exception
			throw new Exception("Sequence of 2 graphs to compare differ, can't compare them.");
		}
		EdgeSet common = new EdgeSet();
		EdgeSet onlythis = new EdgeSet();
		EdgeSet onlyother = new EdgeSet();
		for (Edge cont:this.contacts){
			if (other.contacts.contains(cont)) {
				common.add(cont);
			} else{
				onlythis.add(cont);
			}
		}
		for (Edge cont:other.contacts){
			if (!this.contacts.contains(cont)){
				onlyother.add(cont);
			}
		}
		Graph commongraph = new Graph (common,getNodes(),sequence,cutoff,ct,pdbCode,chainCode,pdbChainCode);
		Graph onlythisgraph = new Graph (onlythis,getNodes(),sequence,cutoff,ct,pdbCode,chainCode,pdbChainCode);
		Graph onlyothergraph = new Graph (onlyother,getNodes(),sequence,cutoff,ct,other.pdbCode,other.chainCode,other.pdbChainCode);
		HashMap<String,Graph> result = new HashMap<String,Graph>();
		result.put("common", commongraph);
		result.put("onlythis", onlythisgraph);
		result.put("onlyother",onlyothergraph);
		return result;
	}
	
	public boolean isModified(){
		return modified;
	}
	
	public boolean isDirected(){
		return directed;
	}
	
	public String getPdbCode() {
		return pdbCode;
	}
	
	public String getPdbChainCode(){
		return pdbChainCode;
	}
	
	public String getChainCode(){
		return chainCode;
	}
	
	public String getSequence(){
		return sequence;
	}
	
	public int getFullLength(){
		return fullLength;
	}
	
	public int getObsLength(){
		return obsLength;
	}
	
	public int getNumContacts(){
		// in theory we could return just numContacts, because we have taken care of updating it every time contacts changed
		// however we call directly contacts.size() as I feel is safer
		return contacts.size(); 
	}
	
	public String getContactType() {
		return ct;
	}
	
	public double getCutoff(){
		return cutoff;
	}
	
	public boolean containsContact(Edge cont){
		// be careful with order, this checks strictly whether the cont.i, cont.j is given, strictly in that order!
		// in undirected case contacts are stored only in 1 direction (j>i) and thus if wrong order given it won't be found 
		return contacts.contains(cont);
	}
	
	public void resetContacts(){
		this.contacts = new EdgeSet();
	}
}


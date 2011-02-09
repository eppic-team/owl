package owl.embed;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.LinkedHashMap;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.uci.ics.jung.graph.util.Pair;
import owl.embed.contactmaps.Individuals;

import owl.core.structure.AminoAcid;
import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbasePdb;
import owl.core.structure.graphs.RIGCommonNbhood;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.Goodies;
import owl.core.util.IntPairComparator;
import owl.core.util.IntPairSet;
import owl.core.util.MySQLConnection;
/**
 * Class to generate minimal subsets of a given contact map according to Sathyapriya et al. 2009. The basic idea
 * is, that all contacts with the greatest common neighborhood size and sequence separation build a minimal subset.
 * Since this class does not provide any public or protected constructors the only way to compute a minimal subset,
 * one has to invoke one of the static methods provided, which return a minimal subset as an RIGraph instance.
 * Example:
 * <p><tt>RIGraph subset = ConePeeler.getMinSubset(...);</tt>//overloaded method returning the minimal subset
 * @author gmueller
 *
 */
public class ConePeeler extends RIGraph {
	private static final long serialVersionUID = 1L;
	/** the default database*/
	private static final String DEF_DB = "pdbase_20090728";
	private static Pdb structure;
	/** the minimal sequence separation*/
	private static int seqSep = 3;
	/** the distance map, mapping all contacts onto their distance*/
	private HashMap<Pair<Integer>,Double> distMap;
	/** mapping all contacts onto their corresponding RIGEdge*/
	private HashMap<Pair<Integer>,RIGEdge> allEdgeMap;
	/** a set containing all edges the are not going to be deleted*/
	private HashSet<Pair<Integer>> protectedEdges;
	/** the minimal subset of all residues*/
	private TreeSet<Integer>       minSetRes;
	/** a set of all initial contacts*/
	private IntPairSet allEdgePairs;
	/**
	 * Zero parameter constructor
	 */
	ConePeeler (){
		super();
	}
	/**
	 * Five parameter constructor: initializes an instance of this class by loading a PDB file from the specified
	 * <tt>database</tt>. Additional parameters as <tt>chainCode</tt>, <tt>contactType</tt>, <tt>cutOff</tt> and
	 * <tt>seqSep</tt> (which is the minimal sequence separation) are required. Note, that the <tt>pdbCode</tt> must
	 * be in lower case PDB code.
	 * @param pdbCode the PBD code of the protein
	 * @param database the database to use
	 * @param chainCode the chain code of the given protein
	 * @param contactType the contact type, i.e. C-alpha = Ca, C-beta = Cb,...
	 * @param cutOff the cutoff distance
	 * @param seqSep the sequence separation, default: 3
	 * @throws ConePeelerException if either an SQLException, or an error due to the pdb code was thrown
	 */
	private ConePeeler (String pdbCode, String database, String chainCode, String contactType, double cutOff, int seqSep) throws ConePeelerException {
		super();
		setFields(pdbCode,database,chainCode,contactType,cutOff,seqSep);
	}
	/**
	 * Two parameter constructor: constructs a RIGraph instance from the <tt>pdb</tt> instance, with the specified
	 * <tt>cutOff</tt>.
	 * @param pdb a Pdb instance representing the protein
	 * @param cutOff the cutoff distance
	 * @throws PdbLoadError
	 */
	protected ConePeeler (Pdb pdb, double cutOff) throws PdbLoadError{
		super();
		structure = new Pdb(pdb.getSequence());
		structure.load(pdb.getChainCode());
		setRIGraph(structure.getRIGraph("Ca", cutOff),3);
	}
	/**
	 * One parameter constructor: uses the RIGraph instance <tt>rig</tt> to 
	 * initialize an object of this class.
	 * @param rig an RIGraph instance
	 */
	private ConePeeler (RIGraph rig){
		super();
		super.setChainCode(rig.getChainCode());
		super.setContactType(rig.getContactType());
		super.setCutoff(rig.getCutoff());
		super.setPdbChainCode(rig.getPdbChainCode());
		super.setPdbCode(rig.getPdbCode());
		String seq = rig.getSequence();
		super.setSequence(seq);
		super.setSerials2NodesMap();
		for(RIGEdge edge : rig.getEdges()){
			Pair<RIGNode> pairNodes = rig.getEndpoints(edge);
			RIGNode node1 = pairNodes.getFirst();
			RIGNode node2 = pairNodes.getSecond();
			addVertex(node1); addVertex(node2);
			addEdge(edge,pairNodes);
		}
		allEdgeMap = new HashMap<Pair<Integer>,RIGEdge>();
		protectedEdges = new HashSet<Pair<Integer>>();
	}
	/**
	 * Setter of the five parameter constructor: initializes an instance of this class by loading a PDB file from the specified
	 * <tt>database</tt>. Additional parameters as <tt>chainCode</tt>, <tt>contactType</tt>, <tt>cutOff</tt> and
	 * <tt>seqSep</tt> (which is the minimal sequence separation) are required. Note, that the <tt>pdbCode</tt> must
	 * be in lower case PDB code.
	 * @param pdbCode the PDB code of the protein
	 * @param database the database to use
	 * @param chainCode the chain code of the protein
	 * @param contactType the contact type
	 * @param cutOff the cutoff distance
	 * @param seqSep the sequence separation
	 * @throws ConePeelerException if either an SQLException, or an error due to the pdb code was thrown
	 */
	private void setFields (String pdbCode, String database, String chainCode, String contactType, double cutOff, int seqSep) throws ConePeelerException {
		try{
		MySQLConnection conn = new MySQLConnection();
		if(structure==null) {
			if(database==null||database.length()==0) structure = new PdbasePdb(pdbCode,DEF_DB,conn);
			else structure = new PdbasePdb(pdbCode,database,conn);
		}
		else if(!pdbCode.matches(structure.getPdbCode()))
			throw new ConePeelerException ("The pdb code of the constant structure field and the denoted pdb code do not match!");
		System.out.println("Loading strucutre of "+pdbCode+" from database...");
		structure.load(chainCode);
		setRIGraph(structure.getRIGraph(contactType, cutOff),3);
		conn.close();
		}
		catch(SQLException e){
			throw new ConePeelerException(e.getMessage());
		}
		catch(PdbCodeNotFoundException e){
			throw new ConePeelerException (e.getMessage());
		}
		catch(PdbLoadError e){
			throw new ConePeelerException (e.getMessage());
		}
	}
	/**
	 * Auxiliary setter: needed to pick all contacts of an RIGraph instance,
	 * below the the specified <tt>cutOff</tt> and to specify the sequence separation. 
	 * @param rig the RIGraph instance
	 * @param seq_Sep the sequence separation
	 */
	private void setRIGraph(RIGraph rig, int seq_Sep){
		if(rig!=null){
			seqSep = seq_Sep;
			allEdgeMap = new HashMap<Pair<Integer>,RIGEdge>();
			protectedEdges = new HashSet<Pair<Integer>>();
			setPdbCode(rig.getPdbCode());
			String seq = rig.getSequence();
			setSequence(seq);
			super.setChainCode(rig.getChainCode());
			super.setPdbChainCode(rig.getPdbChainCode());
			super.setContactType(rig.getContactType());
			double cutOff = rig.getCutoff();
			super.setCutoff(cutOff);
			super.setSerials2NodesMap();
			distMap    = structure.calcDistMatrix(rig.getContactType());
			for(Pair<Integer> pair: distMap.keySet()){
				if(distMap.get(pair).doubleValue()<=cutOff){
					int fVal = pair.getFirst().intValue(), sVal = pair.getSecond().intValue();
					RIGNode node1 = new RIGNode(fVal,AminoAcid.one2three(seq.charAt(fVal-1)));
					RIGNode node2 = new RIGNode(sVal,AminoAcid.one2three(seq.charAt(sVal-1)));
					addVertex(node1);addVertex(node2); addEdgeIJ(fVal,sVal);
				}
			}
		}
	}
	/**
	 * Auxiliary setter: clears all fields of this class.
	 */
	void clear (){
		allEdgePairs = null;
		allEdgeMap = null;
		distMap = null;
		minSetRes = null;
		protectedEdges = null;
	}
	/**
	 * The actual cone peeling algorithm. A LinkedHashMap, containing the number of
	 * neighbors to a given residue, is used to find the contacts with the highest
	 * connectivity. Note, that this method is automatically called by the invocation
	 * of {@link conePeeler()} on a ConePeeler object, so this method does not need to
	 * be directly invoked. Example:
	 * <p><tt>ConePeeler cone = new ConePeeler(...)</tt>//the cone peeler
	 * <p><tt>cone.conePeeler();						//peeling algorithm
	 * <p><tt>RIGraph subset  = cone.getSubset()</tt>	//returns the minimal subset as an RIGraph
	 * @param sortedMap a LinkedHashMap mapping each residue serial to the degree of
	 * neighbors
	 */
	private void conePeeler (LinkedHashMap<Integer,Integer> sortedMap){
		for(Integer i: sortedMap.keySet()){
			RIGNode node = getNodeFromSerial(i.intValue());
			Collection<RIGEdge> edgeSet = getIncidentEdges(node);
			TreeMap<Pair<Integer>,Integer> seqRange    = new TreeMap<Pair<Integer>,Integer>(new IntPairComparator());
			TreeMap<Pair<Integer>,Integer> cnbhSize    = new TreeMap<Pair<Integer>,Integer>(new IntPairComparator());
			TreeMap<Pair<Integer>,RIGCommonNbhood> cnhe = new TreeMap<Pair<Integer>,RIGCommonNbhood> (new IntPairComparator());
			for(RIGEdge edge : edgeSet){
				RIGNode first = getFirstNode(edge);
				RIGNode secon = getSecondNode(edge);
				Pair<Integer> pair = new Pair<Integer> (first.getResidueSerial(),secon.getResidueSerial());
				allEdgeMap.put(pair,edge);
				seqRange.put(pair, getContactRange(edge));
				cnhe.put(pair,getCommonNbhood(first,secon));
				if(!getCommonNbhood(first,secon).equals(null)) cnbhSize.put(pair, new Integer(getCommonNbhood(first,secon).size()));
				else cnbhSize.put(pair,new Integer(0));
			}
			LinkedHashMap<Pair<Integer>,Integer> edgesSorted = Goodies.sortMapByValue(cnbhSize, false);
			for(Pair<Integer> pairs : edgesSorted.keySet()){
				if(cnhe.containsKey(pairs)){
					for(Integer val : cnhe.get(pairs).keySet()){
						int seqRangeI = Math.abs(pairs.getFirst().intValue()-val.intValue());
						int seqRangeJ = Math.abs(pairs.getSecond().intValue()-val.intValue());
						if(pairs.size()>0&&(seqRangeI>seqSep||seqRangeJ>seqSep)) protectedEdges.add(pairs);
						if(seqRangeI<=seqSep||seqRangeJ<=seqSep){
							if(pairs.getFirst()>val){
								Pair<Integer> npair = new Pair<Integer>(val,pairs.getFirst());
								if(!protectedEdges.contains(npair)) allEdgePairs.remove(npair);
							}
							if(pairs.getFirst()<val){
								Pair<Integer> npair = new Pair<Integer>(pairs.getFirst(),val);
								if(!protectedEdges.contains(npair)) allEdgePairs.remove(npair);
							}
							if(pairs.getSecond()>val){
								Pair<Integer> npair = new Pair<Integer>(val,pairs.getSecond());
								if(!protectedEdges.contains(npair)) allEdgePairs.remove(npair);
							}
							if(pairs.getSecond()<val){
								Pair<Integer> npair = new Pair<Integer>(pairs.getSecond(),val);
								if(!protectedEdges.contains(npair)) allEdgePairs.remove(npair);
							}
						}
					}
				}
				if(edgesSorted.get(pairs).intValue()==0||seqRange.get(pairs).intValue()<=5) 
					allEdgePairs.remove(pairs);
			}
		}
	}
	/**
	 * The cone peeler method: generates a map of all residues, according to their number of
	 * next neighbors and sorts them. All residue pairs are used for the actual peeling. Note, that
	 * this method
	 */
	private void conePeeler (){
		TreeMap<Integer,Integer> nodeDeg = getDegree4AllNodes();
		LinkedHashMap<Integer,Integer> sortedNodes = getAllNodesSorted(nodeDeg);
		allEdgePairs = new IntPairSet();
		allEdgePairs.addAll(getAllEdgePairs());
		conePeeler(sortedNodes);
	}
	/**
	 * A method removing all contact pairs of the minimal subset and counts
	 * the number of both deleted and added contacts.
	 * @return the minimal subset as an RIGraph 
	 * @throws ConePeelerException if this method is called before the cone peeling took
	 * place
	 */
	private RIGraph getSubSet () throws ConePeelerException {
		if(allEdgeMap.size() >0){
		int sel = 0;
		int del = 0;
		RIGraph rig = copy();
		rig.removeAllEdges();
		minSetRes = new TreeSet<Integer>();
		for(Pair<Integer> pair: allEdgeMap.keySet()){
			if(allEdgePairs.contains(pair)){
				sel++;
				rig.addEdgeIJ(pair.getFirst().intValue(), pair.getSecond().intValue());
				minSetRes.add(pair.getFirst()); minSetRes.add(pair.getSecond());
			}
			if(!allEdgePairs.contains(pair)) del++;
		}
		System.out.println("The number of edges to be seleted are "+sel+", deleted are "+del+" from "+getEdgeCount());
		System.out.println("Missing Residues "+(rig.getFullLength()-rig.getObsLength()));
		System.out.println("No of residues in the minimal subset is "+minSetRes.size()+" "+minSetRes);
		return rig;
		}
		else throw new ConePeelerException ("This method can only be called after invocation of conePeeler() method!");
	}
	/**
	 * Iterates over all <tt>RIGNode</tt> entries of the superclass <tt>RIGraph</tt>, 
	 * computes the degree of each <tt>RIGNode</tt> and returns a TreeMap of each RIGNode
	 * serial mapping onto their degree.
	 * @return a TreeMap of all RIGNode entries mapping onto their degree
	 */
	private TreeMap<Integer,Integer> getDegree4AllNodes (){
		TreeMap<Integer,Integer> degOfAll = new TreeMap<Integer,Integer>();
		for(RIGNode node: getVertices()){
			degOfAll.put(new Integer(node.getResidueSerial()), new Integer(degree(node)));
		}
		return degOfAll;
	}
	/**
	 * Generates a LinkedHashMap of <tt>map</tt> according to the each value in a
	 * descending fashion.
	 * @param map a Map with Integer keys and values
	 * @return a sorted LinkedHashMap according to the value
	 */
	private LinkedHashMap<Integer,Integer> getAllNodesSorted (Map<Integer,Integer> map){
		return Goodies.sortMapByValue(map, false);
	}
	/**
	 * Iterates over all RIGNodes an adds them to a HashSet iff
	 * this RIGraph instance contains such edge.
	 * @return a HashSet of all RIGNode serials forming an edge
	 */
	private HashSet<Pair<Integer>> getAllEdgePairs (){
		HashSet<Pair<Integer>> edgeSet = new HashSet<Pair<Integer>>();
		for(RIGNode node1 : getVertices()){
			for(RIGNode node2 : getVertices()){
				int serial1 = node1.getResidueSerial(), serial2 = node2.getResidueSerial();
				if(containsEdgeIJ(serial1,serial2)){
					if(serial1<serial2) edgeSet.add(new Pair<Integer>(new Integer(serial1),new Integer(serial2)));
					else edgeSet.add(new Pair<Integer>(new Integer(serial2),new Integer(serial1)));
				}
			}
		}
		return edgeSet;
	}
	/**
	 * Returns the first RIGNode of <tt>edge</tt> of this instance.
	 * @param edge the RIGEdge instance (of which RIGNode is a vertex)
	 * @return the first RIGNode
	 */
	private RIGNode getFirstNode (RIGEdge edge){
		return getEndpoints(edge).getFirst();
	}
	/**
	 * Returns the second RIGNode of <tt>edge</tt> of this instance.
	 * @param edge the RIGEdge instance (of which RIGNode is a vertex)
	 * @return the second RIGNode
	 */
	private RIGNode getSecondNode (RIGEdge edge){
		return getEndpoints(edge).getSecond();
	}
	/**
	 * Static methods computing the minimal subsets according to the cone peeling algorithm. Takes an
	 * RIGraph instance <tt>graph</tt> as parameter and performs computation.
	 * @param graph a residue interaction graph object
	 * @return the minimal set as an RIGraph 
	 * @throws SQLException
	 * @throws PdbCodeNotFoundException
	 * @throws PdbLoadError
	 */
	public static RIGraph getMinSubset (RIGraph graph) throws ConePeelerException {
		ConePeeler cone = new ConePeeler (graph);
		cone.conePeeler();
		return cone.getSubSet();
	}
	/**
	 * Static methods computing the minimal subsets according to the cone peeling algorithm. Takes
	 * String instances as <tt>pdbCode</tt>, <tt>chainCode</tt>, <tt>db</tt> and an integer value
	 * for sequence separation as parameter and performs computation.
	 * @param pdbCode the PDB code
	 * @param chainCode the chain code
	 * @param db the database to use
	 * @param seqSep the sequence separation
	 * @return the minimal set as an RIGraph
	 */
	public static RIGraph getMinSubset (String pdbCode, String chainCode, String db, int seqSep){
		return getMinSubset(pdbCode,chainCode,db,9.0,seqSep);
	}
	/**
	 * Static methods computing the minimal subsets according to the cone peeling algorithm. Takes
	 * String instances as <tt>pdbCode</tt>, <tt>chainCode</tt>, <tt>db</tt>, as well as <tt>cufOff</tt>
	 * and an integer value for sequence separation as parameter and performs computation.
	 * @param pdbCode the PDB code
	 * @param chainCode the chain code
	 * @param db the database to use
	 * @param cutOff the cutoff distance in aangstroem
	 * @param seqSep the sequence separation
	 * @return the minimal set as an RIGraph
	 * @throws ConePeelerException
	 */
	public static RIGraph getMinSubset (String pdbCode, String chainCode, String db, double cutOff, int seqSep) throws ConePeelerException{
		ConePeeler peeler = new ConePeeler (pdbCode, db, chainCode,"Ca",cutOff,seqSep);
		peeler.conePeeler();
		return peeler.getSubSet();
	}
	/**
	 * Static methods computing the minimal subsets according to the cone peeling algorithm. Takes
	 * String instances as <tt>pdbCode</tt>, <tt>chainCode</tt> and <tt>db</tt> as parameter and
	 * performs computation.
	 * @param pdbCode the PDB code
	 * @param chainCode the chain code
	 * @param db the database to use
	 * @return the minimal set as an RIGraph
	 * @throws ConePeelerException
	 */
	public static RIGraph getMinSubset (String pdbCode, String chainCode, String db) throws ConePeelerException{
		return getMinSubset(pdbCode,chainCode,db,3);
	}
	public static void main(String[] args) throws SQLException, PdbCodeNotFoundException, PdbLoadError, IOException{
		RIGraph subset = getMinSubset("1bkr","A","pdbase_20090728",8.0,1);
		String output = "/project/StruPPi/gabriel/Arbeiten/ConePeeler/1bkr1bkr-31.indi";
		Individuals in2 = new Individuals (output);
		@SuppressWarnings("unused")
		RIGraph rig = getMinSubset((RIGraph) in2);
		subset.writeToFile(output);
		Individuals in = new Individuals(output);
		in.printToFile("/project/StruPPi/gabriel/Arbeiten/ConePeeler/","","1bkr");
		Individuals in1 = new Individuals("1bkr",5.0,true);
		in1.printToFile("/project/StruPPi/gabriel/Arbeiten/ConePeeler/","","1bkrRandom");
	}
	
	class ConePeelerException extends RuntimeException {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		public ConePeelerException (String message){
			super(message);
		}
	}

}

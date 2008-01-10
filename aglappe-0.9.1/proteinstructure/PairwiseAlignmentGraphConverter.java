package proteinstructure;

import java.io.IOException;
import java.sql.SQLException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.uci.ics.jung.graph.util.Pair;

import sadp.ContactMap;
import sadp.ContactMapConstructorError;
import sadp.SADP;

/**
 * Instances of this class create two so called aligned graphs given two input 
 * graphs and a non-crossing matching between the nodes of the graphs. We 
 * denote with G1={V1,E1} the first input graph and with G2={V2,E2} the second 
 * input graph where Vx is the set of nodes and Ex the set of edges (x={1,2}) 
 * of the respective graph. The node matching M=V1 x V2. We denote the weight 
 * of an arbitrary edge e with w(e). Each output graph has the following 
 * properties:
 * <ul>
 * <li>The number of <b>nodes</b> is the same in both output graphs.</li>
 * <li>The number of <b>edges</b> of any output graph corresponds to the number 
 *  edges in the respective input graph.</li>
 * <li>Let i,j in V1, i &lt j, and k,l in V2, k &lt l. If (i,j) in E1 and (k,l) in E2 
 *  and (i,l) in M and (j,k) in M then w(i',j')=w(k',l')=w(i,j)+w(k,l) where 
 *  i',j',k',l' denote the new indices of the original nodes in the output 
 *  graphs.</li>
 * </ul>
 * @author Lars Petzold
 *
 */

public class PairwiseAlignmentGraphConverter {

	private RIGraph   alignedG1          = null;
	private RIGraph   alignedG2          = null;
	private HashSet<Pair<Integer>> commonEdges        = new HashSet<Pair<Integer>>();
	private boolean isCommonEdgeFilled = false;

	/**
	 * Constructs two aligned graphs based on the given matching pointed to by
	 * the passed iterator. This constructor requires the graphs to hold 
	 * sequence information.
	 * @param it  iterator pointing to a matching
	 * @param g1  first graph
	 * @param g2  second graph
	 * @param fi  indicates the index of the first node in the matching (note:
	 *             the first node does not necessarily need to be referenced in
	 *             the matching. Do not confuse this with the value of the lowest
	 *             node index in the matching which could be greater, though
	 *             never less!)
	 */
	public PairwiseAlignmentGraphConverter( Iterator<Pair<Integer>> it, RIGraph g1, RIGraph g2, int fi )
	throws AlignmentConstructionError {
		this( new PairwiseAlignmentConverter(it,g1.getSequence(),g2.getSequence(),"1","2",fi).getAlignment(), "1", "2", g1, g2 );
	}

	/**
	 * 
	 * @param a
	 * @param tag1
	 * @param tag2
	 * @param g1
	 * @param g2
	 */
	public PairwiseAlignmentGraphConverter( Alignment a, String tag1, String tag2, RIGraph g1, RIGraph g2 ) {
		alignedG1 = convertGraph(a,tag1,tag2,g1,g2);
		isCommonEdgeFilled = true;
		alignedG2 = convertGraph(a,tag2,tag1,g2,g1);
	}

	/**
	 * Creates aligned graph for g1 given g2.
	 * 
	 * @param a     alignment between the two graphs
	 * @param tag1  tag-name for g1 in the alignment
	 * @param tag2  tag-name for g2 in the alignment
	 * @param g1    graph to be converted
	 * @param g2    graph g1 is aligned to
	 */
	private RIGraph convertGraph( Alignment a, String tag1, String tag2, RIGraph g1, RIGraph g2 ) {
		String                  sequence = a.getAlignedSequence(tag1);

		RIGraph alignedGraph = new RIGraph();
		alignedGraph.setPdbCode(g1.getPdbCode());
		alignedGraph.setPdbChainCode(g1.getPdbChainCode());
		alignedGraph.setChainCode(g1.getChainCode());
		alignedGraph.setModel(g1.getModel());
		alignedGraph.setGroupNum(g1.getGroupNum());
		alignedGraph.setTargetNum(g1.getTargetNum());
		alignedGraph.setCaspModelNum(g1.getCaspModelNum());
		alignedGraph.setContactType(g1.getContactType());
		alignedGraph.setCutoff(g1.getCutoff());
		alignedGraph.setSequence(sequence);
		
		// secondary structure
		SecondaryStructure secStruct = convertSecondaryStruct(a, g1.getSecondaryStructure(), tag1);
		alignedGraph.setSecondaryStructure(secStruct);
		
		// adding nodes
		TreeMap<Integer,RIGNode> serials2nodes = new TreeMap<Integer,RIGNode>();
		for (RIGNode node: g1.getVertices()) {
			// mapped node, +1 as the node numbering in the graph 
			// corresponds to the one found in the pdb file which  
			// starts with 1
			int alignedResser = a.seq2al(tag1,node.getResidueSerial()) + 1;
			RIGNode alignedNode = new RIGNode(alignedResser,
					node.getResidueType(),
					secStruct==null?null:secStruct.getSecStrucElement(alignedResser));
			alignedGraph.addVertex(alignedNode);
			serials2nodes.put(a.seq2al(tag1,node.getResidueSerial())+1,alignedNode);
		}
		// and then gap nodes
		for (int resser=1;resser<=sequence.length();resser++) {
			if (!serials2nodes.containsKey(resser)) {
				RIGNode alignedNode = new RIGNode(resser,AAinfo.getGapCharacterThreeLetter());
				alignedGraph.addVertex(alignedNode);
				serials2nodes.put(resser, alignedNode);
			}
		}
		alignedGraph.setSerials2NodesMap(serials2nodes);
		
		// now edges
		for (RIGEdge edge: g1.getEdges()) { 
			Pair<RIGNode> pair = g1.getEndpoints(edge);
			// positions in the gapped sequence for tag1
			int i1 = a.seq2al(tag1,pair.getFirst().getResidueSerial());
			int j1 = a.seq2al(tag1,pair.getSecond().getResidueSerial());
			double weight = edge.getWeight();
			
			// positions in the ungapped sequence for tag2
			int p2  = a.al2seq(tag2,i1);
			int q2  = a.al2seq(tag2,j1);
			RIGEdge edge2;
			if( !(p2 == -1 || q2 == -1) ) { 
				if ((edge2 = g2.findEdge(g2.getNodeFromSerial(p2),g2.getNodeFromSerial(q2)))!=null) {
					// conserved contacts get the sum of weights of both source contacts
					weight += edge2.getWeight();
				}
			}
			if( !isCommonEdgeFilled ) {
				commonEdges.add(new Pair<Integer>(i1+1,j1+1));
			}
			// mapped Edge, +1 as the node numbering in the graph 
			// corresponds to the one found in the pdb file which 
			// starts with 1
			alignedGraph.addEdge(new RIGEdge(weight), alignedGraph.getNodeFromSerial(i1+1), alignedGraph.getNodeFromSerial(j1+1), g1.getEdgeType(edge));
		}
		
		return alignedGraph;
	}

	/**
	 * Gets the aligned graph corresponding to the first graph passed to the 
	 * constructor.
	 * @return a graph 
	 */
	public RIGraph getFirstGraph() {
		return alignedG1;
	}

	/**
	 * Gets the aligned graph corresponding to the first graph passed to the 
	 * constructor.
	 * @return a graph 
	 */
	public RIGraph getSecondGraph() {
		return alignedG2;
	}

	/**
	 * Gets common edge in the resulting graphs.
	 * @return an edge-set containing common edge only
	 * */
	public HashSet<Pair<Integer>> getCommonEdges() {
		return commonEdges;
	}

	/**
	 * @param args
	 * @throws PdbLoadError
	 * @throws SQLException 
	 * @throws ContactMapConstructorError 
	 * @throws IOException 
	 * @throws PdbCodeNotFoundError 
	 */
	public static void main(String[] args) throws PdbLoadError, SQLException, ContactMapConstructorError, IOException, PdbCodeNotFoundError {

		String pdbcode1="1bxy";
		String chaincode1="A";
		String pdbcode2="1sha";
		String chaincode2="A";

		Pdb pdb1 = new PdbasePdb(pdbcode1);
		pdb1.load(chaincode1);
		Pdb pdb2 = new PdbasePdb(pdbcode2);
		pdb2.load(chaincode2);

		RIGraph g1 = pdb1.get_graph("ALL", 4.2);
		RIGraph g2 = pdb2.get_graph("ALL", 4.2);

		g1.writeToSADPFile("1bxy.sadp");
		g2.writeToSADPFile("1sha.sadp");
		
		ContactMap x = new ContactMap(g1);
		ContactMap y = new ContactMap(g2);
		
		SADP sadp = new SADP(x,y);
		sadp.run();

		TreeSet<Pair<Integer>> matching = sadp.getMatching();

		PairwiseAlignmentGraphConverter pac = null;

		try {
			pac = new PairwiseAlignmentGraphConverter(matching.iterator(),g1,g2,0);
		} catch(Exception e) {
			System.err.println(e.getMessage());
			System.exit(-1);
		}

		try {
			pac.getFirstGraph().write_graph_to_file(pdbcode1+".cm-aglappe");	    
		} catch (IOException e) {
			System.out.println("Error: Writing aligned graph g1 to file failed!");
		}

		try {
			pac.getSecondGraph().write_graph_to_file(pdbcode2+".cm-aglappe");	    
		} catch (IOException e) {
			System.out.println("Error: Writing aligned graph g1 to file failed!");
		}

		//System.out.println("Common edges: "+pac.getCommonEdges().toString());
	}

	/**
	 * Creates a modified copy of a secondary structure object where the residue serials are mapped through an alignment.
	 * @param a the alignment mapping the original serials to aligned serials
	 * @param secStruct the original secondary structure object
	 * @param tag the id in the alignment of the sequence which the secondary structure refers to
	 * @return the converted secondary structure object or null if the original sec. struc. object was null
	 */
	public static SecondaryStructure convertSecondaryStruct(Alignment a, SecondaryStructure secStruct, String tag) {
		SecondaryStructure newSecStruct = null;
		if (secStruct!=null) {
			newSecStruct = new SecondaryStructure();
			Iterator<SecStrucElement> it = secStruct.getIterator();
			while (it.hasNext()) {
				SecStrucElement oldSselem = it.next();
				SecStrucElement sselem = new SecStrucElement(oldSselem.getType(),
						a.seq2al(tag,oldSselem.getInterval().beg) + 1,
						a.seq2al(tag,oldSselem.getInterval().end) + 1,
						oldSselem.getId());
				newSecStruct.add(sselem);
			}
		}
		return newSecStruct;
	}
}

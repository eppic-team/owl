package proteinstructure;

import java.io.IOException;
import java.util.Iterator;
import java.util.TreeMap;

import sadp.ContactMap;
import sadp.IOUtil;
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

    private Graph   alignedG1          = null;
    private Graph   alignedG2          = null;
    private EdgeSet commonEdges        = new EdgeSet();
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
    public PairwiseAlignmentGraphConverter( Iterator<Edge> it, Graph g1, Graph g2, int fi )
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
    public PairwiseAlignmentGraphConverter( Alignment a, String tag1, String tag2, Graph g1, Graph g2 ) {
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
    private Graph convertGraph( Alignment a, String tag1, String tag2, Graph g1, Graph g2 ) {
	TreeMap<Integer,String> nodes    = new TreeMap<Integer, String>();
	EdgeSet                 contacts = new EdgeSet();
	String                  sequence = a.getAlignedSequence(tag1);
	
	Iterator<Edge> it = g1.getContactIterator();
	Edge oE1=new Edge(0,0); // holds 'o'riginal 'E'dge in g1
	int i1=0,j1=0; // positions in the gapped sequence for tag1
	int p2=0,q2=0; // positions in the ungapped sequence for tag2
	
	// make the contacts
	while( it.hasNext() ) {
	    oE1 = it.next();
	    i1  = a.seq2al(tag1,oE1.i); 
	    j1  = a.seq2al(tag1,oE1.j);
	    p2  = a.al2seq(tag2,i1);
	    q2  = a.al2seq(tag2,j1);
	    
	    // 'm'apped 'E'dge, +1 as the node numbering in the graph 
	    // corresponds to the one found in the pdb file which usually 
	    // starts with 1
	    Edge mE1 = new Edge(i1+1,j1+1);
	    
	    // add contact with the contact-anchors being mapped to the aligned sequence
	    contacts.add(mE1);
	    
	    if( !(p2 == -1 || q2 == -1) ) { 
		// 'o'riginal 'E'dge in g'2' 
		Edge oE2 = new Edge(p2,q2);
		
		// assign weight to 'mE' in 'contacts'
		if( g2.hasContact(oE2) ) {
		    // conserved contacts get the sum of weights of both source contacts
		    mE1.weight = oE1.weight + oE2.weight;
		    if( !isCommonEdgeFilled ) {
			commonEdges.add(mE1);
		    }
		} else {
		    // unconserved contact get the weight of original contact in g1 only
		    mE1.weight = oE1.weight;
		}
	    } else {
		mE1.weight = oE1.weight; 
	    }
	}
			
	// make the nodes
	for( int i=0; i<sequence.length(); ++i ) {
	    if( sequence.charAt(i) == '-' ) {
		nodes.put(i+1, AAinfo.getGapCharacterThreeLetter());
	    } else {
		nodes.put( i+1, AAinfo.oneletter2threeletter(sequence.substring(i,i+1)) );
	    }
	    
	}
	
	return new Graph(contacts, nodes, sequence, g1.getCutoff(), g1.getContactType(), g1.getPdbCode(), g1.getChainCode(), g1.getPdbChainCode(), g1.getModel(), null);
    }
    
    /**
     * Gets the aligned graph corresponding to the first graph passed to the 
     * constructor.
     * @return a graph 
     */
    public Graph getFirstGraph() {
	return alignedG1;
    }
    
    /**
     * Gets the aligned graph corresponding to the first graph passed to the 
     * constructor.
     * @return a graph 
     */
    public Graph getSecondGraph() {
	return alignedG2;
    }
    
    /**
     * Gets common edge in the resulting graphs.
     * @return an edge-set containing common edge only
     * */
    public EdgeSet getCommonEdges() {
	return commonEdges;
    }
    
    /**
     * @param args
     */
    public static void main(String[] args) {
	
	ContactMap x = IOUtil.read(args[0]);
	ContactMap y = IOUtil.read(args[1]);
	
	SADP sadp = new SADP(x,y);
	sadp.run();
	
	EdgeSet matching = sadp.getMatching();
	Graph g1 = new Graph(x,null,6.0,AAinfo.getAllContactTypes().iterator().next(),null,null,null,1,null);
	Graph g2 = new Graph(y,null,6.0,AAinfo.getAllContactTypes().iterator().next(),null,null,null,1,null);

	PairwiseAlignmentGraphConverter pac = null;
	
	try {
	    pac = new PairwiseAlignmentGraphConverter(matching.iterator(),g1,g2,0);
	} catch(Exception e) {
	    System.err.println(e.getMessage());
	    System.exit(-1);
	}
	
	try {
	    pac.getFirstGraph().write_graph_to_file(args[0]+".cm-aglappe");	    
	} catch (IOException e) {
	    System.out.println("Error: Writing aligned graph g1 to file failed!");
	}

	try {
	    pac.getSecondGraph().write_graph_to_file(args[1]+".cm-aglappe");	    
	} catch (IOException e) {
	    System.out.println("Error: Writing aligned graph g1 to file failed!");
	}
	
	//System.out.println("Common edges: "+pac.getCommonEdges().toString());
    }

}

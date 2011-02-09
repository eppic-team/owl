package owl.embed.contactmaps;



import java.sql.SQLException;
import java.util.*;

import owl.core.structure.*;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.MySQLConnection;
import owl.embed.SparseMatrix;

import edu.uci.ics.jung.graph.util.Pair;

public class ContactMatrix extends owl.core.structure.graphs.RIGraph {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	SparseMatrix matrix;
	SparseMatrix matrix2;
	HashMap<Integer,HashSet<Pair<Integer>>> conePeeled;
	private int peeledCounter;
	
	public ContactMatrix (String pdb_code) throws SQLException, PdbCodeNotFoundException, PdbLoadError{
		super();
		setFields(pdb_code);
	}
	
	public void setFields (String pdb_code) throws SQLException, PdbCodeNotFoundException, PdbLoadError{
		MySQLConnection conn = new MySQLConnection ();
		Pdb pdb = new PdbasePdb(pdb_code, "pdbase_20090728", conn);
		pdb.load("A");
		RIGraph rig = pdb.getRIGraph("Ca",9.0);
		setPdbCode(rig.getPdbCode());
		setSequence(rig.getSequence());
		setChainCode(rig.getChainCode());
		setCutoff(9.0);
		setContactType(rig.getContactType());
		setSerials2NodesMap();
		setGraph(pdb,9.0);
		setMatrices();
		conn.close();
	}
	/**
	 * Initial setter: initializes the contact matrix and the cone peeler
	 * fields.
	 */
	private void setMatrices(){
		HashMap<Pair<Integer>,Integer> nbhMap = getAllCommonNbhSizes();
		Set<Pair<Integer>> keys = nbhMap.keySet();
		Iterator<Pair<Integer>> it = keys.iterator();
		int size = (int) (Math.log((double) nbhMap.size())/Math.log(2.0));
		int power = (int) (Math.pow(2.0,size+1)), length = getFullLength();
		HashMap<Pair<Integer>,Double> map = new HashMap<Pair<Integer>,Double>(power);
		conePeeled = new HashMap<Integer,HashSet<Pair<Integer>>> ();
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			int fVal = pair.getFirst().intValue()-1, sVal = pair.getSecond().intValue()-1;
			Pair<Integer> npair1 = new Pair<Integer>(new Integer(fVal),new Integer(sVal));
			Pair<Integer> npair2 = new Pair<Integer>(new Integer(sVal),new Integer(fVal));
			Integer val = nbhMap.get(pair);
			Double value = new Double((int) val.intValue());
			map.put(npair1, value);
			map.put(npair2, value);
			if(conePeeled.containsKey(val)){
				HashSet<Pair<Integer>> subset = conePeeled.get(val); 
				subset.add(npair1);
				conePeeled.put(val, subset);
				peeledCounter++;
			}
			else{
				HashSet<Pair<Integer>> subset = new HashSet<Pair<Integer>>(2);
				subset.add(npair1);
				conePeeled.put(val,subset);
				peeledCounter++;
			}
		}
		matrix  = (new SparseMatrix(map,length,length)).add(ContactMap.addBackbone(length));
		matrix2 = matrix.multiply(matrix);
	}
	/**
	 * Initial setter: initializes the vertex and edge map of the RIGraph superclass.
	 * @param pdb the Pdb instance
	 * @param cutOff the cut off distance
	 */
	public void setGraph(Pdb pdb, double cutOff){
		String seq = pdb.getSequence();
		HashMap<Pair<Integer>,Double> distMap = pdb.calcAtomDistMatrix("Ca");
		Set<Pair<Integer>> keys = distMap.keySet();
		Iterator<Pair<Integer>> it = keys.iterator();
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			Double dist = distMap.get(pair);
			if(dist.doubleValue()<=cutOff){
				Integer fVal = pair.getFirst(), sVal = pair.getSecond();
				int first = pdb.getResSerFromAtomSer(fVal.intValue()), secon = pdb.getResSerFromAtomSer(sVal.intValue());
				Character c1  = new Character(seq.charAt(first-1)), c2 = new Character(seq.charAt(secon-1));				
				addVertex(new RIGNode(first,AminoAcid.one2three(c1)));
				addVertex(new RIGNode(secon,AminoAcid.one2three(c2)));
				addEdgeIJ(first,secon);
			}
		}
	}
	/**
	 * Method, to return all keys (integers) of the cone peeler field
	 * in a ascending order. The key values correspond to the number
	 * of common neighbors of a given contact <tt>(i,j)</tt>.
	 * @return an array of Integers representing the number of common neighbors
	 * @throws NullPointerException if the cone peeler field is not initialized before calling 
	 */
	public Integer[] getConePeelerKeys() throws NullPointerException {
		if(conePeeled!=null){
			Set<Integer> keys = conePeeled.keySet();
			Integer[] frqArray = new Integer[conePeeled.size()];
			Iterator<Integer> it = keys.iterator();
			int counter = 0;
			while(it.hasNext()){
				frqArray[counter] = it.next();
				counter++;
			}
			Arrays.sort(frqArray);
			return frqArray;
		}
		else throw new NullPointerException ("Cone peeler must be initialized before calling this method!");
	}
	/**
	 * Method, returning all index pairs as an array of all native contacts present in the given
	 * protein. <b>Note</b>, that the index pairs are ordered in a ascending fashion, according to 
	 * the number of common neighbors. So, the first entry as the smallest common neighborhood, whereas
	 * the last entry has the greatest common neighborhood.
	 * @return an array of pairs of integers <tt>(i,j)</tt> representing a contact between
	 * amino acid <tt>i</tt> and <tt>j</tt>
	 */
	@SuppressWarnings("unchecked")
	public Pair<Integer>[] getOrderedIndexPairs (){
		Pair<Integer>[] indArray = new Pair[peeledCounter];
		Integer[] frq = (conePeeled.keySet()).toArray(new Integer[1]);
		Arrays.sort(frq);
		for(int i = 0; i < frq.length; i++){
			HashSet<Pair<Integer>> sub = conePeeled.get(new Integer((int) frq[i].doubleValue()));
			Iterator<Pair<Integer>> it = sub.iterator();
			int addCounter = 0;
			while(it.hasNext()){
				indArray[i+addCounter] = it.next();
				addCounter++;
			}
		}
		return indArray;
	}
	/**
	 * Method, returning all index pairs above the <tt>threshold</tt> as an array of all native contacts present in the given
	 * protein. <b>Note</b>, that the index pairs are ordered in a ascending fashion, according to 
	 * the number of common neighbors. So, the first entry as the smallest common neighborhood equal to <tt>threshold</tt>, whereas
	 * the last entry has the greatest common neighborhood.
	 * @param threshold
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public Pair<Integer>[] getOrderedIndexPairs (int threshold){
		Integer[] peeled = getConePeelerKeys();
		int counter = 0, length = peeled.length;
		while(peeled[counter].intValue()<threshold){
			counter++;
		}
		int counter1 = counter, counter2 = 0;
		while(counter<length){
			counter2+=conePeeled.get(peeled[counter]).size();
			counter++;
		}
		Pair<Integer>[] array = new Pair[counter2];
		counter = 0;
		for(int i = 0; i < length-counter1; i++){
			HashSet<Pair<Integer>> sub = conePeeled.get(peeled[i+counter1]);
			Iterator<Pair<Integer>> it = sub.iterator();
			while(it.hasNext()){
				array[counter]=it.next();
				counter++;
			}
		}
		return array;
	}

	public static void main (String[] args) throws SQLException, PdbCodeNotFoundException, PdbLoadError{
		ContactMatrix cm = new ContactMatrix("1bkr");
		Pair<Integer>[] ar = cm.getOrderedIndexPairs(10);
		for(int i = 0; i < ar.length; i++){
			System.out.println(ar[i].toString());
		}
	}
}

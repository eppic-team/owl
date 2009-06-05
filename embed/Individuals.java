package embed;

import proteinstructure.*;

import java.io.*;
import java.sql.SQLException;
import java.util.*;

//import edu.uci.ics.jung.*;
import edu.uci.ics.jung.graph.util.Pair;

import Jama.Matrix;

/*import embed.BoundsSmoother.BoundsDigraphNode;
import embed.BoundsSmoother.SimpleEdge;*/
import tools.*;
//import org.jgap.Chromosome;

/**
 * this class generates a random subset to a given number of contacts and provides methods to calculate the two 
 * different error estimation methods: DMError() and CMError()
 */
public class Individuals extends HashSet<Pair<Integer>> {
	
	/**
	 * field
	 */
	private static final long serialVersionUID = 1L;
	
	/**
	 * field: contact type - default: C-alpha
	 */
	private final String ct = "Ca";
	
	/**
	 * field: cutoff distance, set to default 9 angstroem
	 */
	private final double di = 9.0;

	/**
	 * field: array of integers representing contact pairs
	 */
	private int [][] entries;
	
	/**
	 * field
	 */
	protected boolean fifty;
	
	/**
	 * field: same as 'entries', in this case all contact pairs are stored in a HashSet, to avoid redundancy
	 */
	private HashSet<Pair<Integer>> store;
	
	/**
	 * field: number of contacts
	 */
	private int numOfContacts;
	
	/**
	 * String fields: 'name' = PDB code, 'sequence' = protein sequence and 'chainCode' = chain identifier
	 */
	private String name, sequence, chainCode;
	
	private double CMError, DMError;
	
	/*----------------------Constructors--------------------------------------------------*/
	/**
	 * zeror parameter constructor, defines the default values of all non final fields
	 */
	public Individuals () {
		this.chainCode = " ";
		this.CMError = 0.0;
		this.DMError = 0.0;
		this.entries = new int[1][2];
		this.name = " ";
		this.sequence = " ";
		this.storer();
		this.numOfContacts = 0;
	}
	
	/**
	 * one parameter constructor, reads instance of Individuals as input parameter
	 * @param in
	 */
	public Individuals (Individuals in) {
		this.setIndis(in);
	}
	
	/**
	 * one parameter constructor, reads array of integers as input parameter
	 * @param indices
	 */
	public Individuals (int[][] indices) {
		this.setIndis(indices);
	}
	
	/**
	 * one parameter constructor, reads instance of RIGraph as input parameter
	 * @param rig
	 * @throws Exception
	 */
	public Individuals (RIGraph rig) {
		this.setIndis(rig);
	}
	
	/**
	 * three parameter constructor, reads instances of RIGraph, MySQLConnection and a boolean expression
	 * as input parameter, the boolean parameter defines whether (true) or not (false) random sampling is needed
	 * @param rig
	 * @param conn
	 * @param rand
	 * @throws PdbLoadError 
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundError 
	 * @throws NullPointerException 
	 * @throws ArrayIndexOutOfBoundsException 
	 * @throws Exception 
	 * @throws Exception
	 */
	public Individuals (RIGraph rig, MySQLConnection conn, boolean rand) throws ArrayIndexOutOfBoundsException, NullPointerException, PdbCodeNotFoundError, SQLException, PdbLoadError  {
		this.setIndis(rig, conn, rand);
	}
	
	/**
	 * four parameter constructor, can either randomly sample a specified number of contacts or takes a specified number
	 * of contacts, the boolean parameter defines whether (true) or not (false) random sampling is needed
	 * @param rig
	 * @param conn
	 * @param rand
	 * @param val
	 * @throws ArrayIndexOutOfBoundsException
	 * @throws NullPointerException
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public Individuals (RIGraph rig, MySQLConnection conn, boolean rand, int val) throws ArrayIndexOutOfBoundsException, NullPointerException, PdbCodeNotFoundError, SQLException, PdbLoadError  {
		this.setIndis(rig, conn, rand, val);
	}
	
	/**
	 * four parameter constructor, can either randomly sample a specified number of contacts or takes a specified number
	 * of contacts, the boolean parameter defines whether (true) or not (false) random sampling is needed 
	 * @param rig
	 * @param conn
	 * @param rand
	 * @param val
	 * @throws ArrayIndexOutOfBoundsException
	 * @throws NullPointerException
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public Individuals (RIGraph rig, MySQLConnection conn, boolean rand, double val) throws ArrayIndexOutOfBoundsException, NullPointerException, PdbCodeNotFoundError, SQLException, PdbLoadError  {
		this.setIndis(rig, conn, rand, val);
	}
	
	/**
	 * three parameter constructor, reads instances of Bound[][], String and a double array
	 * as input parameter
	 * @param bound
	 * @param i
	 * @param dm
	 */
	public Individuals (Bound[][] bound, String i, double[][] dm) {
		this.setIndis(bound, i, dm);
	}

	/*--------------------setters---------------------------------------------------------*/
	
	/**
	 * setter, copies all instance variables of input parameter 'in' as an instance of Individuals
	 */
	public void setIndis (Individuals in) {
		this.setIndis(in.getEntries());
		this.setName(in.getName());
		this.store = in.getHashSet();
		this.CMError = in.getCM();
		this.DMError = in.getDM();
		this.numOfContacts = in.getNumOfContacts();
		this.sequence = in.getSequence();
		this.chainCode = in.getChainCode();
	}
	
	/**
	 * setter, reads an array of integers and uses them for the 'entries' instance variable
	 * @param indices
	 */
	public void setIndis (int [][] indices) {
		int dim = indices.length;
		int k = 0;
		while(k < dim && (indices[k][0] != 0 || indices[k][1] != 0)){
			k++;
		}
		this.entries = new int[k][2];
		int counter1 = 0;
		int counter2 = 0;
		for(int i = 0; i < dim; i++){
			//for(int j = i + 1; j < dim; j++){
				if(indices[i][0] != indices[i][1] - 1){
					this.entries[counter1][0] = indices[i][0];
					this.entries[counter1][1] = indices[i][1];
					counter1++;
					counter2++;
				
			}
		}
	}
	
	/**
	 * setter, converts a HashSet to an array of integers
	 * @param hash
	 */
	public void setEntries (HashSet<Pair<Integer>> hash) {
		int dim = hash.size();
		Iterator<Pair<Integer>> i = hash.iterator();
		this.entries = new int[dim][2];
		int j = 0;
		while(i.hasNext()) {
			Pair<Integer> pair = new Pair<Integer>(i.next());
			entries[j][0] = pair.getFirst();
			entries[j][1] = pair.getSecond();
			j++;
		}
	}
	
	/**
	 * setter, converts a RIGraph instance to a Individuals instance
	 * @param rig
	 */
	public void setIndis (RIGraph rig) {
		this.setName(rig.getPdbCode());
		Bound[][] bound = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		int lengt = rig.getEdgeCount(), l = 0;		
		this.entries = new int[lengt][2];
		for(int i = 0; i < lengt; i++){
			for(int j = i + 1; j < lengt; j++){
				if((bound[i][j] != null)&&(l < lengt)){
					if(i != j - 1){
						this.entries[l][0] = i;
						this.entries[l][1] = j;
						l++;
					}
				}
			}
		}
		this.setSequence(rig.getSequence());
		
	}
	
	/**
	 * setter, does the same as setIndis (RIGraph rig) setter method, can do random sampling over a set of Individuals
	 * instances
	 * @param rig
	 * @param conn
	 * @param randomize
	 * @throws Exception 
	 * @throws NullPointerException 
	 * @throws ArrayIndexOutOfBoundsException 
	 * @throws PdbLoadError 
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundError 
	 * @throws Exception 
	 * @throws Exception
	 */
	public void setIndis (RIGraph rig, MySQLConnection conn, boolean randomize) throws ArrayIndexOutOfBoundsException, NullPointerException, PdbCodeNotFoundError, SQLException, PdbLoadError  {
		this.setName(rig.getPdbCode());
		Bound[][] bound = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		double cut = rig.getCutoff();
		String ct = rig.getContactType();
		int lengt = bound.length, l = 0;
		this.entries = new int[rig.getEdgeCount()][2];
		Pdb n = getFullProt(rig, conn);
		RIGraph r = n.get_graph(ct, cut);
		double[][] dm = distMap(n);
		if(randomize){
			bound = randomSet(rig,conn);
		}
		this.setCMError(bound, r);
		this.setDMError(bound, dm);
		for(int i = 0; i < lengt; i++){
			for(int j = i + 1; j < lengt; j++){
				if((bound[i][j] != null)&&(l < lengt)){
					if(i != j - 1){
						this.entries[l][0] = i;
						this.entries[l][1] = j;
						l++;
					}
				}
			}
		}
		this.storer();
		this.setNumOfContacts(this.getHashSet());
		this.setSequence(rig.getSequence());
		this.setChainCode(rig.getPdbChainCode());
	}
	
	/**
	 * setter, does the same as setIndis (RIGraph rig) setter method, can do random sampling over a set of Individuals
	 * number of contactes specified by 'NumCont' parameter
	 * @param rig
	 * @param conn
	 * @param randomize
	 * @param NumCont
	 * @throws ArrayIndexOutOfBoundsException
	 * @throws NullPointerException
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public void setIndis (RIGraph rig, MySQLConnection conn, boolean randomize, int NumCont) throws ArrayIndexOutOfBoundsException, NullPointerException, PdbCodeNotFoundError, SQLException, PdbLoadError  {
		this.setName(rig.getPdbCode());
		Bound[][] bound = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		double cut = rig.getCutoff();
		String ct = rig.getContactType();
		int lengt = bound.length, l = 0;
		Pdb n = getFullProt(rig, conn);
		RIGraph r = n.get_graph(ct, cut);
		int edgec = r.getEdgeCount();
		int edgepercent = (int) ((double) NumCont/100.0*(double) edgec);
		double[][] dm = distMap(n);
		if(randomize){
			bound = randomSet(rig,conn,edgepercent);
			this.entries = new int[edgepercent][2];
		}
		else{
			this.entries = new int[rig.getEdgeCount()][2];
		}
		this.setCMError(bound, r);
		this.setDMError(bound, dm);
		for(int i = 0; i < lengt; i++){
			for(int j = i + 1; j < lengt; j++){
				if((bound[i][j] != null)&&(l < edgepercent)){
					if(i != j - 1){
						this.entries[l][0] = i;
						this.entries[l][1] = j;
						l++;
					}
				}
			}
		}
		this.storer();
		this.setNumOfContacts(this.getHashSet());
		this.setSequence(rig.getSequence());
		this.setChainCode(rig.getPdbChainCode());
	}
	
	/**
	 * setter, does the same as setIndis (RIGraph rig) setter method, can do random sampling over a set of Individuals
	 * number of contactes specified by 'NumCont' parameter
	 * @param rig
	 * @param conn
	 * @param randomize
	 * @param NumCont
	 * @throws ArrayIndexOutOfBoundsException
	 * @throws NullPointerException
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public void setIndis (RIGraph rig, MySQLConnection conn, boolean randomize, double NumCont) throws ArrayIndexOutOfBoundsException, NullPointerException, PdbCodeNotFoundError, SQLException, PdbLoadError  {
		this.setName(rig.getPdbCode());
		Bound[][] bound = Reconstructer.convertRIGraphToBoundsMatrix(rig);
		double cut = rig.getCutoff();
		String ct = rig.getContactType();
		int lengt = bound.length, l = 0;
		Pdb n = getFullProt(rig, conn);
		RIGraph r = n.get_graph(ct, cut);
		int edgec = r.getEdgeCount();
		int edgepercent = (int) (NumCont/100.0*(double) edgec);
		double[][] dm = distMap(n);
		if(randomize){
			bound = randomSet(rig,conn,edgepercent);
			this.entries = new int[edgepercent][2];
		}
		else{
			this.entries = new int[rig.getEdgeCount()][2];
		}
		this.setCMError(bound, r);
		this.setDMError(bound, dm);
		for(int i = 0; i < lengt; i++){
			for(int j = i + 1; j < lengt; j++){
				if((bound[i][j] != null)&&(l < edgepercent)){
					if(i != j - 1){
						this.entries[l][0] = i;
						this.entries[l][1] = j;
						l++;
					}
				}
			}
		}
		this.storer();
		this.setNumOfContacts(this.getHashSet());
		this.setSequence(rig.getSequence());
		this.setChainCode(rig.getPdbChainCode());
	}
	
	/**
	 * setter, reads Bound instance, String instance an the distance map as a double array
	 * @param bound
	 * @param n
	 * @param dm
	 */
	public void setIndis (Bound[][] bound, String n, double[][] dm) {
		int lengt = bound.length, l = 0;
		this.entries = new int[lengt][2];
		for(int i = 0; i < lengt; i++){
			for(int j = i + 1; j < lengt; j++){
				if((bound[i][j] != null)&&(l < lengt)){
					if(i != j - 1){
						this.entries[l][0] = i;
						this.entries[l][1] = j;
						l++;
					}
				}
			}
		}
		this.setName(n);
		this.setDMError(bound, dm);
	}
	
	/**
	 * setter, sets the name of this instance using the String 'i'
	 * @param i
	 */
	public void setName (String i) {
		this.name = i;
	}
	
	/**
	 * setter, sets the instance variable CMError using Scorer.getCMError(Bound[][], RIGraph) method
	 * @param bound
	 * @param rig
	 * @throws Exception
	 */
	public void setCMError (Bound[][] bound, RIGraph rig) {
		this.CMError = Scorer.getCMError(bound, rig);
	}
	
	/**
	 * setter, sets the instance variable DMError using Scorer.getDMError(Bound[][], double[][]) method
	 * @param bound
	 * @param dm
	 */
	public void setDMError (Bound[][] bound, double[][] dm) {
		this.DMError = Scorer.getDMError(bound, dm);
	}
	
	/**
	 * setter, sets the instance variable numOfContacts using SparseGraph.getEdgeCount() method
	 * @param rig
	 */
	public void setNumOfContacts (RIGraph rig) {
		this.numOfContacts = rig.getEdgeCount();
	}
	
	/**
	 * setter, sets the instance variable numOfContacts using HashSet.size() method
	 * @param in
	 */
	public void setNumOfContacts (Individuals in) {
		this.numOfContacts = in.getHashSet().size();
	}

	/**
	 * setter, sets the instance variable numOfContacts using HashSet.size() method
	 * @param hash
	 */
	public void setNumOfContacts (HashSet<Pair<Integer>> hash) {
		this.numOfContacts = hash.size();
	}
	
	/**
	 * setter, stores all entries in 'entries' instance variable in a HashSet
	 */
	public void storer () {
		store = new HashSet<Pair<Integer>>();
		int leng = (this.entries).length;
		//Integer[] inds = new Integer[leng];
		
		for(int i = 0; i < leng; i++){
			Integer in1 = new Integer(this.getNumbers(i, 0));
			Integer in2 = new Integer(this.getNumbers(i, 1));
				Pair<Integer> in = new Pair<Integer>(in1, in2);
				store.add(in);
		}
		store.hashCode();
		store.iterator();
	}
	
	/**
	 * setter, sets this protein sequence
	 * @param seq
	 */
	public void setSequence(String seq) {
		this.sequence = seq;
	}
	
	/**
	 * setter, sets this protein chain code
	 * @param chain
	 */
	public void setChainCode (String chain) {
		this.chainCode = chain;
	}
	
	public void setFifty (boolean tr) {
		this.fifty = tr;
	}
	
	/*----------------------getters----------------------------------------------*/
	
	/**
	 * getter, returns contact type, by default: "CA" 
	 * @return 'ct' - this contact type
	 */
	public String getContactT () {
		return this.ct;
	}
	
	/**
	 * getter, returns this contact distance
	 * @return 'di' - this cutoff distance
	 */
	public double getContactDist () {
		return this.di;
	}
	
	/**
	 * getter, returns this chain code
	 * @return 'chainCode'
	 */
	public String getChainCode () {
		return this.chainCode;
	}
	
	/**
	 * reconstructer, reconstructs a RIGraph instance using an Individuals instance
	 * @param in
	 * @return graph
	 */
	public static RIGraph reconstructGraph (Individuals in){
		String aa = in.getSequence();
		RIGraph graph = new RIGraph(aa);
		Iterator<Pair<Integer>> i = in.getHashSet().iterator();
		while(i.hasNext()) {
			Pair<Integer> pair = i.next();
			String aa1 = AAinfo.oneletter2threeletter(Character.toString(aa.charAt(pair.getFirst())));
			String aa2 = AAinfo.oneletter2threeletter(Character.toString(aa.charAt(pair.getSecond())));
			graph.addVertex(new RIGNode(pair.getFirst()+1, aa1));
			graph.addVertex(new RIGNode(pair.getSecond()+1, aa2));
			graph.addEdgeIJ(pair.getFirst()+1, pair.getSecond()+1);
		}
		graph.setContactType(in.getContactT());
		graph.setCutoff(in.getContactDist());
		graph.setPdbChainCode(in.getChainCode());
		graph.setPdbCode(in.getName());
		graph.setChainCode(in.getChainCode());
		return graph;
	}
	
	/**
	 * getter, returns this entries at position (i,j)
	 * @param i
	 * @param j
	 * @return 'entries[i][j]'
	 */
	public int getNumbers (int i, int j) {
		return this.entries[i][j];
	}

	/**
	 * getter, returns instance variable 'entries' as array of integers
	 * @return 'entries'
	 */
	public int[][] getEntries () {
		return this.entries;
	}
	
	/**
	 * getter, returns this name of the protein (using PDB code convention)
	 * @return 'name'
	 */
	public String getName () {
		return this.name;
	}
	
	/**
	 * getter, returns this CMError
	 * @return 'CMError'
	 */
	public double getCM () {
		return this.CMError;
	}
	
	/**
	 * getter, returns this DMError
	 * @return 'DMError
	 */
	public double getDM () {
		return this.DMError;
	}

	/**
	 * getter, returns this CMError
	 * @param pop
	 * @return error
	 */
	public static double[] getCM (Individuals[] pop){
		int dim = pop.length;
		double[] error = new double[dim];
		for(int i = 0; i < dim; i++){
			error[i] = pop[i].getCM();
		}
		return error;
	}
	
	/**
	 * getter, returns this DMError
	 * @param pop
	 * @return error
	 */
	public static double[] getDM (Individuals[] pop){
		int dim = pop.length;
		double[] error = new double[dim];
		for(int i = 0; i < dim; i++){
			error[i] = pop[i].getDM();
		}
		return error;
	}
	
	/**
	 * getter, returns this numOfContacts
	 * @return 'numOfContacts'
	 */
	public int getNumOfContacts () {
		return this.numOfContacts;
	}
	
	/**
	 * getter, returns this protein sequence
	 * @return 'sequence'
	 */
	public String getSequence () {
		return this.sequence;
	}
	
	/**
	 * defines a new hash function
	 * @return 'hasher'
	 */
	public int hashCode() {
		Individuals ind;
		int hasher = 17;
		ind = new Individuals();
		hasher = hasher*37 + ind.hashCode();
		return  hasher;
	}
	
	/**
	 * defines a new equals method for all instances of Individuals
	 */
	public boolean equals(Object o) {
		if(!(o instanceof Individuals)){
			return false;
		}
		
		return (entries.equals(((Individuals) o).entries));
	}
	
	/**
	 * getter, returns this HashSet
	 * @return 'store'
	 */
	public HashSet<Pair<Integer>> getHashSet () {
		return this.store;
	}
		
	/**
	 * getter, returns boolean field 'fifty'
	 * @return 'fifty'
	 */
	public boolean getFifty() {
		return this.fifty;
	}
	
	/**
	 * compares this Individuals with in Individuals retaining only all contacts they have in common
	 * @param in
	 * @return 'hash1' - a Hashset where all entries are the same
	 */
	public boolean comPare(Individuals in) {
		HashSet<Pair<Integer>> hash1 = this.getHashSet();
		HashSet<Pair<Integer>> hash2 = in.getHashSet();
		return hash1.retainAll(hash2);
	}
	
	
	/**
	 * compares this Individuals with in Individuals retaining only all contacts they have in common
	 * returns a HashSet
	 * @param in
	 * @return hash1
	 */
	public HashSet<Pair<Integer>> comPares(Individuals in) {
		HashSet<Pair<Integer>> hash1 = new HashSet<Pair<Integer>>(this.getHashSet());
		HashSet<Pair<Integer>> hash2 = new HashSet<Pair<Integer>>(in.getHashSet());
		hash1.retainAll(hash2);
		return hash1;
	}
	
	
	/**
	 * compares this Individuals with in Individuals retaining only all contacts they have in common
	 * returns a HashSet
	 * @param hash1
	 * @param hash2
	 * @return hash1a
	 */
	public static HashSet<Pair<Integer>> comPares(HashSet<Pair<Integer>> hash1, HashSet<Pair<Integer>> hash2) {
		HashSet<Pair<Integer>> hash1a = new HashSet<Pair<Integer>>(hash1);
		HashSet<Pair<Integer>> hash2a = new HashSet<Pair<Integer>>(hash2);
		hash1a.retainAll(hash2a);
		return hash1a;
	}
	
	/*public Object get(Object in) {
		
	}*/
	
	/*public boolean[] comParer (Individuals ind) {
		//boolean[] out = new boolean[];
		ind.getHashSet().get();
	}*/
	
	/**
	 * compares two instances of Individuals by their HashSets retaining only those contacts both have in common,
	 * if contacts are missing the remaining gaps are filled randomly with contacts from both 'parents'
	 * @param in1
	 * @param in2
	 * @return neu
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundError 
	 * @throws PdbLoadError
	 */
	public static Individuals breedIndis (Individuals in1, Individuals in2) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		Individuals neu = new Individuals();
		HashSet<Pair<Integer>> hash1 = new HashSet<Pair<Integer>>(in1.getHashSet());
		HashSet<Pair<Integer>> hash2 = new HashSet<Pair<Integer>>(in2.getHashSet());
		if(in1.getName().matches(in2.getName())){
			neu.setSequence(in1.getSequence());
			HashSet<Pair<Integer>> hash = comPares(hash1, hash2);
			neu.store = hash;
			neu.name = in1.getName();
			neu.numOfContacts = hash.size();
			neu.setEntries(hash);
			neu.setChainCode(in1.getChainCode());
			Iterator<Pair<Integer>> i = hash1.iterator();
			Iterator<Pair<Integer>> j = hash2.iterator();
			while(neu.getNumOfContacts() != in1.getNumOfContacts()) {
				Random rand = new Random();
				if(rand.nextInt(2) % 2 == 0) {
					boolean test = neu.store.add(i.next());
					while(!test && i.hasNext()){
						test = neu.store.add(i.next());
						}
					i = hash1.iterator();
				}
				else {
					boolean test = neu.store.add(j.next()); 
					while(!test&&j.hasNext()){
						test = neu.store.add(j.next());
					}
					j = hash2.iterator();
				}
				neu.setNumOfContacts(neu.store);
			}
		}
		neu.setEntries(neu.getHashSet());
		MySQLConnection conn = new MySQLConnection();
		RIGraph rig = new RIGraph(neu.getSequence());
		rig = reconstructGraph(neu);
		Pdb prot = new PdbasePdb(neu.getName(), "pdbase", conn);
		prot.load(rig.getPdbChainCode());
		RIGraph full = prot.get_graph(neu.getContactT(), neu.getContactDist());
		neu.setCMError(Reconstructer.convertRIGraphToBoundsMatrix(rig), full);
		neu.setDMError(Reconstructer.convertRIGraphToBoundsMatrix(rig), distMap(prot));
		return neu;
	}
	
	/**
	 * Method, to copy a Individuals instance to a specified Individuals instance
	 * @param in source to be copied
	 * @param out destination
	 */
	public static void copyInd (Individuals in, Individuals out){
		if(out != null){
		out.setIndis(in);
		}
		else{
			out = new Individuals(in);
		}
	}
	
	public static Individuals[] shrinkIndArray (Individuals[] in){
		int dim = in.length, counter1 = 0, counter2 = 0;
		boolean[] test = new boolean[dim];
		for(int i = 0; i < dim; i++){
			if(in[i] != null){
				test[i] = true;
				counter1++;
			}
		}
		Individuals[] shrunk = new Individuals[counter1];
		for(int i = 0; i < dim; i++){
			if(test[i]){
				copyInd(in[i],shrunk[counter2]);
			}
		}
		return shrunk;
	}
	/*-----------------------display methods-------------------------------------------*/
	
	/**
	 * display method, to display a HashSet
	 * @param hash
	 */
	public static void displayHash (HashSet<Pair<Integer>> hash) {
		System.out.println("Hash set: "+hash.toString());
	}
	
	/**
	 * display method, displays Individuals instances
	 */
	public void displayIndis (){
		System.out.print("{");
		int ind = (this.getEntries()).length - 1;
		System.out.print("# of contacts :"+this.getEntries().length+" and contacts of "+this.name+": ");
		for(int i = 0; i < ind; i++){
			System.out.print("{ "+this.getNumbers(i, 0)+", ");
			System.out.print(this.getNumbers(i, 1)+"},");//}
		}
		System.out.print("{ "+this.getNumbers(ind,0)+", "+this.getNumbers(ind,1)+"}");//}
		System.out.println("}, with CMError:"+this.getCM()+" and DMError: "+this.getDM());
	}
	
	/*----------------------------statics-------------------------------------------*/
	/*-------------------------random generators------------------------------------*/
	
	/**
	 * random contact map generator, uses RIGraph MySQLConnection and a percentage input parameter
	 * @param input
	 * @param conn
	 * @param percent
	 * @return in
	 * @throws PdbCodeNotFoundError
	 * @throws NullPointerException
	 * @throws ArrayIndexOutOfBoundsException
	 * @throws SQLException
	 * @throws PdbLoadError
	 */
	public static Individuals randomSets (RIGraph input, MySQLConnection conn, double percent) throws NullPointerException, ArrayIndexOutOfBoundsException, PdbCodeNotFoundError, SQLException, PdbLoadError{
		Distiller dist = new Distiller(input);
		Bound[][] subset = dist.sampleSubset((int) (input.getEdgeCount()*percent));
		Individuals in = new Individuals(subset, input.getPdbCode(), distMap(getFullProt(input, conn)));
		return in;
	}
	
	/**
	 * method, that randomly selects contact of the full contact map
	 * @param input
	 * @param conn
	 * @return subset
	 * @throws Exception
	 * @throws SQLException
	 * @throws PdbLoadError
	 * @throws PdbCodeNotFoundError
	 * @throws NullPointerException
	 * @throws ArrayIndexOutOfBoundsException
	 */
	public static Bound[][] randomSet (RIGraph input, MySQLConnection conn) throws SQLException, PdbLoadError, PdbCodeNotFoundError, NullPointerException, ArrayIndexOutOfBoundsException {
		double cuto = input.getCutoff();
		String ct = input.getContactType();
		Pdb prot = getFullProt(input, conn);
		prot.load(input.getPdbChainCode());
		RIGraph full = prot.get_graph(ct, cuto);
		int numOfConts = input.getEdgeCount();
		Distiller dist = new Distiller(full);
		Bound[][] subset = dist.sampleSubset(numOfConts);
		return subset;
	}
	
	/**
	 * method that randomly selects contacts from the full contact map 
	 * @param input
	 * @param conn
	 * @param NumCont
	 * @return
	 * @throws SQLException
	 * @throws PdbLoadError
	 * @throws PdbCodeNotFoundError
	 * @throws NullPointerException
	 * @throws ArrayIndexOutOfBoundsException
	 */
	public static Bound[][] randomSet (RIGraph input, MySQLConnection conn, int NumCont) throws SQLException, PdbLoadError, PdbCodeNotFoundError, NullPointerException, ArrayIndexOutOfBoundsException {
		double cuto = input.getCutoff();
		String ct = input.getContactType();
		Pdb prot = getFullProt(input, conn);
		prot.load(input.getPdbChainCode());
		RIGraph full = prot.get_graph(ct, cuto);
		int numOfConts = NumCont;
		Distiller dist = new Distiller(full);
		Bound[][] subset = dist.sampleSubset(numOfConts);
		return subset;
	}
	
	/**
	 * method to recover a complete Pdb instance
	 * @param rig
	 * @param conn
	 * @return prot
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundError 
	 * @throws PdbLoadError
	 */
	public static Pdb getFullProt (RIGraph rig, MySQLConnection conn) throws PdbCodeNotFoundError, SQLException, PdbLoadError {
		String pdbCode = rig.getPdbCode();
		String pdbaseDb = "pdbase";
		Pdb prot = new PdbasePdb(pdbCode, pdbaseDb, conn);
		prot.load(rig.getPdbChainCode());
		return prot;
	}
	
	/**
	 * method to calculate the distance map
	 * @param prot
	 * @return fullDistanceMap
	 */
	public static double[][] distMap (Pdb prot) {
		Matrix fullDistanceMap = prot.calculateDistMatrix("CA");
		return fullDistanceMap.getArray();
	}
	
	/*------------------------------main to test-----------------------------------------------*/
	
	public static void main (String[] args) throws Exception {
		File dir = new File("/project/LitNet/Essence/Data/Contactmaps");
		File[] test = dir.listFiles(new RegexFileFilter(".*_all\\.graph"));
		Arrays.sort(test);
		File neu = test[0];
		RIGraph sub = new FileRIGraph(neu.getAbsolutePath());
		MySQLConnection conn = new MySQLConnection();
		/*Individuals indi = new Individuals(sub, conn, false);
		indi.displayIndis();*/
		Individuals[] neus = new Individuals[20];
		for(int i = 0; i < neus.length; i++){
			neus[i] = new Individuals(sub, conn, true);
			neus[i].displayIndis();
		}
		System.out.println("Average DMError: "+Population.getAverageDMError(neus));
		Individuals[] off1 = new Individuals[neus.length];
		Individuals[] off2 = new Individuals[neus.length];
		Population.copyInd(neus,off1);
		Population.copyInd(neus,off2);
		for(int i = 0; i < 5; i++){
			Population pop1 = new Population(off1);
			Population pop2 = new Population(off2);
			Population.copyInd(Population.evolve(false, pop1),off1);
			System.out.print("Average DMError of best 50 %: "+Population.getAverageDMError(off1)+", ");
			System.out.println("standart dev. of DMError of best 50 %: "+Population.getDMStDeviation(off1)+" ");
			Population.copyInd(Population.evolve(false, pop2),off2);
			System.out.print("Average DMError: "+Population.getAverageDMError(off2)+", ");
			System.out.println("standart dev of DMError: "+Population.getDMStDeviation(off2));
		}
		for(int zahl = 0; zahl < off1.length; zahl++){
			off1[zahl].displayIndis();
		}
		System.out.println("Average DMError of best 50 %: "+Population.getAverageDMError(off1));
		for(int zahl = 0; zahl < off1.length; zahl++){
			off2[zahl].displayIndis();
		}
		System.out.println("Average DMError of all below mean value: "+Population.getAverageDMError(off2));
		conn.close();
	}
	
	public static int compareTo (Pair<Double> tree, Pair<Double> input){
		double val1 = tree.getFirst().doubleValue();
		double val2 = input.getFirst().doubleValue();
		return (val1 < val2 ? -1 : (val1 == val2 ? 0 : 1));
	}
}
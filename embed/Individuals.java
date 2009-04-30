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

public class Individuals extends HashSet<Pair<Integer>> {
	
	/**
	 * Fields
	 */
	private static final long serialVersionUID = 1L;
	
	private final String ct = "Ca";
	
	private final double di = 9.0;

	private int [][] entries;
	
	private HashSet<Pair<Integer>> store;
	
	private int numOfContacts;
	
	private String name, sequence, chainCode;
	
	private double CMError, DMError;
	
	/*----------------------Constructors--------------------------------------------------*/
	public Individuals () {
		
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
	public Individuals (RIGraph rig) throws Exception {
		this.setIndis(rig);
	}
	
	/**
	 * three parameter constructor, reads instances of RIGraph, MySQLConnection and a boolean expression
	 * as input parameter
	 * @param rig
	 * @param conn
	 * @param rand
	 * @throws Exception
	 */
	public Individuals (RIGraph rig, MySQLConnection conn, boolean rand) throws Exception {
		this.setIndis(rig, conn, rand);
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
	}
	
	/**
	 * setter, reads an array of integers and uses them for the 'entries' instance variable
	 * @param indices
	 */
	public void setIndis (int [][] indices) {
		int dim = indices.length;
		int k = 0;
		while(indices[k][0] != 0 || indices[k][1] != 0){
			k++;
		}
		this.entries = new int[k][2];
		int counter1 = 0;
		int counter2 = 0;
		for(int i = 0; i < dim; i++){
			for(int j = i + 1; j < dim; j++){
				if(i != j - 1){
					this.entries[counter1][counter2] = indices[i][j];
					counter1++;
					counter2++;
				}
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
	 * @throws Exception
	 */
	public void setIndis (RIGraph rig) throws Exception {
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
	 * @throws Exception
	 */
	public void setIndis (RIGraph rig, MySQLConnection conn, boolean randomize) throws Exception  {
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
	public void setCMError (Bound[][] bound, RIGraph rig) throws Exception {
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
	
	/*----------------------getters----------------------------------------------*/
	
	/**
	 * getter, returns contact type, by default: "CA" 
	 */
	public String getContactT () {
		return this.ct;
	}
	
	/**
	 * getter, returns this contact distance
	 * @return
	 */
	public double getContactDist () {
		return this.di;
	}
	
	/**
	 * getter, returns this chain code
	 * @return
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
			String aa1 = Character.toString(aa.charAt(pair.getFirst() + 1));
			String aa2 = Character.toString(aa.charAt(pair.getSecond() + 1));
			graph.addVertex(new RIGNode(pair.getFirst(), aa1));
			graph.addVertex(new RIGNode(pair.getSecond(), aa2));
			graph.addEdgeIJ(pair.getFirst(), pair.getSecond());
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
	 * @return
	 */
	public int getNumbers (int i, int j) {
		return this.entries[i][j];
	}

	/**
	 * getter, returns instance variable 'entries' as array of integers
	 * @return
	 */
	public int[][] getEntries () {
		return this.entries;
	}
	
	/**
	 * getter, returns this name of the protein (using PDB code convention)
	 * @return
	 */
	public String getName () {
		return this.name;
	}
	
	/**
	 * getter, returns this CMError
	 * @return
	 */
	public double getCM () {
		return this.CMError;
	}
	
	/**
	 * getter, returns this DMError
	 * @return
	 */
	public double getDM () {
		return this.DMError;
	}
	
	/**
	 * getter, returns this numOfContacts
	 * @return
	 */
	public int getNumOfContacts () {
		return this.numOfContacts;
	}
	
	/**
	 * getter, returns this protein sequence
	 * @return
	 */
	public String getSequence () {
		return this.sequence;
	}
	
	/**
	 * defines a new hash function
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
	 * @return
	 */
	public HashSet<Pair<Integer>> getHashSet () {
		return this.store;
	}
	
	/**
	 * compares this Individuals with in Individuals retaining only all contacts they have in common
	 * @param in
	 * @return
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
	 * compares two instances of Individuals by their HashSets retaining only those contacts both have in common
	 * @param in1
	 * @param in2
	 * @return neu
	 * @throws Exception
	 */
	public static Individuals breedIndis (Individuals in1, Individuals in2) throws Exception {
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
				if(rand.nextInt() % 2 == 0) {
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
			//if(this.getNumbers(i, 0) + 1 != this.getNumbers(i, 1)){
				System.out.print("{ "+this.getNumbers(i, 0)+", ");
				System.out.print(this.getNumbers(i, 1)+"},");//}
		}
		//if(this.getNumbers(ind, 0) + 1 != this.getNumbers(ind, 1)){
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
	 * @throws Exception
	 * @throws NullPointerException
	 * @throws ArrayIndexOutOfBoundsException
	 */
	public static Individuals randomSets (RIGraph input, MySQLConnection conn, double percent) throws Exception, NullPointerException, ArrayIndexOutOfBoundsException {
		Distiller dist = new Distiller(input);
		Bound[][] subset = dist.sampleSubset((int) (input.getEdgeCount()*percent));
		Individuals in = new Individuals(subset, input.getPdbCode(), distMap(getFullProt(input, conn)));
		return in;
	}
	
	/**
	 * method, to determine the average error (CM and DM) and breeds with the 50% best fittest to return and new array
	 * of Individuals
	 * @param parents
	 * @return offspring
	 * @throws Exception
	 */
	public static Individuals[] evolve (Individuals[] parents) throws Exception {
		int dim = parents.length;
		double[][] array = new double[dim][2];
		for(int i = 0; i < parents.length; i++){
			array[i][0] = parents[i].getDM();
			array[i][1] = (double) i;
		}
		double CMav = getAverageCMError(parents);
		double DMav = getAverageDMError(parents);
		Individuals[] offspring = new Individuals[dim];
		for(int i = 0; i < dim - 1; i++){
			for(int j = i + 1; j < dim; j ++){
				if(parents[i].getCM() <= CMav && parents[j].getCM() <= CMav){
					if(parents[i].getDM() <= DMav && parents[j].getDM() <= DMav){
						offspring[i] = breedIndis(parents[i], parents[j]);
					}
				}
			}
		}
		return offspring;
	}
	
	/*public static void sortArray(double[][] array, int row, int dim2){
		int dim = array.length;
		int[] sorting = new int[dim];
		for(int k = 1; k < (int) ((double) dim)*0.5; k++){
			for(int i = 0; i < dim - k - 1; i++){
				if(array[sorting[i]][row] < array[sorting[i + k]][row]){
					sorting[i + k] = (int) array[i][1];
					sorting[i + k + 1] = (int) array[i + k][1];
				}
				else{
					sorting[i + k] = (int) array[i + k][1];
					sorting[i + k + 1] = (int) array[i][1];
				}
			}
		}
		double[][] l = new double[dim][dim2];
		for(int j = 0; j < dim; j++){
			for(int m = 0; m < dim2 ; m++){
				l[j][m] = array[j][m];
			}
		}
		for(int i = 0; i < dim ; i++){
			for(int j = 0; j < dim2; j++){
				array[i][j] = l[sorting[i]][j];
			}
		}
	}
	/**muss noch ueberarbeitet werden, z.b. mit methode, die in liste nur nach Individuals instanzen sucht, die unterhalb eines thresholds liegen
	 * 
	 * @param array
	 * @param row
	 */
	/*public static void sortA (double[][] array, int row){
		TreeSet<Pair<Double>> h = new TreeSet<Pair<Double>>();
		for(int i = 0; i < array.length; i++){
			Pair<Double> n = new Pair<Double>(new Double(array[i][0]), new Double (array[i][1]));
			//Comparator<? super Pair<Double>> c = h.comparator();
			h.add(n);
			compareTo(h.last(), n);
		}
		Iterator<Pair<Double>> i = h.iterator();
		int j = 0;
		while(i.hasNext()){
			Pair<Double> pair = new Pair<Double>(i.next());
			array[j][0] = pair.getFirst().doubleValue();
			array[j][1] = pair.getSecond().doubleValue();
			j++;
		}
	}*/

	/**
	 * method to determine average DM error
	 * @param pop
	 * @return avverageDMError
	 */
	public static double getAverageDMError (Individuals[] pop) {
		int dim = pop.length;
		double averageDMError = 0.0;
		for(int i = 0; i < dim; i++){
			averageDMError = averageDMError + pop[i].getDM();
		}
		return averageDMError/(double) dim;
	}

	/**
	 * method to determine average CM error
	 * @param pop
	 * @return avverageCMError
	 */
	public static double getAverageCMError (Individuals[] pop) {
		int dim = pop.length;
		double averageCMError = 0.0;
		for(int i = 0; i < dim; i++){
			averageCMError = averageCMError + pop[i].getCM();
		}
		return averageCMError/(double) dim;
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
	public static Bound[][] randomSet (RIGraph input, MySQLConnection conn) throws Exception, SQLException, PdbLoadError, PdbCodeNotFoundError, NullPointerException, ArrayIndexOutOfBoundsException {
		double cuto = input.getCutoff();
		String ct = input.getContactType();
		Pdb prot = getFullProt(input, conn);
		prot.load(input.getPdbChainCode());
		RIGraph full = prot.get_graph(ct, cuto);
		int numOfConts = input.getEdgeCount();
		Distiller dist = new Distiller(full);
		//double percent = (double) numOfConts/(double) full.getEdgeCount();
		Bound[][] subset = dist.sampleSubset(numOfConts);
		//Bound[][] allbounds = Reconstructer.convertRIGraphToBoundsMatrix(full);
		/*allbounds = Scorer.inferAllBounds(allbounds);
		int allboundsdim = allbounds.length;
		Bound[][] subset = new Bound[allboundsdim][allboundsdim];
		int i = 0;
		while( i < numOfConts){
			Random rand = new Random();
			int random1 = rand.nextInt(allboundsdim);
			int random2 = rand.nextInt(allboundsdim - random1) + random1;/*
			random1 = (int) ((long) random1) % allboundsdim;
			random2 = (int) ((long) random2) % allboundsdim;
			if((random1 < allboundsdim && random1 >= 0)&&(random2 < allboundsdim && random2 >= 0)){
				if(allbounds[random1][random2] != null){
					subset[random1][random2] = new Bound(allbounds[random1][random2].lower, allbounds[random1][random2].upper);
					i++;
				}
			}
		}*/
		/*for (int i=0;i< subset.length;i++) {
			for (int j=i+2;j< subset.length;j++) {
				if (subset[i][j]!=null) {
					System.out.println(i+" "+j+": "+subset[i][j]);
				}
			*/
		
		return subset;
	}
	
	/**
	 * method to recover a complete Pdb instance
	 * @param rig
	 * @param conn
	 * @return prot
	 * @throws Exception
	 */
	public static Pdb getFullProt (RIGraph rig, MySQLConnection conn) throws Exception {
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
		Individuals[] neus = new Individuals[6];
		for(int i = 0; i < neus.length; i++){
			neus[i] = new Individuals(sub, conn, true);
			neus[i].displayIndis();
		}
		Individuals[] off = evolve(neus);
		for(int zahl = 0; zahl < off.length; zahl++){
			off[zahl].displayIndis();	
		}
		conn.close();
	}
	
	public static int compareTo (Pair<Double> tree, Pair<Double> input){
		double val1 = tree.getFirst().doubleValue();
		double val2 = input.getFirst().doubleValue();
		return (val1 < val2 ? -1 : (val1 == val2 ? 0 : 1));
	}
}

/*class CompType implements Comparator<Pair<Double>> {
	public int compare (Pair<Double> tree, Pair<Double> input){
		double val1 = tree.getFirst().doubleValue();
		double val2 = input.getFirst().doubleValue();
		return (val1 < val2 ? -1 : (val1 == val2 ? 0 : 1));
	}
}*/
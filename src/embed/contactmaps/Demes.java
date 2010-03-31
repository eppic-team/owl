package embed.contactmaps;


import java.io.*;
import java.sql.SQLException;
import java.util.*;

import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbasePdb;
import owl.core.structure.RIGraph;
import owl.core.util.MySQLConnection;
import owl.core.util.RegexFileFilter;

//import Jama.Matrix;
import edu.uci.ics.jung.graph.util.Pair;
//import tools.*;

/**
 * <p>
 * This class deals with arrays of <code>{@link Individuals}</code> instances and can evolve a set of <code>{@link Individuals}</code> instances.
 * It compares two array entries, determines the entry with smaller error and can generate offspring by aligning and retaining
 * the contacts present in both <code>{@link Individuals}</code> instances and randomly selects the missing contacts (the number
 * of contacts, which is equivalent to field set by the <code>{@link Individuals#setNumOfContact(int)}</code> method, is kept constant over evolution).
 * In order to avoid unwanted 'shrinkage' of the population size, after each evolution step, the best 50 percent of all <tt>Individuals</tt> (i.e. the parental
 * generation) are kept and ranked again. Only the best ranked <tt>Individuals</tt> are in turn transferred to the filial generation.
 * </p>
 * <p>
 * The first step is always the initialization of an array of <tt>Individuals<tt>, as a random sampling over all contact maps. The random distribution of the
 * initial generation can be surveyed by the method <code>{@link #getMet()}</code>. This method returns an index table, representing a pairwise alignment with each <tt>Individuals</tt>.
 * The metric used is defined as follows:
 * </p>
 * <p>
 * </p>
 * <p>
 * <tt>d(x,y) := sum min(|i_x - i_y| + |j_x - j_y|)</tt>, over all contacts, where <tt>i</tt> is the first contact and <tt>j</tt> is the second contact.
 * <tt>x</tt> and <tt>y</tt> are two <tt>Individuals</tt>. On average in each <tt>Demes</tt> instance, the frequency of a given contact pair <tt>(i,j)</tt> is 10 %,
 * that is, in an <tt>Demes</tt> instance, with 20 <tt>Individuals</tt> this contact is present in two different <tt>Individuals</tt>.
 * </p>
 * <p>
 * </p>
 * <p>
 * Note, that multiple runs of instances with non-matching pdb-codes must call the method <code>{@link Individuals#clearFullCMandDM()}</code> prior to initialization,
 * otherwise severe Exception and/or Errors may occur. 
 * @author gmueller
 *
 */
public class Demes {
	
	/*-----------------------------------fields-----------------------------------------------*/
	
	/**
	 * field: an array of 'Individuals' instances
	 */
	private Individuals[] pop;
	
	/**
	 * field: an array of boolean, used to determine the 50% best of this 'Species' according to the 
	 * selected method (DMError or CMError)
	 */
	private boolean[] fiftyfifty;
	
	/**
	 * statistical fields: defined by the given error estimation method
	 */
	private double CMError, CMstdev, DMError, DMstdev;
	
	/**
	 * field: name - PDB codes used
	 */
	private String name;
	
	/**
	 * field: size - number of 'Individuals' instance in this 'Species'
	 * field: gen - number of generation
	 */
	private int size, gen;
	
	/**
	 * field: contacts that all the Individuals have in common
	 */
	private HashSet<Pair<Integer>> subhash;
	
	//private HashMap<Pair<Integer>, Double> metric;
	
	private HashMap<Pair<Integer>,HashSet<Integer>> indexer;
	
	private Metric metric2;
	
	private HashMap<Pair<Integer>, Integer> weighted;
	
	/**
	 * field: returning an array of indices of the best ranked Individuals
	 */
	private int[] ranked;
	
	
	//private int[][] weighted;
	
	private static final String pdbaseDb = "pdbase_20090728";
	
	/*-------------------------Constructors-----------------------------------------------------*/
	
	/**
	 * zero parameter constructor, sets all fields to their default values
	 */
	public Demes() {
		size = 0;
		pop = new Individuals[1];
		pop[0] = new Individuals();
		CMError = 0.0;
		DMError = 0.0;
		CMstdev = 0.0;
		DMstdev = 0.0;
		fiftyfifty = new boolean[1];
	}
	
	/**
	 * one parameter constructor, creates a Species of 'i' individuals using the zero parameter constructor
	 * of 'Individuals'
	 * @param i an Integer, setting the size of this instance
	 */
	public Demes(int i) {
		setSize(i);
		pop = new Individuals[i];
		subhash = new HashSet<Pair<Integer>>();
		for(int j = 0; j < i; j++){
			pop[i] = new Individuals();
		}
		CMError = 0.0;
		DMError = 0.0;
		CMstdev = 0.0;
		DMstdev = 0.0;
		fiftyfifty = new boolean[i];
	}
	
	/**
	 * one parameter constructor, reads an array of Individuals instances
	 * @param p
	 */
	public Demes (Individuals[] p) {
		setPop(p);
		setSize(p.length);
		setCMErrorStats(p);
		setDMErrorStats(p);
		setPopName(p[0].getName());
		fiftyfifty = new boolean[p.length];
		setSubHash(p);
		//setWHash();
		setWHasher();
		setMetrics();
		setIndexer();
	}
	
	/**]
	 * one parameter constructor, reads an instance of 'Species' - similar to a copy method
	 * @param pop a <tt>Deme</tt>
	 */
	public Demes (Demes pop) {//throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = pop.getSize();
		setSize(dim);
		this.pop = new Individuals[dim];
		setPop(pop.getPop());
		//setFullContactMap(pop);
		setCMErrorStats(pop.getPop());
		setDMErrorStats(pop.getPop());
		fiftyfifty = new boolean[dim];
		setPopName(pop.getName());
		setSubHash(pop.getPop());
		//setWHash();
		setWHasher();
		setMetrics();
		setIndexer();
		}
	
	/**
	 * Two parameter constructor: uses an starter <tt>Individuals</tt> and creates a
	 * random starter sub population (deme).
	 * @param starter
	 * @param size
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public Demes (Individuals starter, int size) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		MySQLConnection conn = new MySQLConnection ();
		Pdb pdb = new PdbasePdb(starter.getName(), pdbaseDb, conn);
		pdb.load(starter.getChainCode());
		RIGraph rig = pdb.getRIGraph("Ca", 9.0);
		pop = new Individuals[size];
		for(int i = 0; i < size; i++){
			pop[i] = new Individuals(rig,conn,true,starter.getNumOfContacts());
		}
		fiftyfifty = new boolean[size];
		setIndexer();
	}

	/**
	 * one parameter constructor, reads an array of Individuals instances
	 * @param p 
	 */
	public Demes (Individuals[] p, int generation) {//throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		setPop(p);
		setSize(p.length);
		//setFullContactMap(p[0]);
		setCMErrorStats(p);
		setDMErrorStats(p);
		setPopName(p[0].getName());
		if(generation >= 8){
			System.out.println("test");
		}
		fiftyfifty = new boolean[p.length];
		setSubHash(p);
		//setWHash();
		setWHasher();
		setGen(generation);
		setMetrics();
		setIndexer();
	}
	
	/**
	 * one parameter constructor, reads an instance of 'Species' - similar to a copy method
	 * @param pop 
	 */
	public Demes (Demes pop, int generation){// throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = pop.getSize();
		setSize(dim);
		this.pop = new Individuals[dim];
		setPop(pop.getPop());
		//setFullContactMap(pop);
		setCMErrorStats(pop.getPop());
		setDMErrorStats(pop.getPop());
		if(generation >= 8){
			System.out.println("test");
		}
		fiftyfifty = new boolean[dim];
		setPopName(pop.getName());
		setSubHash(pop.getPop());
		//setWHash();
		setWHasher();
		setGen(generation);
		setMetrics();
		setIndexer();
		}
	
	
	/**
	 * five parameter constructor: uses the 'sathyasSet()' method.
	 * @param popsize
	 * @param startfiles
	 * @param numofCons
	 * @param generation
	 * @param addname
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public Demes (int popsize, int startfiles, int numofCons, int generation, String addname) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		setSize(popsize);
		sathyasSet(popsize, startfiles, numofCons);
		setPopName(getPop(0).getName() + addname);
		//setFullContactMap(pop[0]);
		setCMErrorStats(pop);
		setDMErrorStats(pop);
		fiftyfifty = new boolean[popsize];
		//setWHash();
		setGen(generation);
		setWHasher();
		setMetrics();
		setSubHash(pop);
		setIndexer();
	}
	
	/**
	 * One parameter constructor: takes a String denoting a contact map file and converts the 'cmap' file to <code>{@link Demes}</code> instance.
	 * Note, that any processed file must have 'cmap' extension, otherwise an <tt>IllegalArgumentException</tt> is thrown. Additionally, any initialization
	 * of multiple instances via file reading, this constructor must be called for every single initialization.  
	 * @param file_name a String denoting the 'cmap' file
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 * @throws IllegalArgumentException if the denoted file is neither a file nor is not in the required format
	 */
	public Demes (String file_name) throws FileNotFoundException, IOException, PdbCodeNotFoundError, SQLException, PdbLoadError, IllegalArgumentException {
		File file = new File (file_name);
		if(file.exists()){
			readFile(file);
		}
		else throw new IllegalArgumentException ("The denoted directory does not exist.");
	}
	
	/**
	 * Three parameter constructor: initializes an instance of this class. The first parameter <tt>pdb_code</tt> is the
	 * standard pdb code, identifying the protein. If no such pdb code is present in the database, an <code>{@link PdbCodeNotFoundError}</code>
	 * is issued. Note, that only the default database is used, in order to change the database, please use the four parameter
	 * constructor <code>{@link #Demes(String, String, int, double)}</code>. Additionally, this constructor does only random
	 * samples.
	 * @param pdb_code a String denoting the pdb code
	 * @param size the number of entries in the Individuals array
	 * @param percent_cont the number of contacts each Individuals has, must hold <tt> 0.0 < percent <= 100.0</tt>
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError if an unknown pdb code is used
	 * @throws PdbLoadError
	 * @throws IllegalArgumentException if 'size <= 0' or 'percent not in (0.0,100.0]' 
	 */
	public Demes (String pdb_code, int size, double percent_cont) throws SQLException, PdbCodeNotFoundError, PdbLoadError, IllegalArgumentException {
		setDemes(pdb_code, pdbaseDb, size, percent_cont);
	}
	
	/**
	 * Four parameter constructor: initializes an instance of this class. The first parameter <tt>pdb_code</tt> is the
	 * standard pdb code, identifying the protein. If no such pdb code is present in the database, an <code>{@link PdbCodeNotFoundError}</code>
	 * is issued. Note, that database used must be specified. If the default database shall be used, please use the three parameter constructor
	 * <code>{@link #Demes(String, int, double)}</code>. Additionally, this constructor does only random
	 * samples.
	 * @param pdb_code a String denoting the pdb code
	 * @param db the database
	 * @param size the number of entries in the Individuals array
	 * @param percent the number of contacts each Individuals has, must hold <tt> 0.0 < percent <= 100.0</tt>
	 * @throws IllegalArgumentException  if 'size <= 0' or 'percent not in (0.0,100.0]'
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError if an unknown pdb code is used
	 * @throws PdbLoadError
	 */
	public Demes (String pdb_code, String db, int size, double percent) throws IllegalArgumentException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		setDemes(pdb_code, db, size, percent);
	}
	
	/*-------------------------------------Setters------------------------------------------*/
	
	/**
	 * one parameter setter, reads an array of 'Individuals'
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 */
	public void setPop(Individuals[] pop){// throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		if(!hasNull(pop)){
			int dim = pop.length;
			this.pop = new Individuals[dim];
			//Individuals.fullcontactmap = Species.fullcontactmap;
			for(int i = 0; i < dim; i++){
				this.pop[i] = new Individuals(pop[i]);
			}
			if(!samePDBCode()){
				System.err.println("The input array of Individuals instances do not have a common PDB code, causing abrupt termination...");
				System.exit(1);
			}
		}
		else{
			String nullp = "At least one instance in this Individuals array was not initialized.";
			throw new NullPointerException (nullp);
		}
	}
	
	/**
	 * setter, sets this Species instance's name, usually the pdb code
	 */
	public void setPopName (String na){
		name = na;
	}
	
	/**
	 * setter, sets this average CMError and DMError
	 */
	public void setCMErrorStats (Individuals[] pop) {
		CMError = getAverageCMError(pop);
		CMstdev = getCMStDeviation(pop);
	}
	
	/**
	 * setter, sets this standard deviation
	 */
	public void setDMErrorStats (Individuals[] pop) {
		DMError = getAverageDMError(pop);
		DMstdev = getDMStDeviation(pop);
	}
	
	/**
	 * setter, sets this Species size, i.e. length of the Individuals array
	 */
	public void setSize(int i){
		//if(i >= 14){
			size = i;
		/*}
		else{
			System.err.println("The population size must always be greater or equal to 14.");
			System.exit(1);
		}*/
	}
	
	public void setSize2(int i){
		if(i >= 14){
			size = i;
		}
		else{
			System.err.println("The population size must always be greater or equal to 14.");
			System.exit(1);
		}
	}
	
	/**
	 * Setter, sets this field {@link #fiftyfifty} by pairwise comparing this array of {@link #Individuals} and
	 * setting 'fiftyfifty' true, if this Individual has a better Error value than the worst approx. 50 %.
	 * @param CMDM - determines which Error is taken into account: CMDM = true -> CMError, false -> DMError
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 */
	public void bestFifty (boolean CMDM) {//throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = getSize();
		double[] error = new double[dim];
		int[] better = new int[dim];
		if(CMDM){
			System.arraycopy(Individuals.getCM(getPop()),0,error,0,dim);
		}
		else{
			System.arraycopy(Individuals.getDM(getPop()),0,error,0,dim);
		}
		for(int i = 0; i < dim - 1; i++){
			for(int j = i + 1; j < dim; j++){
				if(error[i] <= error[j]){
					better[i] = better[i] + 1;
				}
				if(error[i] >= error[j]){
					better[j] = better[j] + 1;
				}
				}
			}
		int k = 0;
		for(int i = 0; i < dim; i++){
			if(better[i] >= (int) (((double) dim) - 1)/2.0){
				fiftyfifty[i] = true;
				k++;
			}
		}
		ranked = new int[k];
		int counter = 0;
		for(int i = 0; i < dim; i++){
			if(fiftyfifty[i]){
				ranked[counter] = i;
				counter++;
			}
		}
	}
	
	/**
	 * Setter, sets this field {@link #fiftyfifty} by pairwise comparing this array of {@link #Individuals} and
	 * setting {@link #fiftyfifty} true, if this Individual has a better Error value than the worst approx. 50 %.
	 * Main difference to the setter method {@link #bestFifty(boolean)} is, that the metric values are used as
	 * weights in order to ensure convergence of the error value as well as the contact maps.
	 * @param dummy - dummy parameter, to distinguish this method from {@link #bestFifty(boolean)}
	 * @param CMDM - determines which Error is taken into account: CMDM = true -> CMError, false -> DMError
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 */
	public void bestFifty (String dummy, boolean CMDM) {//throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = getSize();
		double[] error = new double[dim];
		int[] better = new int[dim];
		double[] avmetrics = new double[dim];
		System.arraycopy(getAvMetrics2(), 0, avmetrics, 0, dim);
		double metrixsum = entrySum(avmetrics);
		if(CMDM){
			System.arraycopy(Individuals.getCM(getPop()),0,error,0,dim);
		}
		else{
			System.arraycopy(Individuals.getDM(getPop()),0,error,0,dim);
		}
		for(int i = 0; i < dim - 1; i++){
			for(int j = i + 1; j < dim; j++){
				if(avmetrics[i]*error[i]/metrixsum <= avmetrics[j]*error[j]/metrixsum){
					better[i] = better[i] + 1;
				}
				if(avmetrics[i]*error[i]/metrixsum >= avmetrics[j]*error[j]/metrixsum){
					better[j] = better[j] + 1;
				}
				}
			}
		int k = 0;
		for(int i = 0; i < dim; i++){
			if(better[i] >= (int) (((double) dim) - 1)/(2.0)){
				fiftyfifty[i] = true;
				k++;
			}
		}
		ranked = new int[k];
		int counter = 0;
		for(int i = 0; i < dim; i++){
			if(fiftyfifty[i]){
				ranked[counter] =  i;
				counter++;
			}
		}
	}
	
	/**
	 * similar method to that of <code>{@link #bestFifty(String,boolean)}</code>.
	 * @param CMDM
	 * @param dummy
	 */
	public void bestFifty (boolean CMDM, String dummy) {//throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = getSize();
		double[] error = new double[dim];
		int[] better = new int[dim];
		double[] avmetrics = new double[dim];
		System.arraycopy(getAvMetrics2(), 0, avmetrics, 0, dim);
		double metrixsum = entrySum(avmetrics);
		if(CMDM){
			System.arraycopy(Individuals.getCM(getPop()),0,error,0,dim);
		}
		else{
			System.arraycopy(Individuals.getDM(getPop()),0,error,0,dim);
		}
		for(int i = 0; i < dim - 1; i++){
			for(int j = i + 1; j < dim; j++){
				double weight_i = avmetrics[i]*metrixsum;
				double weight_j = avmetrics[j]*metrixsum;
				if(weight_i*error[i] <= weight_j*error[j]){
					better[i] = better[i] + 1;
				}
				if(weight_i*error[i] >= weight_j*error[j]){
					better[j] = better[j] + 1;
				}
			}
		}
		int k = 0;
		for(int i = 0; i < dim; i++){
			if(better[i] >= (int) (((double) dim) - 1)/(2.0)){
				fiftyfifty[i] = true;
				k++;
			}
		}
		ranked = new int[k];
		int counter = 0;
		for(int i = 0; i < dim; i++){
			if(fiftyfifty[i]){
				ranked[counter] =  i;
				counter++;
			}
		}
	}
	
	/**
	 * This auxiliary method ranks all instances of the field <code>{@link #pop}</code> according to their DMError, weighted by
	 * the CMError. So the special feature of this method is, that both error values of Individuals instance is used in order
	 *  to rank them.
	 * @param dummy - a dummy parameter, to distinguish this method from similarly named methods
	 */
	public void bestFifty (String dummy){
		if(pop != null){
			int pop_size = size, counter = 0, counter1 = 0;
			int[] comp_array = new int[pop_size];
			for(int i = 0; i < pop_size - 1; i++){
				for(int j = i + 1; j < pop_size; j++){
					double sum_of_error = getAvCMError() * ((double) pop_size);
					double frac1 = getPop(i).getCM()/sum_of_error, frac2 = getPop(j).getCM()/sum_of_error; 
					double error1 = frac1*getPop(i).getDM(), error2 = frac2*getPop(j).getDM();
					if(error1 >= error2){
						comp_array[i] += 1;
					}
					else{
						comp_array[j] += 1;
					}
				}
			}
			for(int i = 0; i < pop_size; i++){
				if(comp_array[i] >= (int)((double) pop_size - 1)/2.0){
					fiftyfifty[i] = true;
					counter++;
				}
			}
			ranked = new int[counter];
			for(int i = 0; i < pop_size; i++){
				if(fiftyfifty[i]){
					ranked[counter1] = i;
					counter1++;
				}
			}
		}
	}
	
	/**
	 * setter, sets this 'gen' field, shall only be used in an evolutionary context
	 * @param m
	 */
	public void setGen(int m) {
		gen = m;
	}
	
	/**
	 * setter, sets 'fiftyfifty' field by checking whether this Species instance's error value, specified by 'CMDM',
	 * is less or equal to this Speciess average error, CMDM = true: CMError, CMDM = false: DMError
	 * @param CMDM
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 */
	public void compareToAverage(boolean CMDM) {// throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = getSize();
		double avError = 0.0;
		double[] Errors = new double[dim];
		if(CMDM){
			avError = getAvCMError();
			System.arraycopy(getErrorArray(true), 0, Errors, 0, dim);
		}
		else{
			avError = getAvDMError();
			System.arraycopy(getErrorArray(false), 0, Errors, 0, dim);
		}
		for(int i = 0; i < dim; i++){
			if(Errors[i] <= avError){
				fiftyfifty[i] = true;
			}
		}
	}
	
	/**
	 * method, to compare all fields 'storer' retaining only contacts all have in common
	 * @param p
	 */
	public void setSubHash(Individuals[] p){
		int dim = p.length;
		subhash = new HashSet<Pair<Integer>>(p[0].getHashSet()); 
		for(int i = 1; i < dim; i++){
			subhash.retainAll(p[i].getHashSet());
		}
	}
	
	/*
	 * setter, counting the frequency of all contacts in this Species by instantiating the field 'weighted'
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 
	public void setWHash (){// throws SQLException, PdbCodeNotFoundError, PdbLoadError {
		int counter = 0, numConts = getNumOfContacts();
		int dim = getSize(), dim2 = dim*numConts;
		int[][] weighter = new int[dim2][3];
		HashSet<Pair<Integer>> hashs = new HashSet<Pair<Integer>>();
		for(int i = 0; i < dim; i++){																//looping through the HashSet field 'store'
			HashSet<Pair<Integer>> hash = new HashSet<Pair<Integer>>(getPop(i).getHashSet());	//a new HashSet to check for contacts already present
			Iterator<Pair<Integer>> it = hash.iterator();
			while(it.hasNext()){
				Pair<Integer> pa = new Pair<Integer>(it.next());
				if(hashs.add(pa)){																	//if the contact not present in the new HashSet
					weighter[counter][0] = 1;
					weighter[counter][1] = pa.getFirst().intValue();
					weighter[counter][2] = pa.getSecond().intValue();
					counter++;
				}
				else{																				//if the contact already present in the HashSet
					for(int k = 0; k < counter; k++){												//a loop checks in 'weighter' where entry was placed 
						if(weighter[k][1] == pa.getFirst() && weighter[k][2] == pa.getSecond()){	//and increments the zeroths poistion of 'weighter'
							weighter[k][0] = 1 + weighter[k][0];
							break;
						}
					}
				}
			}
		}
		weighted = new int[counter][3];
		for(int i = 0; i < counter; i++){															//initializing the field 'weighted'
			weighted[i][0] = weighter[i][0];
			weighted[i][1] = weighter[i][1];
			weighted[i][2] = weighter[i][2];
		}
		//sortWeighted();
	}*/
	
	/**
	 * setter, setting the field 'weighteed'. Since both, 'weighted' and 'weighteed' are doing
	 * the same, this field should be favoured over the 'weighted' one, because looking up
	 * is much quicker...
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 */
	public void setWHasher (){// throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = getSize(), compare = dim * getPop(0).getNumOfContacts();
		weighted = new HashMap<Pair<Integer>, Integer> (compare);					//initializing the field 'weighteed'
		for(int i = 0; i < dim; i++){
			HashSet<Pair<Integer>> hash1 = getPop(i).getHashSet();					//HashSet of the i-th Individual
			
			Iterator<Pair<Integer>> it = hash1.iterator();								//Iterator over this HashSet
			
			while(it.hasNext()){														//iteration over the HashSet
				
				int counter = 1;														//default value, each contact is at least once present
				
				Pair<Integer> pair = it.next();
				if(weighted.containsKey(pair)){									//if the Pair of Integers is already present in the HashMap
																						//the value is incremented according to the times of occurrence
					
					counter = weighted.get(pair).intValue() + 1;					
					weighted.put(pair, new Integer(counter));
				}
				else{
					weighted.put(pair, new Integer(counter));						//if the Pair of Integers is not present in the HashMap the default
																						//occurrence value (1) is used
				}
			}
		}
	}
	
	/**
	 * method, using Sathya's preselected set as input and generating random samples by picking 
	 * random contacts specified by 'numofConts' parameter 
	 * @param popsize
	 * @param startfiles
	 * @param endfiles
	 * @param numofCons
	 * @throws SQLException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public void sathyasSet(int popsize, int startfiles, int numofCons) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		String dir1 = "/project/StruPPi/gabriel/workspace_old/aglappe/embed/teststructures/";
		File dir = new File(dir1);
		File[] test = dir.listFiles(new RegexFileFilter(".*_all\\.graph"));
		Arrays.sort(test);
		//File neu = test[2];
		BufferedReader fileIn = new BufferedReader(new FileReader(dir1+"sathya_set.list"));		//reading the contact files of Sathya
		String line;
		RIGraph fullCM = new RIGraph();
		MySQLConnection conn = new MySQLConnection();
		int counter = 0;
		while((line = fileIn.readLine()) != null) {												//reading loop
			if(counter < startfiles){
				counter++;
			}
			else{
				if(counter == startfiles){			
					String[] tokens = line.split(" ");
					String pdbCode=tokens[0];
					String pdbChainCode=tokens[1];
					Pdb pdb = new PdbasePdb(pdbCode, "pdbase_20090728", conn);
					pdb.load(pdbChainCode);
					fullCM = pdb.getRIGraph("Ca", 9);
					break;//counter++;
				}
			}
		}
		fileIn.close();
		pop = new Individuals[popsize];													//initialize an array of Individuals with 'popsize' entries
		if(startfiles >= 0){																	//tests, if the parameters 'startfiles' and 'endfiles' generate
																								//any exceptions
			Individuals.setFullContactMap(fullCM);
			for(int i = 0; i < popsize; i++){
				pop[i] = new Individuals(Individuals.fullcontactmap, conn, true, numofCons);
			}
		}		
		else{
			if(startfiles < 0){																	//if some wrong parameter are used...
				conn.close();
				String notzero = "Parameter 'startfiles' = " + startfiles + " must never be less than zero!";
				throw new IllegalArgumentException("\n" + notzero);
			}
		}
		conn.close();	
	}
	
	
	/**
	 * setter method introducing a metric d : X -> [0,1], where X denotes the contact space and
	 * d is the function from X to the unity intervall. If two Individuals x and y are having the
	 * same contacts than d(x,y) = 0; whereas x and y have no contact in common d(x,y) = 0. In
	 * this sense X is a subset of the power set P(T) of T, where T is a finite set {1,...,n} x {1,...,n}.
	 * Accordingly, all the elements of X, being subsets of T, are having the same number of elements. 
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 */
	public void setMetrics() {//throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		/*int dim = getSize();
		int compare = (int) ((double) (dim*(dim - 1))*2.0/3.0);
		metric = new HashMap<Pair<Integer>, Double>(compare);
		for(int i = 0; i < dim - 1; i++){
			for(int j = i + 1; j < dim; j++){
				HashSet<Pair<Integer>> hash1 = getPop(i).getHashSet();						//HashSet instance of the i-th Individual
				
				HashSet<Pair<Integer>> hash2 = getPop(j).getHashSet();						//HashSet instance of the j-th Individual
				int counter = 0, hashsize = hash1.size();
				
				if(hash2.retainAll(hash1)){														//all common contacts are retained in 'hash2'
					
					counter = hash2.size();														//the number of common contacts are given by the
																								//size of 'hash2'
					
				}
				Pair<Integer> pair = new Pair<Integer> (new Integer (i), new Integer (j));		//the indices are use for a Pair<Integer> instance
				
				Double metrics = new Double (1.0 - ((double) counter/ ((double) hashsize)));			//the number of the common contacts between i and j
																								//are devided by the number of all contacts and substracted
																								//from one
				
				metric.put(pair, metrics);													//initializing the field 'metric'
			}
		}*/
		metric2 = new Metric(this);
	}
	
	/**
	 * method, to write field 'weighted' to a file named via PDB code plus generation and number of contacts with '.cmap' as file extension,
	 * <p>
	 * the standard output is as follows:
	 * </p>
	 * <p>
	 * the header: with sequence, pdb code, chain code, generation, error values and deviations,
	 * </p>
	 * <p>
	 * an the actual contact map table, where the first and second column correspond to the first and second contact pair,
	 * the third column is the frequency and the fourth is the index of the <tt>Individuals</tt>, where this contact occurs.
	 * </p>
	 * @param path a String denoting the output directory
	 * @param name a String specifying the output file name
	 * @throws IOException
	 * @throws FileNotFoundException
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 */
	public void printToFile(String path, String name) throws IOException, FileNotFoundException {//SQLException, PdbCodeNotFoundError, PdbLoadError {
		File dirtest = new File(path);
		if(!dirtest.exists()){
			dirtest.mkdirs();
		}
		FileOutputStream file = new FileOutputStream(path+name+"-"+getNumOfContacts()+".cmap");
		PrintStream printa = new PrintStream(file);
		printa.print(toString());
		printa.close();
		file.close();
		System.out.println(getName()+" at generation "+gen+" written to file...");
	}
	
	/**
	 * Writes an instance of this class to a 'cmap' file. The major difference of this method compared to <code>{@link #printToFile(String, String)}</code>
	 * is, that the indexing column is obsolete, which is necessary for usage of this file in 'cmview'.
	 * @param path
	 * @param name
	 * @param dummy
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	public void printToFile(String path, String name, String dummy) throws IOException, FileNotFoundException {//SQLException, PdbCodeNotFoundError, PdbLoadError {
		File dirtest = new File(path);
		if(!dirtest.exists()){
			dirtest.mkdirs();
		}
		FileOutputStream file = new FileOutputStream(path+name+"-"+getNumOfContacts()+".cmap");
		PrintStream printa = new PrintStream(file);
		printa.print(toString(dummy));
		printa.close();
		file.close();
		System.out.println(getName()+" at generation "+gen+" written to file...");
	}
	
	/**
	 * additional printer method: writes the field <code>{@link #metric}</code> and <code>{@link #metric2}</code> to a
	 * predefined file.
	 * @param path
	 * @param filename
	 * @param k
	 * @throws IOException
	 * @throws FileNotFoundException
	 */
	public void printMetric (String path, String filename, int k) throws IOException, FileNotFoundException {
		String cont = "#Metric File " + getName() + "\n";
		HashMap<Pair<Integer>, Integer> hashmap2 = getMet();
		HashSet<Pair<Integer>> hashset = getMetKey();
		Iterator<Pair<Integer>> it = hashset.iterator();
		cont = cont + "\n";
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			cont = cont + pair.getFirst().intValue() + "\t" + pair.getSecond().intValue() + "\t" + hashmap2.get(pair) + "\n"; 
		}
		FileOutputStream file = new FileOutputStream(path + filename + "-"+k+getNumOfContacts()+".met");
		PrintStream printa = new PrintStream(file);
		printa.print(cont);
		printa.close();
		file.close();
		System.out.println(getName()+" at generation "+gen+" written to file...");
	}
	
	/**
	 * additional printer method: writes all <tt>Individuals</tt> instances in this <code>{@link #pop}</code> to
	 * a predefined file.
	 * @param dir
	 * @param addname
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public void printAllIndisToTempDir (String dir, String addname) throws FileNotFoundException, IOException{
		int length = getSize();
		for(int i = 0; i < length; i++){
			Integer index = new Integer (i); 
			getPop(i).printToFile(dir, addname, index.toString());
		}
	}
	
	/**
	 * auxiliary setter: sets the field <code>{@link #indexer}</code> by mapping all contacts of each Individuals
	 * in the field <code>{@link #pop}</code> onto their array index. So, any output of an instance of this class,
	 * has a table with all contacts and the corresponding Individuals array index.
	 */
	public void setIndexer (){
		if(pop != null && pop.length > 0){
			indexer = new HashMap<Pair<Integer>,HashSet<Integer>> ();
			for(int i = 0; i < size; i++){
				HashSet<Pair<Integer>> in_set = getPop(i).getHashSet();
				Iterator<Pair<Integer>> it = in_set.iterator();
				Integer index = new Integer (i);				
				while(it.hasNext()){
					Pair<Integer> pair = it.next();
					if(indexer.containsKey(pair)){
						HashSet<Integer> subset = indexer.get(pair);
						subset.add(index);
						indexer.put(pair, subset);
					}
					else{
						HashSet<Integer> subset = new HashSet<Integer> ();
						subset.add(index);
						indexer.put(pair, subset);
					}
				}
			}
		}
	}
	
	/**
	 * Setter: reads a standard output file of this class and instantiates an object
	 * of this class by taking all information contained in the file.
	 * @param file
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws PdbCodeNotFoundError
	 * @throws SQLException
	 * @throws PdbLoadError
	 * @throws IllegalArgumentException
	 */
	public void readFile (File file) throws FileNotFoundException, IOException, PdbCodeNotFoundError, SQLException, PdbLoadError, IllegalArgumentException {
		if(file.exists() && file.getAbsolutePath().contains(".cmap")){
			BufferedReader reader = new BufferedReader (new FileReader (file));
			String linereader = null, chain_code = null, seq = null;
			HashMap<Integer,HashSet<Pair<Integer>>> contact_map = new HashMap<Integer,HashSet<Pair<Integer>>>();
			while((linereader = reader.readLine()) != null){ 
				if(linereader.contains("#")){
					if(linereader.contains("SEQUENCE")){
						seq = linereader.split(": ")[1];
					}
					if(linereader.contains("PDB") && !linereader.contains("CHAIN CODE")){
						setPopName(linereader.split(": ")[1]);
					}
					if(linereader.contains("PDB CHAIN CODE")){
						chain_code = linereader.split(": ")[1];
					}
					if(linereader.contains("GENERATION")){
						int m = (int) Double.parseDouble(linereader.split(": ")[1]);
						setGen(m);
					}
					if(linereader.contains("Species SIZE")){
						int size = (int) Double.parseDouble(linereader.split(": ")[1]);
						pop = new Individuals[size];
						this.size = size;
					}
					if(linereader.contains("CMError")){
						CMError = Double.parseDouble(linereader.split(": ")[1]);
					}
					if(linereader.contains("CMstDev")){
						CMstdev = Double.parseDouble(linereader.split(": ")[1]);
					}
					if(linereader.contains("DMError")){
						DMError = Double.parseDouble(linereader.split(": ")[1]);
					}
					if(linereader.contains("DMstDev")){
						DMstdev = Double.parseDouble(linereader.split(": ")[1]);
					}
				}
				else{
					if(linereader.contains("\t")){
						String[] ar = linereader.split("\t");
						Integer f_val = new Integer ((int) Double.parseDouble(ar[0]));
						Integer s_val = new Integer ((int) Double.parseDouble(ar[1]));
						Pair<Integer> pair = new Pair<Integer> (f_val,s_val);
						Integer index = new Integer((int) Double.parseDouble(ar[3]));
						if(contact_map.containsKey(index)){
							HashSet<Pair<Integer>> subset = contact_map.get(index);
							subset.add(pair);
							contact_map.put(index, subset);
						}
						else{
							HashSet<Pair<Integer>> subset = new HashSet<Pair<Integer>> ();
							subset.add(pair);
							contact_map.put(index, subset);
						}
					}
				}
			}
			Set<Integer> keys = contact_map.keySet();
			Iterator<Integer> it = keys.iterator();
			
			while(it.hasNext()){
				Integer index = it.next();
				HashSet<Pair<Integer>> cm = contact_map.get(index);
				int i = index.intValue();
				pop[i] = new Individuals ();
				pop[i].setName(getName());
				pop[i].storer(cm);
				pop[i].setChainCode(chain_code);
				pop[i].setSequence(seq);
				if(Individuals.fullcontactmap == null){
					Individuals.setFullContactMap(pop[i].reconstructGraph());
				}
				pop[i].setEntries(cm);
				pop[i].setErrorValues();
				pop[i].setFullContact(Individuals.fullcontactmap.getEdgeCount());
				pop[i].setNumOfContacts(cm);
				
			}
			fiftyfifty = new boolean [keys.size()];
			setMetrics();
			setIndexer();
			setSize(pop.length);
			setWHasher();
		}
		else{
			if(!file.exists()){
				throw new FileNotFoundException ("The denoted file does not exist.");
			}
			else{
				throw new IllegalArgumentException ("The only file type, accepted by this method, must have '.cmap' extension.");
			}
		}
	}
	
	/**
	 * <p>
	 * a method checking whether the size of this instance is greater or equal to 14. This is important since any
	 * evolution run will cause a NullPointerException, if the size is below 14. That is because of the following:
	 * </p>
	 * <p>
	 * </p>
	 * <p>
	 * any of the methods <code>{@link #bestFifty(...)}</code> will only use the best 50 % of the Individuals, so the maximal number
	 * of bred Individuals <tt>n</tt> will be <tt>n (n - 1)/2</tt>. Since this number is added, later, one gets: 
	 * </p>
	 * <p>
	 * <tt>n (n - 1)/2 + n = n (n + 1)/2</tt> as intermediate offspring. The intermediate offspring is ranked again and only the best fifty percent
	 * is taken as next generation, so one gets:
	 * </p>
	 * <p>
	 * <tt>
	 */
	public void checkSize (){
		if(size < 14){
			System.err.println("For any evolutionary run, the size of the deme must be greater than 14, but was set to" + size + "!");
			System.exit(1);
		}
	}
	
	/**
	 * 
	 * @param pdb_code
	 * @param pdb_db
	 * @param size
	 * @param percent
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @throws IllegalArgumentException
	 */
	public void setDemes (String pdb_code, String pdb_db, int size, double percent) throws SQLException, PdbCodeNotFoundError, PdbLoadError, IllegalArgumentException {
		if(size > 0 && (percent <= 100.0 && percent > 0.0) ){
			MySQLConnection conn = new MySQLConnection ();
			Pdb pdb = new PdbasePdb(pdb_code, pdb_db, conn);
			pdb.load("A");
			RIGraph rig = pdb.getRIGraph("Ca", 9.0);
			pop = new Individuals [size];
			for(int i = 0; i < size; i++){
				pop[i] = new Individuals (rig, conn, true, percent);
			}
			fiftyfifty = new boolean[size];
			conn.close();
			setCMErrorStats(pop);
			setDMErrorStats(pop);
			setMetrics();
			setSize2(size);
			setWHasher();
			setIndexer();
			setSubHash(pop);
			setPopName(rig.getPdbCode());
		}
		else {
			if(size <= 0) throw new IllegalArgumentException ("The size must never be less than or equal to zero.");
			else if (percent > 100.0 || percent <= 0.0) throw new IllegalArgumentException ("Parameter 'percent' must be greater than zero and less than or equal to 100.0.");
		}
	}
	
	/*--------------------------------------Getters-------------------------------------------*/
	
	/**
	 * getter, returns this Individuals array
	 * @return p - an array of Individuals
	 */
	public Individuals[] getPop () {
		int dim = getSize();
		Individuals[] p = new Individuals[dim];
		for(int i = 0; i < dim; i++){
			p[i] = new Individuals(pop[i]); 
		}
		return p;
	}
	
	/**
	 * 
	 * @param i
	 * @return p - instance of 'Individuals'
	 */
	public Individuals getPop (int i) {
		return new Individuals(pop[i]);
	}
	
	/**
	 * getter returns this name
	 * @return a String representing the pdb code
	 */
	public String getName() {
		return new String (name);
	}
	
	/**
	 * getter, returns this average CMError
	 * @return the CMError
	 */
	public double getAvCMError () {
		return CMError;
	}
	
	/**
	 * getter, returns this average DMError
	 * @return the DMError
	 */
	public double getAvDMError () {
		return DMError;
	}
	
	/**
	 * getter, returns this standard deviation of CMError
	 * @return the CMError standard deviation
	 */
	public double getCMstdev () {
		return CMstdev;
	}
	
	/**
	 * getter, returns this standard deviation of DMError
	 * @return the DMError standard deviation
	 */
	public double getDMstdev () {
		return DMstdev;
	}
	
	/**
	 * getter, returns this Demes size
	 * @return
	 */
	public int getSize(){
		return size;
	}
	
	/**
	 * getter, returns this 'fiftyfifty' field at the i-th position
	 * @param i
	 * @return
	 */
	public boolean getFifty(int i){
		return fiftyfifty[i];
	}
	
	/**
	 * getter, returns this instances field 'fiftyfifty' as an array of boolean
	 * @return fifty
	 */
	public boolean[] getFifty(){
		int dim = getSize();
		boolean[] fifty = new boolean[dim];
		System.arraycopy(fiftyfifty, 0, fifty, 0, dim);
		return fifty;
	}
	
	/**
	 * returns an int array (an array of indices) with all Individuals that are best ranked.
	 * @return
	 */
	public int[] getFiftyArray () {
		int counter = 0, dim1 = getSize();
		for(int i = 0; i < dim1; i++){
			if(getFifty(i)){
				counter++;
			}
		}
		int[] array = new int[counter];
		int j = 0;
		for(int i = 0; i < dim1; i++){
			if(getFifty(i)){
				array[j] = i; 
			}
		}
		return array;
	}
	
	/**
	 * getter, returns this field 'n
	 * @return
	 */
	public int getNumOfContacts(){
		return pop[0].getNumOfContacts();
	}
	
	/**
	 * getter, returns the number of generation
	 * @return 'gen'
	 */
	public int getGen() {
		return gen;
	}
	
	/**
	 * method, returing CMError or DMError, specified by the boolean parameter:
	 * true = CMError() method called, false = DMError() method called (method of the Individuals class)
	 * @param CMDM
	 * @return
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 */
	public double[] getErrorArray(boolean CMDM) {//throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = getSize();
		double[] Errors = new double[dim];
		if(CMDM){
			for(int i = 0; i < dim; i++){
				Errors[i] = getPop(i).getCM();
			}
		}
		else{
			for(int i = 0; i < dim; i++){
				Errors[i] = getPop(i).getDM();
			}
		}
		return Errors;
	}
	
	/**
	 * method returning an array of integers, that are the indices of the best ranked 'Individuals' instances of this 'Species' instance
	 * depending on the compare method applied to the 'fiftyfifty' field (i.e. bestFifty(boolean) or compareToAverage(boolean))
	 * @param CMDM
	 * @return
	 */
	public int[] getCompIndis(boolean CMDM){
		int dim = getSize(), counter = 0;
		int[] comps1 = new int[dim];
		for(int i = 0; i < dim; i++){
			if(getFifty(i)){
				comps1[counter] = i;
				counter++;
			}
		}
		int[] comps2 = new int[counter];
		System.arraycopy(comps1,0,comps2,0,counter);
		return comps2;
	}
	
	/**
	 * method returning the field 'ranked'
	 * @return
	 */
	public int[] getRank(){
		return ranked;
	}
	
	/**
	 * method returning the i-th entry of field 'ranked' 
	 * @param i
	 * @return
	 */
	public int getRank (int i){
		return ranked[i];
	}
	
	/**
	 * canonical <code>toString()</code>. Converts an instance of this class to the standard CMView output
	 * format. That is, a header with all essential information of the protein as pdb code, sequence etc. and
	 * a table with all contact pairs.
	 */
	public String toString (){
		Individuals ne = new Individuals(getPop(0));
		String cont = "#CMVIEW GRAPH FILE ver: 1.0\n#SEQUENCE: "+ne.getSequence()+"\n"+
		"#PDB: "+ne.getName()+ "\n#PDB CHAIN CODE: "+ne.getChainCode()+"\n#CT: "+Individuals.getContactT()+ "\n#CUTOFF: "+Individuals.getContactDist()+"\n"+ 
		"#GENERATION: "+getGen() + "\n#Species SIZE: " + getSize() + "\n#CMError: "+ getAvCMError() + "\n#CMstDev: " + getCMstdev() + "\n#DMError: "+ getAvDMError() + 
		"\n#DMstDev: " + getDMstdev() + "\n#NUMB. CONTACTS: " + ne.getNumOfContacts() +"\n" + "#NUMB. OF ALL CONTACTS: " + ne.getFullContact() + "\n#NUMB. OF CONTS. IN POP: " + weighted.size() + "\n";
		int[][] index_array = sortWeighted2();
		int length = index_array[0].length;
		HashMap<Pair<Integer>, Integer> hashmap = getMap();
		cont = cont + "\n";
		for(int i = 0; i < length; i++){
			int index1 = index_array[0][i], index2 = index_array[1][i];
			Pair<Integer> pair = new Pair<Integer> (new Integer(index1),new Integer (index2));
			Integer[] ar = getIndexer(pair);
			int length2 = ar.length;
			for(int j = 0; j < length2; j++){
				cont += index1 + "\t" + index2 + "\t" + ((double) hashmap.get(new Pair<Integer>(new Integer(index1),new Integer(index2))))/((double) getSize()) + "\t" + ar[j] + "\n"; 
			}
		}
		return cont;
	}

	public String toString (String dummy){
		Individuals ne = new Individuals(getPop(0));
		String cont = "#CMVIEW GRAPH FILE ver: 1.0\n#SEQUENCE: "+ne.getSequence()+"\n"+
		"#PDB: "+ne.getName()+ "\n#PDB CHAIN CODE: "+ne.getChainCode()+"\n#CT: "+Individuals.getContactT()+ "\n#CUTOFF: "+Individuals.getContactDist()+"\n"+ 
		"#GENERATION: "+getGen() + "\n#Species SIZE: " + getSize() + "\n#CMError: "+ getAvCMError() + "\n#CMstDev: " + getCMstdev() + "\n#DMError: "+ getAvDMError() + 
		"\n#DMstDev: " + getDMstdev() + "\n#NUMB. CONTACTS: " + ne.getNumOfContacts() +"\n" + "#NUMB. OF ALL CONTACTS: " + ne.getFullContact() + "\n#NUMB. OF CONTS. IN POP: " + weighted.size() + "\n";
		int[][] index_array = sortWeighted2();
		int length = index_array[0].length;
		HashMap<Pair<Integer>, Integer> hashmap = getMap();
		cont = cont + "\n";
		for(int i = 0; i < length; i++){
			int index1 = index_array[0][i], index2 = index_array[1][i];
			Pair<Integer> pair = new Pair<Integer> (new Integer(index1),new Integer (index2));
			Integer[] ar = getIndexer(pair);
			int length2 = ar.length;
			for(int j = 0; j < length2; j++){
				cont += index1 + "\t" + index2 + "\t" + ((double) hashmap.get(new Pair<Integer>(new Integer(index1),new Integer(index2))))/((double) getSize()) +"\n"; 
			}
		}
		return cont;
	}
	/**
	 * returns a HashMap with all contacts and their corresponding frequency in this deme.
	 * @return a HashMap of all contact pairs mapped onto their frequency
	 */
	public HashMap<Pair<Integer>, Integer> getMap(){
		return new HashMap<Pair<Integer>, Integer> (weighted); 
	}
	
	/**
	 * returning a HashSet of all contact pair present in this <code>{@link Demes}</code>
	 * @return a HashSet of all contact pairs
	 */
	public HashSet<Pair<Integer>> getKey(){
		return new HashSet<Pair<Integer>> (weighted.keySet());
	}
	
	/*
	 * returning the metric instance, which corresponds to a table of all Individuals where
	 * a pairwise comparison is made 
	 * @return a HashMap of all Individuals <tt>i</tt> and <tt>j</tt> mapped onto their distance
	 
	public HashMap<Pair<Integer>, Double> getMet(){
		return new HashMap<Pair<Integer>, Double> (metric);
	}*/
	
	/*
	 * returning a HashSet of all Individuals in the field <code>{@link #pop}</code>, which were
	 * compared in order to generate a metric
	 * @return a HashSet of Pairs of all Individuals indices
	 
	public HashSet<Pair<Integer>> getMetKey(){
		return new HashSet<Pair<Integer>> (metric.keySet());
	}*/
	
	/**
	 * returns the field <code>{@link #metric2}</code> as a HashMap
	 * @return
	 */
	public HashMap<Pair<Integer>, Integer> getMet(){
		return new HashMap<Pair<Integer>, Integer> (metric2.getMetMap());
	}
	
	/**
	 * returns the key of the field <code>{@link #matric2}</code>
	 * @return
	 */
	public HashSet<Pair<Integer>> getMetKey(){
		return new HashSet<Pair<Integer>> (metric2.getMetMap().keySet());
	}
	
	/*
	 * compares the pairwise distances given by the field <code>{@link #metric}</code>
	 * and computes an average distance value for each entry. 
	 * @return a double array representing the average distance
	 
	public double[] getAvMetrics(){
		int dim = getSize(), counter = 0;
		double[] array = new double[dim];
		HashSet<Pair<Integer>> hash = getMetKey();
		Iterator<Pair<Integer>> it1 = hash.iterator();
		Iterator<Pair<Integer>> it2 = hash.iterator();
		while(it1.hasNext()){
			Pair<Integer> pair1 = it1.next();
			while(it2.hasNext()){
				Pair<Integer> pair2 = it2.next();
				if(pair1 != pair2){
					if(pair1.getFirst().intValue() == pair2.getFirst().intValue()){
						counter = pair1.getFirst().intValue();
						array[counter] = array[counter] + getMet().get(pair1).doubleValue() + getMet().get(pair2).doubleValue();
					}
					else{
						if(pair1.getSecond().intValue() == pair2.getSecond().intValue()){
							counter = pair1.getSecond().intValue();
							array[counter] = array[counter] + getMet().get(pair1).doubleValue() + getMet().get(pair2).doubleValue();
						}
						else{
							if(pair1.getFirst().intValue() == pair2.getSecond().intValue()){
								counter = pair1.getFirst().intValue();
								array[counter] = array[counter] + getMet().get(pair1).doubleValue() + getMet().get(pair2).doubleValue();
							}
							else{
								if(pair1.getSecond().intValue() == pair2.getFirst().intValue()){
									counter = pair1.getSecond().intValue();
									array[counter] = array[counter] + getMet().get(pair1).doubleValue() + getMet().get(pair2).doubleValue();
								}
							}
						}
					}
				}
			}
			it2 = hash.iterator();
		}
		double[] array1 = new double[dim];
		for(int i = 0; i < dim; i++){
			array1[i] = array[i]/((double) dim - 1); 
		}
		return array1;
	}*/
	
	/**
	 * compares the pairwise distances given by the field <code>{@link #metric2}</code>
	 * and computes an average distance value for each entry. 
	 * @return a double array representing the average distance
	 */
	public double[] getAvMetrics2(){
		int dim = getSize(), counter = 0;
		double[] array = new double[dim];
		HashSet<Pair<Integer>> hash = getMetKey();
		Iterator<Pair<Integer>> it1 = hash.iterator();
		Iterator<Pair<Integer>> it2 = hash.iterator();
		while(it1.hasNext()){
			Pair<Integer> pair1 = it1.next();
			while(it2.hasNext()){
				Pair<Integer> pair2 = it2.next();
				if(pair1 != pair2){
					if(pair1.getFirst().intValue() == pair2.getFirst().intValue()){
						counter = pair1.getFirst().intValue();
						array[counter] = array[counter] + getMet().get(pair1).intValue() + getMet().get(pair2).intValue();
					}
					else{
						if(pair1.getSecond().intValue() == pair2.getSecond().intValue()){
							counter = pair1.getSecond().intValue();
							array[counter] = array[counter] + getMet().get(pair1).intValue() + getMet().get(pair2).intValue();
						}
						else{
							if(pair1.getFirst().intValue() == pair2.getSecond().intValue()){
								counter = pair1.getFirst().intValue();
								array[counter] = array[counter] + getMet().get(pair1).intValue() + getMet().get(pair2).intValue();
							}
							else{
								if(pair1.getSecond().intValue() == pair2.getFirst().intValue()){
									counter = pair1.getSecond().intValue();
									array[counter] = array[counter] + getMet().get(pair1).intValue() + getMet().get(pair2).intValue();
								}
							}
						}
					}
				}
			}
			it2 = hash.iterator();
		}
		double[] array1 = new double[dim];
		for(int i = 0; i < dim; i++){
			array1[i] = array[i]/((double) dim - 1); 
		}
		return array1;
	}
	
	/**
	 * converts the key set of the field <code>{@link #weighted}</code> into a
	 * sorted integer matrix. The first value corresponds to the first contact index
	 * and the second one to the second contact index.
	 * @return a sorted integer matrix 
	 */
	public int[][] sortWeighted2(){
		HashSet<Pair<Integer>> set = getKey();
		int[][] sorted_array = SortIntArray.converter(set);
		return sorted_array;
	}
	
	/**
	 * this method extracts all contacts with a frequency above <tt>threshold</tt>
	 * @param threshold a double value less than 1 and greater than 0, specifying which contacts will be
	 * considered
	 * @return a HashMap with the contacts and their corresponding frequency
	 */
	public HashMap<Pair<Integer>,Double> contactsWithHighestFrequency (double threshold){
		HashMap<Pair<Integer>,Integer> map = new HashMap<Pair<Integer>,Integer> (weighted);
		HashMap<Pair<Integer>,Double> subhash = new HashMap<Pair<Integer>,Double> ();
		HashSet<Pair<Integer>> keyset = new HashSet<Pair<Integer>> (map.keySet());
		Iterator<Pair<Integer>> it = keyset.iterator();
		double num_of_conts = (double) getPop(0).getNumOfContacts();
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			double percent = ((double) (map.get(pair)).doubleValue())/num_of_conts;
			Double per_val = new Double(percent);
			if(percent >= threshold){
				subhash.put(pair, per_val);
			}
		}
		return subhash;
	}
	
	/**
	 * computes the average distance in this <code>{@link Demes}</code>
	 * @return a double value representing the average distance of all Individuals in this instance
	 */
	public double getAverageMetric2 (){
		HashMap<Pair<Integer>, Integer> met_hash = getMet();
		Set<Pair<Integer>> keyset = met_hash.keySet();
		Iterator<Pair<Integer>> it = keyset.iterator();
		int size = met_hash.size();
		double avmetric = 0.0;
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			avmetric += met_hash.get(pair).doubleValue();
		}
		return avmetric/size;
	}
	
	/**
	 * sanity checker: during each run of this class one must always be sure, that only 
	 * contact maps of the same protein (matching pdb codes) are treated at the same time.
	 * @return true, if all pdb codes match in the field <tt>pop</tt>
	 */
	public boolean samePDBCode (){
		boolean tester = true;
		String pdbcode = "";
		for(int i = 0; i < size; i++){
			pdbcode = new String (getPop(i).getName());
			if(i > 0){
				if(!pdbcode.matches(getPop(i - 1).getName())){
					tester = false;
					break;
				}
			}
		}
		return tester;
	}
	
	/**
	 * returns all common contacts
	 * @return a HashSet of all common contacts
	 */
	public HashSet<Pair<Integer>> getSubHash (){
		return new HashSet<Pair<Integer>> (subhash);
	}
	
	/**
	 * a method returning the field <code>{@link #indexer}</code> as a HashMap, where all contact 
	 * pairs are the key and a HashSet of Integer s is the value. This HashSet contains all indices
	 * of each <code>{@link Individuals}</code>, where the contact appears. 
	 * @return a HashMap, where each contact pair is mapped onto a Set of indices representing
	 * the Individuals
	 */
	public HashMap<Pair<Integer>,HashSet<Integer>> getIndexer (){
		return new HashMap<Pair<Integer>,HashSet<Integer>> (indexer);
	}
	
	/**
	 * a method returning the value of the field <code>{@link #indexer}</code> as a sorted 
	 * Integer array. That is, each contact pair's index set is returned in a descending order.
	 * @param pair
	 * @return
	 */
	public Integer[] getIndexer (Pair<Integer> pair){
		HashSet<Integer> index_set = indexer.get(pair);
		Integer[] ar = new Integer [index_set.size()];
		ar = index_set.toArray(ar);
		Arrays.sort(ar);
		return ar;
	}
	
	/**
	 * <p>
	 * evolution method: firstly, this method takes all entries in the field <code>{@link #pop}</code> and ranks
	 * them according to the <code>{@link #bestFifty(String)}</code> method. This, in turn, initializes the field <code>{@link #ranked}</code>
	 * which is nothing else then an <tt>int</tt> array with all indices pointing to the best ranked Individuals in <tt>pop</tt>.
	 * Those best ranked Individuals are bred according to the <code>{@link Individuals#breedIndis(Individuals, Individuals)}</code> method, yielding
	 * <tt>n (n - 1)/2</tt> new Individuals, where <tt>n</tt> is the number of all best ranked (the length of the field <tt>ranked</tt>). Thereafter,
	 * all parental Individuals (best ranked) are transferred to the intermediate generation and an additional ranking takes place. Only a random 
	 * selection of the best ranked intermediate generation are transferred to actual next generation.
	 * </p>
	 * <p>
	 * Note, that the random selection is checked for redundancies, that is - each best ranked Individuals is only taken once!
	 * </p>
	 * @param CMDM specifies whether CMError or DMError are selection criterion
	 * @param dummy a dummy parameter
	 * @return an array of Individuals
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */			
	public Individuals[] evolve (boolean CMDM, String dummy) throws SQLException, PdbCodeNotFoundError, PdbLoadError {
		checkSize();
		bestFifty(CMDM);
		//ranking all Individuals in the field 'pop
		
		int dim = getRank().length, counter = 0, length = getSize(), square = nSquaredMinusNHalf(dim);
		Individuals[] off = new Individuals[square + dim];
		//intermediate offspring, all bred Individuals plus the best ranked parental Individuals a placed in this 
		//array of Individuals
		
		for(int i = 0; i < dim - 1; i++){
			//looping over all best ranked Individuals in the field 'pop'
			
			for(int j = i + 1; j < dim; j++){
				int index1 = getRank(i), index2 = getRank(j);
				off[counter] = getPop(index1).breedIndis(getPop(index2));
				//breeding the offspring by calling the method 'breedIndis' from the class Individuals
				
				counter++;
			}
		}
		if(hasNull(off)){
			//check, whether the intermediate Individuals array still contains null entries
			
			int i = 0;
			while(counter < square + dim && i < dim){
				off[counter] = getPop(getRank(i));
				//filling the null entries with the best ranked parental Individuals
				
				counter++; i++;
			}
		}
		System.out.println("# of bred Individuals: " + counter);
		Individuals[] off2 = new Individuals[length];
		//the actual offspring Individuals array
		
		Demes p = new Demes(off);
		//creating an intermediate instance of the class 'Demes'
		
		p.bestFifty(CMDM);
		//ranking them
		
		System.out.println("# of best ranked Individuals: " + p.getRank().length);
		int counter1 = 0;
		HashSet<Integer> indexer = new HashSet<Integer> ();
		while(counter1 < length){
			// a loop that only takes the best ranked Individuals in the field 'pop' of 'p'
			//as a random selection
			
			int rank = p.getRank().length;
			Random rand = new Random();
			int index = rand.nextInt(rank);
			//a random number between zero and the maximal number of best ranked Individuals in the
			//field 'pop' of 'p'
			
			Integer in = new Integer (index);
			if(indexer.add(in)){
				//check, whether the index is already present in the variable 'indexer'
				//only if the index is not present 
				off2[counter1] = new Individuals(p.getPop(p.getRank(index)));
				counter1++;
			}
		}
		return off2;
	}
	
	/**
	 * <p>
	 * evolution method: firstly, this method takes all entries in the field <code>{@link #pop}</code> and ranks
	 * them according to the <code>{@link #bestFifty(String)}</code> method. This, in turn, initializes the field <code>{@link #ranked}</code>
	 * which is nothing else then an <tt>int</tt> array with all indices pointing to the best ranked Individuals in <tt>pop</tt>.
	 * Those best ranked Individuals are bred according to the <code>{@link Individuals#breedIndis(Individuals, Individuals)}</code> method, yielding
	 * <tt>n (n - 1)/2</tt> new Individuals, where <tt>n</tt> is the number of all best ranked (the length of the field <tt>ranked</tt>). Thereafter,
	 * all parental Individuals (best ranked) are transferred to the intermediate generation and an additional ranking takes place. Only a random 
	 * selection of the best ranked intermediate generation are transferred to actual next generation.
	 * </p>
	 * <p>
	 * Note, that the random selection is NOT checked for redundancies. That is the major difference to the method <code>{@link #evolve(boolean, boolean)}</code>.
	 * </p>
	 * @param CMDM specifies whether CMError or DMError are selection criterion
	 * @param dummy a dummy parameter
	 * @return an array of Individuals
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public Individuals[] evolve (String dummy, boolean CMDM) throws SQLException, PdbCodeNotFoundError, PdbLoadError {
		checkSize();
		bestFifty(CMDM);
		//ranking all Individuals in the field 'pop
		
		int dim = getRank().length, counter = 0, length = getSize(), square = nSquaredMinusNHalf(dim);
		Individuals[] off = new Individuals[square + dim];
		//intermediate offspring, all bred Individuals plus the best ranked parental Individuals a placed in this 
		//array of Individuals
		
		for(int i = 0; i < dim - 1; i++){
			//looping over all best ranked Individuals in the field 'pop'
			
			for(int j = i + 1; j < dim; j++){
				int index1 = getRank(i), index2 = getRank(j);
				off[counter] = new Individuals(getPop(index1).breedIndis(getPop(index2)));
				//breeding the offspring by calling the method 'breedIndis' from the class Individuals
				
				counter++;
			}
		}
		if(hasNull(off)){
			//check, whether the intermediate Individuals array still contains null entries
			
			int i = 0;
			while(counter < square + dim && i < dim){
				off[counter] = getPop(getRank(i));
				//filling the null entries with the best ranked parental Individuals
				
				counter++; i++;
			}
		}
		System.out.println("# of bred Individuals: " + counter);
		Individuals[] off2 = new Individuals[length];
		//the actual offspring Individuals array
		
		Demes p = new Demes(off);
		//creating an intermediate instance of the class 'Demes'
		
		p.bestFifty(CMDM);
		//ranking them
		
		System.out.println("# of best ranked Individuals: " + p.getRank().length);
		int counter1 = 0;
		while(counter1 < length){
			// a loop that only takes the best ranked Individuals in the field 'pop' of 'p'
			//as a random selection
			
			int rank = p.getRank().length;
			Random rand = new Random();
			int index = rand.nextInt(rank);
			//a random number between zero and the maximal number of best ranked Individuals in the
			//field 'pop' of 'p'
			
			off2[counter1] = new Individuals(p.getPop(p.getRank(index)));
			counter1++;
			
		}
		return off2;
	}
	
	/*--------------------------------------statics--------------------------------------------*/
	
	/**
	 * method, to copy an array of Individuals to a specified array of Individuals
	 * @throws PdbLoadError 
	 * @throws PdbCodeNotFoundError 
	 * @throws SQLException 
	 */
	public static void copyInd(Individuals[] in, Individuals[] out) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = in.length;
		if(out != null){
			if(in.length != out.length){
				for(int i = 0; i < dim; i++){
					out[i] = new Individuals();
				}
			}
		}
		else{
			out = new Individuals[dim];
		}
		for(int i = 0; i < dim; i++){
			if(in[i] != null){
				out[i] = new Individuals(in[i]);
			}
		}
		
	}
	
	/**
	 * method to trim an array of Individuals. This method checks all entries in the array an copies only
	 * non-null entries to new array of Individuals, with <tt>n - i</tt> entries, where <tt>i</tt> is the number of
	 * removed entries.
	 * @param input
	 * @return
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public static Individuals[] trim (Individuals[] input) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		int dim = input.length, counter = 0;
		boolean[] tester = new boolean[dim];
		for(int i = 0; i < dim; i++){
			if(input[i] == null){
				tester[i] = true;
				counter++;
			}
		}
		Individuals[] output = new Individuals[dim - counter];
		counter = 0;
		for(int i = 0; i < dim; i++){
			if(!tester[i]){
			output[counter] = new Individuals(input[i]);
			counter++;
			}
		}
		return output;
	}
	
	/**
	 * method, to compare all Individuals in this array with each other, using either CMError or DMError, specified
	 * by the boolean parameter CMDM, if CMDM = true, CMError is used, otherwise DMError
	 * @param input array of Individuals
	 * @param value
	 * @param CMDM
	 * @return test
	 */
	public static boolean[] compareInd (Individuals[] input, double value, boolean CMDM){
		int dim = input.length;
		boolean[] test = new boolean[dim];
		for(int i = 0; i < dim; i++){
			if(CMDM){
				if(input[i].getCM() <= value) {
					test[i] = true;
				}
			}
			else{
				if(input[i].getDM() <= value){
					test[i] = true;
				}
			}
		}
		return test;
	}
	
	/**
	 * counter method, counts how many instances of an Individuals array have been compared and tested positiv
	 * @param pop
	 * @param CMDM
	 * @param value
	 * @return
	 */
	public static int compareInd(Individuals[] pop, boolean CMDM, double value){
		int dim = pop.length;
		int counter = 0;
		boolean[] test = compareInd(pop, value, CMDM);
		for(int i = 0; i < dim; i++){
			if(test[i]){
				counter++;
			}
		}
		return counter;
	}

	/**
	 * auxiliary: computes the average of the array <tt>array</tt>
	 * @param array
	 * @return
	 */
	public static double entrySum (double[] array){
		int dim = array.length;
		double average = 0.0;
		for(int i = 0; i < dim; i++){
			average = average + array[i];
		}
		return average;
	}		

	/**
	 * method to determine average DM error
	 * @param pop
	 * @return avverageDMError
	 */
	public static double getAverageDMError (Individuals[] pop) {
		int dim = pop.length;
		double averageDMError = 0.0;
		for(int i = 0; i < dim; i++){
			if(pop[i] != null){
				averageDMError = averageDMError + pop[i].getDM();
			}
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
	 * method to evaluate the DMError standard deviation
	 * @param pop
	 * @return stdev
	 */
	public static double getDMStDeviation (Individuals[] pop){
		double average = getAverageDMError(pop);
		double var = 0.0, stdev = 0.0;
		int dim = pop.length;
		for(int i = 0; i < dim; i++){
			var = var + Math.pow(average - pop[i].getDM(), 2);
		}
		stdev = Math.pow(var/((double) dim*(dim - 1)), 0.5);
		return stdev;
	}
	
	/**
	 * method to evaluate the DMError standard deviation
	 * @param pop
	 * @return
	 */
	public static double getCMStDeviation (Individuals[] pop){
		double average = getAverageCMError(pop);
		double var = 0.0, stdev = 0.0;
		int dim = pop.length;
		for(int i = 0; i < dim; i++){
			var = var + Math.pow(average - pop[i].getCM(), 2);
		}
		stdev = Math.pow(var/((double) dim*(dim - 1)), 0.5);
		return stdev;
	}
	
	/**
	 * sanity checker: checks whether the parameter <tt>pop</tt> has 
	 * any null entries
	 * @param pop
	 * @return true, if no null entries were found
	 */
	public static boolean hasNull (Individuals[] pop){
		int length = pop.length;
		boolean nulltester = false;
		for(int i = 0; i < length; i++){
			if(pop[i] == null){
				nulltester = true;
				break;
			}
		}
		return nulltester;
	}
	
	/**
	 * computes <tt>i (i - 1)/2</tt>
	 * @param i
	 * @return <tt>i (i - 1)/2</tt>
	 */
	public static int nSquaredMinusNHalf(int i){
		return (int) ((double) (i*(i - 1))*0.5);
	}
	
	/**
	 * intended to locally optimize a given contact map
	 * TODO might be better in a new class and has to be improved
	 * @param spec
	 * @return
	 */
	public static Demes localOptimizer (Demes spec){
		Individuals[] array = spec.getPop();
		HashMap<Integer,HashSet<Pair<Integer>>> container = new HashMap<Integer,HashSet<Pair<Integer>>> (2*array.length);
		for(int i = 0; i < array.length; i++){
			container.put(new Integer (i), array[i].getHashSet());
		}
		return new Demes (array);
	}
	
	/**
	 * intended to locally optimize a given contact map
	 * TODO might be better in a new class and has to be improved
	 * @param spec
	 * @return
	 */
	public static HashSet<Pair<Integer>> alterHashSet (HashSet<Pair<Integer>> set, int number_of_alteration, int max_entry){
		HashSet<Pair<Integer>> newset = new HashSet<Pair<Integer>> ();
		HashSet<Pair<Integer>> survei = new HashSet<Pair<Integer>> ();
		Iterator<Pair<Integer>> it = set.iterator();
		int counter = 0, size = set.size();
		while(it.hasNext() && counter < size - number_of_alteration){
			Random rand = new Random ();
			int rand2 = rand.nextInt(2);
			if(rand2 == 0){
				Pair<Integer> pair = it.next();
				newset.add(pair);
				counter++;
			}
			else{
				Pair<Integer> pair = it.next();
				survei.add(pair);
			}
		}
		Iterator<Pair<Integer>> it2 = survei.iterator();
		while(it2.hasNext()){
			Pair<Integer> pair = it2.next();
			int f_val = pair.getFirst().intValue(), s_val = pair.getSecond().intValue();
			Random rand = new Random();
			int index1 = rand.nextInt(3);
			int index2 = rand.nextInt(3);
			boolean test1 = f_val + index1 < s_val + index2, test2 = s_val + index2 < max_entry, test3 = f_val - index1 > 0;
			if(test1 && test2){
				Integer ind1 = new Integer (f_val + index1), ind2 = new Integer (s_val + index2);
				Pair<Integer> pair1 = new Pair<Integer> (ind1,ind2);
				newset.add(pair1);
			}
			else{
				if(test3 && test2){
					Integer ind1 = new Integer (f_val - index1), ind2 = new Integer (s_val + index2);
					Pair<Integer> pair1 = new Pair<Integer> (ind1,ind2);
					newset.add(pair1);
				}
			}
		}
		return newset;
	}
	
	public static String getDefaultDatabase (){
		return new String (pdbaseDb);
	}

	public static void main (String[] args) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		/*Demes pop = new Demes(14,0,5,0,"run");
		String path = "/home/gmueller/workspace/aglappe/embed/teststructures/";
		pop.printToFile(path, pop.getName());
		//Individuals*/
		/*String dir = "/project/StruPPi/ga_dan/";
		File dr = new File(dir);
		File[] list =  dr.listFiles(new RegexFileFilter(".*.cm"));
		int length = list.length;
		for(int i = 0; i < length; i++){
			Individuals in = new Individuals(list[i].getAbsolutePath());
			in.printToFile(dir, "dan", in.getName());
			Individuals.clearFullCMandDM();
		}*/
		String pdb = "2jo7";
		Demes dm = new Demes (pdb,20,10.0);
		System.out.println(dm.toString());
	}

	/**
	 * a helper class for the outer class. this class computes the distances of two Individuals
	 * by summing over the minimal absolute distance of the contact pairs.
	 * @author gmueller
	 *
	 */
	private class Metric {

		private static final String path = "/home/gmueller/Arbeiten/ContactMaps/tests/main_run/";

		private HashMap<Pair<Integer>,Integer> metricmap;

		/**
		 * the only constructor of this inner class. Takes an instance of the Demes class and
		 * computes the metric to each entry in the outer class field <code>{@link Demes#pop}</code>.
		 * This class contains only one field <code>{@link #metricmap}</code>, which is a HashMap.
		 * The keys are Pairs of Individuals (indices pointing to the entries in <tt>pop</tt>) and the value
		 * is an Integer, which is the sum over all minimal contact pair differences. 
		 * mapping 
		 * @param pop
		 */
		public Metric (Demes pop) {
			setMetricMap(pop);
		}

		/**
		 * initial setter: takes all common contacts, and computes the minimal difference between this
		 * contact pair and all other contact pairs present in <tt>pop</tt>.
		 * <p>
		 * the metric <tt>d</tt> is defined as <tt>d(x,y) = sum (|i_x - i_y| + |j_x - j_y|)</tt> over all contact pairs, where <tt>x</tt> and
		 * <tt>y</tt> are to distinct Individuals in the field <tt>pop</tt> and <tt>i</tt>, <tt>y</tt> are the first end second index of a
		 * given contact pair of a specific Individuals.
		 * </p> 
		 * @param pop
		 */
		public void setMetricMap(Demes pop) {
			Individuals[] array = pop.getPop();
			int length = array.length, counter1 = 0;//, counter2 = 0;
			metricmap = new HashMap<Pair<Integer>,Integer> ();
			HashMap<Integer,Integer> helpmap = new HashMap<Integer,Integer> ();
			for(int i = 0; i < length - 1; i++){
				for(int j = i + 1; j < length; j++){
					HashSet<Pair<Integer>> set1 = array[i].getHashSet();
					HashSet<Pair<Integer>> set2 = array[j].getHashSet();
					HashSet<Pair<Integer>> commons = new HashSet<Pair<Integer>>(set1);
					commons.retainAll(set2);
					set1.removeAll(commons);
					set2.removeAll(commons);
					Iterator<Pair<Integer>> it1 = set1.iterator();
					while(it1.hasNext()){
						Pair<Integer> pair1 = it1.next();
						Iterator<Pair<Integer>> it2 = set2.iterator();
						while(it2.hasNext()){
							Pair<Integer> pair2 = it2.next();
							int value11 = pair1.getFirst().intValue();
							int value21 = pair1.getSecond().intValue();
							int value12 = pair2.getFirst().intValue();
							int value22 = pair2.getSecond().intValue();
							int difference = Math.abs(value11 - value12) + Math.abs(value21 - value22);
							Integer counter_val = new Integer (counter1);
							if(helpmap.containsKey(counter_val)){
								int map_val = helpmap.get(counter_val).intValue();
								if(difference < helpmap.get(counter_val)){
									helpmap.put(counter_val, new Integer(difference));
								}
								else{
									if(map_val == 0){
										break;
									}
								}
							}
							else{
								helpmap.put(counter_val, new Integer (difference));
							}
							//set2.remove(pair2);
						}
						//counter2++;
						counter1++;
					}
					counter1 = 0; //counter2 = 0;
					Pair<Integer> pair = new Pair<Integer>(new Integer(i), new Integer(j));
					metricmap.put(pair, new Integer(calcSum(helpmap)));
				}
			}
		}

		/**
		 * returns the metric map
		 * @return a HashMap of index pairs and the metric (as value)
		 */
		public HashMap<Pair<Integer>,Integer> getMetMap(){
			return new HashMap<Pair<Integer>,Integer> (metricmap);
		}

		/**
		 * computes the sum of all values of an Integer HashMap.
		 * @param map
		 * @return the sum of all values in the <tt>map</tt>
		 */
		public int calcSum (HashMap<Integer, Integer> map){
			HashSet<Integer> keyset = new HashSet<Integer> (map.keySet());
			Iterator<Integer> it = keyset.iterator();
			int sum = 0;
			while(it.hasNext()){
				Integer keypair = it.next();
				sum = sum + map.get(keypair).intValue();
			}
			return sum;
		}

		/**
		 * returns the field <tt>metricmap</tt> as a String
		 */
		public String toString(){
			HashSet<Pair<Integer>> keyset = new HashSet<Pair<Integer>> (metricmap.keySet());
			Iterator<Pair<Integer>> it = keyset.iterator();
			String output = "";
			while(it.hasNext()){
				Pair<Integer> pair = it.next();
				output = output + pair.getFirst().intValue() + "\t" + pair.getSecond() + "\t" + metricmap.get(pair) + "\n";
			}
			return output;
		}

		/**
		 * method, writing the metric tables to a file 
		 * @throws IOException
		 */
		public void printToFile() throws IOException{
			FileOutputStream output = new FileOutputStream(path + "metrix.met");
			PrintStream printer = new PrintStream(output);
			printer.print(toString());
			output.close();
			printer.close();
		}

	}
}

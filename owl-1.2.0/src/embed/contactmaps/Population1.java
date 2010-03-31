package embed.contactmaps;

import java.sql.SQLException;
import java.util.Random;


import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;

/**
 * This class deals with arrays of 'Individuals' instances. It compares two array entries, determines the entry with smaller error
 * and can generate offspring by aligning and retaining the contacts present in both 'Individuals' instances and randomly selects 
 * the missing contacts (the number of contacts is kept constant over evolution) 
 * @author gmueller
 *
 */
public class Population1 {
	
	/*-----------------------------------fields-----------------------------------------------*/
	
	/**
	 * field: an array of 'Individuals' instances
	 */
	public Individuals[] pop;
	
	/**
	 * field: an array of boolean, used to determine the 50% best of this 'Population' according to the 
	 * selected method (DMError or CMError)
	 */
	public boolean[] fiftyfifty;
	
	/**
	 * statistical fields: defined by the given error estimation method
	 */
	public double CMError, CMstdev, DMError, DMstdev;
	
	/**
	 * field: name - PDB codes used
	 */
	public String name;
	
	/**
	 * field: size - number of 'Individuals' instance in this 'Population'
	 */
	public int size;
	
	/*-------------------------Constructors-----------------------------------------------------*/
	
	/**
	 * zero parameter constructer, sets all fields to their default values
	 */
	public Population1() {
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
	 * one parameter constructor, creates a population of 'i' individuals using the zero parameter constructor
	 * of 'Individuals'
	 * @param i
	 */
	public Population1(int i) {
		setSize(i);
		pop = new Individuals[i];
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
	public Population1 (Individuals[] p){
		setPop(p);
		setSize(p.length);
		setCMErrorStats(p);
		setDMErrorStats(p);
		setPopName(p[0].getName());
		fiftyfifty = new boolean[p.length];
	}
	
	/**]
	 * one parameter constructor, reads an instance of 'Population' - similar to a copy method
	 * @param pop
	 */
	public Population1 (Population1 pop){
		int dim = pop.getSize();
		setSize(dim);
		this.pop = new Individuals[dim];
		setPop(pop.getPop());
		setCMErrorStats(pop.getPop());
		setDMErrorStats(pop.getPop());
		fiftyfifty = new boolean[dim];
		setPopName(pop.getName());
		}
	
	/*-------------------------------------Setters------------------------------------------*/
	
	/**
	 * one parameter setter, reads an array of 'Individuals'
	 */
	public void setPop(Individuals[] pop){
		int dim = pop.length;
		pop = new Individuals[dim];
		for(int i = 0; i < dim; i++){
			pop[i] = new Individuals(pop[i]);
		}
	}
	
	/**
	 * setter, sets this Population instance's name, usually the pdb code
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
	 * setter, sets this Population size, i.e. length of the Individuals array
	 */
	public void setSize(int i){
		size = i;
	}
	
	/**
	 * Setter, sets this field 'fiftyfifty' by pairwise comparing this array of 'Individuals'
	 * @param pop
	 * @param CMDM
	 */
	public void bestFifty (boolean CMDM){
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
				if(error[i] < error[j]){
					better[i] = better[i] + 1;
				}
				if(error[i] > error[j]){
					better[j] = better[j] + 1;
				}
				}
			}
		for(int i = 0; i < dim; i++){
			if(better[i] >= (int) ((double) dim - 1)/2.0){
				fiftyfifty[i] = true;
			}
		}
	}
	
	/*--------------------------------------Getters-------------------------------------------*/
	
	/**
	 * getter, returns this Individuals array
	 * @return p - array of Individuals
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
	public Individuals getPop (int i){
		Individuals p = new Individuals(pop[i]);
		return p;
	}
	
	/**
	 * getter returns this name
	 * @return
	 */
	public String getName() {
		return name;
	}
	
	/**
	 * getter, returns this average CMError
	 * @return
	 */
	public double getAvCMError () {
		return CMError;
	}
	
	/**
	 * getter, returns this average DMError
	 * @return
	 */
	public double getAvDMError () {
		return DMError;
	}
	
	/**
	 * getter, returns this standard deviation of CMError
	 * @return
	 */
	public double getCMstdev () {
		return CMstdev;
	}
	
	/**
	 * getter, returns this standard deviation of DMError
	 * @return
	 */
	public double getDMstdev () {
		return DMstdev;
	}
	
	/**
	 * getter, returns this Population size
	 * @return
	 */
	public int getSize(){
		return size;
	}
	
	/**
	 * getter, returns this 'fiftyfifty' field at the i-th entry
	 * @param i
	 * @return
	 */
	public boolean getFifty(int i){
		return fiftyfifty[i];
	}
	
	/**
	 * getter, returns this field 'n
	 * @return
	 */
	public int getNumOfContacts(){
		return pop[0].getNumOfContacts();
	}
	
	/*--------------------------------------statics--------------------------------------------*/
	
	/**
	 * method, to copy an array of Individuals to a specified array of Individuals
	 */
	public static void copyInd(Individuals[] in, Individuals[] out){
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
	 * method, to get the best fifty percent of a Population evaluated by either CMError or DMError
	 * @param CMDM  true - CMError taken
	 * @param pop
	 */
	public static void bestFifty (boolean CMDM, Population1 pop){
		int dim = pop.getSize();
		double[] error = new double[dim];
		double errorval = 0.0;
		double stmd = 0.0;
		if(CMDM){
			System.arraycopy(Individuals.getCM(pop.getPop()),0,error,0,dim);
			errorval = getAverageCMError(pop.getPop());
			stmd = getCMStDeviation(pop.getPop());
		}
		else{
			System.arraycopy(Individuals.getDM(pop.getPop()),0,error,0,dim);
			errorval = getAverageDMError(pop.getPop());
			stmd = getDMStDeviation(pop.getPop());
		}
		for(int i = 0; i < dim - 1; i++){
				
				if(error[i] <= errorval + 0.5*stmd){
					pop.fiftyfifty[i] = true;
				}
				else{
					pop.fiftyfifty[i] = false;
				}
		}
	}
	
	/**
	 * method, to get the best fifty percent of an array of Individuals evaluated by either CMError or DMError
	 * @param input
	 * @param CMDM
	 * @return
	 */
	public static boolean[] bestFifty (Individuals[] input, boolean CMDM){
		int dim = input.length;
		boolean[] result = new boolean[dim];
		double[] errorarray = new double[dim];
		double error = 0.0, stdev = 0.0;
		if(CMDM){
			error = getAverageCMError(input);
			stdev = getCMStDeviation(input);
			System.arraycopy(Individuals.getCM(input), 0, errorarray, 0, dim);
		}
		else{
			error = getAverageDMError(input);
			stdev = getDMStDeviation(input);
			System.arraycopy(Individuals.getDM(input), 0, errorarray, 0, dim);
		}
		boolean[] compare = compareInd(input, error, CMDM);
		for(int i = 0; i < dim; i++){
			if(!compare[i]){
				double diff = Math.abs(errorarray[i] - error);
				if(diff < stdev){
					result[i] = true;
				}
			}					
		}
		return result;
	}
	
	/**
	 * evolution method, creates offspring by comparing the hashset, retaining common contacts and
	 * fills up the remaining number of contacts randomly by picking contacts from either of the
	 * parents
	 * @param parents
	 * @return offspring
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public static Individuals[] evolve (Individuals[] parents) throws SQLException, PdbCodeNotFoundError, PdbLoadError {
		int dim = parents.length;
		//double CMav = getAverageCMError(parents);
		double DMav = getAverageDMError(parents);
		Individuals[] offspring = new Individuals[dim];
		int counter = 0;
		boolean[] tester1 = new boolean[dim];
		//boolean[] tester2 = new boolean[dim];
		tester1 = compareInd(parents, DMav, false);
		//tester2 = bestFifty(parents, false);
		for(int i = 0; i < dim - 1; i++){
			for(int j = i + 1; j < dim; j++){
				//if(parents[i].getCM() <= CMav && parents[j].getCM() <= CMav){
				if(parents[i] != null && parents[j] != null){
					if(tester1[i] && tester1[j] && counter < dim){
						offspring[counter] = new Individuals(parents[i].breedIndis(parents[j]));
						counter++;
					}
				}
			}
		}
		if(counter < dim - 1 && offspring[counter + 1] == null){
			int i = 0, j = 0;
			while(counter + j < dim && i < dim){
				if(tester1[i]){
					offspring[counter+j] = new Individuals(parents[i]);
					j++;
				}
				i++;
			}
		}
		return offspring;
	}
	
	/**
	 * evolution method, creates offspring by comparing the hashset, retaining common contacts and
	 * fills up the remaining number of contacts randomly by picking contacts from either of the
	 * parents
	 * @param CMDM
	 * @param pop
	 * @return off
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public static Individuals[] evolve (boolean CMDM, Population1 pop) throws SQLException, PdbCodeNotFoundError, PdbLoadError {
		int dim = pop.getSize(), counter = 0;
		pop.bestFifty(CMDM);
		Individuals[] off = pop.getPop();
		for(int i = 0; i < dim - 1; i++){
			for(int j = i + 1; j < dim; j++){
			if(pop.getFifty(i) && pop.getFifty(j) && counter < dim){
				off[counter] = new Individuals(pop.getPop(i).breedIndis(pop.getPop(j)));
				counter++;
			}
		}
		}
		if(off[dim - 1] == null) {//test, if the last entry is not null
			Random rand = new Random();
			int j = 0;
			int[] array = new int[dim];
			while(counter + j < dim){
				int x = rand.nextInt(dim - 1);
				if(off[x].getFifty()){
					for(int i = 0; i <= j; i++){
						if(x == array[i] && x != 0){
							break;
						}
						else{
							off[counter + j] = new Individuals(pop.getPop()[x]);
							array[j] = x;
							j++;
							break;
						}
					}
				}
			}
		}
		return off;
	}
	
	/**
	 * evolution method, creates offspring by comparing the hashset, retaining common contacts and
	 * fills up the remaining number of contacts randomly by picking contacts from either of the
	 * parents
	 * @param pop
	 * @param CMDM
	 * @return
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public static Individuals[] evolve (Individuals[] pop, boolean CMDM) throws SQLException, PdbCodeNotFoundError, PdbLoadError {
		int dim = pop.length;
		int counter = 0;
		double average = 0.0;
		boolean[] tester = new boolean[dim];
		Individuals[] offspring = new Individuals[dim];
		if(CMDM){
			average = getAverageCMError(pop);
			tester = compareInd(pop, average, CMDM);
		}
		else{
			average = getAverageDMError(pop);
			tester = compareInd(pop, average, CMDM);
		}
		while(counter < dim){
			for(int k = 0; k < dim - 1; k++){
				for(int l = k + 1; l < dim; l++){
					if(tester[k] && tester[l] && counter < dim){
						offspring[counter] = new Individuals(pop[k].breedIndis(pop[l]));
						counter++;
					}
				}
			}
		}
		return offspring;
			
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

}

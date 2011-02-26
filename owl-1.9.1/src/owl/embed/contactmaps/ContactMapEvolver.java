package owl.embed.contactmaps;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadError;

//import embed.SparseMatrix;


public class ContactMapEvolver {
	
	private ContactMap[] pop;			//the population of contact maps
	
	private String pdb;					//pdb identifier
	
	private boolean indexer2_used;		//distinguishes which indexing method was used
	
	private int size;					//the population size
	
	private int curr_gen;				//current generation
	
	private int evo_steps;				//number of maximal evolution steps
	
	private HashSet<Integer> indexer2;	//a collection of integers specifying the best ranked contact maps,
										//according to setIndexer2() method
	
	
	private HashSet<Integer> indexer3;	//a collection of integers specifying the best ranked contact maps,
										//according to setIndexer2() method
	
	private double avg_natives2;		//the average of native contacts of the contact map squared
	
	private double non_natives2;		//the average of non-native contacts of the contact map squared
	
	private double avg_natives3;		//the average of native contacts of the contact map 3
	
	private double non_natives3;		//the average of non-native contacts of the contact map 3
	
	private double av_sensit2;			//sensitivity test: true positives/(true positives + false positives)
										//in contact map squared
	
	private double av_sensit3;			//sensitivity test: true positives/(true positives + false positives)
										//in contact map powered 3
	
	public ContactMapEvolver (int numb_of_in, String pdb_code, int steps) throws SQLException, PdbCodeNotFoundException, PdbLoadError{
		setInitialPop(numb_of_in, pdb_code, steps);
	}
	
	public ContactMapEvolver (ContactMap[] ar, int steps, int current){
		setPredecessors(ar);
		setEvoSteps(steps);
		curr_gen = current;
	}
	/**
	 * One parameter constructor: initializes a shallow copy of this instance.
	 * @param cmv
	 */
	public ContactMapEvolver (ContactMapEvolver cmv){
		pop = new ContactMap[cmv.size];
		size = cmv.size();
		for(int i = 0; i < size; i++){
			pop[i] = cmv.getPop(i);
		}
		setPdb(cmv.getPdb());
		curr_gen = cmv.getCurrentGeneration();
		setEvoSteps(cmv.evo_steps);
		avg_natives2 = cmv.getAverageNatives();
		avg_natives3 = cmv.getAverageNatives3();
		non_natives2 = cmv.getNonNatives();
		non_natives3 = cmv.getNonNatives3();
		Random rand = new Random();
		int val = rand.nextInt(2);
		if(val == 1) setIndexer2();
		else setIndexer3();
		setAverageSensitivity();
	}
	
	public void setInitialPop (int numb_of_in, String pdb_code, int steps) throws SQLException, PdbCodeNotFoundException, PdbLoadError{
		pop  = new ContactMap [numb_of_in];
		size = numb_of_in;
		setPdb(pdb_code);
		double av21 = 0.0, av22 = 0.0, av31 = 0.0, av32 = 0.0;
		for(int i = 0; i < numb_of_in; i++){
			pop[i] = new ContactMap(pdb_code,4,0.05);
			av21 += (double) pop[i].getNativeCounter();
			av22 += (double) pop[i].getNonNatives();
			av31 += (double) pop[i].getNativeCounter3();
			av32 += (double) pop[i].getNonNatives3();
		}
		avg_natives2 = av21/((double) size);
		non_natives2 = av22/((double) size);
		avg_natives3 = av31/((double) size);
		non_natives3 = av32/((double) size);
		setEvoSteps(steps);
		Random rand = new Random();
		int val = rand.nextInt(2);
		if(val == 1) setIndexer2();
		else setIndexer3();
		setAverageSensitivity();
	}
	
	public void setPredecessors (ContactMap[] array) throws IllegalArgumentException {
		int length = array.length;
		if(length > 1){
			pop = new ContactMap [length];
			double av21 = 0.0, av22 = 0.0, av31 = 0.0, av32 = 0.0;
			for(int i = 0;  i < length; i++){
				if(i < length - 1){
					if(array[i].getPdbCode().matches(array[i+1].getPdbCode())){
						pop[i] = new ContactMap (array[i]);
						av21 += (double) pop[i].getNativeCounter();
						av22 += (double) pop[i].getNonNatives();
						av31 += (double) pop[i].getNativeCounter3();
						av32 += (double) pop[i].getNonNatives3();
					}
					else throw new IllegalArgumentException ("Never store two contact maps with different pdb codes in one population!"); 
				}
				else {if(array[i].getPdbCode().matches(array[i-1].getPdbCode())) pop[i] = new ContactMap (array[i]);
				else throw new IllegalArgumentException ("Never store two contact maps with different pdb codes in one population!");
				}
			}
			size = length;
			avg_natives2 = av21/((double) size);
			non_natives2 = av22/((double) size);
			avg_natives3 = av31/((double) size);
			non_natives3 = av32/((double) size);
			setPdb(array[0].getPdbCode());
			Random rand = new Random();
			int val = rand.nextInt(2);
			if(val == 1) setIndexer2();
			else setIndexer3();
			setAverageSensitivity();
		}
		else{
			System.err.println("Loss of contact maps...");
			System.exit(1);
		}
	}
	
	public void setPdb (String pdb_code){
		pdb = new String (pdb_code);
	}
	
	public void setEvoSteps (int steps){
		evo_steps = steps;
	}
	
	public void setIndexer2 (){
		if(pop != null){
			int length = size();
			double thresh = ((double) length)/4.0 - 1.0;
			double[] array = new double[length];
			for(int i = 0; i < length - 1; i++){
				for(int j = i + 1; j < length; j++){
					//double sum_nat = pop[i].getNativeCounter() + pop[j].getNativeCounter();
					double sum_non = pop[i].getNonNatives() + pop[j].getNonNatives() + pop[i].getNativeCounter() + pop[j].getNativeCounter();
					double ratio_nat_i = pop[i].getNativeCounter()/sum_non;
					double ratio_nat_j = pop[j].getNativeCounter()/sum_non;
					double ratio_non_i = pop[i].getNonNatives()/sum_non;
					double ratio_non_j = pop[j].getNonNatives()/sum_non;
					if(ratio_nat_i >= ratio_nat_j && ratio_non_i <= ratio_non_j) array[i]++;
					if(ratio_nat_j >= ratio_nat_i && ratio_non_j <= ratio_non_i) array[j]++;
				}
			}
			indexer2 = new HashSet<Integer> ();
			for(int i = 0; i < length; i++){
				if(array[i] >= thresh){
					indexer2.add(new Integer (i));
				}
			}
			indexer2_used = true;
		}
	}
	
	public void setIndexer3 (){
		if(pop != null){
			int length = size();
			double thresh = ((double) length)/4.0 - 1.0;
			double[] array = new double[length];
			for(int i = 0; i < length - 1; i++){
				for(int j = i + 1; j < length; j++){
					//double sum_nat = pop[i].getNativeCounter3() + pop[j].getNativeCounter3();
					double sum_non = pop[i].getNonNatives3() + pop[j].getNonNatives3() + pop[i].getNativeCounter3() + pop[j].getNativeCounter3();
					double ratio_nat_i = pop[i].getNativeCounter3()/sum_non;
					double ratio_nat_j = pop[j].getNativeCounter3()/sum_non;
					double ratio_non_i = pop[i].getNonNatives3()/sum_non;
					double ratio_non_j = pop[j].getNonNatives3()/sum_non;
					if(ratio_nat_i >= ratio_nat_j && ratio_non_i <= ratio_non_j) array[i]++;
					if(ratio_nat_j >= ratio_nat_i && ratio_non_j <= ratio_non_i) array[j]++;
					//if(pop[i].getNativeCounter3() >= pop[j].getNativeCounter3() && pop[i].getNonNatives() <= pop[j].getNonNatives()) array[i]++;
					//if(pop[j].getNativeCounter3() >= pop[i].getNativeCounter3() && pop[j].getNonNatives() <= pop[i].getNonNatives()) array[j]++;
				}
			}
			indexer3 = new HashSet<Integer> ();
			for(int i = 0; i < length; i++){
				if(array[i] >= thresh){
					indexer3.add(new Integer (i));
				}
			}
			indexer2_used = false;
		}
	}
	
	public void emptyAll (){
		pop = null;
		indexer2 = null;
		indexer3 = null;
	}
	
	public void clearAllFields (){
		emptyAll();
		pdb = null;
	}
	
	public void evolve () {
		if((indexer2 != null || indexer3 != null) && pop != null){
			HashSet<Integer> cp1 = getIndexer();
			HashSet<Integer> cp2 = getIndexer();
			HashMap<Integer,ContactMap> maps = new HashMap<Integer,ContactMap> ();
			Iterator<Integer> it1 = cp1.iterator();
			int counter = 0;
			while(it1.hasNext()){
				Integer index1 = it1.next();
				cp2.remove(index1);
				Iterator<Integer> it2 = cp2.iterator();
				while(it2.hasNext()){
					Integer index2 = it2.next();
					ContactMap merged = getPop(index1.intValue()).merge(getPop(index2.intValue()));
					maps.put(new Integer (counter), merged);
					counter++;
				}
			}
			it1 = cp1.iterator();
			while(it1.hasNext()){
				Integer index = it1.next();
				maps.put(new Integer (counter), getPop(index.intValue()));
				counter++;
			}
			int length = maps.size(), counter1 = 0;
			ContactMap[] ar = new ContactMap [length];
			Set<Integer> keys = maps.keySet();
			Iterator<Integer> it = keys.iterator();
			while(it.hasNext()){
				Integer index = it.next();
				ar[index.intValue()] = new ContactMap (maps.get(index));
			}
			ContactMapEvolver ev = new ContactMapEvolver (ar,evoSteps(),getCurrentGeneration());
			HashSet<Integer> indexer_ = ev.getIndexer();
			int length1 = indexer_.size(), cur_size = size(), final_size = 0;
			if(length1 <= cur_size) final_size = length1;
			else final_size = cur_size;
			ContactMap[] off = new ContactMap[final_size];
			it = indexer_.iterator();
			while(it.hasNext() && counter1 < final_size){
				Integer index = it.next();
				off[counter1] = new ContactMap(ev.getPop(index.intValue()));
				counter1++;
			}
			curr_gen++;
			emptyAll();
			ev.clearAllFields();
			setPredecessors(off);
		}
		else throw new NullPointerException ("Population and indexer fields must be initialized before calling this method!");		
	}
	
	public void writeToFile (String dir_name, String file_name) throws SQLException, PdbCodeNotFoundException, PdbLoadError, FileNotFoundException, IOException{
		File file = new File(dir_name);
		if(!file.exists()) file.mkdirs();
		Individuals[] indar = new Individuals[size];
		for(int i = 0; i < size; i++){
			indar[i] = getPop(i).convertToIndividuals();
		}
		Demes deme = new Demes (indar,getCurrentGeneration());
		deme.printToFile(dir_name, file_name);
	}
	
	public void setAverageSensitivity(){
		if(pop != null){
			double av2 = 0.0, av3 = 0.0;
			for(int i = 0; i < size; i++){
				av2 += pop[i].getSensitivity2();
				av3 += pop[i].getSensitivity3();
			}
			av_sensit2 = av2/((double) size);
			av_sensit3 = av3/((double) size);
		}
	}
	
	public int size (){
		return size;
	}
	
	public int evoSteps (){
		return evo_steps;
	}
	
	public double getAverageNatives (){
		return avg_natives2;
	}
	
	public double getNonNatives (){
		return non_natives2;
	}
	
	public double getAverageNatives3 (){
		return avg_natives3;
	}
	
	public double getNonNatives3 (){
		return non_natives3;
	}
	
	public double getAverageSensitivity2(){
		return av_sensit2; 
	}
	
	public double getAverageSensitivity3(){
		return av_sensit3;
	}
	
	public int getCurrentGeneration (){
		return curr_gen;
	}
	
	public HashSet<Integer> getIndexer (){
		if(indexer2_used) return new HashSet<Integer> (indexer2);
		else return new HashSet<Integer> (indexer3);
	}
	
	public String getPdb (){
		return new String (pdb);
	}
	
	public ContactMap[] getPop (){
		int length = size();
		ContactMap[] newcm = new ContactMap [length];
		for(int i = 0; i < length; i++){
			newcm[i] = new ContactMap (pop[i]);
		}
		return newcm;
	}
	
	public ContactMap getPop (int i){
		if(i < size) return new ContactMap(pop[i]);
		else throw new ArrayIndexOutOfBoundsException ("The index exceeded the length of the population array!");
	}
	
	/*public String toString(){
		
		String str = "#CMVIEW GRAPH FILE ver: 1.0\n#SEQUENCE: "+cm1..getSequence()+"\n"+
		"#PDB: "+rig.getPdbCode()+ "\n#PDB CHAIN CODE: "+rig.getChainCode()+"\n#CT: "+rig.getContactType()+ "\n#CUTOFF: "+rig.getCutoff() + "\n";
		SparseMatrix mat = map.get(new Integer(i));
	}*/
	
	public static void main (String[] args) throws SQLException, PdbCodeNotFoundException, PdbLoadError, FileNotFoundException, IOException{
		ContactMapEvolver ev  = new ContactMapEvolver (100,"1bkr",20);
		ContactMapEvolver off = new ContactMapEvolver (ev);
		//String dir_name = "/project/StruPPi/gabriel/Arbeiten/cmEvolver/evo_run/";
		//String file_name = "CMEvolved";
		while(off.getCurrentGeneration() < off.evoSteps()){
			if(off.indexer2_used) System.out.println("generation: "+off.getCurrentGeneration()+"\tindexer: 2\tsize: " + off.size() + "\tnatives2: " + off.getAverageNatives() + "\t" +
					"natives3: "+off.getAverageNatives3()+"\tnon natives: " + off.getNonNatives()+"\tnon natives3: "+off.getNonNatives3()+"\tsensitivity2: "+off.getAverageSensitivity2()+"\tsensitivity3: "+off.getAverageSensitivity3());
			else System.out.println("generation: "+off.getCurrentGeneration()+"\tindexer: 3\tsize: " + off.size() + "\tnatives2: " + off.getAverageNatives() + "\t" +
					"natives3: "+off.getAverageNatives3()+"\tnon natives: " + off.getNonNatives()+"\tnon natives3: "+off.getNonNatives3()+"\tsensitivity2: "+off.getAverageSensitivity2()+"\tsensitivity3: "+off.getAverageSensitivity3());
			//off.writeToFile(dir_name, file_name+off.getCurrentGeneration());
			off.evolve();
		}
	}

}

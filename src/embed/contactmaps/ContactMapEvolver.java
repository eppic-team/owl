package embed.contactmaps;

import java.sql.SQLException;
import java.util.*;

import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;

public class ContactMapEvolver {
	
	private ContactMap[] pop;
	
	private String pdb;
	
	private int size;
	
	private int curr_gen;
	
	private int evo_steps;
	
	private HashSet<Integer> indexer;
	
	private double avg_natives;
	
	private double non_natives;
	
	public ContactMapEvolver (int numb_of_in, String pdb_code, int steps) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		setInitialPop(numb_of_in, pdb_code, steps);
	}
	
	public ContactMapEvolver (ContactMap[] ar, int steps, int current){
		setPredecessors(ar);
		setEvoSteps(steps);
		curr_gen = current;
	}
	
	public void setInitialPop (int numb_of_in, String pdb_code, int steps) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		pop  = new ContactMap [numb_of_in];
		size = numb_of_in;
		setPdb(pdb_code);
		double av1 = 0.0, av2 = 0.0;
		for(int i = 0; i < numb_of_in; i++){
			pop[i] = new ContactMap(pdb_code);
			av1 += (double) pop[i].getNativeCounter();
			av2 += (double) pop[i].getNonNatives();
		}
		avg_natives = av1/((double) size);
		non_natives = av2/((double) size);
		setEvoSteps(steps);
		setIndexer();
	}
	
	public void setPredecessors (ContactMap[] array) throws IllegalArgumentException {
		int length = array.length;
		if(length > 1){
			pop = new ContactMap [length];
			double av1 = 0.0, av2 = 0.0;
			for(int i = 0;  i < length; i++){
				if(i < length - 1){
					if(array[i].getPdbCode().matches(array[i+1].getPdbCode())){
						pop[i] = new ContactMap (array[i]);
						av1 += (double) pop[i].getNativeCounter();
						av2 += (double) pop[i].getNonNatives();
					}
					else throw new IllegalArgumentException ("Never store two contact maps with different pdb codes in one population!"); 
				}
				else {if(array[i].getPdbCode().matches(array[i-1].getPdbCode())) pop[i] = new ContactMap (array[i]);
				else throw new IllegalArgumentException ("Never store two contact maps with different pdb codes in one population!");
				}
			}
			size = length;
			avg_natives = av1/((double) size);
			non_natives = av2/((double) size);
			setPdb(array[0].getPdbCode());
			setIndexer();
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
	
	public void setIndexer (){
		if(pop != null){
			int length = size();
			double thresh = ((double) length)/2.0 - 1.0;
			double[] array = new double[length];
			for(int i = 0; i < length - 1; i++){
				for(int j = i + 1; j < length; j++){
					if(pop[i].getNativeCounter() >= pop[j].getNativeCounter()) array[i]++;// && pop[i].getNonNatives() <= pop[j].getNonNatives()) array[i]++;
					if(pop[j].getNativeCounter() >= pop[i].getNativeCounter()) array[j]++;// && pop[j].getNonNatives() <= pop[i].getNonNatives()) array[j]++;
				}
			}
			indexer = new HashSet<Integer> ();
			for(int i = 0; i < length; i++){
				if(array[i] >= thresh){
					indexer.add(new Integer (i));
				}
			}
		}
	}
	
	public int size (){
		return size;
	}
	
	public int evoSteps (){
		return evo_steps;
	}
	
	public double getAverageNatives (){
		return avg_natives;
	}
	
	public double getNonNatives (){
		return non_natives;
	}
	
	public int getCurrentGeneration (){
		return curr_gen;
	}
	
	public HashSet<Integer> getIndexer (){
		return new HashSet<Integer> (indexer);
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
	
	public ContactMap[] evolve (){
		if(indexer != null && pop != null){
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
			return off;
		}
		else throw new NullPointerException ("Population and indexer fields must be initialized before calling this method!");
		
	}
	
	public static void main (String[] args) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		ContactMapEvolver ev  = new ContactMapEvolver (50,"1bkr",20);
		ContactMapEvolver off = new ContactMapEvolver (ev.getPop(),ev.evoSteps(),ev.getCurrentGeneration());
		while(off.getCurrentGeneration() < off.evoSteps()){
			System.out.println("size: " + off.size() + "\tnatives: " + off.getAverageNatives() + "\tnon natives: " + off.getNonNatives());
			off = new ContactMapEvolver (off.evolve(),off.evoSteps(),off.getCurrentGeneration());
		}
	}

}

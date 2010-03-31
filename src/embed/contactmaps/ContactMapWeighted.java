package embed.contactmaps;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.*;

import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;

import edu.uci.ics.jung.graph.util.Pair;
import embed.SparseMatrix;

public class ContactMapWeighted extends ContactMap {
	
	private SparseMatrix weighted_cm;
	
	private HashSet<Pair<Integer>> negs;
	
	private HashSet<Pair<Integer>> poss;
	
	public ContactMapWeighted (String pdb_code) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		super();
		super.setContactMap(pdb_code);
		setWeightedContactMap();
	}
	
	public void setWeightedContactMap (){
		if(super.getMap() != null && super.getSquareMap() != null){
			HashMap<Pair<Integer>,Double> sq_map = super.getSquareMap().getMatrix();
			SparseMatrix map    = super.getMap();
			Set<Pair<Integer>> keys = sq_map.keySet();
			Iterator<Pair<Integer>> it = keys.iterator();
			negs = new HashSet<Pair<Integer>> ();
			poss = new HashSet<Pair<Integer>> ();
			HashMap<Pair<Integer>,Double> weighted = new HashMap<Pair<Integer>,Double>(); 
			while(it.hasNext()){
				Pair<Integer> pair = it.next();
				Integer f_ind = pair.getFirst(), s_ind = pair.getSecond();
				HashSet<Pair<Integer>> first  = map.getPairSetFromFirstIndex(f_ind.intValue());
				HashSet<Pair<Integer>> second = map.getPairSetFromSecondIndex(s_ind.intValue());
				Iterator<Pair<Integer>> it1 = first.iterator();
				while(it1.hasNext()){
					Pair<Integer> testP1 = it1.next();
					Iterator<Pair<Integer>> it2 = second.iterator();
					while(it2.hasNext()){
						Pair<Integer> testP2 = it2.next();
						if(testP1.getSecond().intValue() == testP2.getFirst().intValue()){
							if(super.getFullMap().containsIndexPair(pair)){
								weighted.put(testP1, new Double (1.0));
								poss.add(testP1);
							}
							else{
								weighted.put(testP1, new Double (-1.0));
								negs.add(testP1);
							}
							if(super.getMap().containsIndexPair(pair)){
								weighted.put(testP2, new Double(1.0));
								poss.add(testP2);
							}
							else{
								weighted.put(testP2, new Double(-1.0));
								negs.add(testP2);
							}
						}
					}
				}
			}
			weighted_cm = new SparseMatrix (weighted,map.getColumnDimension(),map.getColumnDimension());
		}
	}
	
	public SparseMatrix getWeightedContactMap (){
		if(weighted_cm != null) return new SparseMatrix (weighted_cm);
		else throw new NullPointerException ("The weighted contact map was not initialized before calling this method!");
	}
	
	public HashSet<Pair<Integer>> getIndexPairsPresent(){
		return new HashSet<Pair<Integer>> (poss);
	}
	
	public HashSet<Pair<Integer>> getIndexPairsNonPresent(){
		return new HashSet<Pair<Integer>> (negs);
	}
	
	public static void main (String[] args) throws SQLException, PdbCodeNotFoundError, PdbLoadError, FileNotFoundException, IOException{
		String dir      = "/project/StruPPi/gabriel/Arbeiten/";
		String addname1 = "cmEvolver/"; 
		ContactMapWeighted cmh = new ContactMapWeighted ("1bkr");
		Individuals in = cmh.convertToIndividuals();
		System.out.println(cmh.getNativeCounter() +"\n" + in.toString());
		in.printToFile(dir, addname1, "test02"+in.getName());
	}

}

package embed.contactmaps;

import java.io.*;
import java.sql.SQLException;
import java.util.*;

import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import edu.uci.ics.jung.graph.util.*;

public class Metric2 {
	
	private static final String path = "/home/gmueller/Arbeiten/ContactMaps/tests/main_run/";
	
	public HashMap<Pair<Integer>,Integer> metricmap;
	
	public Metric2 (Demes pop) {//throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		setMetricMap(pop);
	}

	public void setMetricMap(Demes pop) {//throws SQLException, PdbCodeNotFoundError, PdbLoadError{
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
	
	public HashMap<Pair<Integer>,Integer> getMetMap(){
		return new HashMap<Pair<Integer>,Integer> (metricmap);
	}
	
	public static int calcSum (HashMap<Integer, Integer> map){
		HashSet<Integer> keyset = new HashSet<Integer> (map.keySet());
		Iterator<Integer> it = keyset.iterator();
		int sum = 0;
		while(it.hasNext()){
			Integer keypair = it.next();
			sum = sum + map.get(keypair).intValue();
		}
		return sum;
	}
	
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
	
	public void printToFile() throws IOException{
		FileOutputStream output = new FileOutputStream(path + "metrix.met");
		PrintStream printer = new PrintStream(output);
		printer.print(toString());
		output.close();
		printer.close();
	}
	
	public static void main (String[] args) throws SQLException, IOException, PdbCodeNotFoundError, PdbLoadError{
		Demes p = new Demes(5,0,5,0,"dummy");
		Metric2 met = new Metric2(p);
		met.printToFile();
	}

}

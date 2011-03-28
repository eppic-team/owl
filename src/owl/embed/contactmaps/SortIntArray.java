package owl.embed.contactmaps;

import java.util.*;
import edu.uci.ics.jung.graph.util.*;

/**
 * A helper class that converts unsorted <code> HashSet<Pair<Integer>> </code> instances into sorted <code> int[][] </code> arrays.
 * This class is intended to help converting the <code> HashSet<Pair<Integer>> </code> field <code> store </code> of the <code> Individuals </code>
 * class to a lexicographical sorted <code> int </code> array, providing only static methods.
 * @author gmueller
 *
 */
public class SortIntArray {
	
	/**
	 * This method converts a HashSet<Pair<Integer>> instance into an one dimensional int array, taking only the first value of each
	 * Pair<Integer> instance. The array is sorted before return.
	 * @param set
	 * @return
	 */
	private static int[] convert_HashSet (Set<Pair<Integer>> set){
		int size = set.size(), counter = 0;
		int[] array = new int[size];
		//the resulting array is initialized
		
		Iterator<Pair<Integer>> it = set.iterator();
		//an iterator over the HashSet is initialized
		
		while(it.hasNext()){
			//looping over the HashSet
			
			array[counter] = it.next().getFirst().intValue();
			//first value taken
			
			counter++;
		}
		Arrays.sort(array);
		//sorting of the array
		
		return array;
	}
	
	/**
	 * 
	 * @param set
	 * @return
	 */
	private static HashMap<Integer,HashSet<Integer>> convertHashSet (Set<Pair<Integer>> set){
		HashMap<Integer,HashSet<Integer>> result = new HashMap<Integer,HashSet<Integer>> ();
		//later output
		
		HashSet<Pair<Integer>> dublicate = new HashSet<Pair<Integer>> (set);
		//a duplicate of the parameter 'set'
		
		HashSet<Integer> firstindex = new HashSet<Integer>();
		//a HashSet containing the first value of each Pair present in the parameter 'set'
		
		Iterator<Pair<Integer>> it = dublicate.iterator();
		while(it.hasNext()){
			//loop over all entries of the parameter 'set'
			
			Pair<Integer> pair = it.next();
			//pair
			
			Integer first = pair.getFirst();
			Integer secon = pair.getSecond();
			if(!firstindex.contains(first)){
				//check, whether 
				firstindex.add(first);
				HashSet<Integer> seconindex = new HashSet<Integer>();
				seconindex.add(secon);
				result.put(first, seconindex);
			}
			else{
				HashSet<Integer> seconindex = new HashSet<Integer> (result.get(first));
				seconindex.add(secon);
				result.put(first, seconindex);
			}
		}
		return result;
	}
	
	private static int[] convertHashSet (HashSet<Integer> indexset){
		int size = indexset.size(), counter = 0;
		int[] array = new int[size];
		Iterator<Integer> it = indexset.iterator();
		while(it.hasNext()){
			array[counter] = it.next().intValue();
			counter++;
		}
		Arrays.sort(array);
		return array;
	}
	
	private static int[][] convertHashMap (HashMap<Integer,HashSet<Integer>> map, int[] array){
		int sizea = array.length;
		int[][] result = new int [2][sizea];
		System.arraycopy(array, 0, result[0], 0, sizea);
		Arrays.sort(result[0]);
		HashSet<Integer> keyset = new HashSet<Integer> (map.keySet());
		Iterator<Integer> it = keyset.iterator();
		while(it.hasNext()){
			Integer index1 = it.next();
			for(int i = 0; i < sizea; i++){
				if(result[0][i] == index1.intValue()){
					int[] subarray = convertHashSet(map.get(index1));
					int length = subarray.length;
					System.arraycopy(subarray, 0, result[1], i, length);
					break;
				}
			}
		}
		return result;
	}
	
	/**
	 * this is the main method to use in order to convert an HashSet of Pairs of Integers to
	 * an sorted integer array. 
	 * @param set
	 * @return
	 */
	public static int[][] converter(Set<Pair<Integer>> set){
		HashMap<Integer,HashSet<Integer>> indexmap = convertHashSet(set);
		int[] array = convert_HashSet(set);
		int[][] newarray = convertHashMap(indexmap,array);
		return newarray;
	}
	
		
	public static void main(String[] args){
		HashSet<Pair<Integer>> set = new HashSet<Pair<Integer>> ();
		for(int i = 0; i < 100; i++){
			Random rand1 = new Random();
			Random rand2 = new Random();
			int index1 = rand1.nextInt(100);
			int index2 = rand2.nextInt(100 - index1) + index1;
			Pair<Integer> pair = new Pair<Integer> (new Integer (index1), new Integer (index2));
			set.add(pair);
		}
		int[][] newarray = converter(set);
		boolean tester = true;
		String method_tester = ""; 
		for(int i = 0; i < newarray[0].length; i++){
			Pair<Integer> pair = new Pair<Integer> (new Integer(newarray[0][i]), new Integer(newarray[1][i]));
			if(set.contains(pair)){
				tester = tester && true;
			}
			else{
				tester = tester && false;
				break;
			}
		}
		if(tester){
			method_tester = "Test was successful!";
		}
		else{
			method_tester = "Test failed!";
		}
		System.out.println(method_tester);
	}

}

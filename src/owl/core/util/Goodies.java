package owl.core.util;
import java.util.*;

/**
 * Collection of small useful helper methods.
 * @author stehr
 *
 */
public class Goodies {
	
	public static final boolean ASCENDING = true;
	public static final boolean DESCENDING = false;	
	
	/**
	 * Sorts a map by the values and returns a sorted map with the right ordering.
	 * @param map the initial map
	 * @param ascending if true, map will be sorted in ascending order of values, otherwise in descending order
	 * @return the sorted map
	 */
	public static <K,V extends Comparable<V>> LinkedHashMap<K,V> sortMapByValue(Map<K,V> map, boolean ascending) {
		List<Map.Entry<K,V>> list = new LinkedList<Map.Entry<K,V>>(map.entrySet());
		if(ascending) {
		Collections.sort(list, new Comparator<Map.Entry<K,V>>() {
			public int compare(Map.Entry<K,V> o1, Map.Entry<K,V> o2) {
				return (o1.getValue()).compareTo(o2.getValue());				
			}
		});
		} else {
			Collections.sort(list, new Comparator<Map.Entry<K,V>>() {
				public int compare(Map.Entry<K,V> o1, Map.Entry<K,V> o2) {
					return -(o1.getValue()).compareTo(o2.getValue());				
				}
			});			
		}

		LinkedHashMap<K,V> result = new LinkedHashMap<K,V>();
		for (Iterator<Map.Entry<K,V>> it = list.iterator(); it.hasNext();) {
			Map.Entry<K,V> entry = it.next();
			result.put(entry.getKey(), entry.getValue());
		}
		return result;
	}
	
		
}

package tools;
import java.util.*;

/**
 * Collection of small useful helper method.
 * @author stehr
 *
 */
public class Goodies {
	
	/**
	 * Sorts a map by the values and returns a sorted map with the right ordering.
	 * @param map the initial map
	 * @return the sorted map
	 */
	public static <K,V extends Comparable<V>> Map<K,V> sortByValue(Map<K,V> map) {
		List<Map.Entry<K,V>> list = new LinkedList<Map.Entry<K,V>>(map.entrySet());
		Collections.sort(list, new Comparator<Map.Entry<K,V>>() {
			public int compare(Map.Entry<K,V> o1, Map.Entry<K,V> o2) {
				return ((Comparable<V>) o1.getValue()).compareTo(o2.getValue());				
			}
		});

		Map<K,V> result = new LinkedHashMap<K,V>();
		for (Iterator<Map.Entry<K,V>> it = list.iterator(); it.hasNext();) {
			Map.Entry<K,V> entry = it.next();
			result.put(entry.getKey(), entry.getValue());
		}
		return result;
	}

}

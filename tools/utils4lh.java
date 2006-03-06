package tools;

import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Iterator;

public class utils4lh {

    private utils4lh() {}

    public static LinkedHashMap reverse(LinkedHashMap lhm) {

	LinkedHashMap rvLhm = new LinkedHashMap();
	
	Object[] a = lhm.keySet().toArray();
	
	for (int i = a.length-1; i > -1; i--) {
	    rvLhm.put(a[i], lhm.get(a[i]));
	    
	}
	
	return rvLhm;
    
    }

    public static LinkedHashMap putKeyFirst(LinkedHashMap lhm, Object key) {
	
        LinkedHashMap oldLhm = (LinkedHashMap)lhm.clone();
	LinkedHashMap newLhm = new LinkedHashMap();

	newLhm.put(key, oldLhm.get(key));
	oldLhm.remove(key);
	newLhm.putAll(oldLhm);

	return newLhm;

    }

    public static LinkedHashSet linkedHashSetDiff(LinkedHashSet a, LinkedHashSet b) {

	Object obj;
	LinkedHashSet a_b = new LinkedHashSet();
	LinkedHashSet ab = new LinkedHashSet();
	LinkedHashSet all = new LinkedHashSet();

	for (Iterator aItr = a.iterator(); aItr.hasNext();) {
	    obj = aItr.next();
	    if (b.contains(obj)) {
		b.remove(obj);
		ab.add(obj);
	    } else {
		a_b.add(obj);
	    }
	}

	all.add(a_b);
	all.add(ab);
	all.add(b);
	
	return all;

    }

}

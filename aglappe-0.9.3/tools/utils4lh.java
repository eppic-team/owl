package tools;

import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Iterator;

public class utils4lh {

    private utils4lh() {}

    public static LinkedHashMap<?,?> reverse(LinkedHashMap<?,?> lhm) {

	LinkedHashMap<Object,Object> rvLhm = new LinkedHashMap<Object, Object>();
	
	Object[] a = lhm.keySet().toArray();
	
	for (int i = a.length-1; i > -1; i--) {
	    rvLhm.put(a[i], lhm.get(a[i]));
	    
	}
	return rvLhm;
    
    }

    public static LinkedHashMap<?,?> putKeyFirst(LinkedHashMap<?,?> lhm, Object key) {
	
    LinkedHashMap<?,?> oldLhm = (LinkedHashMap<?, ?>) lhm.clone();
	LinkedHashMap<Object,Object> newLhm = new LinkedHashMap<Object,Object>();

	newLhm.put(key, oldLhm.get(key));
	oldLhm.remove(key);
	newLhm.putAll(oldLhm);

	return newLhm;

    }

    public static LinkedHashSet<?> linkedHashSetDiff(LinkedHashSet<?> a, LinkedHashSet<?> b) {

	Object obj;
	LinkedHashSet<Object> a_b = new LinkedHashSet<Object>();
	LinkedHashSet<Object> ab = new LinkedHashSet<Object>();
	LinkedHashSet<Object> all = new LinkedHashSet<Object>();

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

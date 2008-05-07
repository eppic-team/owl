package sequence;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * A convenience class for a list of GTGHits
 *
 */
public class GTGHitList implements Iterable<GTGHit> {

	private ArrayList<GTGHit> hits;
	
	public GTGHitList() {
		this.hits = new ArrayList<GTGHit>();
	}
	
	public void add(GTGHit hit) {
		hits.add(hit);
	}
	
	public int size() {
		return this.hits.size();
	}
	
	public void print() {
		for (GTGHit hit:this.hits) {
			hit.print();
		}
	}
	
	/**
	 * Returns and iterator over this GTG hits list
	 * @return
	 */
	public Iterator<GTGHit> iterator() {
		return this.hits.iterator();
	}
}

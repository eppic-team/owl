package sequence;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * A list of blast hits
 *
 */
public class BlastHitList {
	
	private ArrayList<BlastHit> hits;
	
	public BlastHitList() {
		this.hits = new ArrayList<BlastHit>();
	}
	
	/**
	 * Adds a blast hit to this list
	 * @param hit
	 */
	public void add(BlastHit hit) {
		this.hits.add(hit);
	}

	/**
	 * Applies an e-value cutoff trimming out of this list all hits with e-value 
	 * higher than cutoff given
	 * @param eValueCutoff
	 */
	public void applyCutoff(double eValueCutoff) {
		Iterator<BlastHit> it = hits.iterator();
		while (it.hasNext()) {
			if (it.next().eValue>=eValueCutoff) {
				it.remove();
			}
		}
	}
	
	/**
	 * Prints a few selected fields of all blast hits in this list
	 */
	public void print() {
		for (BlastHit hit:hits) {
			hit.print();
		}
	}
	
	/**
	 * Returns the number of blast hits contained in this list
	 * @return
	 */
	public int size() {
		return hits.size();
	}
}

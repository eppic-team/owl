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
	 * Prints a course ascii-art overview of the hits in this hist list.
	 * The length of the query sequence is scaled to the given output length in screen columns.
	 */
	public void printWithOverview(int queryLength, int outputLength) {
		double scaleFactor = 1.0 * outputLength / queryLength;
		BlastHit.printHeadersWithOverview(queryLength, scaleFactor);
		for(BlastHit hit:hits) {
			hit.printWithOverview(scaleFactor);
		}
	}
	
	/**
	 * Returns the number of blast hits contained in this list
	 * @return
	 */
	public int size() {
		return hits.size();
	}

	/**
	 * Return an array of hits contained in this hit list
	 * @return
	 */
	public BlastHit[] getHits() {
		return (BlastHit[]) this.hits.toArray();
	}
	
	/**
	 * Returns the hit with the best e-value or null if this hit list is empty.
	 * @return
	 */
	public BlastHit getBestHit() {
		if(this.size() == 0) return null;
		BlastHit bestHit = this.hits.get(0);
		for(BlastHit hit:hits) {
			if(hit.getEValue() < bestHit.getEValue()) {
				bestHit = hit;
			}
		}
		return bestHit;
	}
	
}

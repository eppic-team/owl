package sequence;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;

/**
 * A list of blast hits
 *
 */
public class BlastHitList {
	
	private ArrayList<BlastHit> hits;
	private int queryLength;				// needed by print() and printSome()
	
	public BlastHitList() {
		this.hits = new ArrayList<BlastHit>();
		this.queryLength = 0;
	}
	
	/**
	 * Adds a blast hit to this list
	 * @param hit
	 */
	public void add(BlastHit hit) {
		this.hits.add(hit);
	}

	/**
	 * Set the query length needed by the print() and printSome() methods.
	 * TODO: This should be parsed from the output file. However, the tabular output we use at the moment
	 * does not contain this information. So this function provides a way to set it externally. Moving
	 * to parsing the XML output would resolve this issue.
	 * @param l
	 */
	public void setQueryLength(int l) {
		this.queryLength = l;
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
	 * Prints a tabular overview of the hits in this list.
	 * Currently, if the query length is set, (i.e. > 0) a simple ascii-art
	 * overview of the matches is printed for each hit.
	 */
	public void print() {
		printSome(this.size());
	}
	
	/**
	 * Prints a tabular overview of the first numHits hits in this list.
	 * Currently, if the query length is set, (i.e. > 0) a simple ascii-art
	 * overview of the matches is printed for each hit.
	 */
	public void printSome(int numHits) {
		
		int outputLength = 80;		// length of graphical output in screen columns
		
		if(queryLength > 0) {
			double scaleFactor = 1.0 * outputLength / queryLength;
			BlastHit.printHeadersWithOverview(queryLength, scaleFactor);
			for (int i = 0; i < Math.min(hits.size(), numHits); i++) {
				BlastHit hit = hits.get(i);
				hit.printWithOverview(scaleFactor);
			}			
		} else {
			// print without graphical overview
			for (int i = 0; i < Math.min(hits.size(), numHits); i++) {
				BlastHit hit = hits.get(i);
				hit.print();
			}			
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
	
	/**
	 * Returns a list of template ids (pdb codes+chain codes)
	 * @return
	 */
	public String[] getTemplateIds() {
		String[] ids = new String[this.size()];
		for (int i=0;i<this.size();i++) {
			ids[i]=hits.get(i).getTemplateId();
		}
		return ids;
	}
	
	/**
	 * 
	 * @param outFile
	 * @throws IOException
	 */
	public void writeTemplateIdsToFile(File outFile) throws IOException {
		PrintWriter Out = new PrintWriter(outFile);
		String[] ids = getTemplateIds();
		for (String id:ids) {
			Out.println(id);
		}
		Out.close();
	}
}

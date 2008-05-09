package sequence;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;

import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.TemplateList;

import tools.MySQLConnection;

/**
 * A convenience class for a list of GTGHits
 *
 */
public class GTGHitList implements Iterable<GTGHit> {

	private ArrayList<GTGHit> hits;
	
	private int queryLength;
	
	public GTGHitList() {
		this.hits = new ArrayList<GTGHit>();
	}
	
	public void add(GTGHit hit) {
		hits.add(hit);
	}
	
	public int size() {
		return this.hits.size();
	}
	
	public void printTabular() {
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
	
	/**
	 * Reassigns the serials of subjectStart and subjectEnd positions based on the 
	 * pdbase full SEQRES sequence, i.e. converts GTG serials to our internal residue
	 * serials (cif serials)
	 * @param conn
	 * @param pdbaseDb
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 * @see GTGHit.reassignSubjectSerials
	 */
	public void reassignSubjectSerials(MySQLConnection conn, String pdbaseDb) throws SQLException, PdbCodeNotFoundError, PdbLoadError {
		for (GTGHit hit: this) {
			hit.reassignSubjectSerials(conn, pdbaseDb);
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
	 * @param numHits
	 */
	public void printSome(int numHits) {
		
		int outputLength = 80;		// length of graphical output in screen columns
		
		if(queryLength > 0) {
			double scaleFactor = 1.0 * outputLength / queryLength;
			GTGHit.printHeaderWithOverview(queryLength, scaleFactor);
			for (int i = 0; i < Math.min(hits.size(), numHits); i++) {
				GTGHit hit = hits.get(i);
				hit.printWithOverview(scaleFactor);
			}			
		} else {
			// print without graphical overview
			for (int i = 0; i < Math.min(hits.size(), numHits); i++) {
				GTGHit hit = hits.get(i);
				hit.print();
			}			
		}
	}
	
	/**
	 * Set the query length needed by the print() and printSome() methods.
	 * @param l
	 */
	public void setQueryLength(int l) {
		this.queryLength = l;
	}
	
	/**
	 * Applies an total score cutoff trimming out of this list all hits with totalScore 
	 * lower than given cutoff
	 * @param totalScoreCutoff
	 */
	public void applyCutoff(int totalScoreCutoff) {
		Iterator<GTGHit> it = this.iterator();
		while (it.hasNext()) {
			if (it.next().totalScore<=totalScoreCutoff) {
				it.remove();
			}
		}
	}
	
	/**
	 * Filters out from this GTGHitList PDB hits released after the given date
	 * (in format yyyymmdd, e.g. 20060425)  
	 * @param conn
	 * @param pdbaseDb
	 * @param date date in format yyyymmdd
	 * @throws SQLException
	 */
	public void filterByMaxReleaseDate(MySQLConnection conn, String pdbaseDb, String date) throws SQLException {
		String[] filteredIds = TemplateList.filterIdsByMaxReleaseDate(conn, pdbaseDb, this.getTemplateIds(), date);
		for (String filteredId:filteredIds) {
			this.hits.remove(this.getHit(filteredId));
		}
	}
	
	/**
	 * Filters out hits that rank below given rank.
	 * @param rank
	 */
	public void filterByMaxRank(int rank) {
		Iterator<GTGHit> it = this.iterator();
		int i=1;
		while (it.hasNext()) {
			it.next();
			if (i>rank) it.remove();
			i++;
		}
	}

	/**
	 * Returns the hit with the best totalScore or null if this hit list is empty.
	 * @return
	 */
	public GTGHit getBestHit() {
		if(this.size() == 0) return null;
		GTGHit bestHit = this.hits.get(0);
		for(GTGHit hit:hits) {
			if(hit.totalScore > bestHit.totalScore) {
				bestHit = hit;
			}
		}
		return bestHit;
	}
	
	/**
	 * Gets a GTGHit given its subjectId or null if GTGHit not present in this list 
	 * @param templateId
	 * @return
	 */
	public GTGHit getHit(String templateId) {
		for (GTGHit hit:this) {
			if (hit.getTemplateId().equals(templateId)) return hit;
		}
		return null;
	}
	
	/**
	 * Gets the GTGHit corresponding to given rank
	 * @param rank
	 * @return the hit at given rank
	 * @throws IndexOutOfBoundsException if given rank out of bounds
	 */
	public GTGHit getByRank(int rank) {
		if (rank>this.hits.size() || rank<1) {
			throw new IndexOutOfBoundsException("Given rank "+rank+" is not whithin the ranks range: 1 to "+this.hits.size());
		}
		return this.hits.get(rank-1);
	}
	
	/**
	 * Returns a list of template ids (subject ids) as concatenated 
	 * pdbCodes+pdbChainCodes e.g. 1abcA 
	 * @return
	 */
	public String[] getTemplateIds() {
		String[] ids = new String[this.size()];
		for (int i=0;i<this.size();i++) {
			ids[i] = this.hits.get(i).getTemplateId();
		}
		return ids;
	}
	
}

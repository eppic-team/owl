package owl.core.sequence;

import java.util.ArrayList;
import java.util.List;

import owl.core.connections.UnirefXMLParser;

/**
 * A minimal Uniref entry representation as parsed from a UniRef xml file 
 * See uniprot archives at ftp://ftp.uniprot.org/pub/databases/uniprot/previous_releases/
 * 
 * @see UnirefXMLParser
 * @author duarte_j
 *
 */
public class UnirefEntry {

	private String id;
	private String uniprotId;
	private String uniparcId;
	private int ncbiTaxId;
	private String sequence;
	
	private List<String> inactiveUniprotIds;
	private List<String> clusterMembers;
	
	private List<String> taxons;
	
	public UnirefEntry() {
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getUniprotId() {
		return uniprotId;
	}

	public void setUniprotId(String uniprotId) {
		this.uniprotId = uniprotId;
	}

	public String getUniparcId() {
		return uniparcId;
	}

	public void setUniparcId(String uniparcId) {
		this.uniparcId = uniparcId;
	}

	public int getNcbiTaxId() {
		return ncbiTaxId;
	}

	public void setNcbiTaxId(int ncbiTaxId) {
		this.ncbiTaxId = ncbiTaxId;
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	
	/**
	 * The taxonomy as returned by UniProt JAPI's UniprotEntry.getTaxonomy(). 
	 * First value is the most general (domain of life), last is the most specific (species).
	 * @return
	 */
	public List<String> getTaxons() {
		return taxons;
	}	
	
	public void setTaxons(List<String> taxons) {
		this.taxons = taxons;
	}
	
	/**
	 * Adds an inactive uniprot id to this uniref entry.
	 * Uniprot ids get obsoleted (inactive) when they are merged, demerged, deleted, etc
	 * The ids are then never reused and kept in an inactive ids database.
	 * e.g. Q31NP8 is an obsoleted entry replaced by P30340
	 * See http://www.uniprot.org/help/query-fields and 
	 * http://www.uniprot.org/uniprot/?query=active:no 
	 * @param uniprotId
	 */
	public void addInactiveUniprotId(String uniprotId) {
		if (inactiveUniprotIds==null) {
			inactiveUniprotIds = new ArrayList<String>();
		}
		this.inactiveUniprotIds.add(uniprotId);
	}
	
	/**
	 * Adds a cluster member UniProt id to this UniRef entry
	 * Uniref entries are clustering several sequences that share 100%, 90% or 50% identity,
	 * thus they are actually clusters for which a representative has been chosen.
	 * e.g. UniRef100_P30340 contains 4 members, the representative P30340 plus Q2EFX9, Q5N5G6, Q6PQB8
	 * We only keep UniProt cluster members and not UniParc cluster members
	 * @param uniprotId
	 */
	public void addClusterMember(String uniprotId) {
		if (clusterMembers==null) {
			clusterMembers = new ArrayList<String>();
		}
		this.clusterMembers.add(uniprotId);
	}
	
	public List<String> getClusterMembers() {
		return this.clusterMembers;
	}
	
	public boolean hasClusterMembers() {
		return clusterMembers!=null;
	}
	
	public boolean hasTaxons() {
		return taxons!=null;
	}
	
	/**
	 * Returns the domain of life (what used to be kingdom) for this entry 
	 * @return
	 */
	public String getFirstTaxon() {
		return this.taxons.get(0);
	}
	
	/**
	 * Returns the most specific taxonomy annotation (species) for this entry
	 * @return
	 */
	public String getLastTaxon() {
		return this.taxons.get(this.taxons.size()-1);
	}
	
	/**
	 * Returns a UniprotEntry with uniprot id, sequence, taxId and taxons as they are in this UnirefEntry
	 * Does not set any information of EMBL CDS or related
	 * @return
	 */
	public UniprotEntry getUniprotEntry() {
		UniprotEntry uniprotEntry = new UniprotEntry(uniprotId);
		uniprotEntry.setUniprotSeq(new Sequence(isUniprot()?uniprotId:uniparcId,sequence));
		uniprotEntry.setTaxId(ncbiTaxId);
		uniprotEntry.setTaxons(taxons);
		return uniprotEntry;
	}
	
	/**
	 * Returns true if this UniRef entry is a UniProt entry, else (i.e. if a UniParc entry) it returns false
	 * @return
	 */
	public boolean isUniprot() {
		return uniprotId!=null;
	}
	
	/**
	 * If the entry contains null except for case uniprotId==null with a uniparcId!=null
	 * or case where uniprotId==null then ncbiTaxId can be ==0
	 * @return
	 */
	public boolean hasNulls() {
		return (id==null || (uniprotId==null && uniparcId==null) || uniparcId==null || (uniprotId!=null && ncbiTaxId==0) || sequence==null);
	}
	
	public String toString() {
		return id+"\t"+uniprotId+"\t"+uniparcId+"\t"+ncbiTaxId+"\t"+sequence;
	}
	
	public String getTabDelimitedMySQLString() {
		return id+"\t"+(uniprotId==null?"\\N":uniprotId)+"\t"+uniparcId+"\t"+ncbiTaxId+"\t"+sequence;
	}
}

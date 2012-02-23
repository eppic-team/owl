package owl.core.sequence;

import owl.core.connections.UnirefXMLParser;

/**
 * A minimal Uniref entry representation as parsed from a uniref xml file 
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

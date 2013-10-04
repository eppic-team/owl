package owl.core.structure.io;

/**
 * A secondary structure line in a PDB file, to be used while parsing PDB files 
 * to hold parsed text in memory before final secondary structure assignment. 
 * 
 * @author duarte_j
 *
 */
public class SecStructureLine {

	public String type; // HELIX or SHEET
	public String begChain; // the beg PDB chain code
	public String endChain; // the end PDB chain code
	public String begPdbChainCode;
	public String endPdbChainCode;
	public int serial;
	public String id; // blank for HELIX, the sheet id for SHEET
	
	public SecStructureLine(String type, String begChain, String endChain, String begPdbChainCode, String endPdbChainCode, int serial, String id) {
		this.type = type;
		this.begChain = begChain;
		this.endChain = endChain;
		this.begPdbChainCode = begPdbChainCode;
		this.endPdbChainCode = endPdbChainCode;
		this.serial = serial;
		this.id = id;
	}
}

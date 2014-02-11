package owl.core.structure.io;

public class PdbxPolySeqLine {
	
	public String asym_id;  
	public int seq_id; 
	public String mon_id; 
	public int pdb_seq_num; 
	public String pdb_strand_id; 
	public String pdb_ins_code; 

	public PdbxPolySeqLine(String asym_id, int seq_id, String mon_id, int pdb_seq_num, String pdb_strand_id, String pdb_ins_code) {
		this.asym_id = asym_id;
		this.seq_id = seq_id;
		this.mon_id = mon_id;
		this.pdb_seq_num = pdb_seq_num;
		this.pdb_strand_id = pdb_strand_id;
		this.pdb_ins_code = pdb_ins_code;
	}
}

package owl.core.structure.io;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeMap;

import owl.core.structure.AminoAcid;
import owl.core.structure.Nucleotide;
import owl.core.structure.PdbLoadException;

public class PdbxPolySeqGroup implements Iterable<PdbxPolySeqLine> {
	
	private String asymId;
	private String pdbChainCode;
	private TreeMap<String,Integer> pdbresser2resser;
	private TreeMap<Integer,String> resser2pdbresser;
	private String sequence;
	
	private ArrayList<PdbxPolySeqLine> list;
	
	private boolean protein;
	private boolean nucleotide;

	public PdbxPolySeqGroup(){
		list = new ArrayList<PdbxPolySeqLine>();
		pdbresser2resser = new TreeMap<String, Integer>();
		resser2pdbresser = new TreeMap<Integer, String>();
		sequence = "";
		protein = false;
		nucleotide = false;
	}
	
	public void add(PdbxPolySeqLine line) {
		// In some rare entries (e.g. 1k55C, 2ci1A, 1h9hA, 1ejgA) there are alt locs 
		// for a residue for which location A is a of a residue type and location B 
		// of another (usually a modified version of residue in A). The CIF file correctly 
		// assigns the same residue number to both. We check that with following line and take
		// the approach of considering only the first residue. As there's no alt loc annotation
		// in the pdbx_poly_seq_scheme field we can't do better than that. Anyway that will
		// probably match our approach of taking the first alphabetical alt loc only.
		if (resser2pdbresser.containsKey(line.seq_id)) return;
		
		list.add(line);
		if (asymId==null) {
			asymId = line.asym_id;
		}
		if (pdbChainCode==null) {
			pdbChainCode = line.pdb_strand_id;
		}
		if (AminoAcid.isStandardAA(line.mon_id)){
			protein = true;
			sequence+=AminoAcid.three2one(line.mon_id);
		} else if (Nucleotide.isStandardNuc(line.mon_id)) {
			nucleotide = true;
			sequence+=Nucleotide.getByCode(line.mon_id).getOneLetterCode();
		} else {
			sequence+=AminoAcid.XXX.getOneLetterCode();
		}
		
		pdbresser2resser.put(line.pdb_seq_num+(line.pdb_ins_code.equals(".")?"":line.pdb_ins_code),line.seq_id);
		resser2pdbresser.put(line.seq_id,line.pdb_seq_num+(line.pdb_ins_code.equals(".")?"":line.pdb_ins_code));
	}
	
	public String getSequence() {
		return sequence;
	}
	
	public String getPdbChainCode() {
		return pdbChainCode;
	}
	
	public String getChainCode() {
		return asymId;
	}
	
	public TreeMap<String,Integer> getPdbresser2resserMap() {
		return pdbresser2resser;
	}
	
	public TreeMap<Integer,String> getResser2pdbresserMap() {
		return resser2pdbresser;
	}
	
	public boolean isProtein() throws PdbLoadException {
		if (protein && nucleotide) {
			throw new PdbLoadException("Mix of protein and nucleotide sequences in chain");
		}
		if (!protein && !nucleotide) {
			// in case where not a single std aa or nucleotide are found then we have a chain of Xs, we assume a protein
			protein = true;
		}
		return protein;
	}
	
	public boolean isEmpty() {
		return list.isEmpty();
	}
	
	@Override
	public Iterator<PdbxPolySeqLine> iterator() {
		return list.iterator();
	}
}

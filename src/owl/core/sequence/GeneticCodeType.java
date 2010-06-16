package owl.core.sequence;

import java.util.HashMap;
import java.util.Map;

/**
 * Genetic code types as defined by NCBI in 
 * http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
 * 
 * @author duarte_j
 *
 */
public enum GeneticCodeType {
	
	STANDARD     (1, true,false,"Standard"),
	VERT_MIT     (2,false, true,"Vertebrate Mitochondrial"),
	YEAST_MIT    (3,false, true,"Yeast Mitochondrial"),
	MOLD_MIT     (4,false, true,"Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma"),
	INV_MIT      (5,false, true,"Invertebrate Mitochondrial"),
	CIL_NUC      (6,true ,false,"Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear"),
	ECHI_MIT     (7,false, true,"Echinoderm Mitochondrial; Flatworm Mitochondrial"),
	EUPLOTID_NUC (8,true, false,"Euplotid Nuclear");
	// there are still more exotic genetic codes defined by ncbi, didn't add them yet

	private static Map<Integer, GeneticCodeType> id2GCT = initId2GCT();

	
	private int id;
	private boolean mitochondrial;
	private boolean nuclear;
	private String name;
	
	private GeneticCodeType(int id, boolean nuclear, boolean mitochondrial, String name) {
		this.id = id;
		this.nuclear = nuclear;
		this.mitochondrial = mitochondrial;
		this.name = name;
	}
	
	private static Map<Integer,GeneticCodeType> initId2GCT() {
		Map<Integer,GeneticCodeType> map = new HashMap<Integer, GeneticCodeType>();
		for (GeneticCodeType gct:GeneticCodeType.values()) {
			map.put(gct.getId(), gct);
		}
		return map;
	}
	
	public boolean isMitochondrial() {
		return mitochondrial;
	}
	
	public boolean isNuclear() {
		return nuclear;
	}
	
	public int getId() {
		return id;
	}
	
	public String getName() {
		return name;
	}
	
	public static GeneticCodeType getById(int id) {
		return id2GCT.get(id);
	}
}

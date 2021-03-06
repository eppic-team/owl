package owl.core.sequence;

import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Genetic code types as defined by NCBI in 
 * http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
 * 
 * @author duarte_j
 *
 */
public enum GeneticCodeType {
	
	STANDARD     (1, 0, true, false,"Standard"),
	VERT_MIT     (2, 4, false, true,"Vertebrate Mitochondrial"),
	YEAST_MIT    (3,10, false, true,"Yeast Mitochondrial"),
	MOLD_MIT     (4, 9, false, true,"Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma"),
	INV_MIT      (5, 8, false, true,"Invertebrate Mitochondrial"),
	CIL_NUC      (6, 2, true ,false,"Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear"),
	// selecton has different ids for flatworm (7) and echinoderm mitochondrial (6), why?
	ECHI_MIT     (7, 6, false, true,"Echinoderm Mitochondrial; Flatworm Mitochondrial"), 
	EUPLOTID_NUC (8, 3, true, false,"Euplotid Nuclear");
	// there are still more exotic genetic codes defined by ncbi, didn't add them yet

	private static Map<Integer, GeneticCodeType> id2GCT = initId2GCT();
	
	private static Pattern MITOCHONDRIAL_REGEX = Pattern.compile(".*mitochond.*",Pattern.CASE_INSENSITIVE);
	
	private int id;
	private int selectonId;
	private boolean mitochondrial;
	private boolean nuclear;
	private String name;
	
	private GeneticCodeType(int id, int selectonId, boolean nuclear, boolean mitochondrial, String name) {
		this.id = id;
		this.selectonId = selectonId;
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
	
	public int getSelectonId() {
		return selectonId;
	}
	
	public String getName() {
		return name;
	}
	
	public static GeneticCodeType getById(int id) {
		return id2GCT.get(id);
	}
	
	public static GeneticCodeType getByOrganelleAndOrganism(String organelle, String organism) {
		if (organelle==null) {
			return GeneticCodeType.STANDARD;
		}
		boolean mitochondrial = false;
		Matcher mitochondMatcher = MITOCHONDRIAL_REGEX.matcher(organelle);
		if (mitochondMatcher.matches()) {
			mitochondrial = true;
		}
		// TODO check organism too and do something sensible when finding mitochondrial
		if (mitochondrial) {
			return null;
		}
		return GeneticCodeType.STANDARD;
	}
}

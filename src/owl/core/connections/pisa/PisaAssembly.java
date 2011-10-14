package owl.core.connections.pisa;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PisaAssembly implements Comparable<PisaAssembly> {

	private int id;
	private int mmsize;
	private String score;
	private double dissEnergy;
	private String formula;
	private String composition;
	private List<Integer> interfaceIds;
	
	public PisaAssembly() {
		interfaceIds = new ArrayList<Integer>();
	}

	public int getId() {
		return id;
	}

	public void setId(int id) {
		this.id = id;
	}

	public int getMmsize() {
		return mmsize;
	}

	public void setMmsize(int mmsize) {
		this.mmsize = mmsize;
	}

	public String getScore() {
		return score;
	}

	public void setScore(String score) {
		score = score.trim();
		score = score.replaceAll("\n", " ");
		score = score.replaceAll("\\s+", " ");
		this.score = score;
	}

	public double getDissEnergy() {
		return dissEnergy;
	}

	public void setDissEnergy(double dissEnergy) {
		this.dissEnergy = dissEnergy;
	}

	public String getFormula() {
		return formula;
	}

	public void setFormula(String formula) {
		this.formula = formula;
	}
	
	public String getComposition() {
		return composition;
	}
	
	public void setComposition(String composition) {
		composition = composition.trim();
		composition = composition.replaceAll("\n", "");
		this.composition = composition;
	}
	
	public List<Integer> getInterfaceIds() {
		return interfaceIds;
	}
	
	public String getInterfaceIdsString() {
		StringBuffer sb = new StringBuffer();
		sb.append("[");
		for (int i=0;i<interfaceIds.size();i++) {
			if (i!=0) sb.append(" ");
			sb.append(interfaceIds.get(i));
		}
		sb.append("]");
		return sb.toString();
	}
	
	public void addInterfaceId(int interfaceId) {
		interfaceIds.add(interfaceId);
	}

	/**
	 * Returns true if this assembly contains at least one macromolecule (either protein or nucleic acid), false otherwise.
	 * It is detected from the case of PISA chain identifiers in the assembly output. If chains are capital we consider 
	 * them to be macromolecules, if lower case we consider they are small molecules.
	 * @return
	 */
	public boolean isMacromolecular() {
		Pattern p = Pattern.compile("[A-Z]");
		Matcher m = p.matcher(formula);
		if (m.find()) {
			return true;
		} else {
			return false;
		}
	}
	
	@Override
	public int compareTo(PisaAssembly o) {
		return Double.compare(this.dissEnergy, o.dissEnergy);
	}
}

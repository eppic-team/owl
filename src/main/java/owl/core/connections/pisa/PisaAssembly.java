package owl.core.connections.pisa;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class PisaAssembly implements Comparable<PisaAssembly> {

	// the fantastically diplomatic PISA score strings ;)
	private static final String SCORE_STABLE = 
			"This assembly appears to be stable in solution.";
	private static final String SCORE_UNSTABLE = 
			"This assembly may be formed from crystallographic considerations, however it does not appear to be stable in solution.";
	// this message seems to be the one appearing on top of the list when there are no stable assemblies, however in some cases it appears
	// as the score of a specific assembly, e.g. 1eer (set 1, assembly 2)
	// we'll considered that it means simply unstable 
	private static final String SCORE_UNSTABLE2 =   
			"Analysis of crystal interfaces has not revealed any strong indications that this assembly may form stable complexes in solution.";
	private static final String SCORE_GRAY = 
			"This assembly falls into a grey region of complex formation criteria. It may or may not be stable in solution.";
	
	public enum PredictionType {
		STABLE, UNSTABLE, GRAY;
	}
	
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
	
	/**
	 * Returns the prediction type (STABLE, UNSTABLE, GRAY) from the score string.
	 * @return
	 */
	public PredictionType getPredictionType() {
		if (score.equals(SCORE_STABLE)) return PredictionType.STABLE;
		if (score.equals(SCORE_UNSTABLE) || score.equals(SCORE_UNSTABLE2)) return PredictionType.UNSTABLE;
		if (score.equals(SCORE_GRAY)) return PredictionType.GRAY;
		return null;
	}
	
	@Override
	public int compareTo(PisaAssembly o) {
		return Double.compare(this.dissEnergy, o.dissEnergy);
	}
	
	public String toString() {
		return "id:"+this.id+"-size:"+this.mmsize;
	}
}

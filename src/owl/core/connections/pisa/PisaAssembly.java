package owl.core.connections.pisa;

import java.util.ArrayList;
import java.util.List;

public class PisaAssembly {

	private int id;
	private int mmsize;
	private String score;
	private double dissEnergy;
	private String formula;
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
	
	public List<Integer> getInterfaceIds() {
		return interfaceIds;
	}
	
	public void addInterfaceId(int interfaceId) {
		interfaceIds.add(interfaceId);
	}
}

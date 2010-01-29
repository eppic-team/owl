package ppi;

import java.util.ArrayList;

public class PPINode {

	/*--------------------------- member variables --------------------------*/
	
	private int	idx;
	private String proteinName;
	private String description;
	private ArrayList<GOAnnotation> goAnnotations;
	
	/*----------------------------- constructors ----------------------------*/
	
	public PPINode(int idx) {
		this.idx = idx;
		this.proteinName = "";
		this.description = "";
		this.goAnnotations = new ArrayList<GOAnnotation>();
	}
	
	public PPINode(int idx, String proteinName, String description) {
		this.idx = idx;
		this.proteinName = proteinName;
		this.description = description;
		this.goAnnotations = new ArrayList<GOAnnotation>();		
	}

	/*---------------------------- public methods ---------------------------*/
	
	public String toString() {
		return String.format("Node:(idx=%d, name=%s, go=%s, desc=%s)", idx, proteinName, goAnnotations, description);
	}
	
	public void addGOAnnotation(GOAnnotation newAnno) {
		this.goAnnotations.add(newAnno);
	}
	
	/**
	 * @return the idx
	 */
	public int getIdx() {
		return idx;
	}

	/**
	 * @return the proteinName
	 */
	public String getProteinName() {
		return proteinName;
	}

	/**
	 * @return the description
	 */
	public String getDescription() {
		return description;
	}

	/**
	 * @return the goAnnotation
	 */
	public ArrayList<GOAnnotation> getGoAnnotations() {
		return goAnnotations;
	}

	/**
	 * @param proteinName the proteinName to set
	 */
	public void setProteinName(String proteinName) {
		this.proteinName = proteinName;
	}

	/**
	 * @param description the description to set
	 */
	public void setDescription(String description) {
		this.description = description;
	}
	

	
}

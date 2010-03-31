package ppi;

public class GOAnnotation {
	
	/*--------------------------- type definitions --------------------------*/
	public enum Domain {F,P,C};	// molecular function, biological process, cellular component
	
	/*--------------------------- member variables --------------------------*/
	int number;
	Domain domain;
	String description;
	
	/*----------------------------- constructors ----------------------------*/
	public GOAnnotation(int number, Domain domain, String description) {
		this.number = number;
		this.domain = domain;
		this.description = description;
	}

	/*---------------------------- public methods ---------------------------*/
	
	public String toString() {
		return String.format("%d", number);
	}
	
	/**
	 * @return the number
	 */
	public int getNumber() {
		return number;
	}

	/**
	 * @return the domain
	 */
	public Domain getDomain() {
		return domain;
	}

	/**
	 * @return the description
	 */
	public String getDescription() {
		return description;
	}
	
}

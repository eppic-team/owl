package owl.litNet;

public class Article {

	/*--------------------------- member variables --------------------------*/
	
	int idx;
	String title;
	String author;
	String journal;
	int year;
	
	/*----------------------------- constructors ----------------------------*/
	
	public Article(int idx, String author, int year, String journal, String title) {
		this.idx = idx;
		this.author = author;
		this.year = year;
		this.journal = journal;
		this.title = title;
	}

	/*---------------------------- public methods ---------------------------*/
	
	public String toString() {
		return author + " " + year;
	}
	
	/**
	 * @return the idx
	 */
	public int getIdx() {
		return idx;
	}

	/**
	 * @return the title
	 */
	public String getTitle() {
		return title;
	}

	/**
	 * @return the author
	 */
	public String getAuthor() {
		return author;
	}

	/**
	 * @return the journal
	 */
	public String getJournal() {
		return journal;
	}

	/**
	 * @return the year
	 */
	public int getYear() {
		return year;
	}
	
}

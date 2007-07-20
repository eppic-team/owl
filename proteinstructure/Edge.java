package proteinstructure;
import java.lang.Comparable;

public class Edge implements Comparable {
	
	public int i;
	public int j;

	public Edge(int i,int j){
		this.i=i;
		this.j=j;
	}
	
	public int compareTo(Object o) {
		Edge other = (Edge) o;
		if (this.i>other.i){
			return 1;
		} 
		else if (this.i<other.i){
			return -1;
		}
		else { // only remains case this.i==other.i
			if (this.j>other.j){
				return 1;
			}
			else if (this.j<other.j){
				return -1;
			}
		}				
		return 0; // if none of the conditions before returned, then both i and j are equal 
	}

	public boolean equals(Object o){
		Edge other = (Edge) o;
		if (this.i==other.i && this.j==other.j){
			return true;
		}
		return false;
	}

	public String toString() {
		return this.i+" "+this.j;
	}
	
	public int hashCode() {
		return i*100000+j; // hash function found after a lot of experimenting!! do not touch!
	}

	/**
	 * Gets range (i.e. sequence separation) of contact
	 * @return
	 */
	public int getRange(){
		return Math.abs(this.i-this.j);
	}
	
}


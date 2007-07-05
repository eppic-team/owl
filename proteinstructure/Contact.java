package proteinstructure;
import java.lang.Comparable;

public class Contact implements Comparable {
	
	public int i;
	public int j;

	public Contact(int i,int j){
		this.i=i;
		this.j=j;
	}
	
	public int compareTo(Object o) {
		int diff=0;
		Contact other = (Contact) o;
		if (this.i>other.i){
			diff=1;
		} 
		else if (this.i<other.i){
			diff=-1;
		}
		else if (this.i==other.i){
			if (this.j>other.j){
				diff=1;
			}
			else if (this.j<other.j){
				diff=-1;
			}
		}				
		return diff;
	}

	public boolean equals(Object o){
		boolean eq = false;
		Contact other = (Contact) o;
		if (this.i==other.i && this.j==other.j){
			eq=true;
		}
		return eq;
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


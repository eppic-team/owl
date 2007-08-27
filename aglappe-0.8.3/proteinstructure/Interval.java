package proteinstructure;
import java.lang.Comparable;
/**
 * Class representing a residue interval with a beginning 
 * and an end serial
 * The compareTo method is at the moment exactly the same as
 * in Contact class. It doesn't make much sense like this but we 
 * keep it as it allows Interval objects to be keys in maps
 * 
 * @author Jose Duarte
 *
 */
public class Interval implements Comparable {
	
	public int beg;
	public int end;

	public Interval(int beg,int end){
		this.beg=beg;
		this.end=end;
	}
	
	public int compareTo(Object o) {
		int diff=0;
		Interval other = (Interval) o;
		if (this.beg>other.beg){
			diff=1;
		} 
		else if (this.beg<other.beg){
			diff=-1;
		}
		else if (this.beg==other.beg){
			if (this.end>other.end){
				diff=1;
			}
			else if (this.end<other.end){
				diff=-1;
			}
		}				
		return diff;
	}
	
	public int getLength(){
		return (end - beg);
	}

	public boolean equals(Object o){
		boolean eq = false;
		Interval other = (Interval) o;
		if (this.beg==other.beg && this.end==other.end){
			eq=true;
		}
		return eq;
	}

	public String toString() {
		return this.beg+" "+this.end;
	}
	
	public int hashCode() {
		return this.toString().hashCode();
	}

}


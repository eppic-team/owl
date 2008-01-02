package proteinstructure;
import java.lang.Comparable;
/**
 * Class representing an integer interval with a beginning 
 * and an end integers
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
		Interval other = (Interval) o;
		if (this.beg>other.beg){
			return 1;
		} 
		else if (this.beg<other.beg){
			return -1;
		}
		else if (this.beg==other.beg){
			if (this.end>other.end){
				return 1;
			}
			else if (this.end<other.end){
				return -1;
			}
		}				
		return 0; // if none of the conditions before returned, then both beg and end are equal
	}
	
	public int getLength(){
		return (end - beg);
	}

	public boolean equals(Object o){
		Interval other = (Interval) o;
		if (this.beg==other.beg && this.end==other.end){
			return true;
		}
		return false;
	}

	public String toString() {
		return this.beg+" "+this.end;
	}
	
}


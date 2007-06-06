package proteinstructure;
import java.lang.Comparable;

public class GenericContact implements Comparable {
	
	public String i_res;
	public String j_res;

	public GenericContact(String i_res,String j_res){
		this.i_res=i_res;
		this.j_res=j_res;
	}
	
	public int compareTo(Object o) {
		int diff=0;
		GenericContact other = (GenericContact) o;
		if (this.i_res.compareTo(other.i_res)>0){
			diff=1;
		} 
		else if (this.i_res.compareTo(other.i_res)<0){
			diff=-1;
		}
		else if (this.i_res.equals(other.i_res)){
			if (this.j_res.compareTo(other.j_res)>0){
				diff=1;
			}
			else if (this.j_res.compareTo(other.j_res)<0){
				diff=-1;
			}
		}				
		return diff;
	}

	public boolean equals(Object o){
		boolean eq = false;
		GenericContact other = (GenericContact) o;
		if (this.i_res==other.i_res && this.j_res==other.j_res){
			eq=true;
		}
		return eq;
	}
	
	public String toString() {
		return this.i_res+"\t"+this.j_res;
	}
}


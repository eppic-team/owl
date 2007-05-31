package proteinstructure;

import java.util.TreeMap;
import java.util.Collections;

public class EdgeNbh extends TreeMap<Integer,String> {
	

	private static final long serialVersionUID = 1L;

	// central edge residues
	public int i_resser;
	public String i_resType;
	public int j_resser;
	public String j_resType;
	
	public EdgeNbh(int i_resser, String i_resType, int j_resser, String j_resType){
		super();
		this.i_resser=i_resser;
		this.i_resType=i_resType;
		this.j_resser=j_resser;
		this.j_resType=j_resType;		
	}
	
	public String getMotif(){
		String motif="";
		int min=Math.min(Math.min(i_resser,j_resser), Collections.min(this.keySet()));
		int max=Math.max(Math.max(i_resser,j_resser), Collections.max(this.keySet()));
		for (int i=min;i<=max;i++){
			if (this.containsKey(i)){
				motif+=AA.threeletter2oneletter(this.get(i));	
			} else if (i==i_resser){
				motif+="x";
			} else if (i==j_resser){
				motif+="y";				
			} else {
				motif+="_";
			}
		}
		return motif;
	}
	
}

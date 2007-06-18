package proteinstructure;

import java.util.TreeMap;
import java.util.HashMap;
import java.util.ArrayList;

public class NbhProbDistribution {

	double entropy;
	HashMap<String,Double> dist;
	TreeMap<String,Integer> ranks;
	
	public NbhProbDistribution(HashMap<String,Double> dist) {
		this.dist=dist;
		this.entropy=calculateEntropy();
		getRanks(); //initialises ranks TreeMap
	}
	
	public double getProb(String res){
		return dist.get(res);
	}
	
	public int getRank(String res) {
		return ranks.get(res);
	}

	public double getEntropy(){
		return entropy;
	}

	private void getRanks(){
		ranks = new TreeMap<String, Integer>();
		ArrayList<String> doneRes = new ArrayList<String>(); 
		for (int rank=1;rank<=20;rank++){
			double max = 0.0;
			String maxres="";
			for (String res:dist.keySet()){
				if (!doneRes.contains(res)) {
					double prob = dist.get(res);
					if (prob>=max){
						max = prob;
						maxres = res;
					}
				}
			}
			ranks.put(maxres,rank);
			doneRes.add(maxres);			
		}
	}
	
	private double calculateEntropy(){
		double sumplogp=0.0;
		for (double prob:dist.values()){
			if (prob!=0){ // plogp is defined to be 0 when p=0 (because of limit). If we let java calculate it, it gives NaN (-infinite) because it tries to compute log(0) 
				sumplogp += prob*(Math.log(prob)/Math.log(2));
			}
		}
		return (double) (-1)*sumplogp;
	}
	
}


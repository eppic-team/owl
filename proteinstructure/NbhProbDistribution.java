package proteinstructure;

import java.util.TreeMap;
import java.util.HashMap;

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
		//TODO implement getRanks method, don't know how to do it!!
	}
	
	private double calculateEntropy(){
		double sumplogp=0.0;
		for (double prob:dist.values()){
			sumplogp += prob*(Math.log(prob)/Math.log(2));
		}
		double entropy = (double) (-1)*sumplogp;
		return entropy;
	}
	
}


package owl.core.structure.graphs;

import java.util.TreeMap;
import java.util.ArrayList;


public class NbhProbDistribution {

	final static int MAXRANK = 21;
	
	double entropy;
	TreeMap<String,Double> dist;
	TreeMap<String,Integer> ranks;
	
	public NbhProbDistribution(TreeMap<String,Double> dist) {
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
	
	public ArrayList<String> getResiduesSortedByRank(){
		ArrayList<String> sortedResidues = new ArrayList<String>();
		for (int i=1;i<=21;i++){
			for (String res:ranks.keySet()){
				if (ranks.get(res)<=i && !sortedResidues.contains(res)){
					sortedResidues.add(res);
				}
			}
		}
		return sortedResidues;
	}

	private void getRanks(){
		ranks = new TreeMap<String, Integer>();
		// first we set the residues with prob=0.0 to MAXRANK
		for (String res:dist.keySet()){
			if (!ranks.containsKey(res)) { // we don't check the ones already assigned
				double prob = dist.get(res);
				if (prob==0.0){
					ranks.put(res,MAXRANK);
				}
			}
		}
		// now we set ranks for the rest
		int lastRank = 0;
		double lastMax = 0.0;
		int numberNonZeroProb = dist.size()-ranks.size();
		for (int rank=1;rank<=numberNonZeroProb;rank++){
			double max = 0.0;
			String maxres="";
			for (String res:dist.keySet()){
				if (!ranks.containsKey(res)) {
					double prob = dist.get(res);
					if (prob>=max){
						max = prob;
						maxres = res;
					}
				}
			}
			if (max==lastMax){
				ranks.put(maxres, lastRank);
			} else {
				ranks.put(maxres,rank);
				lastRank = rank;
			}
			lastMax = max;
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


package owl.gmbp;

import java.util.Vector;

public class NbhString_ClusterAnalysis {
	
	public static final int defaultEpsilon = 5;
	public static final int defaultMinNumNBs = 10;
	public static final double defaultRadius = 100;
	
	// --- input arguments
	private Vector<float[]> nbhsNodes;
	private double epsilon;
	private int minNumNBs;
	// --- private members for cluster analysis
	private boolean[] visistedN;
	private int[] clusterN;      // contains clusterIDs: -1=noClusterAllocatedYet 0=Noise ID>0-->belongsToCluster
	private Vector<int[]> edges; // contains array with index of start-node and end-node for each edge --> int[2]
								 // index=start-node --> outgoing edge; index=end-node --> inbound edge
	private int numFoundClusters;	
	private double[][] clusterProp; // [][0]:minL, [][1]:averageL, [][2]:maxL, [][3]:minP, [][4]:averP, [][5]:maxP
	private Vector<Integer>[] clusters; // clusters.length=numFoundClusters+1; 
									//clusters contains a vector of all IDs for noise (background) and each found cluster
	private Vector<Double> slopesOfE;	

	private double rSphere;
	
	public NbhString_ClusterAnalysis(Vector<float[]> nodes){
		this.nbhsNodes = nodes;
		this.epsilon = defaultEpsilon; // should be chosen with respect to rSphere
		this.minNumNBs = defaultMinNumNBs;
		this.rSphere = defaultRadius;
		
		initialiseArguments();
	}
	
	public NbhString_ClusterAnalysis(Vector<float[]> nodes, double eps, int num, double radius){
		this.nbhsNodes = nodes;
		this.epsilon = eps; // should be chosen with respect to rSphere
		this.minNumNBs = num;
		this.rSphere = radius;
		
		initialiseArguments();
	}
	
	// public methods
	public void runClusterAnalysis(){
//		// test functionality of distance methods
//		int steps = 20;
//		double deltaRad = Math.PI/steps;
//		double lC = Math.PI;
//		double pC = Math.PI/2;
//		System.out.println("Test euclidian distance");
//		for (int i=0; i<2*steps; i++){
//			double lambda = i*deltaRad;
//			for (int j=0; j<steps; j++){
//				double phi = j*deltaRad;
//				double dist = getEuclideanDist(lC, pC, lambda, phi);
//				dist = (double)(Math.round(dist*100))/100;
//				System.out.print(dist+"\t");
//			}
//			System.out.println();
//		}
//		System.out.println("Test geodesic distance");
//		for (int i=0; i<2*steps; i++){
//			double lambda = i*deltaRad;
//			for (int j=0; j<steps; j++){
//				double phi = j*deltaRad;
//				double dist = getGeodesicDist(lC, pC, lambda, phi);
//				dist = (double)(Math.round(dist*100))/100;
//				System.out.print(dist+"\t");
//			}
//			System.out.println();
//		}
//		System.out.println("End Test Distance");
//		// end test
		
		int clusterID = 0;
		int index;
		for(int i=0; i<this.nbhsNodes.size(); i++){			
			Vector<Integer> nbIndices = new Vector<Integer>();
			if (!this.visistedN[i]){ // not visited yet
				this.visistedN[i] = true;
				nbIndices = getAllDirectNeighbours(i);
				int numNB = nbIndices.size();
				float[] node = (float[]) this.nbhsNodes.get(i);
				if (nbIndices.size()<this.minNumNBs){
					this.clusterN[i] = 0;
					// Testoutput
//					System.out.println(0+": "+(node[4]+Math.PI)+" , "+node[3]);
				}
				else{
					clusterID++;
					this.clusterN[i] = clusterID;
//					// Testoutput
//					System.out.println(clusterID+": "+(node[4]+Math.PI)+" , "+node[3]+ "  nbIndices.size()="+numNB);
					
					for (int j=0; j<nbIndices.size(); j++){
						index = (int) nbIndices.get(j);
						if (this.clusterN[index] < 0){ // not yet member of any cluster or noise
							expandCluster(index,clusterID);
						}
					}
				}					
			}
		}
		// at the end each node should be element of a cluster (id>0) or noise (id=0)
		this.numFoundClusters = clusterID;
		System.out.println("number of extracted Clusters= "+this.numFoundClusters);
	}
	
	@SuppressWarnings("unchecked")
	public void analyseClusters(){
		if (this.numFoundClusters>0){
//			int[] clusters = new int[this.numFoundClusters+1]; 
			int clusterID;
			clusters = new Vector[this.numFoundClusters+1];
			for (int i=0; i<clusters.length; i++)
				clusters[i] = new Vector<Integer>();
			for (int nodeID=0; nodeID<this.clusterN.length; nodeID++){
				clusterID = this.clusterN[nodeID];
				clusters[clusterID].add(nodeID);
			}
			
			// Compute min, max and average values of lambda and phi for each cluster (besides background)
			// [0]:minL, [1]:averageL, [2]:maxL, [3]:minP, [4]:averP, [5]:maxP
			clusterProp = new double[this.numFoundClusters][6];
			System.out.println("ClusterProperties: ");
			for (int i=1; i<clusters.length; i++){
				double minL=2*Math.PI, minP=Math.PI, maxL=-minL, maxP=-minP;
				double averL=0, averP=0;
				double lambda, phi;
				for (int j=0; j<clusters[i].size(); j++){
					// clusters[i].get(j) == id of node
					float[] node = (float[]) this.nbhsNodes.get(clusters[i].get(j));
					lambda = node[4]; //+Math.PI;
					phi = node[3];
					if (lambda<minL)
						minL = lambda;
					if (lambda>maxL)
						maxL = lambda;
					if (phi<minP)
						minP = phi;
					if (phi>maxP)
						maxP = phi;
					averL += lambda;
					averP += phi;
				}
				averL = averL/clusters[i].size();
				averP = averP/clusters[i].size();
				clusterProp[i-1][0] = minL;
				clusterProp[i-1][1] = averL; 
				clusterProp[i-1][2] = maxL; 
				clusterProp[i-1][3] = minP; 
				clusterProp[i-1][4] = averP; 
				clusterProp[i-1][5] = maxP;  
				System.out.println("ID="+i+"\t"+"#nodes="+clusters[i].size()+"  lambda: "+minL+" : "+averL+" : "+maxL+"\t"+"\t"+"phi: "+minP+" : "+averP+" : "+maxP);
			}
			
			// Test printout
			System.out.println("Clusters: ");
			for (int i=0; i<clusters.length; i++){
				System.out.print("ID="+i+"\t");
				for (int j=0; j<clusters[i].size(); j++)
					System.out.print(clusters[i].get(j)+", ");
				System.out.println();
			}
		}
		
	}
	
	// private methods
	
	private void initialiseArguments(){
		this.visistedN = new boolean[this.nbhsNodes.size()];
		this.clusterN = new int[this.nbhsNodes.size()];
		this.numFoundClusters = 0;
		for(int i=0; i<this.nbhsNodes.size(); i++){
			this.visistedN[i] = false;
			this.clusterN[i] = -1;
		}
		// --extract edges from nbhsNodes-Vector
		float[] node, nbNode;
		this.edges = new Vector<int[]>();
		this.slopesOfE = new Vector<Double>();
		for(int i=0; i<this.nbhsNodes.size()-1; i++){
			node = (float[]) this.nbhsNodes.get(i);	
			nbNode = (float[]) this.nbhsNodes.get(i+1);
			if (node[0]==nbNode[0] && node[1]==nbNode[1]){
				// same trace
				double slope = getSlope(node, nbNode);
				this.slopesOfE.add(slope);
				int[] indicesEdge = {i, i+1};
				this.edges.add(indicesEdge);
			}
		}
	}
	
	private double getSlope(float[] node1, float[] node2){
		double m;
		m = (node2[3]-node1[3])/(node2[4]-node1[4]);
		return m;
	}
	
	/*  Computes the geodesic great circle distance between two points on a sphere,
	 *  given by (lambda1,phi1) and (lambda2,phi2), whereas both angles need to be >0:
	 *  lambda[0:2Pi] and phi[0:Pi].
	 * */
	private double getGeodesicDist(double l1, double p1, double l2, double p2){
		double dist = 0;
		double spherAng = 0;
		p1 -= Math.PI/2;
		p2 -= Math.PI/2;
//		p1 = Math.PI/2 - p1;
//		p2 = Math.PI/2 - p2;
		spherAng = Math.acos((Math.sin(p1)*Math.sin(p2)) + (Math.cos(p1)*Math.cos(p2)*Math.cos(l1-l2)));
		dist = this.rSphere * spherAng; // just a scaling --> doesn't matter how large rSphere is --> just with respect to threshold epsilon
//		dist = spherAng;
		return dist;
	}
	
	/*  Computes the Euclidian distance between two points on a sphere,
	 *  given by (lambda1,phi1) and (lambda2,phi2), whereas both angles need to be >0:
	 *  lambda[0:2Pi] and phi[0:Pi].
	 * */
	@SuppressWarnings("unused")
	private double getEuclideanDist(double l1, double p1, double l2, double p2){
		double dist = 0;
//		double R = this.g2dSize.getHeight() / 2;
		
		double x1 = this.rSphere*Math.sin(p1)*Math.cos(l1);
		double y1 = this.rSphere*Math.sin(p1)*Math.sin(l1);
		double z1 = this.rSphere*Math.cos(p1);
		double x2 = this.rSphere*Math.sin(p2)*Math.cos(l2);
		double y2 = this.rSphere*Math.sin(p2)*Math.sin(l2);
		double z2 = this.rSphere*Math.cos(p2);
		
		dist = Math.sqrt(Math.pow(x1-x2,2) + Math.pow(y1-y2,2) + Math.pow(z1-z2,2));
		return dist;
	}
	
	private Vector<Integer> getAllDirectNeighbours(int index){
		Vector<Integer> nb = new Vector<Integer>();		
		float[] centralN, node;
		centralN = (float[]) this.nbhsNodes.get(index);
		for(int i=0; i<this.nbhsNodes.size(); i++){
//			if (!this.visistedN[i]) // not visited yet 
			if (i!=index)
			{
				node = (float[]) this.nbhsNodes.get(i);
				double dist;
//				dist = getEuclideanDist(centralN[4]+Math.PI, centralN[3], node[4]+Math.PI, node[3]);
				dist = getGeodesicDist(centralN[4]+Math.PI, centralN[3], node[4]+Math.PI, node[3]);
				if (dist<this.epsilon)
					nb.add(i);
			}
		}
		return nb;
	}
		
	private void expandCluster(int index, int clusterID){
		this.clusterN[index] = clusterID;
		int indexNB;
//		// Testoutput
//		float[] node = (float[]) this.nbhsNodes.get(index);
//		System.out.println(clusterID+": "+(node[4]+Math.PI)+" , "+node[3]);
		
		if (!this.visistedN[index]){  // not visited yet
			this.visistedN[index] = true;
			Vector<Integer> nbIndices = getAllDirectNeighbours(index);
			if (nbIndices.size() >= this.minNumNBs){
				for (int j=0; j<nbIndices.size(); j++){
					indexNB = (int) nbIndices.get(j);
					if (this.clusterN[indexNB] == -1){ // not yet member of any cluster or noise
						expandCluster(indexNB,clusterID);
					}
				}
			}
		}
	}
	
	// --- Getters and Setters ---	
	public Vector<float[]> getNbhsNodes() {
		return nbhsNodes;
	}

	public void setNbhsNodes(Vector<float[]> nbhsNodes) {
		this.nbhsNodes = nbhsNodes;
	}

	public double getEpsilon() {
		return epsilon;
	}

	public void setEpsilon(double epsilon) {
		this.epsilon = epsilon;
	}

	public int getminNumNBs() {
		return minNumNBs;
	}

	public void setminNumNBs(int minNumNBs) {
		this.minNumNBs = minNumNBs;
	}
	
	public int[] getClusterN() {
		return clusterN;
	}
	
	public double[][] getClusterProp() {
		return clusterProp;
	}

	public Vector<Integer>[] getClusters() {
		return clusters;
	}

	// --- End Getters and Setters ---	
}

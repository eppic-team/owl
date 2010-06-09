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
	private int numFoundClusters;	
	private double[][] clusterProp; // [][0]:minL, [][1]:averageL, [][2]:maxL, [][3]:minP, [][4]:averP, [][5]:maxP
	private Vector<Integer>[] clusters; // clusters.length=numFoundClusters+1; 
									//clusters contains a vector of all IDs for noise (background) and each found cluster
	private Vector<int[]> edges; // contains array with index of start-node and end-node for each edge --> int[2]
	 							// index=start-node --> outgoing edge; index=end-node --> inbound edge
	private Vector<Double> slopesOfE; // contains slope of each edge (of this.edges) from start-node to end-node
	private double[][] clusterDirProp; // [][0]:minIncSlope, [][1]:averageIncSlope, [][2]:maxIncSlope, [][3]:minOutSlope, [][4]:averOutSlope, [][5]:maxOutSlope
	private int[] nbCluster;    // holds ID of neighbouring cluster (at the other end of edge); if any neighbouring cluster exists, ID=0

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
	
	// ______public methods
	public void runClusterAnalysis(){		
		int clusterID = 0;
		int index;
		for(int i=0; i<this.nbhsNodes.size(); i++){			
			Vector<Integer> nbIndices = new Vector<Integer>();
			if (!this.visistedN[i]){ // not visited yet
				this.visistedN[i] = true;
				nbIndices = getAllDirectNeighbours(i);
//				int numNB = nbIndices.size();
//				float[] node = (float[]) this.nbhsNodes.get(i);
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
	
	public void analyseEdgeDirection(){
		if (this.numFoundClusters<=0)
			runClusterAnalysis();
		if (this.numFoundClusters>0){
			if (this.clusters==null)
				analyseClusters();
			float[] node, nbNode;
			System.out.println("EdgeDirectionAnalysis");
			// Compute min, max and average values of lambda and phi for each cluster (besides background)
			// [][0]:minIncSlope, [][1]:averageIncSlope, [][2]:maxIncSlope, [][3]:minOutSlope, [][4]:averOutSlope, [][5]:maxOutSlope
			clusterDirProp = new double[this.numFoundClusters][6];
			// for each cluster; start at 1, because 0 hold IDs for noise nodes
			for (int i=1; i<this.clusters.length; i++){ 
				double minIncSlope=100, averIncSlope=0, maxIncSlope=0;
				double minOutSlope=100, averOutSlope=0, maxOutSlope=0;
				int cntInc=0, cntOut=0;
				//for each node of cluster
				for (int j=0; j<clusters[i].size(); j++){
					double incSlope, outSlope;
					int nodeID = clusters[i].get(j); // nodeID[0:nbhsNodes.size()-1]
					// look for incoming and outgoing edges
					// compute average slope and deviation of slope for both (incoming and outgoing)
					node = (float[]) this.nbhsNodes.get(nodeID);	
					if (nodeID>0){
						nbNode = (float[]) this.nbhsNodes.get(nodeID-1);
						if (node[0]==nbNode[0] && node[1]==nbNode[1]){
							incSlope = getSlope(node, nbNode);
							averIncSlope += incSlope;
							if (incSlope<minIncSlope)
								minIncSlope = incSlope;
							if (incSlope>maxIncSlope)
								maxIncSlope = incSlope;
							cntInc++;
						}
					}
					if (nodeID<this.nbhsNodes.size()-1){
						nbNode = (float[]) this.nbhsNodes.get(nodeID+1);
						if (node[0]==nbNode[0] && node[1]==nbNode[1]){
							outSlope = getSlope(node, nbNode); //this.slopesOfE.get(nodeID);
							averOutSlope += outSlope;
							if (outSlope<minOutSlope)
								minOutSlope = outSlope;
							if (outSlope>maxOutSlope)
								maxOutSlope = outSlope;
							cntOut++;
						}
					}
				}
				averIncSlope = averIncSlope/cntInc; // /clusters[i].size();
				averOutSlope = averOutSlope/cntOut; // /clusters[i].size();
				double deltaIncSlope = Math.abs(maxIncSlope-minIncSlope);
				double deltaOutSlope = Math.abs(maxOutSlope-minOutSlope);
				clusterDirProp[i-1][0] = minIncSlope;
				clusterDirProp[i-1][1] = averIncSlope; 
				clusterDirProp[i-1][2] = maxIncSlope; 
				clusterDirProp[i-1][3] = minOutSlope; 
				clusterDirProp[i-1][4] = averOutSlope; 
				clusterDirProp[i-1][5] = maxOutSlope;  
				// --- Test output ---
				System.out.println("C_ID="+i+"\t"+"Inc: "+minIncSlope+":"+averIncSlope+":"+maxIncSlope+"\t"
						+"Out: "+minOutSlope+":"+averOutSlope+":"+maxOutSlope+"\t"+"deltaInc="+deltaIncSlope+" deltaOut="+deltaOutSlope);
			}
			
			// for each cluster: is there another cluster, where the majority of edges runs to?
			// for each node in Cluster C1
			this.nbCluster = new int[this.numFoundClusters];
			for (int i=1; i<this.clusters.length; i++){ 
				// for each node in Cluster C1 find following node (nodeNB) in sequence
				// determine which cluster (or noise) nodeNB belongs to
				// if min% of nodeNBs element of same cluster C2 --> C2 is the following cluster of C1
				// otherwise C1 has no following cluster
				int[] countNBclusters = new int[this.clusters.length+1];
				int cnt=0;
				for (int j=0; j<this.clusters[i].size(); j++){
					int nodeID = clusters[i].get(j); // nodeID[0:nbhsNodes.size()-1]
					int nbNodeID = nodeID+1;
					node = (float[]) this.nbhsNodes.get(nodeID);
					if (nbNodeID<this.nbhsNodes.size()){
						nbNode = (float[]) this.nbhsNodes.get(nbNodeID);
						if (!(node[0]==nbNode[0] && node[1]==nbNode[1]))
							nbNodeID = -1; // has no successor
						if (nbNodeID>-1) { // has successor
							int clusterID = this.clusterN[nbNodeID];
							countNBclusters[clusterID]++;
							cnt++;
							// if centre inbetween
							if ((int) node[1]>(int) node[2] && (int) node[1]<(int) nbNode[2]){
								countNBclusters[countNBclusters.length-1]++;
							}
						}							
					}
				}
				System.out.println(this.clusters[i].size()+"nodes of cluster "+i+" end up in: ");
				int idSuccClus = 0;
				boolean centralRes = false;
				for (int j=0; j<countNBclusters.length; j++){
					System.out.print(countNBclusters[j]+" "+(double)countNBclusters[j]/cnt+"\t");
					if ((double)countNBclusters[j]/cnt > 0.5){
						if (j<=numFoundClusters)
							idSuccClus = j;
						else
							centralRes = true;
					}
				}
//				if (idSuccClus>this.numFoundClusters)
				if (centralRes)
					idSuccClus *= -1;
				System.out.println("Successor cluster="+idSuccClus);
				this.nbCluster[i-1] = idSuccClus;
			}
		}
	}
	
	@SuppressWarnings("unchecked")
	public void analyseClusters(){
		if (this.numFoundClusters<=0)
			runClusterAnalysis();
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
	
	//  ________private methods	
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
//		m = (node2[4]-node1[4])/(node2[3]-node1[3]);
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
	
	public double[][] getClusterDirProp() {
		return clusterDirProp;
	}
	
	public int[] getNbCluster(){
		return this.nbCluster;
	}

	public Vector<Integer>[] getClusters() {
		return clusters;
	}

	// --- End Getters and Setters ---	
}

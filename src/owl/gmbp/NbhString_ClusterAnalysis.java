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
	private Vector<Integer>[] clusters; // clusters.length=numFoundClusters+1;  //Vector<Integer>[] clusters;
									//clusters contains a vector of all IDs for noise (background) and each found cluster
	private Vector<int[]> edges; // contains array with index of start-node and end-node for each edge --> int[2]
	 							// index=start-node --> outgoing edge; index=end-node --> inbound edge
	private Vector<Double> slopesOfE; // contains slope of each edge (of this.edges) from start-node to end-node
	private double[][] clusterDirProp; // [][0]:minIncSlope, [][1]:averageIncSlope, [][2]:maxIncSlope, [][3]:minOutSlope, [][4]:averOutSlope, [][5]:maxOutSlope
	private int[] nbCluster;    // holds ID of neighbouring cluster (at the other end of edge); if any neighbouring cluster exists, ID=0
	private Vector<double[]> clusterAverDirec;   // holds direction (vector saved as array) for each cluster

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
	
	/**
	 Performs DBScan-algorithm
	 */		
	public void runClusterAnalysis(){		
		int clusterID = 0;
		int index;
		for(int i=0; i<this.nbhsNodes.size(); i++){			
			Vector<Integer> nbIndices = new Vector<Integer>();
			if (!this.visistedN[i]){ // not visited yet
				this.visistedN[i] = true;
				nbIndices = getAllDirectNeighbours(i);
				if (nbIndices.size()<this.minNumNBs){
					this.clusterN[i] = 0;
				}
				else{
					clusterID++;
					this.clusterN[i] = clusterID;					
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
	
	/**
	 Computes the set of all direct neighbours, for which 
	 the distance is smaller than the given threshold epsilon.
	 @param index of node
	 @return vector with all direct neighbours
	 */	
	private Vector<Integer> getAllDirectNeighbours(int index){
		Vector<Integer> nb = new Vector<Integer>();		
		float[] centralN, node;
		centralN = (float[]) this.nbhsNodes.get(index);
		for(int i=0; i<this.nbhsNodes.size(); i++){
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
	
	
	/**
	 Computes all properties of edge bundles outgoing from each found cluster
	 that might be important for displaying purposes.
	 */		
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
							incSlope = getSlope(nbNode, node);
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
//				minIncSlope+=Math.PI; averIncSlope+=Math.PI; maxIncSlope+=Math.PI;
//				minOutSlope+=Math.PI; averOutSlope+=Math.PI; maxOutSlope+=Math.PI;
				System.out.println("C_ID="+i+"\t"+"Inc: "+minIncSlope+":"+averIncSlope+":"+maxIncSlope+"\t"
						+"Out: "+minOutSlope+":"+averOutSlope+":"+maxOutSlope+"\t"+"deltaInc="+deltaIncSlope+" deltaOut="+deltaOutSlope);
			}
			
			// for each cluster: average direction (vectors)
			// use average lambda and phi as centroid
			this.clusterAverDirec = new Vector<double[]>();
			for (int i=1; i<this.clusters.length; i++){ 
//				System.out.println("C_ID="+i);
				double[] averDir = new double[2];
				double[] dir = new double[2];
				// --- histogram of quadrants --> use quadrant with most outgoing edges for average (remove outliers)
				int[] histQuad = new int[4];
				for (int j=0; j<clusters[i].size(); j++){
					int nodeID = clusters[i].get(j);
					node = (float[]) this.nbhsNodes.get(nodeID);
					nodeID = nodeID+1;
					if (nodeID<this.nbhsNodes.size()){
						nbNode = (float[]) this.nbhsNodes.get(nodeID);
						if (node[0]==nbNode[0] && node[1]==nbNode[1]){
							// node[4]:lambda node[3]:phi
							dir[0] = nbNode[4]-node[4];
							dir[1] = nbNode[3]-node[3];
							dir = normaliseVector(dir);
							int quad = getQuadrant4Vector(dir);
							histQuad[quad]++;
						}
					}
				}
				int quad2use = 0;
				for (int k=1; k<4; k++)
					if (histQuad[k]>histQuad[quad2use])
						quad2use = k;
				System.out.println("C_ID="+i+"  Quadrants occupied: "+histQuad[0]+" - "+histQuad[1]+" - "+histQuad[2]+" - "+histQuad[3]+" --> Quadrant2Use = "+quad2use);
				//for each node of cluster
				int startNode = 0;
				for (int j=0; j<clusters[i].size(); j++){
					int nodeID = clusters[i].get(j);
					node = (float[]) this.nbhsNodes.get(nodeID);
					nodeID = nodeID+1;
					if (nodeID<this.nbhsNodes.size()){
						nbNode = (float[]) this.nbhsNodes.get(nodeID);
						if (node[0]==nbNode[0] && node[1]==nbNode[1]){
							// node[4]:lambda node[3]:phi
							dir[0] = nbNode[4]-node[4];
							dir[1] = nbNode[3]-node[3];
							dir = normaliseVector(dir);
							int quad = getQuadrant4Vector(dir);
							if (quad==quad2use){
								if (j==startNode)
									averDir = dir;
								else{
									averDir[0] = averDir[0]+dir[0];
									averDir[1] = averDir[1]+dir[1];
									averDir = normaliseVector(averDir);
								}								
							}
							else{
								if (j==startNode)
									startNode++;
							}
//							System.out.println(quad+" dir="+String.valueOf(dir[0])+" , "+String.valueOf(dir[1])+" --> averDir="+String.valueOf(averDir[0])+" , "+String.valueOf(averDir[1]));													
						}						
					}
				}
				this.clusterAverDirec.add(averDir);
				System.out.println("C_ID="+i+"  averDir="+String.valueOf(averDir[0])+" , "+String.valueOf(averDir[1]));
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
//				System.out.println(this.clusters[i].size()+"nodes of cluster "+i+" end up in: ");
				int idSuccClus = 0;
				boolean centralRes = false;
				for (int j=0; j<countNBclusters.length; j++){
//					System.out.print(countNBclusters[j]+" "+(double)countNBclusters[j]/cnt+"\t");
					if ((double)countNBclusters[j]/cnt > 0.5){
						if (j<=numFoundClusters)
							idSuccClus = j;
						else
							centralRes = true;
					}
				}
				if (centralRes)
					idSuccClus *= -1;
//				System.out.println("Successor cluster="+idSuccClus);
				System.out.println(this.clusters[i].size()+"nodes of cluster "+i+" end up in: "+idSuccClus);
				this.nbCluster[i-1] = idSuccClus;
			}
		}
	}
	
	/**
	 returns the quadrant for vector
	 0:NorthEast, 1:SouthEast, 2:SouthWest, 3:NorthWest
	 @param vector of direction
	 @return quadrant id
	 */		  
	private int getQuadrant4Vector(double[] vec){
		int quad = 0;
		if (vec[0]>=0 && vec[1]<0)
			quad = 0;
		else if (vec[0]>=0 && vec[1]>0)
			quad = 1;
		else if (vec[0]<0 && vec[1]>0)
			quad = 2;
		else //if (vec[0]<0 && vec[1]<0)
			quad = 3;
		return quad;
	}
	
	private double[] normaliseVector(double[] vec){
		double[] nVec = new double[2];
		double absVec = Math.sqrt((vec[0]*vec[0])+(vec[1]*vec[1]));
		nVec[0] = vec[0]/absVec;
		nVec[1] = vec[1]/absVec;
		return nVec;
	}
	
	/**
	 Computes Cluster properties of contained nodes:
	 minimum, maximum and average values of node positions (lambda and phi).
	 */		
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
				minL+=Math.PI; averL+=Math.PI; maxL+=Math.PI;
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
	
	/**
	 Computes distribution of nodes within certain cluster (clusterID>0)
	 @param ID that defines the cluster 
	 */		
	public double[] analyseNodeDistributionInCluster(int clusterID){
		double[] nodeDistr = null;
		if (clusterID>0 && clusterID<=this.numFoundClusters){			

//			if (this.clusterProp==null)
				analyseClusters();
			// cluster centre
			double averL = this.clusterProp[clusterID-1][1];
		    double averP = this.clusterProp[clusterID-1][4];
		    double maxRad = 0;
		    if (Math.abs(this.clusterProp[clusterID-1][0]-averL)>maxRad)
		    	maxRad = Math.abs(this.clusterProp[clusterID-1][0]-averL);
		    if (Math.abs(this.clusterProp[clusterID-1][2]-averL)>maxRad)
		    	maxRad = Math.abs(this.clusterProp[clusterID-1][2]-averL);
		    if (Math.abs(this.clusterProp[clusterID-1][3]-averP)>maxRad)
		    	maxRad = Math.abs(this.clusterProp[clusterID-1][3]-averP);
		    if (Math.abs(this.clusterProp[clusterID-1][5]-averP)>maxRad)
		    	maxRad = Math.abs(this.clusterProp[clusterID-1][5]-averP);
		    double deltaRad = 0.03;
		    int numRanges = (int)(maxRad/deltaRad) + 1;
		    nodeDistr = new double[numRanges];
		    
			// ID=0 -> Background, ID>0 -> ClusterNr
			Vector<Integer> nodeIDs = this.clusters[clusterID];
			int numNodes = nodeIDs.size();
			double lambda, phi;
			double thres = getGeodesicDist(averL-maxRad, averP-maxRad, averL, averP);
		    System.out.println("Node Distribution within cluster "+clusterID+" numNodes="+numNodes+"  maxRange="+String.valueOf(maxRad)+" thres="+thres+"_"+(thres/this.rSphere));
			
			double rad = deltaRad;
			int i = 0;
			int cnt=0;
			double sum = 0;
//			for (int nodeID:nodeIDs){
//				float[] node = (float[]) this.nbhsNodes.get(nodeID);
//				lambda = node[4]; //+Math.PI;
//				phi = node[3];
////				double dist = getGeodesicDist(lambda, phi, averL, averP);		
//				double dist = getEuclideanDist2d(lambda, phi, averL, averP);
////				System.out.println("node"+nodeID+" l="+(lambda+Math.PI)+" p="+phi+" dist="+dist);				
//			}
			while(rad-deltaRad < maxRad){	
//				double innerT = getGeodesicDist(averL-(rad-deltaRad), averP-(rad-deltaRad), averL, averP);
//				double outerT = getGeodesicDist(averL-rad, averP-rad, averL, averP);
//				System.out.println(innerT+"<dist<"+outerT);
//				innerT = getGeodesicDist(averL, averP-(rad-deltaRad), averL, averP);
//				outerT = getGeodesicDist(averL, averP-rad, averL, averP);
//				System.out.println(innerT+"<dist<"+outerT);
//				innerT = getGeodesicDist(averL-(rad-deltaRad), averP, averL, averP);
//				outerT = getGeodesicDist(averL-rad, averP, averL, averP);
//				System.out.println(innerT+"<dist<"+outerT);
//				double innerT = getEuclideanDist2d(averL-(rad-deltaRad), averP-(rad-deltaRad), averL, averP);
//				double outerT = getEuclideanDist2d(averL-rad, averP-rad, averL, averP);
//				System.out.println(innerT+"<dist<"+outerT);
//				innerT = getEuclideanDist2d(averL, averP-(rad-deltaRad), averL, averP);
//				outerT = getEuclideanDist2d(averL, averP-rad, averL, averP);
//				System.out.println(innerT+"<dist<"+outerT);
				double innerT = getEuclideanDist2d(averL-(rad-deltaRad), averP, averL, averP);
				double outerT = getEuclideanDist2d(averL-rad, averP, averL, averP);
//				System.out.println(innerT+"<dist<"+outerT);


				for (int nodeID:nodeIDs){
					float[] node = (float[]) this.nbhsNodes.get(nodeID);
					lambda = node[4]; //+Math.PI;
					phi = node[3];
//					double dist = getGeodesicDist(lambda, phi, averL, averP);				
//					if (dist>=innerT && dist<outerT)
//						cnt++;
					
					double dist = getEuclideanDist2d(lambda, phi, averL, averP);
					if (dist>=innerT && dist<outerT)
						cnt++;

//					dist = dist/rSphere;
//					if (dist>=rad-deltaRad && dist<rad)
//						cnt++;
				}
				
//				double distr = (double)cnt/(double)numNodes;
				nodeDistr[i] = (double)cnt/(double)numNodes;
				sum += nodeDistr[i];
//				System.out.printf("i=%d \t rad<%5.3f \t cnt=%d \t distr=%5.3f sum=%5.3\n", i, rad, cnt, nodeDistr[i], sum); 
				System.out.println("i="+i+"  rad<"+rad+"   "+innerT+"<dist<"+outerT+"   cnt="+cnt+"   distr="+nodeDistr[i]+ "   sum="+sum);
				i++;
				rad+=deltaRad;
				cnt=0;
			}			
			
		}
		else
			System.out.println("Distribution can't be computed of background (noise) cluster!");
		
		return nodeDistr;
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
	
	/**
	 returns the slope of a line through two nodes
	 @param position one
	 @param position two
	 @return slope
	 */		 
	private double getSlope(float[] node1, float[] node2){
		double m;
		m = (node2[3]-node1[3])/(node2[4]-node1[4]);
//		m = (node2[4]-node1[4])/(node2[3]-node1[3]);
		return m;
	}
	
	/**
	 Computes the geodesic great circle distance between two points on a sphere, 
	 whereas lambda[0:2Pi] and phi[0:Pi].
	 @param lambda value of first point
	 @param phi value of first point
	 @param lambda value of second point
	 @param phi value of second point
	 @return geodesic distance 
	 */		 
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
	
	/**
	 Computes the Euclidian distance between two points on a sphere, 
	 whereas lambda[0:2Pi] and phi[0:Pi].
	 @param lambda value of first point
	 @param phi value of first point
	 @param lambda value of second point
	 @param phi value of second point
	 @return Euclidian distance 
	 */	
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
	
	private double getEuclideanDist2d(double l1, double p1, double l2, double p2){
		double dist = Math.sqrt(Math.pow(l1-l2, 2) + Math.pow(p1-p2, 2));
		return dist;
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

	public int getMinNumNBs() {
		return minNumNBs;
	}

	public void setMinNumNBs(int minNumNBs) {
		this.minNumNBs = minNumNBs;
	}
	
	public void setRadiusSphere(double rad){
		this.rSphere = rad;
	}
	
	public int[] getClusterN() {
		return clusterN;
	}
	
	public double[][] getClusterProp() {
		return clusterProp;
	}
	
	public Vector<double[]> getClusterAverDirec(){
		return clusterAverDirec;
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

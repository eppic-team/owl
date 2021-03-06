package owl.core.util;


import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;



/**
 * An implementation of a single linkage clusterer
 * 
 * See http://en.wikipedia.org/wiki/Single-linkage_clustering
 * 
 * @author jose
 *
 */
public class SingleLinkageClusterer {
	
	private class LinkedPair {
		
		private int first;
		private int second;
		private double closestDistance;
		
		public LinkedPair(int first, int second, double minDistance) {
			this.first = first;
			this.second = second;
			this.closestDistance = minDistance;
		}
		
		public int getFirst() {
			return first;
		}
		
		public int getSecond() {
			return second;
		}
		
		public double getClosestDistance() {
			return closestDistance;
		}
		
		public String toString() {
			
			String closestDistStr = null;
			if (closestDistance==Double.MAX_VALUE) {
				closestDistStr = String.format("%6s", "inf");
			} else {
				closestDistStr = String.format("%6.2f",closestDistance);
			}
			
			return "["+first+","+second+"-"+closestDistStr+"]";
		}
		
	}
	
	private boolean debug;

	private double[][] matrix;
	
	private boolean isScoreMatrix;
	
	private int numItems;
	
	private LinkedPair[] dendrogram;
	
	//private Set<Integer> toSkip;
	
	private ArrayList<Integer> indicesToCheck;
	
	
	/**
	 * Constructs a new SingleLinkageClusterer
	 * Subsequently use {@link #getDendrogram()} to get the full tree  
	 * or {@link #getClusters(double)} to get the clusters at a certain cutoff in the tree
	 * Please note that the matrix will be altered during the clustering procedure. A copy must be 
	 * made before by the user if needing to use the original matrix further.
	 * @param matrix the distance matrix with distance values in j>i half, all other values will be ignored
	 * @param isScoreMatrix if false the matrix will be considered a distance matrix: lower values (distances) mean closer objects,
	 * if true the matrix will be considered a score matrix: larger values (scores) mean closer objects
	 * @throws IllegalArgumentException if matrix not square 
	 */
	public SingleLinkageClusterer(double[][] matrix, boolean isScoreMatrix) {
		this.matrix = matrix;
		this.debug = false;
		this.isScoreMatrix = isScoreMatrix;
		
		if (matrix.length!=matrix[0].length) {
			throw new IllegalArgumentException("Distance matrix for clustering must be a square matrix");
		}

		this.numItems = matrix.length;
				
	}
	
	public void setDebug() {
		this.debug = true;
	}
	
	/**
	 * Get the full dendrogram (size n-1) result of the hierarchical clustering
	 * @return
	 */
	public LinkedPair[] getDendrogram() {
		if (dendrogram==null) {
			clusterIt();
		}
		
		return dendrogram;
	}
	
	/**
	 * Calculate the hierarchical clustering and store it in dendrogram array
	 * This is the naive implementation (o(n3)) of single linkage clustering as outlined in wikipedia:
	 * http://en.wikipedia.org/wiki/Single-linkage_clustering
	 */
	private void clusterIt() {
		
		dendrogram = new LinkedPair[numItems-1];
		
		if (debug) {
			System.out.println("Initial matrix");
			printMatrix();
		}
		
		for (int m=0;m<numItems-1;m++) {
			
			updateIndicesToCheck(m);
			LinkedPair pair = getClosestPair();
			merge(pair);
			dendrogram[m] = pair;
			
			//if (debug) {
			//	System.out.println("Matrix after iteration "+m+" (merged "+pair.getFirst()+","+pair.getSecond()+")");
			//	printMatrix();
			//}
		}
		
	}
	
	/**
	 * Merge 2 rows/columns of the matrix by the linkage function (see {@link #link(double, double)}
	 * @param closestPair
	 */
	private void merge(LinkedPair closestPair) {
		
		
		int first = closestPair.getFirst();
		int second = closestPair.getSecond();
		
		for (int other=0;other<numItems;other++) {
			matrix[Math.min(first,other)][Math.max(first, other)] = link(getDistance(first, other), getDistance(second, other));
		}
		
	}
	
	/**
	 * The linkage function: minimum of the 2 distances (i.e. single linkage clustering)
	 * @param d1
	 * @param d2
	 * @return
	 */
	private double link(double d1, double d2) {
		if (isScoreMatrix) {
			return Math.max(d1,d2);
		} else {
			return Math.min(d1,d2);
		}
	}
	
	private double getDistance(int first, int second) {
		return matrix[Math.min(first, second)][Math.max(first, second)];
	}
	
	private void updateIndicesToCheck(int m) {
		
		if (indicesToCheck==null) {
			indicesToCheck = new ArrayList<Integer>(numItems);

			for (int i=0;i<numItems;i++) {
				indicesToCheck.add(i);
			}			
		}
		
		if (m==0) return;
		
		indicesToCheck.remove(new Integer(dendrogram[m-1].getFirst()));
	}
	
	private LinkedPair getClosestPair() {		
		
		LinkedPair closestPair = null;
		
		if (isScoreMatrix) {
			double max = 0.0;
			for (int i:indicesToCheck) {
				
				for (int j:indicesToCheck) {					
					if (j<=i) continue;
					
					if (matrix[i][j]>=max) {
						max = matrix[i][j];
						closestPair = new LinkedPair(i,j,max);
					}

				}
			}
		} else {
			double min = Double.MAX_VALUE;
			for (int i:indicesToCheck) {
				
				for (int j:indicesToCheck) {
					if (j<=i) continue;
					
					if (matrix[i][j]<=min) {
						min = matrix[i][j];
						closestPair = new LinkedPair(i,j,min);
					}

				}
			}
		}
		
		return closestPair;
	}
	
	/**
	 * Get the clusters by cutting the dendrogram at given cutoff
	 * @param cutoff
	 * @return
	 */
	public Map<Integer, Set<Integer>> getClusters(double cutoff) {
		
		if (dendrogram==null) {
			clusterIt();
		}
		
		Map<Integer, Set<Integer>> clusters = new TreeMap<Integer, Set<Integer>>();
		
		int clusterId = 1;
		
		for (int i=0;i<numItems-1;i++) {
			
			if (isWithinCutoff(i, cutoff)) {
				
				//int containingClusterId = getContainingCluster(clusters, dendrogram[i]);
				
				int firstClusterId = -1;
				int secondClusterId = -1;
				for (int cId:clusters.keySet()) {
					Set<Integer> members = clusters.get(cId);
					
					if (members.contains(dendrogram[i].getFirst())) {
						firstClusterId = cId;
					}
					if (members.contains(dendrogram[i].getSecond())) { 
						secondClusterId = cId;
					}
				}
				
				
				if (firstClusterId==-1 && secondClusterId==-1) {
					// neither member is in a cluster yet, let's assign a new cluster and put them both in
					Set<Integer> members = new TreeSet<Integer>();
					members.add(dendrogram[i].getFirst());
					members.add(dendrogram[i].getSecond());
					clusters.put(clusterId, members);
					clusterId++;
				} else if (firstClusterId!=-1 && secondClusterId==-1) {
					// first member was in firstClusterId already, we add second										
					clusters.get(firstClusterId).add(dendrogram[i].getSecond());
				} else if (secondClusterId!=-1 && firstClusterId==-1) {
					// second member was in secondClusterId already, we add first										
					clusters.get(secondClusterId).add(dendrogram[i].getFirst());					
				} else {
					// both were in different clusters already
					// we need to join them: necessarily one must be of size 1 and the other of size>=1
					Set<Integer> firstCluster = clusters.get(firstClusterId);
					Set<Integer> secondCluster = clusters.get(secondClusterId);
					if (firstCluster.size()<secondCluster.size()) {
						if (debug) System.out.println("Joining cluster "+firstClusterId+" to cluster "+secondClusterId);
						// we join first onto second
						for (int member : firstCluster) {
							secondCluster.add(member);
						}
						clusters.remove(firstClusterId);						
					} else {
						if (debug) System.out.println("Joining cluster "+secondClusterId+" to cluster "+firstClusterId);
						// we join second onto first
						for (int member : secondCluster) {
							firstCluster.add(member);
						}
						clusters.remove(secondClusterId);
					}
				}
				
				if (debug) System.out.println("Within cutoff:     "+dendrogram[i]);
			
			} else {
				
				if (debug) System.out.println("Not within cutoff: "+dendrogram[i]);
				
			}
		}
		
		// reassigning cluster numbers by creating a new map (there can be gaps in the numbering if cluster-joining happened)
		Map<Integer,Set<Integer>> finalClusters = new TreeMap<Integer, Set<Integer>>();
		int newClusterId = 1;
		for (int oldClusterId:clusters.keySet()) {			
			finalClusters.put(newClusterId, clusters.get(oldClusterId));
			newClusterId++;
		}
		
		// anything not clustered is assigned to a singleton cluster (cluster with one member)
		for (int i=0;i<numItems;i++) {
			boolean isAlreadyClustered = false;
			for (Set<Integer> cluster:finalClusters.values()) { 
				if (cluster.contains(i)) {
					isAlreadyClustered = true;
					break;
				}
			}
			if (!isAlreadyClustered) {				
				Set<Integer> members = new TreeSet<Integer>();
				members.add(i);
				finalClusters.put(newClusterId, members);
				newClusterId++;
			}

		}

		// print all clusters if in debug mode
		if (debug) {
			for (int cId:finalClusters.keySet()) {
				System.out.print(cId+": ");
				for (int member:finalClusters.get(cId)) {
					System.out.print(member+" ");
				}
				System.out.println();
			}
		}

		return finalClusters;
	}
	
	private boolean isWithinCutoff(int i, double cutoff) {
		if (isScoreMatrix) {
			if (dendrogram[i].getClosestDistance()>cutoff) {
				return true;
			} else {
				return false;
			}
		} else {
			if (dendrogram[i].getClosestDistance()<cutoff) {
				return true;
			} else {
				return false;
			}			
		}
	}
	
	private void printMatrix() {
		for (int i=0;i<numItems;i++) {
			for (int j=0;j<numItems;j++) {
				if (i==j) {
					System.out.printf("%6s ","x");
				}
				else if (i<j) {
					if (matrix[i][j]==Double.MAX_VALUE) System.out.printf("%6s ","inf");
					else System.out.printf("%6.2f ",matrix[i][j]);
				}
				else {
					if (matrix[j][i]==Double.MAX_VALUE) System.out.printf("%6s ","inf");
					else System.out.printf("%6.2f ",matrix[j][i]);
				}
			}
			System.out.println();
		}
		System.out.println();
	}

}

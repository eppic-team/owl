import proteinstructure.*;
import java.util.*;

/**
 * Perform a connected-component clustering for contact graphs (like the one implemented in matlab)
 * @author stehr
 *
 */
public class testGraphClustering {
	
	// TODO: This whole thing does currently not work
	
	HashMap<Edge,Integer> e2c;
	TreeSet<Integer> clusters;
	HashMap<Integer,EdgeSet> members;
	EdgeSet edges;
	int seqLen;
	
	/* constructor */
	public testGraphClustering() {
		try {
			e2c = new HashMap<Edge,Integer>();					// map from edge to cluster
			clusters = new TreeSet<Integer>();					// set of all non-empty clusters
			members = new HashMap<Integer,EdgeSet>();			// cluster members
			
			Pdb pdb = new PdbasePdb("3eca");
			pdb.load("A");
			Graph g = pdb.get_graph("Ca", 8.0);
			edges = g.getContacts();
			seqLen = g.getFullLength();
			
			System.out.println("Sequence length: " + seqLen);
			System.out.println("Number of contacts: " + edges.size());
			
		} catch(Exception e)  {
			e.printStackTrace();
		}
	}
		
	private int dist(Edge e, Edge f) {
		return Math.abs(e.i - f.i) + Math.abs(e.j - f.j);
	}
	
	private void merge(Edge e1, Edge e2) {
		int cid1 = e2c.get(e1);
		int cid2 = e2c.get(e2);
		EdgeSet cluster1 = members.get(cid1);
		EdgeSet cluster2 = members.get(cid2);
		if(cluster1.size() > cluster1.size()) {			// swap clusters
			EdgeSet temp = cluster1;
			cluster1 = cluster2;
			cluster2 = temp;
			int tmp = cid1;
			cid1 = cid2;
			cid2 = tmp;
		}
		//System.out.println("Merging clusters " + cluster1 + " and " + cluster2);
		for(Edge e3:cluster2) {
			e2c.put(e3, cid1);						// point all edges to cluster1
		}
		cluster1.addAll(cluster2);						// put all edges to cluster1
		cluster2.removeAll(cluster2);					// empty cluster2
		if(!clusters.remove(cid2)) {
			System.err.println("Error!");						// trash cluster2
		}
		//System.out.println("Result: " + cluster1 + " and " + cluster2);
	}
	
	public static void main(String[] args) {
		
		testGraphClustering testGC = new testGraphClustering();
		testGC.doGridClustering(4,3,15);
		
	}
	
	public void doGridClustering(int maxManhDist, int minClustSize, int minSeqSep) {
					
			int size = (seqLen / maxManhDist) + 1;
			System.out.println("Grid size: " + size);
			
			EdgeSet[][] grid = new EdgeSet[size][size];						// grid of clusters
			
			int cid = 1;
			for(Edge e:edges) {
				EdgeSet newCluster = new EdgeSet();
				newCluster.add(e);
				clusters.add(cid);									// add cluster to clusters
				members.put(cid, newCluster);
				e2c.put(e, cid);
				cid++;
			}
			
			// put edges in grid cells
			for(Edge e:edges) {
				int x = e.i / maxManhDist;
				int y = e.j / maxManhDist;
				if(grid[x][y] == null) {
					EdgeSet gridCluster = new EdgeSet();
					gridCluster.add(e);
					grid[x][y] = gridCluster;
				} else {
					grid[x][y].add(e);
				}
			}
			
			System.out.println("Before clustering: Number of clusters:" + clusters.size());
			
			boolean change = true;
			while(change) {
				change = false;
				System.out.println("While clustering: Number of clusters:" + clusters.size());
				for(Edge e:edges) {
					int x = e.i / maxManhDist;
					int y = e.j / maxManhDist;
					// search my own bin
					for(Edge e2:grid[x][y]) {
						if(dist(e,e2) <= maxManhDist) {
							int cid1 = e2c.get(e);
							int cid2 = e2c.get(e2);
							if(cid1 != cid2) {
								merge(e,e2);
								change=true;
							}
						}
					}
					// search bin to the right
					if(x+1 < size && grid[x+1][y] != null) {
						for(Edge e2:grid[x+1][y]) {
							if(dist(e,e2) <= maxManhDist) {
								int cid1 = e2c.get(e);
								int cid2 = e2c.get(e2);
								if(cid1 != cid2) {
									merge(e,e2);
									change=true;
								}
							}
						}
					}
					// search bin below
					if(y+1 < size && grid[x][y+1] != null) {
						for(Edge e2:grid[x][y+1]) {
							if(dist(e,e2) <= maxManhDist) {
								int cid1 = e2c.get(e);
								int cid2 = e2c.get(e2);
								if(cid1 != cid2) {
									merge(e,e2);
									change=true;
								}
							}
						}
					}
					// search bin diagonal
					if(y+1 < size && x+1 < size && grid[x+1][y+1] != null) {
						for(Edge e2:grid[x+1][y+1]) {
							if(dist(e,e2) <= maxManhDist) {
								int cid1 = e2c.get(e);
								int cid2 = e2c.get(e2);
								if(cid1 != cid2) {
									merge(e,e2);
									change=true;
								}
							}
						}
					}				
				}
			}
			
			// Verify clustering
			for(Integer cId1:clusters) {
				for(Integer cId2:clusters) {
					if(cId1 != cId2) { 
						for(Edge e1: members.get(cId1)) {
							for(Edge e2: members.get(cId2)) {
								if(dist(e1,e2) <=  maxManhDist) System.err.println("Error: " + e1 + "," + e2 + " in " + cId1 + " and " + cId2);
							}
						}
					}
				}
			}
			
			System.out.println("Finished clustering. Number of clusters:" + clusters.size());
			int numEdges = 0;
			for(Integer cId:clusters) {
				numEdges += members.get(cId).size();
			}
			System.out.println("Number of edges in clusters: " + numEdges);
			
			TreeSet<Integer> filteredClusters = new TreeSet<Integer>();
			
			// delete small clusters
			for(Integer c:clusters) {
				EdgeSet clust = members.get(c);
				if(clust.size() >= minClustSize) filteredClusters.add(c);
			}
			// delete short range clusters
			for(Edge e:edges) {
				if(Math.abs(e.i - e.j) < minSeqSep) {
					filteredClusters.remove(e2c.get(e));
				}
			}
					
			System.out.println("Finished filtering. Number of clusters:" + filteredClusters.size());
			numEdges = 0;
			for(Integer cId:filteredClusters) {
				numEdges += members.get(cId).size();
			}
			System.out.println("Number of edges in clusters: " + numEdges);
			
			int i = 1;
			for(Integer cId:filteredClusters) {
				EdgeSet c = members.get(cId);
				System.out.println(i++ + " (" + c.size() + "):" + c);
			}
			

	}
}

package owl.mutanom;

import java.io.*;
import java.util.*;

import javax.vecmath.Point3d;

import owl.core.structure.PdbChain;
import owl.core.structure.PdbChainSet;


/**
 * Simulate sets of mutations as a background model for statistical analysis 
 * @author stehr
 *
 */
public class MutationSimulator {

	/*--------------------------- member variables --------------------------*/
	
	Collection<String> pdbCodes = null;
	Random rand = new Random();
	
	/*----------------------------- constructors ----------------------------*/
	public MutationSimulator(Collection<String> pdbCodes) {
		this.pdbCodes = pdbCodes;
	}
	
	/*---------------------------- public methods ---------------------------*/
	/**
	 * Returns a random k-subset of indices between 1 and n.
	 */
	public Set<Integer> getRandomSubset(int n, int k) {
		Set<Integer> subSet = new HashSet<Integer>();
		int c = 0;
		while(c < k) {
			int r = rand.nextInt(n);		
			if(!subSet.contains(r+1)) {
				subSet.add(r+1);
				c++;
			}
		}
		return subSet;
	}
	
//	public double getAveragePathLength(RIGraph rig, Set<Integer> nodeIndices) {
//		HashSet<RIGNode> nodeSet = new HashSet<RIGNode>();
//		for(int idx:nodeIndices) {
//			
//		}
//		return 0;
//	}
	
	/**
	 * Returns the average euclidian distance of the given C-alpha positions
	 * or NaN if some of the coordinates are not available.
	 * @param pdb
	 * @param nodeIndices
	 * @return
	 */
	public double getAverageDistance(PdbChain pdb, Set<Integer> nodeIndices) {
		Vector<Integer> indices = new Vector<Integer>(nodeIndices);
		double sumDist = 0;
		int numDist = 0;
		int n = nodeIndices.size();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < i; j++) {
				int ix = indices.get(i);
				int iy = indices.get(j);
				if(pdb.hasCoordinates(ix, "CA") && pdb.hasCoordinates(iy, "CA")) {
					Point3d px = pdb.getAtomCoord(ix, "CA"); // TODO: use other CT?
					Point3d py = pdb.getAtomCoord(iy, "CA");
					double dist = px.distance(py);
					sumDist += dist;
					numDist ++;
				} else {
					System.err.println("No coordinates for " + ix + " and " + iy);
				}
			}
		}
		if(numDist > 0) {
			return sumDist / numDist;		
		} else {
			return Double.NaN;
		}
	}
		
	public TreeSet<Double> generateDistanceStatistics(int nProt, int nRandSets, int nRes) {
		Vector<String> pdbCodes = new Vector<String>(this.pdbCodes);
		TreeSet<Double> stats = new TreeSet<Double>();
		for (int i = 0; i < Math.min(nProt, pdbCodes.size()); i++) {
			// load pdb
			PdbChain pdb = PdbChain.readStructureOrExit(pdbCodes.get(i));
			Vector<Integer> observedResidues = new Vector<Integer>(pdb.getAllStdAaResSerials());
			// generate random subsets
			for (int j = 0; j < nRandSets; j++) {
				Set<Integer> randIndices = getRandomSubset(observedResidues.size(), nRes); // TODO: what about unobserved residues?
				//System.out.println(randIndices); // DEBUG
				// calculate average distance
				Set<Integer> randResSers = new TreeSet<Integer>();
				for(int idx:randIndices) {
					randResSers.add(observedResidues.get(idx-1));
				}
				double avgDist = getAverageDistance(pdb, randResSers);
				// store in treeset				
				if(avgDist > 0) stats.add(avgDist);
			}
		}
		return stats;
	}
	
	public void generateShortestPathStatistics() {
		// for a set of proteins
		// pick a random subset of nodes
		// determine the averagePathLength
		// store in a treeset
		// save tree to file (for testing: store in member variable)
	}
	
	public void getSignificance(double averageShortestPath) {
		// load statistics from file
		// store in appropriate data structure (treemap<double, integer>)
		// determine significance for given asp value
	}
	/*---------------------------- static methods ---------------------------*/
	public static void testGetRandomSubset() {
		int n = 100;
		int k = 20;
		int trials = 500;
		MutationSimulator ms = new MutationSimulator(null);
		for (int i = 0; i < trials; i++) {
			Set<Integer> rSet = ms.getRandomSubset(n, k);
			System.out.println(rSet);
			if(rSet.size() != k) System.err.printf("Wrong size (%d instead of %d)\n", rSet.size(), k);
			for(int j:rSet) {
				if(j < 1 || j > n) System.err.printf("Value %d out of bounds\n", j);
			}
		}
	}
	
	public static void testGetProbability() {
		//File outFile = new File("test.txt");		
		try {
			//PrintWriter out = new PrintWriter(outFile);
			TreeMap<Double, Integer> cdf = loadSetFromFile(new File("dists_2.txt"));
			for (int i = 0; i < 50; i++) {
				double x = Math.random() * 150;
				//out.printf("%5.1f\t%4.2f\n", x, getProbability(cdf, x));
				System.out.printf("P(X < %5.1f) = %4.2f\n", x, getProbability(cdf, x));
			}
			//out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void writeSetToFile(Set<Double> set, File file) throws FileNotFoundException {
		PrintWriter out = new PrintWriter(file);
		for(double num:set) {
			out.printf("%f\n", num);
		}
		out.close();
	}
	
	public static TreeMap<Double, Integer> loadSetFromFile(File file) throws IOException {
		TreeMap<Double, Integer> map = new TreeMap<Double, Integer>();
		BufferedReader in = new BufferedReader(new FileReader(file));
		String line;
		int c = 1;
		while((line=in.readLine()) != null) {
			map.put(Double.parseDouble(line.trim()), c);
			c++;
		}
		in.close();
		return map;
	}
	
	/**
	 * Given an empirical distribution function as a treemap, where
	 * values are mapped to their index in an ordered set of the values,
	 * evaluates for a given value x the probability to draw values
	 * from the distribution that are smaller than x.
	 * Formally, this means evaluating the approximated
	 * cumulative distribution function.
	 * TODO: Change from TreeMap to TreeSet
	 * @return the probability to draw a value smaller than x from the given distribution
	 */
	public static double getProbability(TreeMap<Double, Integer> cdf, double x) {
		int idx = cdf.headMap(x).size();
		return 1.0 * idx / cdf.size();
	}
	
	public static void makeDistanceStatistics() {
		System.out.println("Reading cullpdb_20...");
		Collection<String> pdbCodes = null;
		try {
			pdbCodes = PdbChainSet.readCullPdb20List();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int maxProt = 500;
		int nRandSets = 2000;
		MutationSimulator ms = new MutationSimulator(pdbCodes);		
		TreeSet<Double> stats = null;
		
		for(int setSize = 1; setSize <= 7; setSize++) {
			stats = ms.generateDistanceStatistics(maxProt, nRandSets, setSize);
			System.out.println(stats.size());
			// write to file
			String fileName = "dists_" + setSize + ".txt";
			try {
				writeSetToFile(stats, new File(fileName));
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	/*--------------------------------- main --------------------------------*/
	
	public static void main(String[] args) {
		//testGetRandomSubset();
		testGetProbability();
		System.out.println("done.");
	}
	
}

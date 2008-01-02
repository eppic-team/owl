import proteinstructure.*;
import java.util.*;

/**
 * Calculate a matrix containing the element-wise differences between the distance maps of two proteins
 * @author stehr
 *
 */
public class testDeltaDistanceMap {
	
	/**
	 * Run tests
	 */
	public static void main(String[] args) {
		try {
			String pdbCode1 = "12as";
			String chainCode1 = "A";
			String pdbCode2 = "12as";
			String chainCode2 = "B";
			
			System.out.println("Loading pdb objects...");
			Pdb pdb1 = new PdbasePdb("12as");
			pdb1.load("B");
			assert(pdb1 != null);
			Pdb pdb2 = new PdbasePdb("12as");
			pdb2.load("A");
			assert(pdb2 != null);
			
			System.out.println("Calculating distance maps...");
			HashMap<Edge,Double> distMap1 = pdb1.calculate_dist_matrix("Ca");
			assert(distMap1 != null);
			HashMap<Edge,Double> distMap2 = pdb2.calculate_dist_matrix("Ca");
			assert(distMap2 != null);
			
			System.out.println("Calculating difference distance map...");
			HashMap<Edge,Double> diffDistMap = pdb1.getDiffDistMap("Ca", pdb2, "Ca");
			assert(diffDistMap != null);
			assert(diffDistMap.size() == distMap1.size());
			assert(diffDistMap.size() == distMap2.size());
			
			System.out.print("Missing matrix values (dim=" + pdb1.getFullLength() + "): ");
			int mis = 0;
			for(int i=1; i<=pdb1.getFullLength();i++) {
				for(int j=1; j<i;j++) {
					if(!diffDistMap.containsKey(new Edge(j,i))) {
						//System.out.print("(" + i + "," + j + ") ");
						mis++;
					}
				}
			}
			System.out.println(mis);
			
			Edge e = new Edge(27,30);

			double dist1 = distMap1.get(e);
			double dist2 = distMap2.get(e);
			double diffDist = diffDistMap.get(e);
			assert(diffDist == Math.abs(dist1 - dist2));
			double min = Collections.min(diffDistMap.values());
			double max = Collections.max(diffDistMap.values());
			
			System.out.println("Checking difference distance map for " + pdbCode1 + chainCode1 + " and " + pdbCode2 + chainCode2);
			System.out.println("size=" + diffDistMap.size() + " min=" + min + " max= " + max);
		}
		catch(Exception e) {
			e.printStackTrace();

		}

	}

}

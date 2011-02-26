package owl.core.structure.graphs;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

import owl.core.structure.AminoAcid;
import owl.core.util.MySQLConnection;


public class NbhoodBackgrnd {

	private ArrayList<RIGNbhood> allNbhoods;
	private HashMap<String, ArrayList<RIGNbhood>> hash;
	
	private String db;
	private Integer[] graphids;
	private MySQLConnection conn;
	
	public NbhoodBackgrnd(MySQLConnection conn, String db, Integer[] graphids) {
		this.conn = conn;
		this.db = db;
		this.graphids = graphids;
	}

	public void load() {
		System.out.println(graphids.length+" graphs to read from database "+db);
		
		allNbhoods = new ArrayList<RIGNbhood>();

		int numGraphs = 0;
		int numNodes = 0;
		
		System.out.println("Reading graphs...");
		
		for (int i = 0;i<graphids.length;i++){
			try {
				if (i%1000==0) {
					System.out.println();
				}
				if (i%100==0) {
					System.out.printf("%5d ",i);
				}

				// get graph
				RIGraph graph = new DbRIGraph(db,conn,graphids[i]);	

				for (RIGNode node:graph.getVertices()) {
					numNodes++;
					allNbhoods.add(graph.getNbhood(node));
				}

				numGraphs++;

			} catch (GraphIdNotFoundError e) {
				System.err.println(e.getMessage());
			} catch (SQLException e) {
				System.err.println(e.getMessage());
			} 
		}
		System.out.println("\nDone reading neighbourhoods. Read "+numGraphs+" graphs and "+numNodes+" neighborhoods");

		// hashing
		hash();
	}
	
	private void hash() {
		System.out.println("Hashing...");
		hash = new HashMap<String, ArrayList<RIGNbhood>>();
		for (RIGNbhood bckgrdNbh:allNbhoods) {
			for (AminoAcid aa1:AminoAcid.getAllStandardAAs()) {
				for (AminoAcid aa2:AminoAcid.getAllStandardAAs()) {
					if (aa2.getThreeLetterCode().compareTo(aa1.getThreeLetterCode())>=0 && bckgrdNbh.containsResType(aa1.getThreeLetterCode()) && bckgrdNbh.containsResType(aa2.getThreeLetterCode())) {
						//ResPair pair = new ResPair(aa1,aa2);
						if (hash.containsKey(aa1.getThreeLetterCode()+aa2.getThreeLetterCode())) {
							hash.get(aa1.getThreeLetterCode()+aa2.getThreeLetterCode()).add(bckgrdNbh);
						} else {
							ArrayList<RIGNbhood> al = new ArrayList<RIGNbhood>();
							al.add(bckgrdNbh);
							hash.put(aa1.getThreeLetterCode()+aa2.getThreeLetterCode(), al);
						}
					}
				}
			}
		}
		System.out.println("Done hashing. Number of hash groups: " + hash.size());
	}

	public NbhProbDistribution getNbhProbDist(RIGNbhood nbh) {

		// we take 1st and last elements of the neighborhood as our hash keys, just because there is a convenient way to grab them
		String aa1 = nbh.get(nbh.firstKey()).getResidueType();
		String aa2 = nbh.get(nbh.lastKey()).getResidueType();
		if (aa2.compareTo(aa1)<0) {
			String tmp = aa1;
			aa1 = aa2;
			aa2 = tmp;
		}

		int countTotal = 0;
		int[] countPerAA = new int[20]; // order is as in AAinfo.getAAs()
		//for (RIGNbhood bckgrdNbh:allNbhoods) {      // version without hashing (a lot slower)
		for (RIGNbhood bckgrdNbh:hash.get(aa1+aa2)) { // we get only the set of neighborhoods that contain aa1 and aa2
			if (nbh.match(bckgrdNbh)) {
				countTotal++;
				int i = 0;
				for (AminoAcid aa:AminoAcid.getAllStandardAAs()) {
					if (aa.getThreeLetterCode().equals(bckgrdNbh.getCentralResidue().getResidueType())) {
						countPerAA[i]++;
					}
					i++;
				}
			}
		}
		TreeMap<String,Double> dist = new TreeMap<String, Double>();
		int i =0 ;
		for (AminoAcid aa:AminoAcid.getAllStandardAAs()) {
			if (countTotal!=0) {
				dist.put(aa.getThreeLetterCode(),(double) countPerAA[i]/countTotal);
			} else {
				dist.put(aa.getThreeLetterCode(),0.0); // we assign 0 as the frequency when there's no observations at all (countTotal = 0)
			}
			i++;
		}



		return new NbhProbDistribution(dist);
	}
	
}

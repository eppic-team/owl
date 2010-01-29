package embed;

import java.io.File;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;

import edu.uci.ics.jung.graph.util.Pair;
import graphAveraging.GraphAverager;
import graphAveraging.GraphAveragerError;
import proteinstructure.IntPairSet;
import proteinstructure.Pdb;
import proteinstructure.PdbasePdb;
import proteinstructure.RIGEnsemble;
import proteinstructure.RIGraph;
import tools.MySQLConnection;

/**
 * Class to store a IntPairSet (set of contacts) together with its score.
 * Can be sorted based on the score order. 
 * @author duarte
 *
 */
public class SetScore implements Comparable<SetScore>{
	
	public double score;
	public IntPairSet set;
	
	public SetScore(double score, IntPairSet set) {
		this.score = score;
		this.set = set;
	}
	
	public int compareTo(SetScore o) {
		if (this.score<o.score) return -1;
		if (this.score>o.score) return 1;
		return 0;
	}
	
	/*------------------------- statics  -----------------------------*/
	
	private static ArrayList<SetScore> readSetScoresFromDB(MySQLConnection conn, String db, String scoresTable, String subsetsTable, ArrayList<Integer> ids) throws SQLException{
		ArrayList<SetScore> setScores = new ArrayList<SetScore>();
		Statement st = conn.createStatement();
		ResultSet rs;
		for (int id:ids) {
			double score = -1;
			IntPairSet set = new IntPairSet();
			String sql = "SELECT error FROM "+db+"."+scoresTable+" WHERE ssid="+id;
			rs = st.executeQuery(sql);
			while(rs.next()) {
				score = rs.getDouble(1);
			}
			rs.close();
			sql = "SELECT i,j FROM "+db+"."+subsetsTable+" WHERE ssid="+id;
			rs = st.executeQuery(sql);
			while(rs.next()) {
				set.add(new Pair<Integer>(rs.getInt(1),rs.getInt(2)));
			}
			rs.close();
			setScores.add(new SetScore(score,set));
		}
		st.close();
		return setScores;
	}
	
	/**
	 * Gets a RIGraph consisting of the union of contacts of the best numBestSets subsets 
	 * read from database with edges weighted by the fraction of ocurrence in the subsets.
	 * @param db the database from which to read
	 * @param scoresTable the table containing the scores
	 * @param subsetsTable the table containing the subsets
	 * @param numBestSets the number of best subsets we want the ensemble for
	 * @param sequence
	 * @param contactType
	 * @param cutoff
	 * @return
	 * @throws SQLException
	 */
	public static RIGraph readEnsembleFromDB(String db, String scoresTable, String subsetsTable, int numBestSets, String sequence, String contactType, double cutoff) throws SQLException {
		MySQLConnection conn = new MySQLConnection();
		String sql = "SELECT ssid FROM "+db+"."+scoresTable+" ORDER BY error ASC LIMIT "+numBestSets;
		Statement st = conn.createStatement();
		ResultSet rs = st.executeQuery(sql);
		ArrayList<Integer> ids = new ArrayList<Integer>();
		while (rs.next()) {
			ids.add(rs.getInt(1));
		}
		rs.close();
		st.close();
		ArrayList<SetScore> setScores = readSetScoresFromDB(conn, db, scoresTable, subsetsTable, ids);
		conn.close();
		
 		RIGEnsemble rigs = new RIGEnsemble(contactType, cutoff);
 		double sumScore = 0;
		for (SetScore setScore:setScores) {
			sumScore+=setScore.score;
			rigs.addRIG(Distiller.createRIGraphFromIntPairSet(sequence, setScore.set, contactType, cutoff));
		}
		System.out.printf("Average error value for best sets: "+Distiller.SCORE_PRINT_FORMAT+"\n",(sumScore/(double)numBestSets));
		GraphAverager ga = null;
		try { 
			ga = new GraphAverager(rigs);
		} catch (GraphAveragerError e) {
			// this shouldn't happen
			System.err.println("Unexpected error while creating the average RIG: "+e.getMessage());
		}
		return ga.getAverageGraph();
	}

	/*-------------------------- main  -------------------------------*/
	
	/**
	 * Reads from db some subsets and makes an average ensemble out of them writing 
	 * the final weighted contact map to a file. 
	 */
	public static void main(String[] args) throws Exception {
		if (args.length<5) {
			System.err.println("Usage: SetScore <db> <scoresTable> <subsetsTable> <numBestSets> <outdir>");
			System.exit(1);
		}
		String db = args[0];
		String scoresTable = args[1];
		String subsetsTable = args[2];
		int numBestSets = Integer.parseInt(args[3]);
		String outDir = args[4]; 
		
		String pdbCode = "1sha";
		String pdbChainCode = "A";
		String contactType = "Ca";
		double cutoff = 8;
		Pdb pdb = new PdbasePdb(pdbCode);
		pdb.load(pdbChainCode);
		String sequence = pdb.getSequence();
		String fileName = pdbCode+pdbChainCode+"_ensemble_"+numBestSets+".cm";
		
		RIGraph graph = readEnsembleFromDB(db, scoresTable, subsetsTable, numBestSets, sequence, contactType, cutoff);
		graph.write_graph_to_file(new File(outDir,fileName).getAbsolutePath());
	}
}

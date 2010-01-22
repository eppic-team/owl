package deltaRank;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collection;

import edu.uci.ics.jung.graph.util.Pair;


import proteinstructure.AAinfo;
import proteinstructure.RIGNbhood;
import proteinstructure.RIGNode;
import proteinstructure.RIGraph;
import tools.MySQLConnection;

public class DeltaRank {

	private MySQLConnection conn;
	private RIGraph graph;
	private double[][] matrix;
	private double score = 0;
	private int scoringResiduesCount;
	private String[] vectors;
	
	public DeltaRank(MySQLConnection myConn, RIGraph riGraph) {
		conn= myConn;
		graph = riGraph;
		score = 0.0;
		vectors = new String[graph.getFullLength()];
		matrix = new double[graph.getFullLength()][graph.getFullLength()];
		for (int i = 1; i <= graph.getFullLength();i++) {
			for (int j = 1; j <= graph.getFullLength();j++) {
				matrix[i-1][j-1] = calculateDeltaRank(i,j);
			}
		}
		updateScore();
		updateVectors();
	}
	
	/**
	 * Calculates the delta rank for one cell of a contact map matrix
	 * @param i
	 * @param j
	 * @return delta rank, -100 if not enough data
	 */
	
	private double calculateDeltaRank(int i, int j) {
		RIGNbhood nbhoodj,nbhoodi, nbhoodjAfter, nbhoodiAfter;
		Statement stm;
		ResultSet res;
		double ret = -100;
		Collection<RIGNode> nbj, nbi;
		RIGNode nodeI, nodeJ;
		nodeJ = graph.getNodeFromSerial(j);
		nodeI = graph.getNodeFromSerial(i);

		
		if (nodeJ == null || nodeI == null) {
			return -100;
		}
		
		nbhoodj = graph.getNbhood(nodeJ);
		nbhoodi = graph.getNbhood(nodeI);
		nbj = nbhoodj.getNeighbors();
		nbj.add(nodeI);
		nbhoodjAfter = new RIGNbhood(nodeJ, nbj);
		
		nbi = nbhoodi.getNeighbors();
		nbi.add(nodeJ);
		nbhoodiAfter = new RIGNbhood(nodeI, nbi);
		String sql = "SELECT IFNULL(((SELECT LOCATE('"+AAinfo.threeletter2oneletter(nodeI.getResidueType())+"',rvector) from mw.vectors where nbstring='"+nbhoodi.getNbString()+"') -	" +
					 	"(SELECT LOCATE('"+AAinfo.threeletter2oneletter(nodeI.getResidueType())+"',rvector) from mw.vectors where nbstring='"+nbhoodiAfter.getNbString()+"')) +" +
					 	"((SELECT LOCATE('"+AAinfo.threeletter2oneletter(nodeJ.getResidueType())+"',rvector) from mw.vectors where nbstring='"+nbhoodj.getNbString()+"') -	" +
							 	"(SELECT LOCATE('"+AAinfo.threeletter2oneletter(nodeJ.getResidueType())+"',rvector) from mw.vectors where nbstring='"+nbhoodjAfter.getNbString()+"')),-100);";
		try {
			stm = conn.createStatement();
			res = stm.executeQuery(sql);
			if (res.next()) {
				ret =  (double)res.getInt(1);
			} 
			
			res.close();
			stm.close();
		} catch (SQLException e) {
			e.printStackTrace();
			e.getMessage();
		}
		return ret;
	}
	
	public double[][] getMatrix() {
		return matrix;
	}
	
	public void setGraph(RIGraph graph2) {
		graph = graph2;
	}

	public double getScore() {
		return (double)score/(double)scoringResiduesCount;
	}
	
	
	/**
	 * Updates the delta rank vectors
	 */
	
	private void updateVectors() {
		RIGNbhood nbhood;
		RIGNode node;
		Statement stm;
		ResultSet res;
		
		for (int i = 1; i <= graph.getFullLength(); i++) {
			node = graph.getNodeFromSerial(i);
			if (node == null) { continue; }
			nbhood = graph.getNbhood(node);
			String sql = "SELECT rvector from mw.vectors where nbstring='"+nbhood.getNbString()+"';";
			try {
				stm = conn.createStatement();
				res = stm.executeQuery(sql);
				if (res.next()) {
					vectors[i-1] = res.getString(1);
				} else {
					vectors[i-1] = "";
				}
				
				res.close();
				stm.close();
			} catch (SQLException e) {
				e.printStackTrace();
				e.getMessage();
			} 
		}
	}
	
	/**
	 * The delta rank score is defined as the net sum of rank changes compared with the background distribution 
	 */
	
	private void updateScore() {
		RIGNbhood nbhood;
		RIGNode node;
		Statement stm;
		ResultSet res;
		int ret;
		score = 0;
		scoringResiduesCount = 0;
		
		for (int i = 1; i <= graph.getFullLength(); i++) {
			node = graph.getNodeFromSerial(i);
			if (node == null) { continue; }
			nbhood = graph.getNbhood(node);
			String sql = "SELECT IFNULL(((SELECT LOCATE('"+AAinfo.threeletter2oneletter(node.getResidueType())+"',rvector) from mw.vectors where nbstring='x') -	" +
		 	"(SELECT LOCATE('"+AAinfo.threeletter2oneletter(node.getResidueType())+"',rvector) from mw.vectors where nbstring='"+nbhood.getNbString()+"')),-100);";
			try {
				stm = conn.createStatement();
				res = stm.executeQuery(sql);
				if (res.next()) {
					ret =  res.getInt(1);
					if (ret > -50) {
						scoringResiduesCount++;
						score += ret;
					}
				} 
				
				res.close();
				stm.close();
			} catch (SQLException e) {
				e.printStackTrace();
				e.getMessage();
			} 
		}
	}
	
	public void updateMap(Pair<Integer> cont) {
		for (int i = 1; i <= graph.getFullLength();i++) {
			matrix[i-1][cont.getSecond()-1] = calculateDeltaRank(i,cont.getSecond());
		}
		for (int j = 1; j <= graph.getFullLength();j++) {
			matrix[cont.getFirst()-1][j-1] = calculateDeltaRank(cont.getFirst(),j);				
		}
		updateScore();
		updateVectors();
	}

	public String[] getVectors() {
		return vectors;
	}
}

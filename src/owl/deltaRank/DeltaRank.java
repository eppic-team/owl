package owl.deltaRank;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collection;

import owl.core.structure.AAinfo;
import owl.core.structure.graphs.RIGNbhood;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.MySQLConnection;

import edu.uci.ics.jung.graph.util.Pair;



public class DeltaRank {

	private MySQLConnection conn;
	private String db;
	private RIGraph graph;
	private double[][] matrix;
	private double score = 0;
	private int scoringResiduesCount;
	private String[] vectors;
	private double[][] probabilities;
	
	public DeltaRank(MySQLConnection myConn, RIGraph riGraph, String db) {
		conn= myConn;
		this.db = db;
		graph = riGraph;
		score = 0.0;
		vectors = new String[graph.getFullLength()];
		matrix = new double[graph.getFullLength()][graph.getFullLength()];
		probabilities = new double[graph.getFullLength()][20];
		for (int i = 1; i <= graph.getFullLength();i++) {
			for (int j = 1; j <= graph.getFullLength();j++) {
				matrix[i-1][j-1] = calculateDeltaRank(i,j);
				
			}
			for (int j=1; j <=20; j++) {
				probabilities[i-1][j-1] = 0.05;
			}
		}
		updateScore();
		updateVectors();
		updateProbabilities();
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
		String sql = "SELECT IFNULL(((SELECT LOCATE('"+AAinfo.threeletter2oneletter(nodeI.getResidueType())+"',rvector) from "+db+".vectors where nbstring='"+nbhoodi.getNbString()+"') -	" +
					 	"(SELECT LOCATE('"+AAinfo.threeletter2oneletter(nodeI.getResidueType())+"',rvector) from "+db+".vectors where nbstring='"+nbhoodiAfter.getNbString()+"')) +" +
					 	"((SELECT LOCATE('"+AAinfo.threeletter2oneletter(nodeJ.getResidueType())+"',rvector) from "+db+".vectors where nbstring='"+nbhoodj.getNbString()+"') -	" +
							 	"(SELECT LOCATE('"+AAinfo.threeletter2oneletter(nodeJ.getResidueType())+"',rvector) from "+db+".vectors where nbstring='"+nbhoodjAfter.getNbString()+"')),-100);";
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
		String nbstring;
		
		for (int i = 1; i <= graph.getFullLength(); i++) {
			node = graph.getNodeFromSerial(i);
			if (node == null) { 
				nbstring = "x";
			} else {
				nbhood = graph.getNbhood(node);
				nbstring = nbhood.getNbString();
			}
			String sql = "SELECT rvector from mw.vectors where nbstring='"+nbstring+"';";
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
	
	/* Updates the probabilities for each position in the sequence */
	
	private void updateProbabilities() {
		RIGNbhood nbhood;
		RIGNode node;
		Statement stm;
		ResultSet res;
		for (int i = 0; i < graph.getFullLength(); i++) {
			try {
				node = graph.getNodeFromSerial(i+1);
				if (node != null) {
				
					nbhood = graph.getNbhood(node);
					String sql = "SELECT (L+A+G+V+E+D+S+K+T+I+R+P+N+F+Q+Y+H+M+W+C),L,A,G,V,E,D,S,K,T,I,R,P,N,F,Q,Y,H,M,W,C from mw.avectors where str='"+nbhood.getNbString()+"';";
					stm = conn.createStatement();
					res = stm.executeQuery(sql);
					if (res.next()) {
						double sum = res.getInt(1);
						probabilities[i][0] = res.getDouble(2)/sum;
						probabilities[i][1] = res.getDouble(3)/sum;
						probabilities[i][2] = res.getDouble(4)/sum;
						probabilities[i][3] = res.getDouble(5)/sum;
						probabilities[i][4] = res.getDouble(6)/sum;
						probabilities[i][5] = res.getDouble(7)/sum;
						probabilities[i][6] = res.getDouble(8)/sum;
						probabilities[i][7] = res.getDouble(9)/sum;
						probabilities[i][8] = res.getDouble(10)/sum;
						probabilities[i][9] = res.getDouble(11)/sum;
						probabilities[i][10] = res.getDouble(12)/sum;
						probabilities[i][11] = res.getDouble(13)/sum;
						probabilities[i][12] = res.getDouble(14)/sum;
						probabilities[i][13] = res.getDouble(15)/sum;
						probabilities[i][14] = res.getDouble(16)/sum;
						probabilities[i][15] = res.getDouble(17)/sum;
						probabilities[i][16] = res.getDouble(18)/sum;
						probabilities[i][17] = res.getDouble(19)/sum;
						probabilities[i][18] = res.getDouble(20)/sum;
						probabilities[i][19] = res.getDouble(21)/sum;
						
					}
				}
				
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
		updateProbabilities();
	}

	public String[] getVectors() {
		return vectors;
	}

	public double[][] getProbabilities() {
		return probabilities;
	}
}

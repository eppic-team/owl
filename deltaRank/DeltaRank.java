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
	private int scoringCellsCount = 0;
	
	public DeltaRank(MySQLConnection myConn, RIGraph riGraph) {
		conn= myConn;
		graph = riGraph;
		score = 0.0;
		matrix = new double[graph.getFullLength()][graph.getFullLength()];
		for (int i = 1; i <= graph.getFullLength();i++) {
			for (int j = 1; j <= graph.getFullLength();j++) {
				matrix[i-1][j-1] = calculateDeltaRank(i,j);
				if (matrix[i-1][j-1] > -80) {
					scoringCellsCount += 1;
					score += matrix[i-1][j-1];
				}
			}
		}
	}
	
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
		return (score/(double)scoringCellsCount)*100;
	}
	
	public void updateMap(Pair<Integer> cont) {
		for (int i = 1; i <= graph.getFullLength();i++) {
			if (matrix[i-1][cont.getSecond()-1] > -90) {
				scoringCellsCount--;
				score -= matrix[i-1][cont.getSecond()-1];
			}
			matrix[i-1][cont.getSecond()-1] = calculateDeltaRank(i,cont.getSecond());
			if (matrix[i-1][cont.getSecond()-1] > -90) {
				score += matrix[i-1][cont.getSecond()-1];
				scoringCellsCount++;
			}
		}
		for (int j = 1; j <= graph.getFullLength();j++) {
			if (matrix[cont.getFirst()-1][j-1] > -90) {
				score -= matrix[cont.getFirst()-1][j-1];
				scoringCellsCount--;
			}
			matrix[cont.getFirst()-1][j-1] = calculateDeltaRank(cont.getFirst(),j);				
			if (matrix[cont.getFirst()-1][j-1] > -90) {
				scoringCellsCount++;
				score += matrix[cont.getFirst()-1][j-1];
			}
		}
	}
}

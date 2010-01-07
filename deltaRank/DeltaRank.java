package deltaRank;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Collection;


import proteinstructure.AAinfo;
import proteinstructure.RIGNbhood;
import proteinstructure.RIGNode;
import proteinstructure.RIGraph;
import tools.MySQLConnection;

public class DeltaRank {

	public static double[][] computeMatrix(MySQLConnection conn, RIGraph graph) {
	
		double[][] matrix = new double[graph.getFullLength()][graph.getFullLength()];
		RIGNbhood nbhoodj,nbhoodi, nbhoodjAfter, nbhoodiAfter;
		Statement stm;
		ResultSet res;
		Collection<RIGNode> nbj, nbi;
		for (int i = 1; i <= graph.getFullLength();i++) {
			nbhoodi = graph.getNbhood(graph.getNodeFromSerial(i));

			for (int j = 1; j <= graph.getFullLength();j++) {
				nbhoodj = graph.getNbhood(graph.getNodeFromSerial(j));
				
				nbj = nbhoodj.getNeighbors();
				nbj.add(graph.getNodeFromSerial(i));
				nbhoodjAfter = new RIGNbhood(graph.getNodeFromSerial(j), nbj);
				
				nbi = nbhoodi.getNeighbors();
				nbi.add(graph.getNodeFromSerial(j));
				nbhoodiAfter = new RIGNbhood(graph.getNodeFromSerial(i), nbi);
				
				String sql = "SELECT IFNULL(((SELECT LOCATE('"+AAinfo.threeletter2oneletter(graph.getNodeFromSerial(i).getResidueType())+"',rvector) from mw.vectors where nbstring='"+nbhoodi.getNbString()+"') -	" +
							 	"(SELECT LOCATE('"+AAinfo.threeletter2oneletter(graph.getNodeFromSerial(i).getResidueType())+"',rvector) from mw.vectors where nbstring='"+nbhoodiAfter.getNbString()+"')) +" +
							 	"((SELECT LOCATE('"+AAinfo.threeletter2oneletter(graph.getNodeFromSerial(j).getResidueType())+"',rvector) from mw.vectors where nbstring='"+nbhoodj.getNbString()+"') -	" +
									 	"(SELECT LOCATE('"+AAinfo.threeletter2oneletter(graph.getNodeFromSerial(j).getResidueType())+"',rvector) from mw.vectors where nbstring='"+nbhoodjAfter.getNbString()+"')),-100);";
				
				System.out.println(sql);
				try {
					stm = conn.createStatement();
					res = stm.executeQuery(sql);
					if (res.next()) {
						matrix[i-1][j-1] = (double)res.getInt(1);
						System.out.println(res.getInt(1));
					} else {
						matrix[i-1][j-1] = -100;
					}
					//matrix[i-1][j-1] = 0.2;
					
				} catch (SQLException e) {
					e.printStackTrace();
					e.getMessage();
				}

				

			}
		}
		return matrix;
	}
	
	public double[][] getMatrix() {
		
		return null;
	}

	
}

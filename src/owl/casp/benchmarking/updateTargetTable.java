package owl.casp.benchmarking;

import java.io.File;
import java.sql.SQLException;

import owl.core.structure.PdbChain;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.MySQLConnection;

/**
 * Quick-and-dirty script to update the table scores_ranks_targets with the edge count for the native structures
 * @author stehr
 */
public class updateTargetTable {
	
	public static void main(String[] args) {
		
		File answerDir = new File("/project/StruPPi/CASP8/answers/");
		String dbTable = "casp8.scores_ranks_targets";
		
		MySQLConnection conn = null;
		try {
			conn = new MySQLConnection();
		} catch (SQLException e) {
			System.err.println("Error. Could not connect to database: " + e.getMessage() + ". Exiting.");
			System.exit(1);
		}
		
		String[] fileNames = answerDir.list();
		
		for(String s: fileNames) {
			File f = new File(answerDir, s);
			String n = f.getName();
			//System.out.println(n);
			if(n.startsWith("T0") && n.endsWith(".pdb") && n.length()==9) {
	
				String target = n.substring(0,5);
				System.out.println(target);
				
				PdbChain pdb = PdbChain.readStructureOrNull(f.getPath());
				if(pdb == null) {
					System.err.println("Error. Could not read file " + s);
				} else {
					RIGraph rig = pdb.getRIGraph("Cb", 8.0);
					int cts = 0, cts12=0, cts24=0;
					for(RIGEdge e:rig.getEdges()) {
						int seqSep = rig.getContactRange(e);
						if(seqSep >= 24) cts24++;
						if(seqSep >= 12) cts12++;
						if(seqSep >= 1) cts++;
					}
					
					// write to database
					String sql = String.format("UPDATE %s SET cts=%d, cts12=%d, cts24=%d WHERE target='%s'", dbTable, cts, cts12, cts24, target);
					try {
						System.out.println(sql);
						conn.executeSql(sql);
					} catch (SQLException e1) {
						System.err.println("Error writing to database: " + e1.getMessage());
						//System.err.println(sql);
					}
				}
			}
		}	
		System.out.println("done.");
	}
}

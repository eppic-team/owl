import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;

import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbasePdb;
import owl.core.util.MySQLConnection;



public class findDisorderMatchingOrder {

	private static final int TAIL_THRESHOLD = 15;
	
	/**
	 * @param args
	 * @throws SQLException 
	 * @throws IOException 
	 */
	public static void main(String[] args) throws SQLException, IOException {
		
		File nrFile = new File("/home/duarte_j/tmp/nr.clusters");
		MySQLConnection conn = new MySQLConnection();
		
		BufferedReader br = new BufferedReader(new FileReader(nrFile));
		String line;
		int count = 0;
		while ((line=br.readLine())!=null) {
			String[] tokens = line.split(" ");
			ArrayList<Pdb> pdbs = new ArrayList<Pdb>();
			for (String token:tokens) {
				String pdbCode = token.substring(0, 4);
				String pdbChainCode = token.substring(4, 5);
				Pdb pdb = null;
				try {
					pdb = new PdbasePdb(pdbCode,"pdbase",conn);
					pdb.load(pdbChainCode);
				} catch (PdbLoadError e){
					System.err.println(e.getMessage());
					continue;
				} catch (PdbCodeNotFoundError e) {
					System.err.println(e.getMessage());
					continue;
				}
				pdbs.add(pdb);
			}
			int i = 0;
			int fullLength = 0;
			boolean hasDisordered = false;
			for (Pdb pdb: pdbs) {
				if (i==0) {
					fullLength = pdb.getFullLength();	
				} 
				if (pdb.getObsLength()<fullLength-TAIL_THRESHOLD) {
					hasDisordered = true;
				}
			}
			if (hasDisordered) {
				count++;
			}
		}
		br.close();

		System.out.println(count + " clusters contain disordered tails");
	}

}

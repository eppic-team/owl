package scripts;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import owl.core.connections.pisa.PisaAsmSet;
import owl.core.connections.pisa.PisaAsmSetList;
import owl.core.connections.pisa.PisaAssembly;
import owl.core.connections.pisa.PisaConnection;

/**
 * A script that given a file with a list of PDB codes, downloads PISA xml assembly 
 * predictions for them and prints out a list of the oligomeric predictions
 * 
 * @author duarte_j
 *
 */
public class getPisaAssemblyPreds {

	private static List<String> readListFile(File file) throws IOException {
		List<String> pdbCodes = new ArrayList<String>();
		BufferedReader flist = new BufferedReader(new FileReader(file));
		String line;
		while ((line = flist.readLine() ) != null ) {
			if (line.startsWith("#")) continue;
			if (line.isEmpty()) break;
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			pdbCodes.add(pdbCode);
		}
		flist.close();
		return pdbCodes;
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {

		File file = new File(args[0]);
		
		PisaConnection pc = new PisaConnection();
		
		List<String> pdbCodes = readListFile(file);
		
		Map<String, PisaAsmSetList> asms = pc.getAssembliesDescription(pdbCodes);
		
		for (String pdbCode:pdbCodes) {
			PisaAsmSetList pal = asms.get(pdbCode);
			if (pal==null) continue; // pisa didn't return the pdb code, probably obsoleted
			if (pal.size()==0) {
				System.out.printf("%4s\t%2d\n",pdbCode,1);
			} else {
				PisaAsmSet pas = pal.get(0);
				int i = 0;
				for (PisaAssembly pa:pas) {
					if (!pa.isMacromolecular()) continue;
					if (i==0) System.out.print(pdbCode);
					else System.out.printf("%4s","");
					
					int mmsizePred = 1;
					if (pa.getDissEnergy()>0) mmsizePred = pa.getMmsize();

					System.out.printf("\t%2d\t%2d\t%5.1f\t%20s\t%s\t%s\n",
							mmsizePred,pa.getMmsize(),
							pa.getDissEnergy(),
							pa.getFormula(),
							pa.getInterfaceIdsString(),
							pa.getScore());
					i++;
				}
			}
		}
		
		
	}

}

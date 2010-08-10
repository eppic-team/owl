package owl.core.structure;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.TreeMap;

/**
 * A class containing static methods to parse the symop.lib file from the 
 * CCP4 package. The file contains the transformations belonging to all 
 * protein crystallography space groups.
 * 
 * @author duarte_j
 *
 */
public class SymoplibParser {
	
	// the symop library from the CCP4 package
	private static final String SYMOPFILE = "/owl/core/structure/symop.lib";
	
	private static final InputStream symoplibIS = SymoplibParser.class.getResourceAsStream(SYMOPFILE);
	private static final TreeMap<Integer, SpaceGroup> sgs = parseSymopLib();
	
	/**
	 * Gets the space group for the given standard identifier.
	 * See for example http://en.wikipedia.org/wiki/Space_group
	 * @param id
	 * @return
	 */
	public static SpaceGroup getSpaceGroup(int id) {
		return sgs.get(id);
	}
	
	public static TreeMap<Integer,SpaceGroup> getAllSpaceGroups() {
		return sgs;
	}
	
	private static TreeMap<Integer,SpaceGroup> parseSymopLib() {
		TreeMap<Integer, SpaceGroup> map = new TreeMap<Integer, SpaceGroup>();
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(symoplibIS));
			String line;
			SpaceGroup currentSG = null;
			while ((line=br.readLine())!=null) {
				if (!line.startsWith(" ")) {
					if (currentSG!=null) {
						map.put(currentSG.getId(),currentSG);
					}
					String[] tokens = line.split("\\s+");
					int id = Integer.parseInt(tokens[0]);
					String shortSymbol = tokens[4];
					currentSG = new SpaceGroup(id, shortSymbol);
				} else {
					currentSG.addTransformation(line.trim());
				}
			}
			br.close();
			// and we add the last SG
			map.put(currentSG.getId(), currentSG);
		} catch (IOException e) {
			System.err.println("Fatal error! Can't read resource file "+SYMOPFILE+". Error: "+e.getMessage()+". Exiting.");
			System.exit(1);
		}
		return map;
	}
	
	

}

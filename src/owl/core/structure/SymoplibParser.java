package owl.core.structure;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
	
	private static final Pattern namePat = Pattern.compile(".*\\s([A-Z]+)(?:\\s'.+')?\\s+'(.+)'.*");
	
	private static final InputStream symoplibIS = SymoplibParser.class.getResourceAsStream(SYMOPFILE);
	private static final TreeMap<Integer, SpaceGroup> sgs = parseSymopLib();
	
	private static HashMap<String, SpaceGroup> name2sgs; // map for lookups based on short names
	
	/**
	 * Gets the space group for the given standard identifier.
	 * See for example http://en.wikipedia.org/wiki/Space_group
	 * @param id
	 * @return
	 */
	public static SpaceGroup getSpaceGroup(int id) {
		return sgs.get(id);
	}
	
	/**
	 * Gets the space group for the given international short name, using
	 * the PDB format, e.g. 'P 21 21 21' or 'C 1 c 1'
	 * @param shortName
	 * @return
	 */
	public static SpaceGroup getSpaceGroup(String shortName) {
		return name2sgs.get(shortName);
	}
	
	public static TreeMap<Integer,SpaceGroup> getAllSpaceGroups() {
		return sgs;
	}
	
	private static TreeMap<Integer,SpaceGroup> parseSymopLib() {
		TreeMap<Integer, SpaceGroup> map = new TreeMap<Integer, SpaceGroup>();
		name2sgs = new HashMap<String, SpaceGroup>();
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(symoplibIS));
			String line;
			SpaceGroup currentSG = null;
			while ((line=br.readLine())!=null) {
				if (!line.startsWith(" ")) {
					if (currentSG!=null) {
						map.put(currentSG.getId(),currentSG);
						name2sgs.put(currentSG.getShortSymbol(), currentSG);
					}
					
					int id = Integer.parseInt(line.substring(0, line.indexOf(' ')));
					Matcher m = namePat.matcher(line);
					String shortSymbol = null;
					String brav = null;
					if (m.matches()) {
						brav = m.group(1);
						shortSymbol = m.group(2);
					}
					currentSG = new SpaceGroup(id, shortSymbol,SpaceGroup.BravaisLattice.getByName(brav));
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

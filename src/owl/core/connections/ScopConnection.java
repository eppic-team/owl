package owl.core.connections;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.features.Scop;
import owl.core.structure.features.ScopRegion;

/**
 * 
 * Connection class to get Scop annotation data from the SCOP website or
 * from a local copy of the SCOP data.
 *
 */
public class ScopConnection {
	
	private static final String SCOP_URL_PREFIX = "http://scop.mrc-lmb.cam.ac.uk/scop/parse/";
	private static final String SCOP_DIR = "/project/StruPPi/Databases/SCOP";


	/**
	 * Parses SCOP annotation for the PDB code of the given PdbChain object, setting its
	 * Scop member with SCOP annotation.
	 * @param pdb the PdbChain object for which we want Scop annotation
	 * @param version the SCOP version that we want to parse
	 * @param online if true SCOP annotation will be taken from web, if false 
	 * from local file
	 * @throws IOException 
	 */
	public static void parseScop(PdbChain pdb, String version, boolean online) throws IOException {
		String pdbCode = pdb.getPdbCode();
		// if this is not a pdb entry with a pdb code there's no SCOP id to get
		if (pdbCode==null || pdbCode.equals(PdbAsymUnit.NO_PDB_CODE)) return; 
		
		String pdbChainCode = pdb.getPdbChainCode();
		
		Scop scop = new Scop();	
		ScopRegion sr = null;
		String startPdbResSer = "", endPdbResSer = "";
		BufferedReader in;
		String inputLine;
	
		if (online) {
			URL scop_cla = new URL(SCOP_URL_PREFIX+"dir.cla.scop.txt_"+version);
			URLConnection sc = scop_cla.openConnection();
			in = new BufferedReader(new InputStreamReader(sc.getInputStream()));
		} else {
			File scop_cla = new File(SCOP_DIR,"dir.cla.scop.txt_"+version);
			in = new BufferedReader(new FileReader(scop_cla));
		}

		while ((inputLine = in.readLine()) != null) 
			if (inputLine.startsWith(pdbCode,1)) {
				String[] fields = inputLine.split("\\s+");
				String[] regions = fields[2].split(",");
				for (int j=0; j < regions.length; j++) {
					Pattern p = Pattern.compile("^(-)|([a-zA-Z\\d]):(-?\\d+[a-zA-Z]*)-(-?\\d+[a-zA-Z]*)|(-?\\d+[a-zA-Z]*)-(-?\\d+[a-zA-Z]*)|([a-zA-Z\\d]):");
					Matcher m = p.matcher(regions[j]);
					if (m.find()) {
						if (((pdbChainCode.equals(PdbAsymUnit.NULL_CHAIN_CODE) && ((m.group(1) != null && m.group(1).equals("-")) || m.group(5) != null))) || 
								(m.group(2) != null && m.group(2).equals(pdbChainCode)) || 
								(m.group(7) != null && m.group(7).equals(pdbChainCode)) ||
								(m.group(2) != null && m.group(2).equals("A") && pdbChainCode.equals(PdbAsymUnit.NULL_CHAIN_CODE)) ||
								(m.group(7) != null && m.group(7).equals("A") && pdbChainCode.equals(PdbAsymUnit.NULL_CHAIN_CODE))) {
							if (m.group(3) != null) {
								startPdbResSer = m.group(3);
								endPdbResSer = m.group(4);
							} else if (m.group(5) != null) {
								startPdbResSer = m.group(5);
								endPdbResSer = m.group(6);								
							} else {
								startPdbResSer = pdb.getPdbResSerFromResSer(pdb.getMinObsResSerial());
								endPdbResSer = pdb.getPdbResSerFromResSer(pdb.getMaxObsResSerial());
							}
							sr = new ScopRegion(fields[0], fields[3], Integer.parseInt(fields[4]), j, regions.length, startPdbResSer, endPdbResSer, pdb.getResSerFromPdbResSer(startPdbResSer), pdb.getResSerFromPdbResSer(endPdbResSer));
							scop.add(sr);
						}
					}
				}
				//break; // we can't break: pdbCodes are not necessarily ordered in the scop annotation file: we need to parse to the end of the file
			}
		in.close();

		scop.setVersion(version);
		
		pdb.setScop(scop);
	}	

}

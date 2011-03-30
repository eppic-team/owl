package owl.core.connections;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.zip.GZIPInputStream;

import owl.core.structure.AminoAcid;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.Residue;
import owl.core.structure.features.CatalSiteSet;
import owl.core.structure.features.CatalyticSite;

/**
 * Connection class to get catalytic site data from the EBI CSA web resource 
 * or from local copies of the CSA data file. 
 * 
 *
 */
public class CSAConnection {

	private static final String CSA_DIR = "/project/StruPPi/Databases/CSA";
	private static final String CSA_URL_PREFIX = "http://www.ebi.ac.uk/thornton-srv/databases/CSA/archive/";

	/**
	 * Queries the Catalytic Site Atlas (Porter et al. NAR 32: D129-D133, 2004) for the 
	 * PDB code of the given PdbChain object and sets its member variable <code>catalSiteSet</code> with the results.
	 * If no entry is found,<code>catalSiteSet</code> is set to null.
	 * @param pdb
	 * @param version the CSA version to use (e.g. {@link CatalSiteSet.LATEST_VERSION})
	 * @param online whether to access online version of CSA or use a local copy
	 * @return number of errors (mismatching residues between this structure and the CSA annotation)
	 * @throws IOException if local file or online database could not be read
	 */
	public static int parseCSA(PdbChain pdb, String version, boolean online) throws IOException {
		BufferedReader in;
		String inputLine;

		if (online) {
			URL csa = new URL(CSA_URL_PREFIX+"CSA_"+version.replaceAll("\\.","_")+".dat.gz");
			URLConnection conn = csa.openConnection();
			InputStream inURL = conn.getInputStream();
			OutputStream out = new BufferedOutputStream(new FileOutputStream("CSA_"+version.replaceAll("\\.","_")+".dat.gz"));
			byte[] buffer = new byte[1024];
			int numRead;
			long numWritten = 0;
			while ((numRead = inURL.read(buffer)) != -1) {
				out.write(buffer, 0, numRead);
				numWritten += numRead;
			}
			inURL.close();
			out.close();
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream("CSA_"+version.replaceAll(".","_")+".dat.gz"))));
		} else {
			File csaFile = new File(CSA_DIR,"CSA_"+version.replaceAll("\\.","_")+".dat");
			in = new BufferedReader(new FileReader(csaFile));			
		}

		int csaMistakes = 0;
		CatalSiteSet catalSiteSet = new CatalSiteSet();
		String curPdbCode, curPdbChainCode, curPdbResSerial, curResType;
		String prevPdbCode = "";
		int curSiteId;
		int prevSiteId = -1;
		CatalyticSite cs = null;
		while ((inputLine = in.readLine()) != null) { 
			String[] fields = inputLine.split(",");
			curPdbCode = fields[0]; 
			curPdbChainCode = (fields[3].equals(""))?PdbAsymUnit.NULL_CHAIN_CODE:fields[3];
			if (curPdbCode.equals(pdb.getPdbCode())) {
				if (curPdbChainCode.equals(pdb.getPdbChainCode())) {
					curPdbResSerial = fields[4];
					curResType = fields[2];
					curSiteId = Integer.valueOf(fields[1]);
					// only if the pdb residue type is a standard AA
					if (AminoAcid.isStandardAA(curResType)) {
						// only if the pdb residue type agrees with our residue type
						Residue residue = pdb.getResidue(pdb.getResSerFromPdbResSer(curPdbResSerial));
						if (residue.getAaType().getThreeLetterCode().equals(curResType)) {
							// each time the site changes except for the first site of a chain,
							// add the site to the set
							if ((curSiteId != prevSiteId) & (prevSiteId != -1)) {
								catalSiteSet.add(cs);
							}
							// each time the site changes or it is the first site of a chain,
							// create a new site
							if ((curSiteId != prevSiteId) | (prevSiteId == -1)) {
								String littEntryPdbCode = fields[7].substring(0,4);
								String littEntryPdbChainCode = fields[7].substring(4).equals("")?PdbAsymUnit.NULL_CHAIN_CODE:fields[7].substring(4);
								cs = new CatalyticSite(curSiteId, fields[6], littEntryPdbCode, littEntryPdbChainCode); 
							}
							// add the res to the site
							cs.addRes(pdb.getResSerFromPdbResSer(curPdbResSerial), fields[5]);
							prevSiteId = curSiteId;
						} else {
							csaMistakes++;
						}
					}
				}
			} else if (prevPdbCode.equals(pdb.getPdbCode())) {
				if (prevSiteId != -1) {
					catalSiteSet.add(cs);
				}
				break;
			}
			prevPdbCode = curPdbCode;
		}

		in.close();
		if (online) {
			File csaFile = new File("CSA_"+version.replaceAll(".","_")+".dat.gz");
			csaFile.delete();
		}

		if ((csaMistakes > 0) | (prevSiteId == -1)) {
			catalSiteSet = null;
		}
		
		pdb.setCatalSiteSet(catalSiteSet);

		return csaMistakes;
	}
}

package owl.core.connections.pisa;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.xml.sax.SAXException;

import owl.core.structure.ChainInterface;
import owl.core.structure.ChainInterfaceList;
import owl.core.structure.PdbAsymUnit;

/**
 * Connection class to get interface data from the PISA server.
 * 
 * See http://www.ebi.ac.uk/msd-srv/prot_int/pi_download.html
 *  
 * @author duarte_j
 *
 */
public class PisaConnection {

	public static final String PISA_INTERFACES_URL = "http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/interfaces.pisa?"; //pdbcodelist
	public static final String PISA_ASSEMBLIES_URL = "http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/multimers.pisa?";  //pdbcodelist
	public static final String PISA_PDB_ASSEMBLIES_URL = "http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/multimer.pdb?";//pdbcode:n,m
	
	private static final String LOCAL_CIF_DIR = "/nfs/data/dbs/pdb/data/structures/all/mmCIF";
	
	private static final int MAX_ENTRIES_PER_REQUEST = 10; // pisa doesn't specify a limit but recommends 20-50 per request
														   // first we used 50 and seemed to work, lately (14Oct 2011) we had to change to 10 
														   // because the server would stall at downloading
	
	private String interfacesUrl;
	private String assembliesUrl;
	//private String pdbAssembliesUrl;
	
	public PisaConnection()	{
		this.interfacesUrl = PISA_INTERFACES_URL;
		this.assembliesUrl = PISA_ASSEMBLIES_URL;
		//this.pdbAssembliesUrl = PISA_PDB_ASSEMBLIES_URL;
	}
	
	/**
	 * Retrieves the XML PISA interface description from the PISA web server dividing the 
	 * query into chunks of {@value #MAX_ENTRIES_PER_REQUEST}, parses it and 
	 * returns the result as a map of pdb codes to lists of PISA interfaces 
	 * @param pdbCodesList pdb codes list for which we want to retrieve pisa interfaces
	 * @return
	 * @throws IOException
	 * @throws SAXException
	 */
	public Map<String,PisaInterfaceList> getInterfacesDescription(List<String> pdbCodesList) throws IOException, SAXException {
		Map<String,PisaInterfaceList> allPisaInterfaces = new HashMap<String,PisaInterfaceList>();
		// we do batches of MAX_ENTRIES_PER_REQUEST
		for (int i=0;i<pdbCodesList.size();i+=MAX_ENTRIES_PER_REQUEST) {
			String commaSepList = "";
			for (int c=i;c<i+MAX_ENTRIES_PER_REQUEST && c<pdbCodesList.size();c++) {
				if (c!=i) commaSepList+=",";
				commaSepList+=pdbCodesList.get(c);
			}
			allPisaInterfaces.putAll(getInterfacesDescription(commaSepList));
		}
		return allPisaInterfaces;
	}
	
	/**
	 * Retrieves the XML PISA interface description from the PISA web server, parses it and 
	 * returns the result as a map of pdb codes to lists of PISA interfaces 
	 * @param commaSepList
	 * @return
	 * @throws IOException
	 * @throws SAXException
	 */
	private Map<String,PisaInterfaceList> getInterfacesDescription(String commaSepList) throws IOException, SAXException {
		URL interfacesURL = new URL(interfacesUrl+commaSepList);
		URLConnection conn = interfacesURL.openConnection();
		
		PisaInterfaceXMLParser pxmlParser = new PisaInterfaceXMLParser(conn.getInputStream());
		Map<String,PisaInterfaceList> map = pxmlParser.getAllInterfaces();
		if (map.isEmpty()) throw new IOException("The PISA server returned no interface data for URL "+interfacesURL.toString());
		return map;
	}
	
	/**
	 * Retrieves the XML PISA assembly description from the PISA web server dividing the 
	 * query into chunks of {@value #MAX_ENTRIES_PER_REQUEST}, parses it and 
	 * returns the result as a map of pdb codes to lists of PISA assemblies 
	 * @param pdbCodesList pdb codes list for which we want to retrieve pisa assemblies
	 * @return
	 * @throws IOException
	 * @throws SAXException
	 */
	public Map<String,PisaAsmSetList> getAssembliesDescription(List<String> pdbCodesList) throws IOException, SAXException {
		Map<String,PisaAsmSetList> allPisaAsmSets = new HashMap<String,PisaAsmSetList>();
		// we do batches of MAX_ENTRIES_PER_REQUEST
		for (int i=0;i<pdbCodesList.size();i+=MAX_ENTRIES_PER_REQUEST) {
			String commaSepList = "";
			for (int c=i;c<i+MAX_ENTRIES_PER_REQUEST && c<pdbCodesList.size();c++) {
				if (c!=i) commaSepList+=",";
				commaSepList+=pdbCodesList.get(c);
			}
			allPisaAsmSets.putAll(getAssembliesDescription(commaSepList));
		}
		return allPisaAsmSets;
	}
	
	/**
	 * Retrieves the XML PISA assembly description from the PISA web server, parses it and 
	 * returns the result as a map of pdb codes to lists of PISA assemblies 
	 * @param commaSepList
	 * @return
	 * @throws IOException
	 * @throws SAXException
	 */
	private Map<String,PisaAsmSetList> getAssembliesDescription(String commaSepList) throws IOException, SAXException {
		URL assembliesURL = new URL(assembliesUrl+commaSepList);
		URLConnection conn = assembliesURL.openConnection();
		
		PisaAssembliesXMLParser pxmlParser = new PisaAssembliesXMLParser(conn.getInputStream());
		Map<String,PisaAsmSetList> map = pxmlParser.getAllAssemblies();
		if (map.isEmpty()) throw new IOException("The PISA server returned no assembly data for URL "+assembliesURL.toString());
		return map;
	}
	
	public static void main(String[] args) throws Exception {
		PisaConnection pc = new PisaConnection();
		String pdbCode1 = "1aor";
		String pdbCode2 = "1bxy";
		
		File cifFile1 = new File(System.getProperty("java.io.tmpdir"),"pisa_conn"+"_"+pdbCode1+".cif");
		File cifFile2 = new File(System.getProperty("java.io.tmpdir"),"pisa_conn"+"_"+pdbCode2+".cif");
		PdbAsymUnit.grabCifFile(LOCAL_CIF_DIR, null, pdbCode1, cifFile1, false);
		PdbAsymUnit.grabCifFile(LOCAL_CIF_DIR, null, pdbCode2, cifFile2, false);
		PdbAsymUnit pdb1 = new PdbAsymUnit(cifFile1);
		PdbAsymUnit pdb2 = new PdbAsymUnit(cifFile2);

		List<String> pdbCodesList = new ArrayList<String>();
		pdbCodesList.add(pdbCode1);
		pdbCodesList.add(pdbCode2);
		
		Map<String,PisaInterfaceList> allPisaInterfaces = pc.getInterfacesDescription(pdbCodesList);
		Map<String,ChainInterfaceList> all = new HashMap<String,ChainInterfaceList>();
		all.put(pdbCode1, allPisaInterfaces.get(pdbCode1).convertToChainInterfaceList(pdb1));
		all.put(pdbCode2, allPisaInterfaces.get(pdbCode2).convertToChainInterfaceList(pdb2));
		

		for (String pdbCode:all.keySet()) {
			System.out.println("#####");
			System.out.println("# "+pdbCode);
			System.out.println("#####");
			ChainInterfaceList interfaces = all.get(pdbCode);
			for (ChainInterface interf:interfaces) {
				
				interf.printTabular(System.out, true);
			}
		}
		
		Map<String,PisaAsmSetList> allPisaAsmSets = pc.getAssembliesDescription(pdbCodesList);
		for (String pdbCode:allPisaAsmSets.keySet()) {
			System.out.println(pdbCode);
			for (PisaAsmSet asmSet:allPisaAsmSets.get(pdbCode)) {
				for (PisaAssembly ass:asmSet) {
					System.out.println(pdbCode+"\t"+ass.getId()+"\t"+ass.getMmsize()+"\t"+ass.getFormula()+"\t"+String.format("%6.2f",ass.getDissEnergy())+"\t"+ass.getInterfaceIdsString());
				}
			}
		}
		

	}
}

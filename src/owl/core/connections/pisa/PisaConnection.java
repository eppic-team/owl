package owl.core.connections.pisa;

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
	
	private static final int MAX_ENTRIES_PER_REQUEST = 50; // pisa doesn't specify a limit but recommends 20-50 per request
	
	private String interfacesUrl;
	//private String assembliesUrl;
	//private String pdbAssembliesUrl;
	
	public PisaConnection(String interfacesUrl, String assembliesUrl, String pdbAssembliesUrl)	{
		this.interfacesUrl = interfacesUrl;
		//this.assembliesUrl = assembliesUrl;
		//this.pdbAssembliesUrl = pdbAssembliesUrl;
	}
	
	/**
	 * Retrieves the XML PISA interface description from the PISA web server dividing the 
	 * query into chunks of {@value #MAX_ENTRIES_PER_REQUEST}, parses it and 
	 * returns the result as a map of pdb codes to lists of PISA interfaces 
	 * @param pdbCodesList
	 * @return
	 * @throws IOException
	 * @throws SAXException
	 */
	public Map<String,ChainInterfaceList> getInterfacesDescription(List<String> pdbCodesList) throws IOException, SAXException {
		Map<String,ChainInterfaceList> allInterfaces = new HashMap<String,ChainInterfaceList>();
		// we do batches of MAX_ENTRIES_PER_REQUEST
		for (int i=0;i<pdbCodesList.size();i+=MAX_ENTRIES_PER_REQUEST) {
			String commaSepList = "";
			for (int c=i;c<i+MAX_ENTRIES_PER_REQUEST && c<pdbCodesList.size();c++) {
				if (c!=i) commaSepList+=",";
				commaSepList+=pdbCodesList.get(c);
			}
			allInterfaces.putAll(getInterfacesDescription(commaSepList));
		}		
		return allInterfaces;
	}
	
	/**
	 * Retrieves the XML PISA interface description from the PISA web server, parses it and 
	 * returns the result as a map of pdb codes to lists of PISA interfaces 
	 * @param commaSepList
	 * @return
	 * @throws IOException
	 * @throws SAXException
	 */
	private Map<String,ChainInterfaceList> getInterfacesDescription(String commaSepList) throws IOException, SAXException {
		URL interfacesURL = new URL(interfacesUrl+commaSepList);
		URLConnection conn = interfacesURL.openConnection();
		
		PisaInterfaceXMLParser pxmlParser = new PisaInterfaceXMLParser(conn.getInputStream());
		return pxmlParser.getAllInterfaces();
	}
	
	public static void main(String[] args) throws Exception {
		PisaConnection pc = new PisaConnection(PISA_INTERFACES_URL, null, null);
		List<String> pdbCodes = new ArrayList<String>();
		pdbCodes.add("1aor");
		pdbCodes.add("1bxy");
		Map<String,ChainInterfaceList> all = pc.getInterfacesDescription(pdbCodes);
		for (String pdbCode:all.keySet()) {
			System.out.println("#####");
			System.out.println("# "+pdbCode);
			System.out.println("#####");
			ChainInterfaceList interfaces = all.get(pdbCode);
			for (ChainInterface interf:interfaces) {
				
				interf.printTabular(System.out);
			}
		}
	}
}

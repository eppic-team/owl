package owl.core.connections;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Collection;
import java.util.TreeMap;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.features.StructuralDomainFeature;
import owl.core.features.StructuralDomainType;
import owl.core.util.Interval;
import owl.core.util.IntervalSet;



/**
 * Connector class for querying the PDomains online service for Structural Domain Assignment at http://pdomains.sdsc.edu/v2/.
 * @author stehr
 */
public class PDomainsConnection {
	
	/*------------------------------- constants -----------------------------*/
	//private static final String queryURL = "http://pdomains.sdsc.edu/v2/compareassignments.php?pdbId=%s&chain=%s&CATH=%b&SCOP=%b&pdp=%b&dp=%b&DHcL=%b&DDomain=%b&PUU=%b&NCBI=%b&Dodis=%b";
	private static final String QUERY_URL = "http://pdomains.sdsc.edu/v2/compareassignments.php?pdbId=%s&chain=%s&CATH=true&SCOP=true&pdp=true&dp=true&DHcL=true&DDomain=true&PUU=true&NCBI=true&Dodis=true";
	//private static final String REGEXP_DOMAINS = "<td width=\"18%\" align=center>dp</td>";
	private static final String REGEXP_PARAM = "<PARAM NAME=\"([a-z]+\\d+)\" VALUE=\"(\\w+)\">";

	/*--------------------------------- main --------------------------------*/
	
	/**
	 * Queries PDomains service for domain assignments with the given method for the given pdb and chain id
	 * @return a possibly empty collection of StructuralDomainFeature objects representing the domain assignments
	 */
	public Collection<StructuralDomainFeature> getDomainsAsFeatures(String pdbId, String chainId, StructuralDomainType method) {
		Collection<StructuralDomainFeature> features = new LinkedList<StructuralDomainFeature>();
		TreeMap<String, IntervalSet> domainMap = getDomains(pdbId, chainId, method);
		for(String domName:domainMap.keySet()) {
			IntervalSet intervals = domainMap.get(domName);
			features.add(new StructuralDomainFeature(method,domName,intervals));
		}
		return features;
	}
	
	/**
	 * Queries the PDomains service for domain assignments with the given method for the given pdb and chain id
	 * @return a possibly empty collection of StructuralDomainFeatures containing the domain assignments
	 */
	public TreeMap<String, IntervalSet> getDomains(String pdbId, String chainId, StructuralDomainType method) {
		TreeMap<String, IntervalSet> domainMap = new TreeMap<String, IntervalSet>();		
		IntervalSet fragments = new IntervalSet();
		URL searchUrl = null;
		Pattern p = Pattern.compile(REGEXP_PARAM);
		try {
			searchUrl = new URL(String.format(QUERY_URL, pdbId, chainId));
			BufferedReader in = new BufferedReader(new InputStreamReader(searchUrl.openStream()));
			String line;
			StringBuilder sb = new StringBuilder();
			while((line = in.readLine()) != null) {
				sb.append(line);
			}
			Matcher m = p.matcher(sb);
		    int numDomains = 0;
		    int numFragments = 0;
	    	int start = -1;
	    	int end = -1;
	    	String domainName = null;
		    while(m.find()) {
		    	String param = m.group(1);
		    	String val = m.group(2);
		    	if(param.startsWith("domain")) { domainName = val; }
		    	if(param.startsWith("start")) { start = Integer.parseInt(val); }	// pDomains reports cif positions!
		    	if(param.startsWith("end")) { end = Integer.parseInt(val); }		// pDomains reports cif positions!
		    	if(param.startsWith("method")) {
		    		if(val.equals(method.toString())) {
		    			// add new fragment
		    			numFragments++;
		    			fragments.addInterval(start, end);
		    			if(!domainMap.containsKey(domainName)) {
		    				// create new domain
		    				IntervalSet newSet = IntervalSet.createFromInterval(new Interval(start, end));
		    				domainMap.put(domainName, newSet);
		    				numDomains++;
		    			} else {
		    				IntervalSet existingSet = domainMap.get(domainName);
		    				existingSet.add(new Interval(start, end));
		    			}
		    		}
		    	}
		    }
		} catch (MalformedURLException e) {
			System.err.println("Malformed Url: " + searchUrl);
		} catch (IOException e) {
			System.err.println("Error reading from Url " + searchUrl + ": " + e.getMessage());
		}
		return domainMap;
	}
	
	public static void main(String[] args) {
		PDomainsConnection pDomains = new PDomainsConnection();
		StructuralDomainType myMethod = StructuralDomainType.DP;
		String pdbId = "2rd0";
		String chainId = "A";
		System.out.println("Searching for " + pdbId + chainId);
		Collection<StructuralDomainFeature> domains = pDomains.getDomainsAsFeatures(pdbId, chainId, myMethod);
		if(domains.size() == 0) {
			System.out.println("No data found for " + pdbId + chainId);
		} else {
			for(StructuralDomainFeature f: domains) {
				System.out.println(f);
			}
		}
		System.out.println("done.");
	}
}

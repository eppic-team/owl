package owl.core.connections;

import java.io.*;
import java.net.*;
import java.util.*;
import java.util.regex.*;

import owl.core.features.PhosphoSitePlusFeature;
import owl.core.features.ProteinModificationType;
import owl.core.structure.AminoAcid;


/**
 * Provides methods to query the PhosphoSitePlus database. Note: the online connection to the database is not implemented yet, so
 * currently this class only provides methods for parsing html files which were manually downloaded from the site http://www.phosphosite.org.
 * @author stehr
 */
public class PhosphoSitePlusConnection {
	
	/**
	 * A posttranslational modification record from the PhosphoSitePlus database.
	 * @author stehr
	 */
	public class Modification {
				
		public int modId;				// the accession number of the modification site in PhoshoSitePlus
		public ProteinModificationType type;			// the type of modication (e.g. phosphorilation)
		public int pos;					// the position of the modification site (in uniprot sequence coordinates)
		public AminoAcid aa;			// the base type of the modified amino acid
		
		Modification(int modId, ProteinModificationType type, int pos, AminoAcid aa) {
			this.modId = modId;
			this.type = type;
			this.pos = pos;
			this.aa = aa;
		}
		
		public String toString() {
			return String.format("%s%d-%s",this.aa.getOneLetterCode(), this.pos, this.type.getSymbol());
		}
	}
	
	/*------------------------------ constants ------------------------------*/
	public static final String SEARCH_URL = "http://www.phosphosite.org/simpleSearchSubmitAction.do?queryId=-1&from=0&searchStr=%s&x=0&y=0";
	public static final String SHOW_PROTEIN_URL = "http://www.phosphosite.org/proteinAction.do;jsessionid=%s?id=%d&showAllSites=true";
	public static final String SHOW_SITE_URL = "http://www.phosphosite.org/siteAction.do?id=%d";
	
	public static final String REGEXP_PROT = "/proteinAction.do;jsessionid=(\\w+)\\?id=(\\d+)&showAllSites=true\" class=\"link13HoverRed\">human";
	public static final String REGEXP_SITE = "/siteAction.do\\?id=(\\d+)\" class=\"linkSite\">([A-Z])(\\d+)-(\\w+)</a>";
	public static final String REGEXP_SPECIES = "/proteinAction.do\\?id=\\d+&showAllSites=true\" class=\"link12UnderBlue\">(.+)</a>";
	
	/*---------------------------- public methods ---------------------------*/
	
	public static void testCookies() {
	    try {
		String urlString = "http://www.phosphosite.org/proteinAction.do;jsessionid=77F3444971272A0D490DA94937560A1F?id=14709&showAllSites=true";
	    CookieManager manager = new CookieManager();
	    manager.setCookiePolicy(CookiePolicy.ACCEPT_ALL);
	    CookieHandler.setDefault(manager);
	    URL url = new URL(urlString);
	    URLConnection connection = url.openConnection();
	    @SuppressWarnings("unused")
		Object obj = connection.getContent();
	    url = new URL(urlString);
	    connection = url.openConnection();
	    obj = connection.getContent();
	    CookieStore cookieJar = manager.getCookieStore();
	    List<HttpCookie> cookies = cookieJar.getCookies();
	    for (HttpCookie cookie: cookies) {
	      System.out.printf("Cookie: %s%n", cookie);
	    }
	    } catch(IOException e) {
	    	e.printStackTrace();
	    }
	}
	
	/**
	 * Parses modification sites from an HTML file downloaded from PhosphoSitePlus. This method processes
	 * files which were retrieved from PhosphoSitePlus using the URL
	 * http://www.phosphosite.org/proteinAction.do?id=14709&showAllSites=true. Currently, this cannot be
	 * done automatically beacause the website does not allow access to robots.
	 * @param htmlFile the HTML file downloaded from PhoshoSitePlus.
	 * @return the modifications found in the given HTML file.
	 */
	public static Collection<PhosphoSitePlusFeature> getModifications(File htmlFile) {
		LinkedList<PhosphoSitePlusFeature> sites = new LinkedList<PhosphoSitePlusFeature>();
		Pattern p1 = Pattern.compile(REGEXP_SPECIES);
		Pattern p2 = Pattern.compile(REGEXP_SITE);
		try {
			BufferedReader in = new BufferedReader(new FileReader(htmlFile));
			String line;
			boolean human = false;
			while((line = in.readLine()) != null) {
				//System.out.println(line);
				Matcher m1 = p1.matcher(line);
				if(m1.find()) {
					String species = m1.group(1);
					if(species.equals("human")) {
						human = true;
					} else {
						if(human && species.indexOf("iso") >= 0) {
							System.out.print("(skipping alternative isoforms)");
						}
						human = false;
					}
				}
				if(human) {
					Matcher m2 = p2.matcher(line);
					if(m2.find()) {
						int siteId = Integer.parseInt(m2.group(1));
						AminoAcid aa = AminoAcid.getByOneLetterCode(m2.group(2).charAt(0));
						int pos = Integer.parseInt(m2.group(3));
						ProteinModificationType type = ProteinModificationType.getBySymbol(m2.group(4).charAt(0));
						sites.add(new PhosphoSitePlusFeature(siteId,type,pos,aa));
					}
				}
			}
		} catch (IOException e) {
			System.err.println("Error reading from file " + htmlFile + ": " + e.getMessage());
		}
		return sites;
	}
	
	/**
	 * Returns a collection of modifications for the given protein retrieved from PhosphoSitePlus.
	 * TODO: This currently does not work due to cookie-issues. Workaround: Parse from HTML file.
	 */
	public Collection<Modification> getModifications(String protName) {
		LinkedList<Modification> sites = new LinkedList<Modification>();
		
	    CookieManager manager = new CookieManager();
	    manager.setCookiePolicy(CookiePolicy.ACCEPT_ALL);
	    CookieHandler.setDefault(manager);
		
		// search for protein id
		URL searchUrl = null;
		Pattern p = Pattern.compile(REGEXP_PROT);
		String sessionId = null;
		int protId = -1;
		try {
			searchUrl = new URL(String.format(SEARCH_URL,protName));
			BufferedReader in = new BufferedReader(new InputStreamReader(searchUrl.openStream()));
			String line;
			while((line = in.readLine()) != null) {
				//System.out.println(line);
				Matcher m = p.matcher(line);
				if(m.find()) {
					sessionId = m.group(1);
					String g = m.group(2);
					protId = Integer.parseInt(g);
				}
			}
		} catch (MalformedURLException e) {
			System.err.println("Malformed Url: " + searchUrl);
		} catch (IOException e) {
			System.err.println("Error reading from Url " + searchUrl);
		}
		if(sessionId == null) {
			System.err.println("No entry found for " + protName);
		} else {
			System.out.println("SessionID: " + sessionId);
			System.out.println("ProteinID: " + protId);
		}
		
		// retrieve cookies
	    CookieStore cookieJar = manager.getCookieStore();
	    List<HttpCookie> cookies = cookieJar.getCookies();
	    String cookiesStr = "";
	    for (HttpCookie cookie: cookies) {
	      System.out.printf("Cookie: %s%n", cookie);
	      cookiesStr += cookie + "; ";
	    }
	    cookiesStr += "testcookie=yes; ";
	    cookiesStr += "cookieEnabled=yes";
	    //cookiesStr = cookiesStr.substring(0, cookiesStr.length()-1);
	    
		// find modification sites
		URL protUrl = null;
		Pattern p2 = Pattern.compile(REGEXP_SITE);
		try {
			protUrl = new URL(String.format(SHOW_PROTEIN_URL,sessionId, protId));
			URLConnection conn = protUrl.openConnection();
		    conn.setRequestProperty("Cookie", cookiesStr);
		    conn.setRequestProperty("User-Agent", "Mozilla/5.0 (X11; U; Linux i686; en-GB; rv:1.9.0.14) Gecko/2009082707 Firefox/3.0.14");
		    conn.setRequestProperty("cookieEnabled","true");
		    conn.setRequestProperty("javaEnabled","true");
		    System.out.println(protUrl);
			//protUrl = new URL("http://www.phosphosite.org/proteinAction.do;jsessionid=77F3444971272A0D490DA94937560A1F?id=14709&showAllSites=true");
			BufferedReader in = new BufferedReader(new InputStreamReader(conn.getInputStream()));
			String line;
			while((line = in.readLine()) != null) {
				System.out.println(line);
				Matcher m = p2.matcher(line);
				if(m.find()) {
					int siteId = Integer.parseInt(m.group(1));
					AminoAcid aa = AminoAcid.getByOneLetterCode(m.group(2).charAt(0));
					int pos = Integer.parseInt(m.group(3));
					ProteinModificationType type = ProteinModificationType.getBySymbol(m.group(4).charAt(0));
					Modification mod = new Modification(siteId,type,pos, aa);
					sites.add(mod);
				}
			}
		} catch (MalformedURLException e) {
			System.err.println("Malformed Url: " + searchUrl);
		} catch (IOException e) {
			System.err.println("Error reading from Url " + searchUrl);
		}
		return sites;
	}
	
	/*--------------------------------- main --------------------------------*/
	public static void main(String[] args) {
		File phophoSiteDir = new File("/project/StruPPi/henning/projects/mutanom/analysis/data/PhosphoSitePlus");
		String protName = "MLH1";
		File htmlFile = new File(phophoSiteDir, protName + ".html");
		//int protId = 14709;
		//int siteId = 78796;
		
		//testCookies();
		Collection<PhosphoSitePlusFeature> sites = PhosphoSitePlusConnection.getModifications(htmlFile);
		if(sites.size() == 0) {
			System.err.println("No modification sites founds");
		} else {
			System.out.println("Modification sites in " + protName + ":");
			for(PhosphoSitePlusFeature m:sites) {
				System.out.println(m);
			}
		}
		System.out.println("done.");
	}
	
}

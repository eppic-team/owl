package owl.mutanom.core;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import owl.core.structure.AminoAcid;

/**
 * Connection class to retrieve SNP entries from NCBI's dbSNP database.
 * @author stehr
 */
public class DbSnpConnection  extends DefaultHandler {
	
	/*------------------------------ constants ------------------------------*/
	// add constants for dbSNP online access
	
	/*--------------------------- member variables --------------------------*/
	// private variables used for xml parsing
	private int snpId;					// fields of current SNP record
	private AminoAcid wtAA;
	private AminoAcid mutAA;
	private String mutStr;
	private int mutPos;
	private String mutDnaStr;
	private double freq;
	
	private boolean npMatches; 
	private String npToMatch;
	private String nmToMatch;	
	
	private LinkedList<SNP> snps;		// snps to be returned
	private HashSet<Integer> rsIds;		// remember ids to avoid duplicates
	private String tempVal;				// temporary variable for xml parsing
	/*----------------------------- constructors ----------------------------*/
	
	/**
	 * Initialize a new dbSNP connection.
	 */
	public DbSnpConnection() {
	}

	/*---------------------------- public methods ---------------------------*/
	
	/**
	 * Returns a collection of SNPs parsed from a dbSNP XML file.
	 * @param snpFile the xml file to parse
	 * @param np a refseq np id, the SNP will only be reported if the NP id matches
	 * @param nm a refseq nm id, the DNA mutations will only be reported if the NM id matches
	 */
	public Collection<SNP> parseSnpsFromXmlFile(File snpFile, String npStr, String nmStr) {
		
		snps = new LinkedList<SNP>();
		rsIds = new HashSet<Integer>();
		
		if(npStr == null || nmStr == null) {
			return snps;
		}
		
		npToMatch = npStr;
		nmToMatch = nmStr;		
		
		//get a factory
		SAXParserFactory spf = SAXParserFactory.newInstance();
		try {

			//get a new instance of parser
			SAXParser sp = spf.newSAXParser();

			//parse the file and also register this class for call backs
			sp.parse(snpFile, this);

		}catch(SAXException se) {
			se.printStackTrace();
		}catch(ParserConfigurationException pce) {
			pce.printStackTrace();
		}catch (IOException ie) {
			ie.printStackTrace();
		}

		return snps;		
	}
	
	/*---------------------------- private methods --------------------------*/
	
	/**
	 * Given two refseq nm or np ids, check whether their base id is equal.
	 * Base id is the number without the subvsersion, e.g. NM_000546.4 -> 000546
	 * @return true if the two ids are of the same type and their base number matches
	 */
	private boolean idsMatch(String id1, String id2) {
		if(!id1.substring(0,2).equals(id2.substring(0,2))) return false;
		String[] fields1 = id1.split("\\.");
		String[] fields2 = id2.split("\\.");
		int n1 = Integer.parseInt(fields1[0].substring(3));
		int n2 = Integer.parseInt(fields2[0].substring(3));
		return n1==n2;
	}
	
	/*---------------------------- event handlers ---------------------------*/
	
	public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
		tempVal = "";
		if(qName.equalsIgnoreCase("Rs")) {
			// new snp element: reset temp values
			snpId = -1;
			wtAA = null;
			mutAA = null;
			mutPos = -1;
			mutStr = null;
			mutDnaStr = "c.?";
			freq = Double.NaN;
			npMatches = false;
			
			// set new snpId
			snpId = Integer.parseInt(attributes.getValue("rsId"));	
		
		} else if (qName.equalsIgnoreCase("RsStruct")) {
			// <RsStruct protAcc="NP_000537" protGi="120407068" protLoc="282" protResidue="R" rsResidue="W" structGi="1310770" structLoc="189" structResidue="R"/>
			String np = attributes.getValue("protAcc");
			if(idsMatch(np, npToMatch)) {
				npMatches = true;
				wtAA = AminoAcid.getByOneLetterCode(attributes.getValue("protResidue").charAt(0));
				mutAA = AminoAcid.getByOneLetterCode(attributes.getValue("rsResidue").charAt(0));
				mutPos = Integer.parseInt(attributes.getValue("protLoc"));
			}
		}	
	}
	
	public void characters(char[] ch, int start, int length) throws SAXException {
		tempVal = new String(ch,start,length);
	}

	public void endElement(String uri, String localName, String qName) throws SAXException {

			if(qName.equalsIgnoreCase("Rs")) {
				// store/print mutation
				if(npMatches) {	// the correct NP id was found in the file
					if(wtAA != null && mutAA != null && mutPos > 0 && snpId > 0) {
						if(!rsIds.contains(snpId)) {
							SNP newSnp = new SNP(wtAA, mutAA, mutPos, snpId, freq);
							newSnp.setDnaMutStr(mutDnaStr);
							snps.add(newSnp);
							rsIds.add(snpId);
						}
					} else {
						//System.err.printf("SNP: %d\t%s\t%s\t%s%d%s\t%s\n", snpId, nmId, npId, wtAA, mutPos, mutAA, mutDnaStr);
					}

				}
			} else if(qName.equalsIgnoreCase("hgvs")) {
				//System.out.println(tempVal);
				if(tempVal.startsWith("NM")) {
					String[] fields = tempVal.split(":");
					if(idsMatch(fields[0],nmToMatch)) {
						mutDnaStr = fields[1];
					}
				} else if(tempVal.startsWith("NP")) {
					String[] fields = tempVal.split(":");
					if(idsMatch(fields[0],npToMatch)) {
						npMatches = true;
						mutStr = fields[1];
						if(mutStr!=null && mutStr.startsWith("p.")) {
							Pattern p = Pattern.compile("^p.(\\w{3})(\\d+)(\\w{3})$");
							Matcher m = p.matcher(mutStr);
							if(m.find()) {
								wtAA = AminoAcid.getByThreeLetterCode(m.group(1));
								mutPos = Integer.parseInt(m.group(2));
								mutAA = AminoAcid.getByThreeLetterCode(m.group(3));
							}
						}
					}
				} 
			} 
		}
	
	/*--------------------------------- main --------------------------------*/

	/**
	 * Main method for testing or command line use.
	 */
	public static void main(String[] args) {
		File snpFile = new File("/project/StruPPi/henning/projects/mutanom/analysis/data/snp/ncbi/xml/s_TP53.xml");
		String nm = "NM_000546";
		String np = "NP_000537";
		
		DbSnpConnection dbSnp = new DbSnpConnection();
		Collection<SNP> snps = dbSnp.parseSnpsFromXmlFile(snpFile,np,nm);
		for(SNP snp:snps) {
			System.out.println(snp + "\t" + snp.getDnaMutStr() + "\t" + snp.rsId);
		}
		System.out.println(snps.size() + " SNPs found.");
		System.out.println("done.");
		
	}
	
}

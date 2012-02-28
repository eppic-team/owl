package owl.core.connections;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.zip.GZIPInputStream;

import org.xml.sax.Attributes;
import org.xml.sax.ContentHandler;
import org.xml.sax.InputSource;
import org.xml.sax.Locator;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.XMLReaderFactory;

import owl.core.sequence.UnirefEntry;


public class UnirefXMLParser implements ContentHandler {

	// xml tags
	private static final String ENTRY_TAG = "entry";
	private static final String REFERENCE_SEQUENCE_TAG = "referenceSequence";
	private static final String DB_REFERENCE_TAG = "dbReference";
	private static final String PROPERTY_TAG = "property";
	private static final String ORGANISM_TAG = "organism";
	private static final String SEQUENCE_TAG = "sequence";
	private static final String REPRESENTATIVE_MEMBER_TAG = "representativeMember";
	private static final String MEMBER_TAG = "member";
	
	private File unirefXMLFile;
	private InputSource input;
		
	private StringBuffer buffer;
	
	private boolean inValue;
	
	private UnirefEntry currentEntry;
	
	private boolean inReferenceSequence;
	private boolean inDbReferenceUniProt;
	private boolean inDbReferenceUniParc;
	private boolean inOrganism;
	private boolean inRepresentativeMember;
	private boolean inMember;
	
	private int uniprotAccessionCounter;
	
	private PrintWriter out;
	private PrintWriter clustersOut;
	
	private boolean isGzipped;
	
	private String uniprotName; // in earlier versions it is "UniProt", later on it became "UniProtKB"
	
	/**
	 * Constructs new UnirefXMLParser, use subsequently {@link #parseIntoTabDelimited(File)}
	 * @param unirefXMLFile either an xml file or a gzipped xml file (with extension gz)
	 */
	public UnirefXMLParser(File unirefXMLFile) {
		this.isGzipped = false;
		this.unirefXMLFile = unirefXMLFile;
		if (this.unirefXMLFile.getName().endsWith(".gz")) isGzipped = true;
		this.uniprotName = "UniProt";
	}

	/**
	 * Parses the UniRef XML file outputting to outFile tab delimited
	 * with columns: uniref_id, uniprot_id, uniparc_id, tax_id, sequence    
	 * @param outFile
	 * @param clustersFile
	 * @throws SAXException
	 * @throws IOException
	 */
	public void parseIntoTabDelimitedFiles(File outFile, File clustersFile) throws SAXException, IOException {
		out = new PrintWriter(outFile);
		clustersOut = new PrintWriter(clustersFile);

		InputStream is = null;
		if (isGzipped) {
			is = new GZIPInputStream(new FileInputStream(this.unirefXMLFile));
		} else {
			is = new FileInputStream(this.unirefXMLFile);
		}
		
		this.input = new InputSource(is); 
		
		XMLReader parser = XMLReaderFactory.createXMLReader();
	 		
		parser.setContentHandler(this);
		
		parser.parse(input);
		
		out.close();
		clustersOut.close();
	}
	
	/**
	 * Initialises the buffer used for value reading.
	 * Subsequently {@link #characters(char[], int, int)} will write to the buffer
	 */
	private void initValueReading() {
		inValue = true;
		buffer = new StringBuffer();
	}
	
	/**
	 * Flushes the buffer returning the String accumulated.
	 * @return
	 */
	private String flushValue() {
		String readValue = buffer.toString();
		inValue = false;
		return readValue;
	}
	
	public void startDocument() throws SAXException {
		
		inValue = false;
		inOrganism = false;
		inReferenceSequence = false;
		inDbReferenceUniProt = false;
		inDbReferenceUniParc = false;
		inRepresentativeMember = false;
		inMember = false;
		uniprotAccessionCounter = 0;
		buffer = null;
	}

	public void endDocument() throws SAXException {

	}

	public void startElement(String uri, String localName, String name,
			Attributes atts) throws SAXException {

		if (name.equals(ENTRY_TAG)) {
			currentEntry = new UnirefEntry();
			currentEntry.setId(atts.getValue("id")); 
		}
		else if (name.equals(REPRESENTATIVE_MEMBER_TAG)) {
			inRepresentativeMember = true;
			// representativeMember tag is used in newer versions of the XML format (at least since 2010_01)
			// where uniprot is called UniProtKB
			uniprotName = "UniProtKB";
		}
		else if (name.equals(ORGANISM_TAG)) {
			inOrganism = true;
		}
		else if (name.equals(REFERENCE_SEQUENCE_TAG)) {
			inReferenceSequence = true;
		}
		else if (name.equals(MEMBER_TAG)) {
			inMember = true;
		}
		else if (name.equals(SEQUENCE_TAG)) {
			initValueReading();
		}
		if (inOrganism) {
			if (name.equals(DB_REFERENCE_TAG)){
				if (atts.getValue("type").equals("NCBI Taxonomy")) {
					currentEntry.setNcbiTaxId(Integer.parseInt(atts.getValue("id")));
				}
			} 
		}
		if (inReferenceSequence) {
			if (name.equals(DB_REFERENCE_TAG)) {
				if (atts.getValue("type").equals(uniprotName+" ID")) {
					inDbReferenceUniProt = true;	
					uniprotAccessionCounter = 0;
				} else if (atts.getValue("type").equals("UniParc ID")) {
					inDbReferenceUniParc = true;
				}
			}
		}
		if (inReferenceSequence && inDbReferenceUniProt) {
			if (name.equals(PROPERTY_TAG)) {
				if (atts.getValue("type").equals(uniprotName+" accession")) {
					if (uniprotAccessionCounter==0) {
						currentEntry.setUniprotId(atts.getValue("value"));
					} else {
						currentEntry.addInactiveUniprotId(atts.getValue("value"));
					}
					uniprotAccessionCounter++;
				} else if (atts.getValue("type").equals("UniParc ID")) {
					currentEntry.setUniparcId(atts.getValue("value"));
				}
			}
		}

		if (inRepresentativeMember) {
			if (name.equals(DB_REFERENCE_TAG)) {
				if (atts.getValue("type").equals(uniprotName+" ID")) {
					inDbReferenceUniProt = true;	
					uniprotAccessionCounter = 0;
				} else if (atts.getValue("type").equals("UniParc ID")) {
					currentEntry.setUniparcId(atts.getValue("id"));
					inDbReferenceUniParc = true;
				}
			}
		}
		if (inRepresentativeMember && inDbReferenceUniProt) {
			if (name.equals(PROPERTY_TAG)) {
				if (atts.getValue("type").equals(uniprotName+" accession")) {
					if (uniprotAccessionCounter==0) {
						currentEntry.setUniprotId(atts.getValue("value"));
					} else {
						currentEntry.addInactiveUniprotId(atts.getValue("value"));
					}
					uniprotAccessionCounter++;
				} else if (atts.getValue("type").equals("UniParc ID")) {
					currentEntry.setUniparcId(atts.getValue("value"));
				} else if (atts.getValue("type").equals("NCBI taxonomy")) {
					currentEntry.setNcbiTaxId(Integer.parseInt(atts.getValue("value")));
				}
			}
		}
		if (inRepresentativeMember && inDbReferenceUniParc) {
			if (name.equals(PROPERTY_TAG)) {
				if (atts.getValue("type").equals("NCBI taxonomy")) {
					currentEntry.setNcbiTaxId(Integer.parseInt(atts.getValue("value")));
				}
			}
		}
		
		if (inMember) {
			if (name.equals(DB_REFERENCE_TAG)) {
				if (atts.getValue("type").equals(uniprotName+" ID")) {
					inDbReferenceUniProt = true;
					uniprotAccessionCounter = 0;
				}
			}
		}
		if (inMember && inDbReferenceUniProt) {
			if (name.equals(PROPERTY_TAG)) {
				if (atts.getValue("type").equals(uniprotName+" accession")) {
					if (uniprotAccessionCounter==0) {
						currentEntry.addClusterMember(atts.getValue("value"));
					}
					uniprotAccessionCounter++;
				}
			}
		}
		
	}
	
	private void writeTabDelimited() {
		out.println(currentEntry.getTabDelimitedMySQLString()); 
		if (currentEntry.isUniprot()) {
			clustersOut.println(currentEntry.getUniprotId()+"\t"+currentEntry.getUniprotId());
			if (currentEntry.hasClusterMembers()) {
				 for (String memberuniprotid:currentEntry.getClusterMembers()){
					 clustersOut.println(currentEntry.getUniprotId()+"\t"+memberuniprotid);
				 }
			}
		}
	}

	public void endElement(String uri, String localName, String name)
			throws SAXException {

		if (name.equals(ENTRY_TAG)){
			if (currentEntry.hasNulls()) {
				System.err.println("Warning! nulls in Uniref entry "+currentEntry.getId());
			}
			writeTabDelimited();			
		}
		else if (name.equals(REPRESENTATIVE_MEMBER_TAG)) {
			inRepresentativeMember = false;
		}
		else if (name.equals(ORGANISM_TAG)){
			inOrganism = false;
		}
		else if (name.equals(REFERENCE_SEQUENCE_TAG)) {
			inReferenceSequence = false;
		}
		else if (name.equals(MEMBER_TAG)) {
			inMember = false;
		}
		else if (name.equals(SEQUENCE_TAG)) {
			String sequence = flushValue().replaceAll("\n", "");
			currentEntry.setSequence(sequence.trim());
		}
		if (inReferenceSequence) {
			if (name.equals(DB_REFERENCE_TAG)) {
				inDbReferenceUniProt = false;
				inDbReferenceUniParc = false;
			}
		}
		if (inRepresentativeMember) {
			if (name.equals(DB_REFERENCE_TAG)) {
				inDbReferenceUniProt = false;
				inDbReferenceUniParc = false;
			}
		}
		if (inMember) {
			if (name.equals(DB_REFERENCE_TAG)) {
				inDbReferenceUniProt = false;
			}
		}
		
	}

	public void characters(char[] ch, int start, int length)
	throws SAXException {
		if (inValue) {
			buffer.append(ch, start, length);
		}
	}

	
	
	
	/*--------------------  empty methods --------------------------*/
	
	public void startPrefixMapping(String prefix, String uri)
	throws SAXException {
	}

	public void endPrefixMapping(String prefix) throws SAXException {
	}

	public void ignorableWhitespace(char[] ch, int start, int length)
			throws SAXException {
	}

	public void processingInstruction(String target, String data)
			throws SAXException {
	}

	public void setDocumentLocator(Locator locator) {
	}

	public void skippedEntity(String name) throws SAXException {
	}
	
	/*-----------------------  main -------------------------------*/
	
	/**
	 * Testing it
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		if (args.length<3) {
			System.err.println("Usage: UnirefXMLParser <input xml or xml.gz file> <out data file> <out cluster members file>");
			System.exit(1);
		}
		UnirefXMLParser urxp = new UnirefXMLParser(new File(args[0]));
		urxp.parseIntoTabDelimitedFiles(new File(args[1]), new File(args[2]));
		
	}

}


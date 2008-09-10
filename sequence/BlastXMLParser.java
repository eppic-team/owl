package sequence;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.xml.sax.Attributes;
import org.xml.sax.ContentHandler;
import org.xml.sax.InputSource;
import org.xml.sax.Locator;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.XMLReaderFactory;

/**
 * Blast XML output parser using SAX
 * @author duarte
 *
 */
public class BlastXMLParser implements ContentHandler {

	// xml tags
	private static final String QUERY_LENGTH_TAG = "BlastOutput_query-len";
	private static final String QUERY_DEF_TAG = "BlastOutput_query-def";
	private static final String ITERATION_ITER_NUM_TAG = "Iteration_iter-num";
	private static final String ITERATIONS_HITS_TAG = "Iteration_hits";
	private static final String HIT_TAG = "Hit";
	private static final String HIT_DEF_TAG = "Hit_def";
	private static final String HIT_LEN_TAG = "Hit_len";
	private static final String HSP_TAG = "Hsp";
	private static final String HSP_BIT_SCORE_TAG = "Hsp_bit-score";
	private static final String HSP_EVALUE_TAG = "Hsp_evalue";
	private static final String HSP_QUERY_FROM_TAG = "Hsp_query-from";
	private static final String HSP_QUERY_TO_TAG = "Hsp_query-to";
	private static final String HSP_HIT_FROM_TAG = "Hsp_hit-from";
	private static final String HSP_HIT_TO_TAG = "Hsp_hit-to";
	private static final String HSP_IDENTITY = "Hsp_identity";
	private static final String HSP_ALIGN_LEN = "Hsp_align-len";
	private static final String HSP_QSEQ = "Hsp_qseq";
	private static final String HSP_HSEQ = "Hsp_hseq";
	
	private static final String ID_REGEX = "^\\s*(\\S+).*";

	private File blastXMLFile;
	private InputSource input;
	
	private BlastHitList hitList;
	
	private String queryId;
	private BlastHit currentHit;
	private String currentSubjectId;
	private int currentSubjectLength;
	
	private String currentQSeq;
	private String currentHSeq;
	private int currentIdentity;
	private int currentAliLength;
	
	private StringBuffer buffer;
	private boolean inValue;
	
	private boolean inIterationHits;
	private boolean inHit;
	private boolean inHsp;
	
	/**
	 * Constructs new BlastXMLParser parsing the given Blast XML file
	 * Get the hits calling {@link #getHits()}
	 * @param blastXMLFile
	 * @throws SAXException
	 * @throws IOException
	 */
	public BlastXMLParser(File blastXMLFile) throws SAXException, IOException{
		this.blastXMLFile = blastXMLFile;
		
		InputStream is = new FileInputStream(this.blastXMLFile);
		this.input = new InputSource(is); 
			
		XMLReader parser = XMLReaderFactory.createXMLReader();
	 		
		//BlastXMLHandler handler = new BlastXMLHandler();
		parser.setContentHandler(this);
		
		parser.parse(input);

	}
	
	/**
	 * Returns the parsed BlastHitList
	 * @return
	 */
	public BlastHitList getHits() {
		return hitList;
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
		hitList = new BlastHitList();
		inValue = false;
		inIterationHits = false;
		inHit = false;
		inHsp = false;
		buffer = null;
	}

	public void endDocument() throws SAXException {

	}

	public void startElement(String uri, String localName, String name,
			Attributes atts) throws SAXException {

		if (name.equals(QUERY_LENGTH_TAG)) {
			initValueReading();
		}
		else if (name.equals(QUERY_DEF_TAG)) {
			initValueReading();
		}
		else if (name.equals(ITERATION_ITER_NUM_TAG)) {
			initValueReading();
		}
		else if (name.equals(ITERATIONS_HITS_TAG)) {
			inIterationHits = true;
		}
		else if (name.equals(HIT_TAG)) {
			inHit = true;
		}
		if (inIterationHits && inHit) {
			if (name.equals(HIT_DEF_TAG)){
				initValueReading();
			} 
			else if (name.equals(HIT_LEN_TAG)) {
				initValueReading();
			}
			else if (name.equals(HSP_TAG)) {
				inHsp = true;
				this.currentHit = new BlastHit();
				this.currentHit.setQueryId(queryId);
				this.currentHit.setSubjectId(currentSubjectId);
				this.currentHit.setSubjectLength(currentSubjectLength);
			}
			if (inHsp) {
				if (name.equals(HSP_BIT_SCORE_TAG)) {
					initValueReading();
				}
				else if (name.equals(HSP_EVALUE_TAG)) {
					initValueReading();
				}
				else if (name.equals(HSP_QUERY_FROM_TAG)) {
					initValueReading();
				}
				else if (name.equals(HSP_QUERY_TO_TAG)) {
					initValueReading();
				}
				else if (name.equals(HSP_HIT_FROM_TAG)) {
					initValueReading();
				}
				else if (name.equals(HSP_HIT_TO_TAG)) {
					initValueReading();
				}
				else if (name.equals(HSP_IDENTITY)) {
					initValueReading();
				}
				else if (name.equals(HSP_ALIGN_LEN)) {
					initValueReading();
				}
				else if (name.equals(HSP_QSEQ)) {
					initValueReading();
				}
				else if (name.equals(HSP_HSEQ)) {
					initValueReading();
				}
			}
		}
		
		
	}

	public void endElement(String uri, String localName, String name)
			throws SAXException {

		if (name.equals(QUERY_LENGTH_TAG)){
			this.hitList.setQueryLength(Integer.parseInt(flushValue()));
		}
		else if (name.equals(QUERY_DEF_TAG)) {
			String val = flushValue();
			Pattern p = Pattern.compile(ID_REGEX);
			Matcher m = p.matcher(val);
			if (m.matches()) {
				this.queryId = m.group(1);
				this.hitList.setQueryId(queryId);
			}
		}
		else if (name.equals(ITERATION_ITER_NUM_TAG)) {
			if (Integer.parseInt(flushValue())>1) {
				System.err.println("WARNING: this BLAST XML file contains more than one iteration. Multiple iterations parsing not supported yet!");
			}
		}
		else if (name.equals(ITERATIONS_HITS_TAG)) {
			inIterationHits = false;
		}
		else if (name.equals(HIT_TAG)) {
			inHit = false;
		}
		if (inIterationHits && inHit) {
			if (name.equals(HIT_DEF_TAG)) {
				String val = flushValue();
				Pattern p = Pattern.compile(ID_REGEX);
				Matcher m = p.matcher(val);
				if (m.matches()) {
					this.currentSubjectId = m.group(1);
				}
			}
			else if (name.equals(HIT_LEN_TAG)){				
				this.currentSubjectLength = Integer.parseInt(flushValue());
			}
			else if (name.equals(HSP_TAG)) {
				this.currentHit.setAlignment(currentQSeq, currentHSeq);
				this.currentHit.setPercentIdentity((double)(100*currentIdentity)/((double)currentAliLength));
				this.hitList.add(currentHit);
				currentHit = null;
				inHsp = false;
			}
			if (inHsp) {
				if (name.equals(HSP_BIT_SCORE_TAG)) {
					this.currentHit.setScore(Double.parseDouble(flushValue()));
				}
				else if (name.equals(HSP_EVALUE_TAG)) {
					this.currentHit.setEValue(Double.parseDouble(flushValue()));
				}
				else if (name.equals(HSP_QUERY_FROM_TAG)) {
					this.currentHit.setQueryStart(Integer.parseInt(flushValue()));
				}
				else if (name.equals(HSP_QUERY_TO_TAG)) {
					this.currentHit.setQueryEnd(Integer.parseInt(flushValue()));
				}
				else if (name.equals(HSP_HIT_FROM_TAG)) {
					this.currentHit.setSubjectStart(Integer.parseInt(flushValue()));
				}
				else if (name.equals(HSP_HIT_TO_TAG)) {
					this.currentHit.setSubjectEnd(Integer.parseInt(flushValue()));
				}
				else if (name.equals(HSP_IDENTITY)) {
					this.currentIdentity = Integer.parseInt(flushValue());
				}
				else if (name.equals(HSP_ALIGN_LEN)) {
					this.currentAliLength = Integer.parseInt(flushValue());
					this.currentHit.setAliLength(currentAliLength);
				}
				else if (name.equals(HSP_QSEQ)) {
					this.currentQSeq = flushValue();
				}
				else if (name.equals(HSP_HSEQ)) {
					this.currentHSeq = flushValue();
				}
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
	
	
	/*--------------------------  main ----------------------------*/
	
	/**
	 * to test the class. See class BlastParsersTest for a complete junit4 test of this class
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {

		//File xmlfile = new File("/scratch/local/temp/blast.out.xml");
		File xmlfile = new File(args[0]);
		
		BlastXMLParser bxp = new BlastXMLParser(xmlfile);
		
		BlastHitList hitListXML = bxp.getHits();
		System.out.println("query length= "+hitListXML.getQueryLength());
		
		hitListXML.printSome(20);
		System.out.println("Total number of hits: "+hitListXML.size());

		// printing alignments
		for (BlastHit hit: hitListXML) {
			hit.getAlignment().printFasta();
			System.out.println();
		}
		
		
		
	}
}

package owl.core.connections.pisa;

import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;

import org.xml.sax.Attributes;
import org.xml.sax.ContentHandler;
import org.xml.sax.InputSource;
import org.xml.sax.Locator;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.XMLReaderFactory;



/**
 * PISA XML output parser using SAX
 * @author duarte
 *
 */
public class PisaAssembliesXMLParser implements ContentHandler {

	// xml tags
	private static final String PDB_ENTRY_TAG = "pdb_entry";
	private static final String ASM_SET_TAG = "asm_set";
	private static final String ASSEMBLY_TAG = "assembly";
	private static final String INTERFACES_TAG = "interfaces";

	private static final String PDB_CODE_TAG = "pdb_code";
	private static final String STATUS_TAG = "status";
	
	private static final String ID_TAG = "id";
	private static final String MMSIZE_TAG = "mmsize";
	private static final String SCORE_TAG = "score";
	private static final String	DISS_ENERGY_TAG = "diss_energy";
	private static final String	FORMULA_TAG = "formula";

	private static final String NOCC_TAG = "nocc";
		
	// members
	private InputSource input;
	
	private Map<String,PisaAsmSetList> allAssemblies;
	
	private PisaAsmSetList currentAsmSetList;
	private String currentPdbCode;
	private PisaAsmSet currentAsmSet;
	private PisaAssembly currentPisaAssembly;
	
	private StringBuffer buffer;
	private boolean inValue;
	
	private boolean inEntry;
	private boolean inAsmSet;
	private boolean inAssembly;
	private boolean inInterfaces;
	
	private int currentInterfaceId;
	
	/**
	 * Constructs new PisaAssembliesXMLParser and parses the given XML stream
	 * Get the list calling {@link #getAllAssemblies()}
	 * @param is the stream with XML data for assemblies description
	 */
	public PisaAssembliesXMLParser(InputStream is) throws SAXException, IOException {
		this.input = new InputSource(is);
		XMLReader parser = XMLReaderFactory.createXMLReader();
 		
		parser.setContentHandler(this);
		
		parser.parse(input);
	}
	
	/**
	 * Returns the parsed assemblies data
	 * @return a map of pdb codes to lists of assemblies
	 */
	public Map<String,PisaAsmSetList> getAllAssemblies() {
		return allAssemblies;
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
		allAssemblies = new HashMap<String,PisaAsmSetList>();
		inEntry = false;
		inValue = false;
		inAsmSet = false;
		inAssembly = false;
		inInterfaces = false;
		buffer = null;
	}

	public void endDocument() throws SAXException {

	}

	public void startElement(String uri, String localName, String name,
			Attributes atts) throws SAXException {
		if (name.equals(PDB_ENTRY_TAG)){
			inEntry = true;
			currentAsmSetList = new PisaAsmSetList();
		}
		if (inEntry) {
			if (name.equals(PDB_CODE_TAG)) {
				initValueReading();
			} else if (name.equals(STATUS_TAG)) {
				initValueReading();
			} else if (name.equals(ASM_SET_TAG)) {
				inAsmSet = true;
				currentAsmSet = new PisaAsmSet();
			}
			if (inAsmSet && !inAssembly) {
				if (name.equals(ASSEMBLY_TAG)) {
					inAssembly = true;
					currentPisaAssembly = new PisaAssembly();
				} 
			}
			
			if (inAsmSet && inAssembly && !inInterfaces) {
				if (name.equals(ID_TAG)) {
					initValueReading();
				} else if (name.equals(MMSIZE_TAG)){
					initValueReading();
				} else if (name.equals(SCORE_TAG)) {
					initValueReading();
				} else if (name.equals(DISS_ENERGY_TAG)) {
					initValueReading();
				} else if (name.equals(FORMULA_TAG)) {
					initValueReading();
				} else if (name.equals(INTERFACES_TAG)) {
					inInterfaces = true;
				}
			}
			if (inAsmSet && inAssembly && inInterfaces) {
				if (name.equals(ID_TAG)) {
					initValueReading();
				} else if (name.equals(NOCC_TAG)) {
					initValueReading();
				}
			}


		} 
	}

	public void endElement(String uri, String localName, String name)
			throws SAXException {

		if (name.equals(PDB_ENTRY_TAG)) {
			inEntry = false;
			allAssemblies.put(currentPdbCode, currentAsmSetList);
		}
		if (inEntry) {
			if (name.equals(PDB_CODE_TAG)) {
				currentPdbCode = flushValue().toLowerCase();
				currentAsmSetList.setPdbCode(currentPdbCode);
			} else if (name.equals(STATUS_TAG)) {
				currentAsmSetList.setStatus(flushValue());
			} else if (name.equals(ASM_SET_TAG)){
				inAsmSet = false;
				currentAsmSetList.add(currentAsmSet);
			}
			if (inAsmSet && inAssembly && !inInterfaces) {
				if (name.equals(ID_TAG)) {
					currentPisaAssembly.setId(Integer.parseInt(flushValue()));
				} else if (name.equals(MMSIZE_TAG)){
					currentPisaAssembly.setMmsize(Integer.parseInt(flushValue()));
				} else if (name.equals(SCORE_TAG)) {
					currentPisaAssembly.setScore(flushValue());
				} else if (name.equals(DISS_ENERGY_TAG)) {
					currentPisaAssembly.setDissEnergy(Double.parseDouble(flushValue()));
				} else if (name.equals(FORMULA_TAG)) {
					currentPisaAssembly.setFormula(flushValue());
				} else if (name.equals(ASSEMBLY_TAG)) {
					inAssembly = false;
					currentAsmSet.add(currentPisaAssembly);
				}	
			}
			if (inAsmSet && inAssembly && inInterfaces) {
				if (name.equals(ID_TAG)) {
					currentInterfaceId = Integer.parseInt(flushValue());
				} else if (name.equals(NOCC_TAG)) {
					int nocc = Integer.parseInt(flushValue());
					if (nocc>0) {
						currentPisaAssembly.addInterfaceId(currentInterfaceId);
					}
				} else if (name.equals(INTERFACES_TAG)) {
					inInterfaces = false;
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
	
	
}

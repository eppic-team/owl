package owl.core.connections.pisa;

import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;

import javax.vecmath.Matrix4d;

import org.xml.sax.Attributes;
import org.xml.sax.ContentHandler;
import org.xml.sax.InputSource;
import org.xml.sax.Locator;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.XMLReaderFactory;

import owl.core.structure.AminoAcid;
import owl.core.structure.ChainInterface;
import owl.core.structure.ChainInterfaceList;
import owl.core.structure.Pdb;
import owl.core.structure.Residue;
import owl.core.structure.SpaceGroup;


/**
 * PISA XML output parser using SAX
 * @author duarte
 *
 */
public class PisaInterfaceXMLParser implements ContentHandler {

	// xml tags
	private static final String PDB_ENTRY_TAG = "pdb_entry";
	private static final String INTERFACE_TAG = "interface";
	private static final String MOLECULE_TAG = "molecule";
	private static final String RESIDUE_TAG = "residue";

	private static final String PDB_CODE_TAG = "pdb_code";
	
	private static final String ID_TAG = "id";
	private static final String TYPE_TAG = "type";
	private static final String INT_AREA_TAG = "int_area";
	private static final String	INT_SOLV_EN_TAG = "int_solv_en";
	//private static final String	PVALUE_TAG = "pvalue";
	private static final String CHAIN_ID_TAG = "chain_id";
	private static final String CLASS_TAG = "class";
	
	private static final String SYMOP_TAG = "symop";
	private static final String RXX_TAG = "rxx";
	private static final String RXY_TAG = "rxy";
	private static final String RXZ_TAG = "rxz";
	private static final String RYX_TAG = "ryx";
	private static final String RYY_TAG = "ryy";
	private static final String RYZ_TAG = "ryz";
	private static final String RZX_TAG = "rzx";
	private static final String RZY_TAG = "rzy";
	private static final String RZZ_TAG = "rzz";
	private static final String TX_TAG  = "tx";
	private static final String TY_TAG  = "ty";
	private static final String TZ_TAG  = "tz";

	private static final String SER_NO_TAG = "ser_no"; 
	private static final String SEQ_NUM_TAG = "seq_num";
	private static final String NAME_TAG = "name";
	private static final String ASA_TAG = "asa";
	private static final String BSA_TAG = "bsa";
	//private static final String SOLV_EN_TAG = "solv_en";
	
	// members
	private InputSource input;
	
	private Map<String,ChainInterfaceList> allInterfaces;
	
	private ChainInterfaceList currentInterfaces;
	private String currentPdbCode;
	private ChainInterface currentPisaInterface;
	private Pdb currentPisaMolecule;
	private Residue currentResidue;
	private int currentMolecId;
	private Matrix4d currentTransfOrth;
	private Matrix4d currentTransf;
	private String currentMolType;
	
	private StringBuffer buffer;
	private boolean inValue;
	
	private boolean inEntry;
	private boolean inInterface;
	private boolean inMolecule;
	private boolean inResidue;
	
	/**
	 * Constructs new PisaXMLParser and parses the given XML stream
	 * Get the list calling {@link #getAllInterfaces()}
	 * @param is the stream with XML data for interfaces description
	 */
	public PisaInterfaceXMLParser(InputStream is) throws SAXException, IOException {
		this.input = new InputSource(is);
		XMLReader parser = XMLReaderFactory.createXMLReader();
 		
		parser.setContentHandler(this);
		
		parser.parse(input);
	}
	
	/**
	 * Returns the parsed interfaces data
	 * @return a map of pdb codes to lists of interfaces
	 */
	public Map<String,ChainInterfaceList> getAllInterfaces() {
		return allInterfaces;
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
		allInterfaces = new HashMap<String,ChainInterfaceList>();
		inEntry = false;
		inValue = false;
		inInterface = false;
		inMolecule = false;
		inResidue = false;
		buffer = null;
	}

	public void endDocument() throws SAXException {

	}

	public void startElement(String uri, String localName, String name,
			Attributes atts) throws SAXException {
		if (name.equals(PDB_ENTRY_TAG)){
			inEntry = true;
			currentInterfaces = new ChainInterfaceList(ChainInterfaceList.AsaCalcMethod.PISA);
		}
		if (inEntry) {
			if (name.equals(PDB_CODE_TAG)) {
				initValueReading();
			} else if (name.equals(INTERFACE_TAG)) {
				inInterface = true;
				currentPisaInterface = new ChainInterface();
			}
			if (inInterface && !inMolecule) {
				if (name.equals(ID_TAG)) {
					initValueReading();
				} else if (name.equals(TYPE_TAG)) {
					initValueReading();
				} else if (name.equals(INT_AREA_TAG)){
					initValueReading();
				} else if (name.equals(INT_SOLV_EN_TAG)){
					initValueReading();
				//} else if (name.equals(PVALUE_TAG)){
				//	initValueReading();
				} else if (name.equals(MOLECULE_TAG)) {
					inMolecule = true;
					currentPisaMolecule = new Pdb();
					currentMolecId = 0;
					currentTransfOrth = new Matrix4d();
					currentMolType = null;
				}
				
			}
			if (inInterface && inMolecule) {
				if (name.equals(ID_TAG)) {
					initValueReading();
				} else if (name.equals(CHAIN_ID_TAG)){
					initValueReading();
				} else if (name.equals(CLASS_TAG)) {
					initValueReading();
				} else if (name.equals(SYMOP_TAG)) {
					initValueReading();
				} else if (name.equals(RXX_TAG)) {
					initValueReading();
				} else if (name.equals(RXY_TAG)) {
					initValueReading();
				} else if (name.equals(RXZ_TAG)) {
					initValueReading();
				} else if (name.equals(RYX_TAG)) {
					initValueReading();
				} else if (name.equals(RYY_TAG)) {
					initValueReading();
				} else if (name.equals(RYZ_TAG)) {
					initValueReading();
				} else if (name.equals(RZX_TAG)) {
					initValueReading();
				} else if (name.equals(RZY_TAG)) {
					initValueReading();
				} else if (name.equals(RZZ_TAG)) {
					initValueReading();
				} else if (name.equals(TX_TAG)) {
					initValueReading();
				} else if (name.equals(TY_TAG)) {
					initValueReading();
				} else if (name.equals(TZ_TAG)) {
					initValueReading();
				} else if (name.equals(RESIDUE_TAG)){
					inResidue = true;
					currentResidue = new Residue();
				} 
				if (inResidue) {
					if (name.equals(SER_NO_TAG)) {
						initValueReading();
					} else if (name.equals(SEQ_NUM_TAG)) {
						initValueReading();
					} else if (name.equals(NAME_TAG)) {
						initValueReading();
					} else if (name.equals(ASA_TAG)) {
						initValueReading();
					} else if (name.equals(BSA_TAG)) {
						initValueReading();
					//} else if (name.equals(SOLV_EN_TAG)) {
					//	initValueReading();
					}
				}
			}

		} 
	}

	public void endElement(String uri, String localName, String name)
			throws SAXException {

		if (name.equals(PDB_ENTRY_TAG)) {
			inEntry = false;
			allInterfaces.put(currentPdbCode, currentInterfaces);
		}
		if (inEntry) {
			if (name.equals(PDB_CODE_TAG)) {
				currentPdbCode = flushValue().toLowerCase();
			} else if (name.equals(INTERFACE_TAG)){
				inInterface = false;
				currentInterfaces.addInterface(currentPisaInterface);
			}
			if (inInterface && !inMolecule) {
				if (name.equals(ID_TAG)) {
					currentPisaInterface.setId(Integer.parseInt(flushValue()));
				} else if (name.equals(TYPE_TAG)) {
					//currentPisaInterface.setType(Integer.parseInt(flushValue()));
				} else if (name.equals(INT_AREA_TAG)){
					currentPisaInterface.setInterfaceArea(Double.parseDouble(flushValue()));
				} else if (name.equals(INT_SOLV_EN_TAG)){
					currentPisaInterface.setScore(Double.parseDouble(flushValue()));
				//} else if (name.equals(PVALUE_TAG)){
				//	currentPisaInterface.setSolvEnergyPvalue(Double.parseDouble(flushValue()));
				} 
				
			}

			if (inInterface && inMolecule) {
				if (name.equals(MOLECULE_TAG)){
					inMolecule = false;
					if (currentMolecId==1){
						currentPisaInterface.setFirstMolecule(currentPisaMolecule);
						currentPisaInterface.setFirstMolType(currentMolType); 
						currentPisaInterface.setFirstTransfOrth(currentTransfOrth);
						currentPisaInterface.setFirstTransf(currentTransf);
					} else if (currentMolecId==2){
						currentPisaInterface.setSecondMolecule(currentPisaMolecule);
						currentPisaInterface.setSecondMolType(currentMolType);
						currentPisaInterface.setSecondTransfOrth(currentTransfOrth);
						currentPisaInterface.setSecondTransf(currentTransf);
					} else {
						System.err.println("Warning: molecule with id other than 1 or 2 in PISA XML interfaces description");
					}
				}

				if (name.equals(ID_TAG)) {
					currentMolecId = Integer.parseInt(flushValue());
				} else if (name.equals(CHAIN_ID_TAG)){
					currentPisaMolecule.setPdbChainCode(flushValue());
				} else if (name.equals(CLASS_TAG)) {
					currentMolType = flushValue();
				} else if (name.equals(SYMOP_TAG)) {
					currentTransf = SpaceGroup.getMatrixFromAlgebraic(flushValue());
				} else if (name.equals(RXX_TAG)) {
					currentTransfOrth.m00=Double.parseDouble(flushValue());
				} else if (name.equals(RXY_TAG)) {
					currentTransfOrth.m01=Double.parseDouble(flushValue());
				} else if (name.equals(RXZ_TAG)) {
					currentTransfOrth.m02=Double.parseDouble(flushValue());
				} else if (name.equals(RYX_TAG)) {
					currentTransfOrth.m10=Double.parseDouble(flushValue());
				} else if (name.equals(RYY_TAG)) {
					currentTransfOrth.m11=Double.parseDouble(flushValue());
				} else if (name.equals(RYZ_TAG)) {
					currentTransfOrth.m12=Double.parseDouble(flushValue());
				} else if (name.equals(RZX_TAG)) {
					currentTransfOrth.m20=Double.parseDouble(flushValue());
				} else if (name.equals(RZY_TAG)) {
					currentTransfOrth.m21=Double.parseDouble(flushValue());
				} else if (name.equals(RZZ_TAG)) {
					currentTransfOrth.m22=Double.parseDouble(flushValue());
				} else if (name.equals(TX_TAG)) {
					currentTransfOrth.m03=Double.parseDouble(flushValue());
				} else if (name.equals(TY_TAG)) {
					currentTransfOrth.m13=Double.parseDouble(flushValue());
				} else if (name.equals(TZ_TAG)) {
					currentTransfOrth.m23=Double.parseDouble(flushValue());
				} else if (name.equals(RESIDUE_TAG)) {
					inResidue = false;
					currentPisaMolecule.addResidue(currentResidue);
				}
				if (inResidue) {
					if (name.equals(SER_NO_TAG)) {
						currentResidue.setSerial(Integer.parseInt(flushValue()));
					} else if (name.equals(SEQ_NUM_TAG)) {
						currentResidue.setPdbSerial(flushValue());
					} else if (name.equals(NAME_TAG)) {
						currentResidue.setAaType(AminoAcid.getByThreeLetterCode(flushValue()));
					} else if (name.equals(ASA_TAG)) {
						currentResidue.setAsa(Double.parseDouble(flushValue()));
					} else if (name.equals(BSA_TAG)) {
						currentResidue.setBsa(Double.parseDouble(flushValue()));
					//} else if (name.equals(SOLV_EN_TAG)) {
					//	currentResidue.setSolvEnergy(Double.parseDouble(flushValue()));
					}
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

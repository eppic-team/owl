package owl.core.structure.graphs;

import java.io.*;

import owl.core.structure.AminoAcid;
import owl.core.structure.Pdb;
import owl.core.structure.PdbAsymUnit;
import owl.core.util.FileFormatException;

import edu.uci.ics.jung.graph.util.EdgeType;

/**
 * A RIG loaded from a Casp contact prediction file.
 * @author stehr
 * @data 2007-12-19
 */
public class CaspRRFileRIGraph extends RIGraph {

	private static final long serialVersionUID = 1L;

	/**
	 * Constructs a RIGraph object based on a Casp contact prediction file.
	 * @param fileName the contact prediction file
	 * @throws IOException
	 * @throws FileFormatException
	 */
	public CaspRRFileRIGraph (String fileName) throws IOException, FileFormatException{
		super();
		this.contactType=ProtStructGraph.NO_CONTACT_TYPE;
		this.distCutoff=ProtStructGraph.NO_CUTOFF;
		// we initialise pdbCode, chainCode and pdbChainCode to corresponding constants (empty strings at the moment) since file doesn't specify then
		this.pdbCode=PdbAsymUnit.NO_PDB_CODE;
		this.chainCode=Pdb.NO_CHAIN_CODE;
		this.pdbChainCode=Pdb.NO_PDB_CHAIN_CODE;
		
		readFromCaspRRFile(fileName);  
		
	}

	private void readFromCaspRRFile(String fileName) throws FileNotFoundException, IOException, FileFormatException {
		CaspRRFileData rrData = new CaspRRFileData();
		File inFile = new File(fileName);
		if(!CaspRRFileData.isValidCaspRRFile(inFile)) throw new FileFormatException(fileName + " is not a valid casp contact prediction file.");
		rrData.readFromFile(new File(fileName));
		
		// convert meta data
		this.sequence = rrData.getSequence();
		this.contactType = CaspRRFileData.DEFAULT_CONT_TYPE;	// i.e. Cb
		this.groupNum = rrData.groupNum;
		this.targetNum = rrData.targetNum;
		this.caspModelNum = rrData.modelNum;
		
		// add vertices
		if (sequence == null || sequence.equals("")) throw new FileFormatException("No sequence found in CASP RR file " + fileName);
		this.fullLength = sequence.length();
		for (int i=0;i<sequence.length();i++){
			RIGNode node = new RIGNode(i+1,AminoAcid.one2three(sequence.charAt(i)));
			this.addVertex(node); // this takes care of updating the serials2nodes map
		}

		// set distCutoff
		CaspRRFileData.RRContact firstContact = rrData.getContacts().get(0);
		this.distCutoff = firstContact.maxDist;
		
		// add edges
		for(CaspRRFileData.RRContact cont:rrData.getContacts()) {
			if(cont.minDist != 0) 
				throw new FileFormatException("Non-zero minimum distance value in CASP RR file " + fileName + ".");		
			if(this.distCutoff != cont.maxDist) 
				throw new FileFormatException("Distance cutoffs in CASP RR file " + fileName + " are not equal.");
			if (cont.i>=cont.j) 
				throw new FileFormatException("Contact "+cont.i+"-"+cont.j+" specified with i>j in CASP RR file " + fileName+". Only j>i contacts are allowed in CASP RR files.");
			this.addEdge(new RIGEdge(cont.weight), getNodeFromSerial(cont.i), getNodeFromSerial(cont.j), EdgeType.UNDIRECTED);
		}
	}
	
}

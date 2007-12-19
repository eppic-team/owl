package proteinstructure;

import java.io.*;
import java.util.*;

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
	 * @throws GraphFileFormatError
	 */
	public CaspRRFileRIGraph (String fileName) throws IOException, GraphFileFormatError{
		super();
		// we set the sequence to blank when we read from file as we don't have the full sequence
		// if sequence is present in RR file then is read from there
		this.sequence="";
		this.contactType="";
		this.distCutoff=0.0;
		// we initialise pdbCode, chainCode and pdbChainCode to empty strings since file doesn't specify then
		// TODO: should these be set to null instead?
		this.pdbCode="";
		this.chainCode="";
		this.pdbChainCode="";
		
		readFromCaspRRFile(fileName);  
		
	}

	private void readFromCaspRRFile(String fileName) throws FileNotFoundException, IOException, GraphFileFormatError {
		CaspRRFileData rrData = new CaspRRFileData();
		File inFile = new File(fileName);
		if(!CaspRRFileData.isValidCaspRRFile(inFile)) throw new GraphFileFormatError(fileName + " is not a valid casp contact prediction file.");
		rrData.readFromFile(new File(fileName));
		
		// convert meta data
		this.sequence = rrData.getSequence();
		this.contactType = CaspRRFileData.DEFAULT_CONT_TYPE;	// i.e. Cb
		this.groupNum = rrData.groupNum;
		this.targetNum = rrData.targetNum;
		this.caspModelNum = rrData.modelNum;
		
		// add vertices
		serials2nodes = new TreeMap<Integer,RIGNode>();
		if (sequence == null || sequence.equals("")) throw new GraphFileFormatError("No sequence found in " + fileName);
		this.fullLength = sequence.length();
		for (int i=0;i<sequence.length();i++){
			String letter = String.valueOf(sequence.charAt(i));
			RIGNode node = new RIGNode(i+1,AAinfo.oneletter2threeletter(letter));
			serials2nodes.put(i+1, node);
			this.addVertex(node);
		}

		// set distCutoff
		CaspRRFileData.RRContact firstContact = rrData.getContacts().get(0);
		this.distCutoff = firstContact.maxDist;
		
		// add edges
		for(CaspRRFileData.RRContact cont:rrData.getContacts()) {
			if(cont.minDist != 0) throw new GraphFileFormatError("Non-zero minimum distance value in " + fileName + ".");		
			if(this.distCutoff != cont.maxDist) throw new GraphFileFormatError("Distance cutoffs in " + fileName + " are not equal.");
			this.addEdge(new RIGEdge(cont.weight), serials2nodes.get(cont.i), serials2nodes.get(cont.j), EdgeType.UNDIRECTED);
		}
	}
	
}

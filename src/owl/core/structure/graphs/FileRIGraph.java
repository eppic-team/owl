package owl.core.structure.graphs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.AminoAcid;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.util.FileFormatException;

import edu.uci.ics.jung.graph.util.EdgeType;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * A RIGraph derived from a single chain pdb protein structure loaded from a 
 * graph file in OWL format
 * 
 */
public class FileRIGraph extends RIGraph {
	
	private static final long serialVersionUID = 1L;

	private static double DEFAULT_WEIGHT = 1.0;
	
	private boolean simple;
	
	/**
	 * Constructs RIGraph object by reading a file with contacts in OWL format
	 * If the contacts file doesn't have the sequence then the RIGraph object won't have sequence or residue types in RIGNodes
	 * @param contactsfile
	 * @throws IOException
	 * @throws FileFormatException
	 */
	public FileRIGraph (String contactsfile) throws IOException, FileFormatException{
		this(contactsfile, false);
	}

	/**
	 * Constructs RIGraph object by reading a file with contacts
	 * If the contacts file doesn't have the sequence then the RIGraph object won't have sequence or residue types in RIGNodes
	 * @param contactsfile
	 * @param simple true if we want to read file in "simple" format, i.e. only a list of edges, false if file in OWL format
	 * @throws IOException
	 * @throws FileFormatException
	 */
	public FileRIGraph (String contactsfile, boolean simple) throws IOException, FileFormatException{
		super();
		this.simple = simple;
		this.contactType=ProtStructGraph.NO_CONTACT_TYPE;
		this.distCutoff=ProtStructGraph.NO_CUTOFF;
		// we initialise pdbCode, chainCode and pdbChainCode to corresponding constants (empty strings at the moment) in case the file doesn't specify then
		this.pdbCode=PdbAsymUnit.NO_PDB_CODE;
		this.chainCode=PdbChain.NO_CHAIN_CODE;
		this.pdbChainCode=PdbChain.NO_PDB_CHAIN_CODE;
		
		readFromFile(contactsfile);  
		
	}
	
	/**
	 * Parses the graph file reading identifiers, sequence and edges
	 * @param contactsfile
	 * @throws IOException
	 * @throws FileFormatException if file does not start with #AGLAPPE or #CMVIEW or #OWL 
	 * if file format not the right version, if sequence is not present
	 */
	private void readFromFile (String contactsfile) throws IOException, FileFormatException {
		HashMap<Pair<Integer>,Double> contacts2weights = new HashMap<Pair<Integer>,Double>();
		HashSet<Integer> allserials = new HashSet<Integer>();
		BufferedReader fcont = new BufferedReader(new FileReader(new File(contactsfile)));
		int linecount=0;
		String line;
		while ((line = fcont.readLine() ) != null ) {
			linecount++;
			if (!simple) {
				Pattern p = Pattern.compile("^#(?:AGLAPPE|CMVIEW|OWL).*ver: (\\d\\.\\d)");
				Matcher m = p.matcher(line);
				if (m.find()){
					if (!m.group(1).equals(GRAPHFILEFORMATVERSION)){
						fcont.close();
						throw new FileFormatException("Contact map file "+contactsfile+" has a wrong file format version. Supported version is "+GRAPHFILEFORMATVERSION+" and found version was "+m.group(1));
					}
				} else if (linecount==1){ // #AGLAPPE/#CMVIEW/#OWL not found in first line
					fcont.close();
					throw new FileFormatException(contactsfile+" is not a valid contact map file");
				}
				Pattern ps = Pattern.compile("^#SEQUENCE:\\s*(\\w+)$");
				Matcher ms = ps.matcher(line);
				if (ms.find()){
					sequence=ms.group(1);
				}
				ps = Pattern.compile("^#PDB:\\s*(\\w+)");
				ms = ps.matcher(line);
				if (ms.find()){
					pdbCode=ms.group(1);
				}
				ps = Pattern.compile("^#PDB CHAIN CODE:\\s*(\\w+)");
				ms = ps.matcher(line);
				if (ms.find()){
					pdbChainCode=ms.group(1);
				}
				ps = Pattern.compile("^#CHAIN:\\s*(\\w)");
				ms = ps.matcher(line);
				if (ms.find()){
					chainCode=ms.group(1);
				}
				ps = Pattern.compile("^#MODEL:\\s*(\\d+)");
				ms = ps.matcher(line);
				if (ms.find()){
					model=Integer.parseInt(ms.group(1));
				}				
				ps = Pattern.compile("^#CT:\\s*([a-zA-Z/]+)");
				ms = ps.matcher(line);
				if (ms.find()){
					contactType=ms.group(1);
				}												
				ps = Pattern.compile("^#CUTOFF:\\s*(\\d+\\.\\d+)");
				ms = ps.matcher(line);
				if (ms.find()){
					distCutoff=Double.parseDouble(ms.group(1));
				}								
			}
			Pattern pcontact = Pattern.compile("^\\s*(\\d+)\\s+(\\d+)(?:\\s+(\\d+(?:\\.\\d+)?))?\\s*$");
			Matcher mcontact = pcontact.matcher(line);
			if (mcontact.find()){
				int i = Integer.valueOf(mcontact.group(1));
				int j = Integer.valueOf(mcontact.group(2));
				allserials.add(i);
				allserials.add(j);
				double weight = DEFAULT_WEIGHT;
				if (mcontact.group(3)!=null) {
					weight = Double.valueOf(mcontact.group(3));
				}
				contacts2weights.put(new Pair<Integer>(i,j),weight);
			}

		}
		fcont.close();
		
		if (simple) {
			int maxResSer = Collections.max(allserials);
			sequence = "";
			for (int resser=1;resser<=maxResSer;resser++) {
				sequence+=AminoAcid.XXX.getOneLetterCode();
			}
		}
		
		// if we didn't find a sequence we throw format exception
		if (sequence==null) {
			throw new FileFormatException("No sequence present in contact map file "+contactsfile);
		}

		// populating this RIGraph with nodes and setting fullLength
		this.fullLength = sequence.length();
		for (int serial:allserials){
			if (serial>sequence.length()) 
				throw new FileFormatException("AaResidue serial "+serial+" found in edges list of contact map file "+contactsfile+" is bigger than length of sequence");
			RIGNode node = new RIGNode(serial,AminoAcid.one2three(sequence.charAt(serial-1)));
			this.addVertex(node);
		}

		//TODO we still use here DIRECTED as default for "/", eventually this should change by taking another parameter "boolean directed", so "/" could have DIRECTED/UNDIRECTED versions
		EdgeType et = EdgeType.UNDIRECTED;
		if (contactType.contains("/")){
			et = EdgeType.DIRECTED;
		}
		// populating this RIGraph with  RIGEdges
		for (Pair<Integer> resPair:contacts2weights.keySet()){
			//TODO we are reading the 3rd column as weights (not as atom weights), we might need to change this or read a 4th column or something else
			this.addEdge(new RIGEdge(contacts2weights.get(resPair)), getNodeFromSerial(resPair.getFirst()), getNodeFromSerial(resPair.getSecond()), et);
		}
		
		
	}

}

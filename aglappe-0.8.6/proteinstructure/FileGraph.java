package proteinstructure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A residue interaction graph derived from a single chain pdb protein structure loaded from a graph file in aglappe's format
 * 
 * @author 		Jose Duarte
 * Class:		FileGraph
 * Package:		proteinstructure
 */
public class FileGraph extends Graph {
	
	private static double DEFAULT_WEIGHT = 1.0;
	
	/**
	 * Constructs Graph object by reading a file with contacts
	 * If the contacts file doesn't have the sequence then the graph object won't have sequence or nodes
	 * That means it won't be possible to get a ContactMap from it using getCM because CM needs both sequence and nodes
	 * @param contactsfile
	 * @throws IOException
	 * @throws FileNotFoundException
	 * @throws GraphFileFormatError
	 */
	public FileGraph (String contactsfile) throws IOException, FileNotFoundException, GraphFileFormatError{
		// we set the sequence to blank when we read from file as we don't have the full sequence
		// if sequence is present in contactsfile then is read from there
		this.sequence="";
		this.ct="";
		this.cutoff=0.0;
		// we initialise pdbCode, chainCode and pdbChainCode to empty strings in case the file doesn't specify then
		this.pdbCode="";
		this.chainCode="";
		this.pdbChainCode="";
		this.directed=false;

		read_graph_from_file(contactsfile); // initialises contacts, and nodes (only if sequence is given)
		
		if (ct.contains("/")){
			directed=true;
		}

		if (!sequence.equals("")){
			this.fullLength=sequence.length();
			this.obsLength=nodes.size(); 
		} else { 
			// if contacts have correct residue numbering then this should get the right full length up to the maximum node that makes a contact,
			// we will miss: nodes without contacts at the end of sequence and gaps (unobserved residues) at the end of the sequence.
			// We don't know more without nodes and sequence
			this.fullLength=contacts.getMaxNode();
			// in this case nodes has not been initialised so we set obsLength=fullLength as we don't have the information
			this.obsLength=fullLength;  
		}
		this.numContacts=contacts.size();
		this.modified=false;
	}

	private void read_graph_from_file (String contactsfile) throws FileNotFoundException, IOException, GraphFileFormatError {
		contacts = new EdgeSet();
		//System.out.println("Reading contacts from file "+contactsfile);
		BufferedReader fcont = new BufferedReader(new FileReader(new File(contactsfile)));
		int linecount=0;
		String line;
		while ((line = fcont.readLine() ) != null ) {
			linecount++;
			Pattern p = Pattern.compile("^#AGLAPPE.*ver: (\\d\\.\\d)");
			Matcher m = p.matcher(line);
			if (m.find()){
				if (!m.group(1).equals(GRAPHFILEFORMATVERSION)){
					throw new GraphFileFormatError("The graph file "+contactsfile+" can't be read, wrong file format version. Supported version is "+GRAPHFILEFORMATVERSION+" and found version was "+m.group(1));
				}
			} else if (linecount==1){ // #AGLAPPE not found and in first line
				throw new GraphFileFormatError("The graph file "+contactsfile+" can't be read, wrong file format");
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
			ps = Pattern.compile("^#CT:\\s*([a-zA-Z/]+)");
			ms = ps.matcher(line);
			if (ms.find()){
				ct=ms.group(1);
			}												
			ps = Pattern.compile("^#CUTOFF:\\s*(\\d+\\.\\d+)");
			ms = ps.matcher(line);
			if (ms.find()){
				cutoff=Double.parseDouble(ms.group(1));
			}								

			Pattern pcontact = Pattern.compile("^\\s*(\\d+)\\s+(\\d+)(?:\\s+(\\d+\\.\\d+))?\\s*$");
			Matcher mcontact = pcontact.matcher(line);
			if (mcontact.find()){
				int i = Integer.valueOf(mcontact.group(1));
				int j = Integer.valueOf(mcontact.group(2));
				double weight = DEFAULT_WEIGHT;
				if (mcontact.group(3)!=null) {
					weight = Double.valueOf(mcontact.group(3));
				}
				Edge cont = new Edge(i,j,weight);
				contacts.add(cont);
			}

		}
		fcont.close();
		// if sequence was given we take nodes from it
		nodes = new TreeMap<Integer, String>();
		for (int i=0;i<sequence.length();i++){
			String letter = String.valueOf(sequence.charAt(i));
			nodes.put(i+1, AAinfo.oneletter2threeletter(letter));
		}		

	}

}

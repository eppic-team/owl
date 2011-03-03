import gnu.getopt.Getopt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbasePdb;
import owl.core.util.MySQLConnection;



public class dumpseq {
	/*------------------------------ constants ------------------------------*/
	
	private static final String			PDB_DB = "pdbase";

	private static final String IDS_REGEX1 = "^(\\d\\w\\w\\w)(\\w)";
	private static final String IDS_REGEX2 = "^(\\d\\w\\w\\w)_(\\w)";
	private static final String IDS_REGEX3 = "^(\\d\\w\\w\\w)\\s+(\\w)";
	private static final String IDS_REGEX4 = "^(\\d\\w\\w\\w)";

	
	public static void main(String[] args) throws IOException {
		
		String progName = "dumpseq";
		
		String help = "Usage, 2 options:\n" +
				"1)  "+progName+" -i <listfile> \n" +
				"2)  "+progName+" -p <pdbCode+chainCode> \n\n" +
				" -i <file>       : file with list of pdbCodes+chainCodes\n"+
				" -p <string>     : comma separated list of pdbCodes+chainCodes, e.g. -p 1bxyA,1josA\n" +
				"                   If only pdbCode (no chainCode) specified, e.g. 1bxy then first chain will be taken\n"+
				" [-o] <dir>      : outputs one file per sequence in given directory. Default: current dir\n" +
				" [-f] <file>     : one output file for all sequences\n"+
				" [-s]            : outputs to stdout instead of file(s)\n" +
				" [-N]            : don't print FASTA header, just raw sequence\n"+
				" [-D] <database> : pdbase database name. Default: "+PDB_DB+"\n\n";				

		String listfile = "";
		String[] pdbIds = null;
		String pdbaseDb = PDB_DB;
		String outputDir = ".";
		File oneOutputFile = null;
		boolean stdout = false;
		boolean fastaHeader = true;
		
		Getopt g = new Getopt(progName, args, "i:p:o:f:D:sNh?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'i':
				listfile = g.getOptarg();
				break;
			case 'p':
				pdbIds = g.getOptarg().split(",");
				break;
			case 'o':
				outputDir = g.getOptarg();
				break;
			case 'f':
				oneOutputFile = new File(g.getOptarg());
				break;				
			case 'D':
				pdbaseDb = g.getOptarg();
				break;
			case 's':
				stdout = true;
				break;
			case 'N':
				fastaHeader = false;
				break;								
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		if (listfile.equals("") && pdbIds==null){
			System.err.println("Either a listfile or some pdb codes/chain codes must be given");
			System.err.println(help);
			System.exit(1);
		}
		if (!listfile.equals("") && pdbIds!=null) {
			System.err.println("Options -p and -i are exclusive. Use only one of them");
			System.err.println(help);
			System.exit(1);			
		}
		
		MySQLConnection conn = null;		

		try{
			conn = new MySQLConnection();
		} catch (Exception e) {
			System.err.println("Error opening database connection. Exiting");
			System.exit(1);
		}


		if (!listfile.equals("")) {	
			pdbIds = readIdsListFile(new File(listfile));
		}

		int numPdbs = 0;

		PrintStream Out = null;
		if (stdout) {
			Out = System.out;
		} else if (oneOutputFile!=null) {
			Out = new PrintStream(new FileOutputStream(oneOutputFile));
		} 

		for (int i=0;i<pdbIds.length;i++) {
			String pdbCode = null;
			String[] pdbChainCodes = null;
			if (pdbIds[i].length()==4) {
				pdbCode = pdbIds[i];
			} else if (pdbIds[i].length()==5){
				pdbCode = pdbIds[i].substring(0, 4);
				pdbChainCodes = new String[1];
				pdbChainCodes[0] = pdbIds[i].substring(4);
			} else {
				System.err.println("The string "+pdbIds[i]+" doesn't look like a PDB id. Skipping");
				continue;
			}
			

			Pdb pdb = null;
			try {
				pdb = new PdbasePdb(pdbCode, pdbaseDb, conn);
				if (pdbChainCodes==null) {
					pdbChainCodes = pdb.getChains();
				}
			} catch (PdbCodeNotFoundException e) {
				System.err.println("Couldn't find pdb code "+pdbCode);
				continue;
			} catch (SQLException e) {
				System.err.println("SQL error for structure "+pdbCode+", error: "+e.getMessage());
				continue;
			} catch (PdbLoadException e) {
				System.err.println("Error loading pdb data for " + pdbCode +", specific error: "+e.getMessage());
				continue;
			}

			for (String pdbChainCode:pdbChainCodes) {

				try {
					pdb.load(pdbChainCode);

					String sequence = pdb.getSequence();

					File outputFile = new File(outputDir,pdbCode+pdbChainCode+".fasta");

					if (!stdout && oneOutputFile==null) {
						Out = new PrintStream(new FileOutputStream(outputFile.getAbsolutePath()));
					}

					if (fastaHeader) { 
						Out.println(">"+pdbCode+pdbChainCode);
					}

					Out.println(sequence);

					if (!stdout && oneOutputFile==null) {
						Out.close();
					}

					if (!stdout) { // if output of sequence is stdout, then we don't want to print anything else to stdout
						System.out.println("Wrote "+pdbCode+pdbChainCode+".fasta");
					}

					numPdbs++;
				} catch (PdbLoadException e) {
					System.err.println("Error loading pdb data for " + pdbCode + pdbChainCode+", specific error: "+e.getMessage());
				}

			}




		}

		if (!stdout && oneOutputFile!=null) {
			Out.close();
		}
		
		// output results
		if (!stdout) { // if output of sequence is stdout, then we don't want to print anything else to stdout
			System.out.println("Number of dumped sequences: " + numPdbs);
		}


	} 
		

	/**
	 * Reads a list file containing a list of pdb codes and chain codes in 3 possible formats:
	 * - 1 column pdbCodes+chainCodes, e.g. 1bxyA
	 * - 1 column underscore-separated pdbCodes and chainCodes, e.g. 1bxy_A
	 * - 2 colums tab/spaces-separated pdbCodes and chainCodes, e.g. 1bxy A or 1bxy   A
	 * See the IDS_REGEX constants of this class for the regex that we are using.
	 * A mix of the formats is also tolerated.
	 * Chain codes can only be a 1 letter code (so we must use an "A" for NULL codes)
	 * @param listFile
	 * @return an array of pdbCodes(lower case)+chainCodes(conserving case) in the format 1bxyA or pdbCodes only e.g. 1bxy
	 * @throws IOException
	 */
	public static String[] readIdsListFile(File listFile) throws IOException {
		ArrayList<String> codesAL = new ArrayList<String>(); 

		BufferedReader fileIn = new BufferedReader(new FileReader(listFile));
		String line;
		int lineCount=0;
		while((line = fileIn.readLine()) != null) {
			lineCount++;
			if (line.length()!=0 && !line.startsWith("#")) {
				Pattern p1 = Pattern.compile(IDS_REGEX1);
				Matcher m1 = p1.matcher(line);
				Pattern p2 = Pattern.compile(IDS_REGEX2);
				Matcher m2 = p2.matcher(line);
				Pattern p3 = Pattern.compile(IDS_REGEX3);
				Matcher m3 = p3.matcher(line);				
				Pattern p4 = Pattern.compile(IDS_REGEX4);
				Matcher m4 = p4.matcher(line);
				
				if (m1.matches()) {
					codesAL.add(m1.group(1).toLowerCase()+m1.group(2));
				} 
				else if (m2.matches()) {
					codesAL.add(m2.group(1).toLowerCase()+m2.group(2));
				} 
				else if (m3.matches()){
					codesAL.add(m3.group(1).toLowerCase()+m3.group(2));
				}
				else if (m4.matches()) {
					codesAL.add(m4.group(1));
				}
				else {
					System.err.println("Line "+lineCount+" in list file "+listFile+" is not in any of the recognised pdbCode+chainCode formats");
				}
			
			}
		}
		String[] codes = new String[codesAL.size()];
		codesAL.toArray(codes);
		return codes;
	}

}

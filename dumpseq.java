import gnu.getopt.Getopt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;

import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbasePdb;
import tools.MySQLConnection;


public class dumpseq {
	/*------------------------------ constants ------------------------------*/
	
	public static final String			PDB_DB = "pdbase";
	public static final String			DB_HOST = "white";								
	public static final String			DB_USER = MySQLConnection.getUserName();
	public static final String			DB_PWD = "nieve";

	public static final String			GAP_CHARACTER = "-";
	
	public static void main(String[] args) throws IOException {
		
		String progName = "dumpseq";
		
		String help = "Usage, 3 options:\n" +
				"1)  "+progName+" -i <listfile> [-o <output_dir> | -f <one_output_file> | -s] [-N] [-D <pdbase_db>] \n" +
				"2)  "+progName+" -p <pdb_code> -c <chain_pdb_code> [-o <output_dir> | -f <one_output_file> | -s] [-N] [-D <pdbase_db>] \n" +
				"Output options: -o one file per sequence, -f one file for all sequences, -s standard output. With -N no FASTA headers will be written\n"+
				"In case 2) also a list of comma separated pdb codes and chain codes can be specified, e.g. -p 1bxy,1jos -c A,A\n" +
				"If pdbase_db not specified, the default pdbase will be used\n"; 

		String listfile = "";
		String[] pdbCodes = null;
		String[] pdbChainCodes = null;
		String pdbaseDb = PDB_DB;
		String outputDir = "";
		File oneOutputFile = null;
		boolean stdout = false;
		boolean fastaHeader = true;
		
		Getopt g = new Getopt(progName, args, "i:p:c:o:f:D:sNh?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'i':
				listfile = g.getOptarg();
				break;
			case 'p':
				pdbCodes = g.getOptarg().split(",");
				break;
			case 'c':
				pdbChainCodes = g.getOptarg().split(",");
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

		if (listfile.equals("") && pdbCodes==null){
			System.err.println("Either a listfile or some pdb codes/chain codes must be given");
			System.err.println(help);
			System.exit(1);
		}
		if (!listfile.equals("") && pdbCodes!=null) {
			System.err.println("Options -p/-c and -i are exclusive. Use only one of them");
			System.err.println(help);
			System.exit(1);			
		}

		
		MySQLConnection conn = null;		

		try{
			conn = new MySQLConnection(DB_HOST, DB_USER, DB_PWD);
		} catch (Exception e) {
			System.err.println("Error opening database connection. Exiting");
			System.exit(1);
		}




		if (!listfile.equals("")) {			
			BufferedReader fpdb = new BufferedReader(new FileReader(listfile));
			String line = "";
			int numLines = 0;
			fpdb.mark(100000);
			while ((line = fpdb.readLine() ) != null ) {
				if (line.length()>0) numLines++;
			}
			fpdb.reset();
			pdbCodes = new String[numLines];
			pdbChainCodes = new String[numLines];
			numLines = 0;
			while ((line = fpdb.readLine() ) != null ) {
				pdbCodes[numLines] = line.split("\\s+")[0].toLowerCase();
				pdbChainCodes[numLines] = line.split("\\s+")[1];
				numLines++;
			}
		}

		int numPdbs = 0;

		PrintStream Out = null;
		if (stdout) {
			Out = System.out;
		} else if (oneOutputFile!=null) {
			Out = new PrintStream(new FileOutputStream(oneOutputFile));
		} 

		for (int i=0;i<pdbCodes.length;i++) {
			String pdbCode = pdbCodes[i];
			String pdbChainCode = pdbChainCodes[i];

			try {

				Pdb pdb = new PdbasePdb(pdbCode, pdbaseDb, conn);
				pdb.load(pdbChainCode);
				
				String sequence = pdb.getSequence();

				File outputFile = new File(outputDir,pdbCode+"_"+pdbChainCode+".fasta");
				
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
					System.out.println("Wrote "+pdbCode+"_"+pdbChainCode+".fasta");
				}

				numPdbs++;

			} catch (PdbLoadError e) {
				System.err.println("Error loading pdb data for " + pdbCode + pdbChainCode+", specific error: "+e.getMessage());
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Couldn't find pdb code "+pdbCode);
			} catch (SQLException e) {
				System.err.println("SQL error for structure "+pdbCode+"_"+pdbChainCode+", error: "+e.getMessage());
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
		

}

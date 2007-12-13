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
		
		
		String help = "Usage, 3 options:\n" +
				"1)  genGraph -i <listfile> -o <output_dir> [-D <pdbase_db>] \n" +
				"2)  genGraph -p <pdb_code> -c <chain_pdb_code> -o <output_dir> [-D <pdbase_db>] \n" +
				"In case 2) also a list of comma separated pdb codes and chain codes can be specified, e.g. -p 1bxy,1jos -c A,A\n" +
				"If pdbase_db not specified, the default pdbase will be used\n"; 

		String listfile = "";
		String[] pdbCodes = null;
		String[] pdbChainCodes = null;
		String pdbaseDb = PDB_DB;
		String outputDir = "";
		
		Getopt g = new Getopt("genGraph", args, "i:p:c:o:D:h?");
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
			case 'D':
				pdbaseDb = g.getOptarg();
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

		for (int i=0;i<pdbCodes.length;i++) {
			String pdbCode = pdbCodes[i];
			String pdbChainCode = pdbChainCodes[i];

			try {

				Pdb pdb = new PdbasePdb(pdbCode, pdbaseDb, conn);
				pdb.load(pdbChainCode);
				
				String sequence = pdb.getSequence();

				File outputFile = new File(outputDir,pdbCode+"_"+pdbChainCode+".fasta");
				
				PrintStream Out = new PrintStream(new FileOutputStream(outputFile.getAbsolutePath()));
				Out.println(">"+pdbCode+"_"+pdbChainCode);
//				for (int pos=1;pos<=sequence.length();pos++) {
//					if (pos%10==0){
//						Out.printf("%10d",pos/10);
//					} 
//				}
//				Out.println();
//				for (int pos=1;pos<=sequence.length();pos++) {
//					Out.print(pos%10);
//				}
//				Out.println();
				Out.println(sequence);
//				for (int pos=1;pos<=sequence.length();pos++) {
//					if (pdb.hasCoordinates(pos)) {
//						Out.print(sequence.charAt(pos-1));
//					} else {
//						Out.print(GAP_CHARACTER);
//					}
//				}
				Out.close();
				
				System.out.println("Wrote "+outputFile.getAbsolutePath());

				numPdbs++;

			} catch (PdbLoadError e) {
				System.err.println("Error loading pdb data for " + pdbCode + pdbChainCode+", specific error: "+e.getMessage());
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Couldn't find pdb code "+pdbCode);
			} catch (SQLException e) {
				System.err.println("SQL error for structure "+pdbCode+"_"+pdbChainCode+", error: "+e.getMessage());
			}

		}

		// output results
		System.out.println("Number of dumped sequences: " + numPdbs);


	} 
		

}

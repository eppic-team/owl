import gnu.getopt.Getopt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;

import proteinstructure.Pdb;
import proteinstructure.PdbChainCodeNotFoundError;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbaseInconsistencyError;
import proteinstructure.PdbasePdb;
import proteinstructure.PdbfileFormatError;
import proteinstructure.PdbfilePdb;
import proteinstructure.RIGraph;
import tools.MySQLConnection;


public class genGraph {
	/*------------------------------ constants ------------------------------*/
	
	public static final String			PDB_DB = "pdbase";
	public static final String			DB_HOST = "white";								
	public static final String			DB_USER = getUserName();
	public static final String			DB_PWD = "nieve";
	public static final String			DSSP_EXE = "/project/StruPPi/bin/dssp";
	public static final String			DSSP_PARAMS = "--";
	public static final String			NACCESS_EXE = "/project/StruPPi/bin/naccess";
	public static final String			NACCESS_PARAMS = "";

	//public static double			cutoff = 4.2;
	//public static String			edgeType = "ALL";
	
	/*---------------------------- private methods --------------------------*/
	/** 
	 * Get user name from operating system (for use as database username). 
	 * */
	private static String getUserName() {
		String user = null;
		user = System.getProperty("user.name");
		if(user == null) {
			System.err.println("Could not get user name from operating system.");
		}
		return user;
	}
	
	public static void main(String[] args) throws IOException {
		
		
		String help = "Usage, 3 options:\n" +
				"1)  genGraph -i <listfile> -d <distance_cutoff> -t <contact_type> -o <output_dir> [-D <pdbase_db>] \n" +
				"2)  genGraph -p <pdb_code> -c <chain_pdb_code> -d <distance_cutoff> -t <contact_type> -o <output_dir> [-D <pdbase_db>] \n" +
				"3)  genGraph -f <pdbfile> -c <chain_pdb_code> -d <distance_cutoff> -t <contact_type> -o <output_dir> \n" +
				"In case 2) also a list of comma separated pdb codes and chain codes can be specified, e.g. -p 1bxy,1jos -c A,A\n" +
				"If pdbase_db not specified, the default pdbase will be used\n"; 

		String listfile = "";
		String[] pdbCodes = null;
		String[] pdbChainCodes = null;
		String pdbfile = "";
		String pdbaseDb = PDB_DB;
		String edgeType = "";
		double cutoff = 0.0;
		String outputDir = "";
		
		Getopt g = new Getopt("genGraph", args, "i:p:c:f:d:t:o:D:h?");
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
			case 'f':
				pdbfile = g.getOptarg();
				break;
			case 'd':
				cutoff = Double.valueOf(g.getOptarg());
				break;
			case 't':
				edgeType = g.getOptarg();
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

		if (outputDir.equals("") || edgeType.equals("") || cutoff==0.0) {
			System.err.println("Some missing option");
			System.err.println(help);
			System.exit(1);
		}
		if (listfile.equals("") && pdbCodes==null && pdbfile.equals("")){
			System.err.println("Either a listfile, some pdb codes/chain codes or a pdbfile must be given");
			System.err.println(help);
			System.exit(1);
		}
		if ((!listfile.equals("") && pdbCodes!=null) || (!listfile.equals("") && !pdbfile.equals("")) || (pdbCodes!=null && !pdbfile.equals(""))) {
			System.err.println("Options -p/-c, -i and -f/-c are exclusive. Use only one of them");
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
		
		
		if (pdbfile.equals("")){
			
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
					
					long start = System.currentTimeMillis();
					
					Pdb pdb = new PdbasePdb(pdbCode, pdbChainCode, pdbaseDb, conn);

					// get graph
					RIGraph graph = pdb.get_graph(edgeType, cutoff);

					File outputFile = new File(outputDir,pdbCode+"_"+pdbChainCode+"_"+edgeType+"_"+cutoff+".graph");
					graph.write_graph_to_file(outputFile.getAbsolutePath());

					long end = System.currentTimeMillis();
					double time = (double) (end -start)/1000;

					System.out.println("Wrote "+outputFile.getAbsolutePath());
					System.out.printf("%5.3f s\n",time);
					
					numPdbs++;

				} catch (PdbaseInconsistencyError e) {
					System.err.println("Inconsistency in " + pdbCode + pdbChainCode);
				} catch (PdbCodeNotFoundError e) {
					System.err.println("Couldn't find pdb code "+pdbCode);
				} catch (SQLException e) {
					System.err.println("SQL error for structure "+pdbCode+"_"+pdbChainCode+", error: "+e.getMessage());
				} catch (PdbChainCodeNotFoundError e) {
					System.err.println("Couldn't find pdb chain code "+pdbChainCode+" for pdb code "+pdbCode);
				}

			}

			// output results
			System.out.println("Number of structures done successfully: " + numPdbs);


		} else {
			String pdbChainCode = pdbChainCodes[0];
			try {
				Pdb pdb = new PdbfilePdb(pdbfile,pdbChainCode);
				if (!pdb.hasSecondaryStructure()) {
					pdb.runDssp(DSSP_EXE, DSSP_PARAMS);
				}
				RIGraph graph = pdb.get_graph(edgeType, cutoff);

				File outputFile = new File(outputDir,pdb.getPdbCode()+"_"+pdbChainCode+"_"+edgeType+"_"+cutoff+".graph");
				graph.write_graph_to_file(outputFile.getAbsolutePath());
				System.out.println("Wrote graph file "+outputFile.getAbsolutePath()+" from pdb file "+pdbfile);
				
			} catch (PdbfileFormatError e) {
				System.err.println("pdb file "+pdbfile+" doesn't have right format");
			} catch (PdbChainCodeNotFoundError e) {
				System.err.println("chain code "+pdbChainCode+" wasn't found in file "+pdbfile);	
			}
		}
	}

}

import gnu.getopt.Getopt;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbasePdb;
import owl.core.structure.PdbfilePdb;
import owl.core.structure.TemplateList;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.MySQLConnection;



public class genGraph {
	/*------------------------------ constants ------------------------------*/
	
	public static final String			PDB_DB = "pdbase";
	private static final int NUM_RANDOM_SUBSETS = 10;
	
	/*--------------------------- type definitions --------------------------*/
	public enum OutputFormat {CMVIEW, PAUL, SADP};
	
	
	// main

	public static void main(String[] args) throws IOException {
		
		String progName = "genGraph";
		
		String help = "Usage, 3 options:\n" +
				"1)  "+progName+" -i <listfile> -d <distance_cutoff> -t <contact_type> [-o <output_dir>] [-D <pdbase_db>] [-C|-P|-S] [-r <percent>] \n" +
				"2)  "+progName+" -p <pdb_code> -d <distance_cutoff> -t <contact_type> [-o <output_dir>] [-D <pdbase_db>] [-C|-P|-S] [-r <percent>] \n" +
				"3)  "+progName+" -f <pdbfile> [-c <chain_pdb_code>] -d <distance_cutoff> -t <contact_type> [-o <output_dir>] [-C|-P|-S] [-r <percent>] \n\n" +
				"In case 2) also a list of comma separated pdb codes+chain codes can be specified, e.g. -p 1bxyA,1josA \n" +
				"In case 3) if -c not specified then the first chain code in the pdb file will be taken\n" +
				"If pdbase_db not specified, the default pdbase will be used\n" +
				"If output dir not specified default is current\n" +
				"\n" +
				"Use -r <percent> to generate random subsets (of size percent) for each of the contact maps calculated.\n" +
				"\n" +
				"In all cases, the output format can be specified with one of the following options:\n" +
				" -C CMView graph file format (default)\n" +
				" -P PAUL format\n" + 
				" -S SADP format\n"; 
				
		String listfile = "";
		String[] pdbIds = null;
		String pdbChainCode4file = null;
		String pdbfile = "";
		String pdbaseDb = PDB_DB;
		String edgeType = "";
		double cutoff = 0.0;
		String outputDir = "."; // we set default to current directory
		OutputFormat outFormat = OutputFormat.CMVIEW;	// output to cmview graph format by default
		double randomSubset = 1;
		
		Getopt g = new Getopt(progName, args, "i:p:c:f:d:t:o:D:CPSr:h?");
		int c;
		while ((c = g.getopt()) != -1) {
			switch(c){
			case 'i':
				listfile = g.getOptarg();
				break;
			case 'p':
				pdbIds = g.getOptarg().split(",");
				break;
			case 'c':
				pdbChainCode4file = g.getOptarg();
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
			case 'C':
				outFormat = OutputFormat.CMVIEW;
				break;	
			case 'P':
				outFormat = OutputFormat.PAUL;				
				break;
			case 'S':
				outFormat = OutputFormat.SADP;
				break;
			case 'r':
				randomSubset = Double.valueOf(g.getOptarg());
				break;
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		if (edgeType.equals("") || cutoff==0.0) {
			System.err.println("Some missing option");
			System.err.println(help);
			System.exit(1);
		}
		if (listfile.equals("") && pdbIds==null && pdbfile.equals("")){
			System.err.println("Either a listfile, some pdb codes/chain codes or a pdbfile must be given");
			System.err.println(help);
			System.exit(1);
		}
		if ((!listfile.equals("") && pdbIds!=null) || (!listfile.equals("") && !pdbfile.equals("")) || (pdbIds!=null && !pdbfile.equals(""))) {
			System.err.println("Options -p, -i and -f/-c are exclusive. Use only one of them");
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
		
		
		if (pdbfile.equals("")){
			
			if (!listfile.equals("")) {		
				pdbIds = TemplateList.readIdsListFile(new File(listfile));
			}

			int numPdbs = 0;

			for (int i=0;i<pdbIds.length;i++) {
				String pdbCode = pdbIds[i].substring(0,4);
				String pdbChainCode = pdbIds[i].substring(4);

				try {
					
					//long start = System.currentTimeMillis();
					
					Pdb pdb = new PdbasePdb(pdbCode, pdbaseDb, conn);
					pdb.load(pdbChainCode);

					// get graph
					RIGraph graph = pdb.getRIGraph(edgeType, cutoff);
					
					
					String edgeTypeStr = edgeType.replaceAll("/", ":");
					
					File outputFile = new File(outputDir,pdbCode+pdbChainCode+"_"+edgeTypeStr+"_"+cutoff+".cm");
					switch(outFormat) {
					case CMVIEW: graph.writeToFile(outputFile.getAbsolutePath()); break;
					case PAUL: graph.writeToPaulFile(outputFile.getAbsolutePath()); break;
					case SADP: graph.writeToSADPFile(outputFile.getAbsolutePath()); break;
					}
				
					//long end = System.currentTimeMillis();
					//double time = (double) (end -start)/1000;

					System.out.println("Wrote "+outputFile.getAbsolutePath());
					//System.out.printf("%5.3f s\n",time);
					
					if (randomSubset<1) {
						for (int n=1;n<=NUM_RANDOM_SUBSETS;n++) {
							RIGraph subset = graph.getRandomSubset(randomSubset);
							subset.writeToFile(new File(outputDir,
									pdbCode+pdbChainCode+"_"+edgeTypeStr+"_"+cutoff+
									"_"+(int)(randomSubset*100.0)+
									"_"+String.format("%03d", n)+".cm").getAbsolutePath());
						}
						System.out.println("Wrote "+NUM_RANDOM_SUBSETS+" subset graph files of "+(int)(randomSubset*100.0)+" percent size to "+outputDir);
					} else if (randomSubset>1) {
						for (int n=1;n<=NUM_RANDOM_SUBSETS;n++) {
							RIGraph noisy = graph.getRandomNoiseMap(randomSubset-1);
							noisy.writeToFile(new File(outputDir,
									pdbCode+pdbChainCode+"_"+edgeTypeStr+"_"+cutoff+
									"_"+(int)((randomSubset-1)*100.0)+
									"_"+String.format("%03d", n)+".cm").getAbsolutePath());							
						}
						System.out.println("Wrote "+NUM_RANDOM_SUBSETS+" cm files with "+(int)((randomSubset-1)*100.0)+" percent noise contacts added to "+outputDir);
					}
					
					numPdbs++;

				} catch (PdbLoadError e) {
					System.err.println("Error loading pdb data for " + pdbCode + pdbChainCode+", specific error: "+e.getMessage());
				} catch (PdbCodeNotFoundError e) {
					System.err.println("Couldn't find pdb code "+pdbCode);
				} catch (SQLException e) {
					System.err.println("SQL error for structure "+pdbCode+pdbChainCode+", error: "+e.getMessage());
				}

			}

			// output results
			System.out.println("Number of structures done successfully: " + numPdbs);


		} else {
			File pdbFile = new File(pdbfile);
			try {
				Pdb pdb = new PdbfilePdb(pdbfile);
				if (pdbChainCode4file==null) {
					pdbChainCode4file = pdb.getChains()[0];
				}

				pdb.load(pdbChainCode4file);
				RIGraph graph = pdb.getRIGraph(edgeType, cutoff);

				String edgeTypeStr = edgeType.replaceAll("/", ":");
				
				String filename = pdbFile.getName().substring(0, pdbFile.getName().lastIndexOf("."));
				File outputFile = new File(outputDir,filename+"_"+edgeTypeStr+"_"+cutoff+".cm");
				switch(outFormat) {
				case CMVIEW: graph.writeToFile(outputFile.getAbsolutePath()); break;
				case PAUL: graph.writeToPaulFile(outputFile.getAbsolutePath()); break;
				case SADP: graph.writeToSADPFile(outputFile.getAbsolutePath()); break;
				}
				System.out.println("Wrote graph file "+outputFile.getAbsolutePath()+" from pdb file "+pdbFile);
				
				if (randomSubset<1) {
					for (int n=1;n<=NUM_RANDOM_SUBSETS;n++) {
						RIGraph subset = graph.getRandomSubset(randomSubset);
						subset.writeToFile(new File(outputDir,
								filename+"_"+edgeTypeStr+"_"+cutoff+
								"_"+(int)(randomSubset*100.0)+
								"_"+String.format("%03d", n)+".cm").getAbsolutePath());
					}
					System.out.println("Wrote "+NUM_RANDOM_SUBSETS+" subset graph files of "+(int)(randomSubset*100.0)+" percent size to "+outputDir);
				} else if (randomSubset>1) {
					for (int n=1;n<=NUM_RANDOM_SUBSETS;n++) {
						RIGraph noisy = graph.getRandomNoiseMap(randomSubset-1);
						noisy.writeToFile(new File(outputDir,
								filename+"_"+edgeTypeStr+"_"+cutoff+
								"_"+(int)((randomSubset-1)*100.0)+
								"_"+String.format("%03d", n)+".cm").getAbsolutePath());							
					}
					System.out.println("Wrote "+NUM_RANDOM_SUBSETS+" cm files with "+(int)((randomSubset-1)*100.0)+" percent noise contacts added to "+outputDir);
				}
				
			} catch (PdbLoadError e) {
				System.err.println("Error loading from pdb file "+pdbFile+", specific error: "+e.getMessage());
			}
		}
	}

}

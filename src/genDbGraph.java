import gnu.getopt.Getopt;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;

import owl.core.connections.CSAConnection;
import owl.core.connections.ConsurfConnection;
import owl.core.runners.DsspRunner;
import owl.core.runners.NaccessRunner;
import owl.core.structure.Pdb;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.structure.PdbasePdb;
import owl.core.structure.PdbfilePdb;
import owl.core.structure.features.SecStrucElement;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.MySQLConnection;



//import proteinstructure.CiffilePdb;


public class genDbGraph {
	/*------------------------------ constants ------------------------------*/
	
	public static final String			PDB_DB = "pdbase";
	public static final String			DSSP_EXE = "/project/StruPPi/bin/dssp";
	public static final String			DSSP_PARAMS = "--";
	public static final String			NACCESS_EXE = "/project/StruPPi/bin/naccess";
	public static final String			NACCESS_PARAMS = "";

	//public static double			cutoff = 4.2;
	//public static String			edgeType = "ALL";
	

	// main
	
	public static void main(String[] args) throws IOException {
		
		
		String help = "Usage, 3 options:\n" +
				"1)  genDbGraph -i <listfile> -d <distance_cutoff> -t <contact_type> [-r directed] -s <seq_sep> -o <output_db> [-D <pdbase_db>] [-m <mode>] \n" +
				"2)  genDbGraph -p <pdb_code> -c <chain_pdb_code> -d <distance_cutoff> -t <contact_type> [-r directed] -s <seq_sep> -o <output_db> [-D <pdbase_db>] [-m <mode>] \n" +
				"3)  genDbGraph -f <pdbfile> -c <chain_pdb_code> -d <distance_cutoff> -t <contact_type> [-r directed] -s <seq_sep> -o <output_db> [-m <mode>] \n" +
				"\nA comma separated list of contact types and distance cutoffs can be given instead of just 1, e.g. -d 8.0,8.5 -t Ca,Cb will generate the graphs for Ca at 8.0 and for Cb at 8.5\n" +
				"If only 1 contact type given and multiple cutoffs, graphs will be generated at all the cutoffs for the one contact type\n"+
				"\nIn case 2) also a list of comma separated pdb codes and chain codes can be specified, e.g. -p 1bxy,1jos -c A,A\n" +
				"\nIf pdbase_db not specified, the default pdbase will be used\n" +
				"\nSecondary structure will be taken from pdbase database. If reading from pdb file and the pdb file is missing the secondary structure, then it will be assigned using dssp\n"; 

		String listfile = "";
		String[] pdbCodes = null;
		String[] pdbChainCodes = null;
		String pdbfile = "";
		String pdbaseDb = PDB_DB;
		String[] edgeTypes = null;
		double[] cutoffs = null;
		int[] seqseps = null;
		boolean[] directed = null; 
		String outputDb = "";
		String mode = "GRAPH";
		
		Getopt g = new Getopt("genDbGraph", args, "i:p:c:f:d:t:r:s:o:D:m:h?");
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
				String[] cutoffsStr = g.getOptarg().split(",");
				cutoffs = new double[cutoffsStr.length];
				for (int i =0;i<cutoffs.length;i++) {
					cutoffs[i] = Double.valueOf(cutoffsStr[i]);
				}
				break;
			case 't':
				edgeTypes = g.getOptarg().split(",");
				break;
			case 'r':
				String[] directedStr = g.getOptarg().split(",");
				directed = new boolean[directedStr.length];
				for (int i =0;i<directed.length;i++) {
					directed[i] = Boolean.valueOf(directedStr[i]);
				}
				break;				
			case 's':
				String[] seqsepsStr = g.getOptarg().split(",");
				seqseps = new int[seqsepsStr.length];
				for (int i=0;i<seqseps.length;i++) {
					seqseps[i] = Integer.valueOf(seqsepsStr[i]);
				}
				break;
			case 'o':
				outputDb = g.getOptarg();
				break;
			case 'D':
				pdbaseDb = g.getOptarg();
				break;
			case 'm':
				mode = g.getOptarg();
				break;
			case 'h':
			case '?':
				System.out.println(help);
				System.exit(0);
				break; // getopt() already printed an error
			}
		}

		if (directed==null) {
			// we set by default all directed to false
			directed = new boolean[edgeTypes.length];
		}
		
		if (outputDb.equals("") || edgeTypes==null || cutoffs==null) {
			System.err.println("Some missing option\n");
			System.err.println(help);
			System.exit(1);
		}
		if (edgeTypes.length!=cutoffs.length && edgeTypes.length!=1) {
			System.err.println("Not same number of contact types as cutoffs given\n");
			System.err.println(help);
			System.exit(1);
		}
		if (seqseps != null && edgeTypes.length!=seqseps.length) {
			System.err.println("Not same number of contact types as sequence separations given\n");
			System.err.println(help);
			System.exit(1);
		}
		if (directed != null && edgeTypes.length!=directed.length) {
			System.err.println("Not same number of contact types as directionalities given\n");
			System.err.println(help);
			System.exit(1);
		}
		if (listfile.equals("") && pdbCodes==null && pdbfile.equals("")){
			System.err.println("Either a listfile, some pdb codes/chain codes or a pdbfile must be given\n");
			System.err.println(help);
			System.exit(1);
		}
		if ((!listfile.equals("") && pdbCodes!=null) || (!listfile.equals("") && !pdbfile.equals("")) || (pdbCodes!=null && !pdbfile.equals(""))) {
			System.err.println("Options -p/-c, -i and -f/-c are exclusive. Use only one of them\n");
			System.err.println(help);
			System.exit(1);			
		}
		if (!(mode.equals("GRAPH") || mode.equals("PDB") || mode.equals("BOTH"))) {
			System.err.println("Allowed values for mode:GRAPH,PDB,BOTH.");
			System.err.println(help);
			System.exit(1);			
		}
		
		// setting edgeTypes in case only 1 was given with multiple cutoffs
		if (edgeTypes.length==1 && cutoffs.length>1) {
			String edgeType = edgeTypes[0];
			edgeTypes = new String[cutoffs.length];
			for (int i=0;i<cutoffs.length;i++){
				edgeTypes[i] = edgeType;
			}
		}

		
		MySQLConnection conn = null;		

		try{
			conn = new MySQLConnection();
			conn.setSqlMode("NO_UNSIGNED_SUBTRACTION,TRADITIONAL");
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
				fpdb.close();
			}

			int numPdbs = 0;

			for (int i=0;i<pdbCodes.length;i++) {
				String pdbCode = pdbCodes[i];
				String pdbChainCode = pdbChainCodes[i];
				
				boolean dssp = false, scop = false, naccess = false, consurf = false, ec = false, csa = false;
				int numGraphs = 0;
				
				try {
					
					System.out.println("Getting pdb data for "+pdbCode+"_"+pdbChainCode);
					
					Pdb pdb = new PdbasePdb(pdbCode, pdbaseDb, conn);	
					//Pdb pdb = new CiffilePdb(new File("/project/StruPPi/BiO/DBd/PDB-REMEDIATED/data/structures/unzipped/all/mmCIF/"+pdbCode+".cif"), pdbChainCode);	
					pdb.load(pdbChainCode);					
					try {
						pdb.setSecondaryStructure(DsspRunner.runDssp(pdb,DSSP_EXE, DSSP_PARAMS, SecStrucElement.ReducedState.THREESTATE, SecStrucElement.ReducedState.THREESTATE));
						dssp = true;
					} catch (Exception e) {
						System.err.println(e.getMessage());
					}
					if (!mode.equals("GRAPH")) {
						try {
							pdb.checkScop("1.73", false);
							scop = true;
						} catch (Exception e) {
							System.err.println(e.getMessage());
						}
						try {
							ConsurfConnection consurfConn = new ConsurfConnection();
							int mistakes = consurfConn.getConsurfDataLocal(pdb, null);
							System.out.println("ConsurfHssp Mistakes:"+mistakes);
							if (mistakes == 0) consurf = true;
						} catch (Exception e) {
							System.err.println(e.getMessage());
						}	
						try {
							pdb.checkEC(false);
							ec = true;
						} catch (Exception e) {
							System.err.println(e.getMessage());
						}
						try {
							int mistakes = CSAConnection.parseCSA(pdb,"2.2.7", false);
							System.out.println("CSA Mistakes:"+mistakes);
							if (mistakes == 0) csa = true;
						} catch (Exception e) {
							System.err.println(e.getMessage());
						}/*
						pdb.restrictToScopDomain("d1pjua2");
						try {
							pdb.runDssp(DSSP_EXE, DSSP_PARAMS, SecStrucElement.ReducedState.THREESTATE, SecStrucElement.ReducedState.THREESTATE);
							//pdb.runDssp(DSSP_EXE, DSSP_PARAMS);
							dssp = true;
						} catch (Exception e) {
							System.err.println(e.getMessage());
						}*/					
						try {
							NaccessRunner naccRunner = new NaccessRunner(new File(NACCESS_EXE), NACCESS_PARAMS);
							naccRunner.runNaccess(pdb);
							naccess = true;
						} catch (Exception e) {
							System.err.println(e.getMessage());
						}
						
						//pdb.writeToDb(conn,outputDb);
						pdb.writeToDbFast(conn, outputDb);
					}
				
					// get graphs
					if (!mode.equals("PDB")) {
						for (int j = 0; j<edgeTypes.length; j++) {
							System.out.print("--> "+(directed[j]?"directed":"")+" graph "+edgeTypes[j]+" for cutoff "+cutoffs[j]);
							
							RIGraph graph = pdb.getRIGraph(edgeTypes[j], cutoffs[j], directed[j]);
							//graph.restrictContactsBetweenSs();
							if (seqseps != null) {
								if (seqseps[j] > 1) {
									System.out.print(" and sequence separation >= "+seqseps[j]);
									graph.restrictContactsToMinRange(seqseps[j]);
								}
							}
							//graph.writeToDb(conn,outputDb);
							graph.writeToDbFast(conn,outputDb);
							
							System.out.println();							
							numPdbs++;
							numGraphs++;
						}
					}
					
				} catch (PdbLoadError e) {
					System.err.println("Error loading pdb data for " + pdbCode + pdbChainCode+", specific error: "+e.getMessage());
				} catch (PdbCodeNotFoundError e) {
					System.err.println("Couldn't find pdb code "+pdbCode);
				} catch (SQLException e) {
					System.err.println("SQL error for structure "+pdbCode+"_"+pdbChainCode+", error: "+e.getMessage());
				} 
				/* catch (CiffileFormatError e) {
					System.err.println(e.getMessage());
				}*/
				
				System.out.println("SUMMARY:"+pdbCode+"_"+pdbChainCode+" dssp:"+dssp+" scop:"+scop+" naccess:"+naccess+" consurf:"+consurf+" ec:"+ec+" csa:"+csa+ " graphs:"+numGraphs);

			}

			// output results
			System.out.println("Number of graphs loaded successfully: " + numPdbs);


		} else {
			String pdbChainCode = pdbChainCodes[0];
			boolean dssp = false, scop = false, naccess = false, consurf = false, ec = false, csa = false;
			int numGraphs = 0;
			
			try {
				
				System.out.println("Getting chain "+pdbChainCode+" from pdb file "+pdbfile);
				
				Pdb pdb = new PdbfilePdb(pdbfile);
				pdb.load(pdbChainCode);
				if (!pdb.hasSecondaryStructure()) {
					pdb.setSecondaryStructure(DsspRunner.runDssp(pdb,DSSP_EXE, DSSP_PARAMS, SecStrucElement.ReducedState.THREESTATE, SecStrucElement.ReducedState.THREESTATE));
				}
				if (!mode.equals("GRAPH")) {
					try {
						pdb.setSecondaryStructure(DsspRunner.runDssp(pdb,DSSP_EXE, DSSP_PARAMS, SecStrucElement.ReducedState.THREESTATE, SecStrucElement.ReducedState.THREESTATE));
						dssp = true;
					} catch (Exception e) {
						System.err.println(e.getMessage());
					}
					try {
						pdb.checkScop("1.73", false);
						scop = true;
					} catch (Exception e) {
						System.err.println(e.getMessage());
					}
					try {
						ConsurfConnection consurfConn = new ConsurfConnection();
						int mistakes = consurfConn.getConsurfDataLocal(pdb, null);
						System.out.println("ConsurfHssp Mistakes:"+mistakes);
						if (mistakes == 0) consurf = true;
					} catch (Exception e) {
						System.err.println(e.getMessage());
					}
					try {
						pdb.checkEC(false);
						ec = true;
					} catch (Exception e) {
						System.err.println(e.getMessage());
					}
					try {
						int mistakes = CSAConnection.parseCSA(pdb,"2.2.7", false);
						System.out.println("CSA Mistakes:"+mistakes);
						if (mistakes == 0) csa = true;
					} catch (Exception e) {
						System.err.println(e.getMessage());
					}
					/*pdb.restrictToScopDomain("d1eaka3");
					try {
						pdb.runDssp(DSSP_EXE, DSSP_PARAMS, SecStrucElement.ReducedState.THREESTATE, SecStrucElement.ReducedState.THREESTATE);
						//pdb.runDssp(DSSP_EXE, DSSP_PARAMS);
						dssp = true;
					} catch (Exception e) {
						System.err.println(e.getMessage());
					}*/					
					try {
						NaccessRunner naccRunner = new NaccessRunner(new File(NACCESS_EXE), NACCESS_PARAMS);
						naccRunner.runNaccess(pdb);
						naccess = true;
					} catch (Exception e) {
						System.err.println(e.getMessage());
					}
					
					//pdb.writeToDb(conn,outputDb);
					pdb.writeToDbFast(conn, outputDb);
				}
				
				// get graphs
				if (!mode.equals("PDB")) {
					for (int j = 0; j<edgeTypes.length; j++) {
						System.out.print("--> "+(directed[j]?"directed":"")+" graph "+edgeTypes[j]+" for cutoff "+cutoffs[j]);
						
						RIGraph graph = pdb.getRIGraph(edgeTypes[j], cutoffs[j], directed[j]);
						//graph.restrictContactsBetweenSs();
						if (seqseps != null) {
							if (seqseps[j] > 1) {
								System.out.print(" and sequence separation >= "+seqseps[j]);
								graph.restrictContactsToMinRange(seqseps[j]);
							}
						}
						//graph.writeToDb(conn,outputDb);
						graph.writeToDbFast(conn,outputDb);
						
						System.out.println();						
						numGraphs++;
					}
				}
				
			} catch (SQLException e) {
				System.err.println("Couldn't write graph to db, error: "+e.getMessage());
			} catch (PdbLoadError e) {
				System.err.println("Error loading from pdb file "+pdbfile+", specific error: "+e.getMessage());
			} 
			
			System.out.println("SUMMARY:"+pdbfile+"_"+pdbChainCode+" dssp:"+dssp+" scop:"+scop+" naccess:"+naccess+" consurf:"+consurf+" ec:"+ec+" csa:"+csa+ " graphs:"+numGraphs);

		}
		
		// closing db connection
		conn.close();
	}

}

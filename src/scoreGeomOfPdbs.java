import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.TreeMap;
import java.util.Vector;
import java.util.Map.Entry;

import edu.uci.ics.jung.graph.util.Pair;

import owl.casp.benchmarking.Benchmarking;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.util.FileFormatException;
import owl.decoyScoring.GeomScorer;

public class scoreGeomOfPdbs {

	public static final String maxClusterExecutable = "/Volumes/StruPPi/bin/maxcluster"; //"/project/StruPPi/bin/maxcluster";
	public static final String tempFileName = "/Users/vehlow/tmp/casp.benchmark.tmp";
	
	private static final String folderCorrectedFiles = "SortedAtoms";
		
	private static final double cutOff=8.0;
	private static final String pdbChainCode=" ";
	private static final String contactType="SC";		
	private static final int minSeqSep = 1;
//	private static final int maxSeqSep = 50;		

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args){
		

	    String zipP = "/Users/vehlow/Documents/workspace/CMView/src/resources/sphoxelBG/";
        String archiveFN = zipP + "SphoxelBGs.zip";

		boolean all = false;
		String target = "T0471"; // only of interest if all=false --> single target
		int scoreFctType = GeomScorer.FCT_SUM_LOS;

		boolean normaliseScore = false;
		boolean allServerMod = false;
		int directNbThres = 1;
		boolean usePreCalcSphoxelMaps = false;
		
		String scoreType = "";
		switch (scoreFctType){
		case GeomScorer.FCT_EC:
			scoreType += "FCT_EC_"; break;
		case GeomScorer.FCT_SUM_LOS:
			scoreType += "FCT_SUM_LOS_"; break;
		case GeomScorer.FCT_EC_LOSg0:
			scoreType += "FCT_EC_LOSg0_"; break;
		case GeomScorer.FCT_EC_LOSs0:
			scoreType += "FCT_EC_LOSs0_"; break;
		case GeomScorer.FCT_SUM_LOSnorth:
			scoreType += "FCT_SUM_LOSnorth_"; break;
		case GeomScorer.FCT_SUM_LOSsouth:
			scoreType += "FCT_SUM_LOSsouth_"; break;
		default :
			scoreType += "";
		}
		if (!usePreCalcSphoxelMaps)
			scoreType += "new_";
		scoreType += "dResNrG";
		scoreType += String.valueOf(directNbThres);
		if (normaliseScore)
			scoreType += "Norm";
		if (allServerMod)
			scoreType += "_all";
		else
			scoreType += "_TS1";
		
		GeomScorer geomScorer;
		if (usePreCalcSphoxelMaps)
			geomScorer = new GeomScorer(archiveFN, contactType, cutOff, minSeqSep);
		else
			geomScorer = new GeomScorer(contactType, cutOff, minSeqSep);
		geomScorer.setNormaliseScore(normaliseScore);
		geomScorer.setDirectNbThres(directNbThres);
		geomScorer.setScoringFctType(scoreFctType);		
		
		if (all){
			// score all server models within casp8 folder
			try {
				scoreAllServerModels(geomScorer, scoreType, allServerMod);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}			
		}
		else{     
			// score server models for one target
			try {
				scoreServerModelsOfTarget(target, geomScorer, scoreType, allServerMod);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}			
		}
		   

	}
	

	private static void scoreServerModelsOfTarget(String target, GeomScorer geomScorer, String scoreType, boolean allServerMod) throws IOException{
		String inputDir = "/Volumes/StruPPi/CASP8/server_models/";
		String dir = inputDir+target+"/"+folderCorrectedFiles+"/";
	    String targetListFN = "scoreTargetList_"+scoreType;
    	String inputFile = inputDir+target+"_"+targetListFN+".txt";
	    String outputFN = "scoreList_"+scoreType; //Nb";
    	String outputFile = inputDir+target+"_"+outputFN+".csv";
    	String rankFile = inputDir+target+"_Rank"+".txt";
    	
    	String[] paths = new String[]{dir, inputFile, outputFile, rankFile, target};
    	scoreServerModels(paths, geomScorer, scoreType, allServerMod);
	}
	
	private static void scoreServerModels(String[] paths, GeomScorer geomScorer, String scoreType, boolean allServerMod) throws IOException{
		// paths = new String[]{dir, inputFile, outputFile, rankFile, target};
		
		File predDir = new File("/Volumes/StruPPi/CASP8/submitted");
		String groupSuffix = "TS183";
		boolean eval3D = true;
		
		String dir = paths[0];
		String inputFile = paths[1];
		String outputFile = paths[2];
		String rankFile = paths[3];
		String target = paths[4];

		File folder = new File(dir);
	
		String logFileFN = "/Volumes/StruPPi/CASP8/server_models/"+target+"_ScoringLogFile_"+scoreType+".txt";
		PrintWriter logOut = new PrintWriter(new FileWriter(logFileFN));
		
		File[] listOfFiles = folder.listFiles();
	    PrintWriter out = new PrintWriter(new FileWriter(inputFile));
//	    int cnt = 0;
	    for (int i=0;i<listOfFiles.length; i++)	{
//	    	String name =listOfFiles[i].getName();
	    	if (!allServerMod){
	    		if(listOfFiles[i].getName().contains("TS1"))
	    			out.println(listOfFiles[i].getName());
	    	}
	    	else
	    		out.println(listOfFiles[i].getName());
    	}
	    out.close();	    	
    	System.out.println("Element: " + target+" "+dir);
    	System.out.println(inputFile);
    	System.out.println(outputFile);	
    	logOut.println("Element: " + target+" "+dir);
    	logOut.println(inputFile);
    	logOut.println(outputFile);

		PrintStream out3 = new PrintStream(new File(rankFile));
    	boolean rankingAvailable = false;
    	TreeMap<String, Pair<Double>> ranking = null; //new TreeMap<String, Pair<Double>>();
    	String baseDir = "/Volumes/StruPPi/CASP8/";
		Benchmarking bm = new Benchmarking(null, new File(baseDir));
		int ret = bm.getTargetResult(predDir, groupSuffix, target, eval3D);
		if(ret == Benchmarking.NO_ERROR) {
			bm.printLastResultTable(out3);
			out3.println();
			bm.printLastResultSummary(out3);
			
			ranking = bm.getLastResult();
		}
		if (ranking != null && ranking.size()>0){
			rankingAvailable = true;
			// iterate over translated coordinates
			for(Entry<String,Pair<Double>> entry : ranking.entrySet()){
				out3.println(entry.getKey()+": "+entry.getValue().getFirst()+"  "+entry.getValue().getSecond());
			}
		}
//		else
//			System.out.println("Error");
    	
		if (rankingAvailable){
	    	BufferedReader in = new BufferedReader(new FileReader(inputFile)); 
			String proteinId; 
			PrintWriter out2 = new PrintWriter(new FileWriter(outputFile));
			//int numStructures = 0;
			out2.println("target"+","+"serverModel"+","+"Rank"+","+"GDT"+","+"GeomScore");
			while ((proteinId = in.readLine()) != null){ 
//				if (proteinId.startsWith("3Dpro"))
				{
					
					int modNo=0;
					char c = proteinId.charAt(proteinId.length()-1);
					if (c=='1'||c=='2'||c=='3'||c=='4'||c=='5')
						modNo = Integer.valueOf(String.valueOf(c)); //Integer.valueOf(proteinId.charAt(proteinId.length()-1));
					if (modNo !=0 ){
						PdbChain pdb = null;
						try {
							PdbAsymUnit fullpdb = new PdbAsymUnit(new File(dir+proteinId),modNo);
							pdb = fullpdb.getChain(pdbChainCode);
						} catch (PdbLoadException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
							continue;
						} catch (FileFormatException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
							continue;
						}
						
//						if (!Scorer.isValidPdb(pdb) || 
						if (pdb==null || pdb.getRIGraph(contactType, cutOff)==null) {
							continue;
						}
						double score = geomScorer.scoreIt(pdb);
						double rank = 0;
						double gdtScore = 0;
						if (ranking.containsKey(proteinId)){
							rank = ranking.get(proteinId).getFirst();
							gdtScore = ranking.get(proteinId).getSecond();
						}
						else{
//							System.out.println("No ranking and score for "+proteinId);
						}

						System.out.println("Final logOddsScore for protein "+proteinId+" is : "+score); //+"  for "+i+" chosen contacts");
						logOut.println("Final logOddsScore for protein "+proteinId+" is : "+score); //+"  for "+i+" chosen contacts");
						
						out2.println(target+","+proteinId+","+rank+","+gdtScore+","+score); // add rmsd etc.
						
						//numStructures++;
					}
					
				}					
			}
			out2.close();
			in.close();
		}
		
		logOut.close();	    		
		
	}
	
	private static void scoreAllServerModels(GeomScorer geomScorer, String scoreType, boolean allServerMod) throws IOException{
		
		/*
		 * Type the path of the directory where you have the models you want to score in the string inputDir.
		 */
	    String inputDir = "/Volumes/StruPPi/CASP8/server_models/";
	    Vector<String> inputDirs = new Vector<String>();
	    Vector<String> targetNames = new Vector<String>();
		/*
		 * This is the input file containing a list of all the files in the folder you have 
		 * given in the inputDir path below. 
		 */		
	    String inputFile = "/Volumes/StruPPi/CASP8/server_models/scoreTargetList.txt";
	    Vector<String> inputFiles = new Vector<String>();
	    String targetListFN = "scoreTargetList_"+scoreType;
	    Vector<String> outputRankFiles = new Vector<String>();
	    /*
	     * This is the file where the output of the scores will be generated.
	     */
	    String outputFile = "/Volumes/StruPPi/CASP8/server_models/scoresList.txt";	
	    Vector<String> outputFiles = new Vector<String>();	
	    String outputFN = "scoreList_"+scoreType; //Nb";
			    	    
	    File folder = new File(inputDir);
	    String[] listOfFolders = folder.list();
	    for (int i=0;i<listOfFolders.length; i++){
	    	File file = new File(inputDir+listOfFolders[i]);
	    	if (file.isDirectory()){
	    		targetNames.add(listOfFolders[i]);
	    		inputDirs.add(inputDir+listOfFolders[i]+"/"+folderCorrectedFiles+"/");
	    		inputFiles.add(inputDir+listOfFolders[i]+"_"+targetListFN+".txt");
	    		outputRankFiles.add(inputDir+listOfFolders[i]+"_Rank"+".txt");
	    		outputFiles.add(inputDir+listOfFolders[i]+"_"+outputFN+".csv");
	    	}
	    }
	    int iF=0;
	    for (String dir : inputDirs){
//	    	if(targetNames.get(iF).equals("T0471"))
	    	{
		    	inputFile = inputFiles.get(iF);
		    	outputFile = outputFiles.get(iF);
		    	String rankFile = outputRankFiles.get(iF);
		    	String target = targetNames.get(iF); // "T0389";
		    	
		    	String[] paths = new String[]{dir, inputFile, outputFile, rankFile, target};
		    	scoreServerModels(paths, geomScorer, scoreType, allServerMod);

//		    	folder = new File(dir);
		    	
//			    File[] listOfFiles = folder.listFiles();
//			    PrintWriter out = new PrintWriter(new FileWriter(inputFile));
////			    int cnt = 0;
//			    for (int i=0;i<listOfFiles.length; i++)	{
////			    	String name =listOfFiles[i].getName();
//			    	if (!allServerMod){
//			    		if(listOfFiles[i].getName().contains("TS1"))
//			    			out.println(listOfFiles[i].getName());
//			    	}
//			    	else
//			    		out.println(listOfFiles[i].getName());
//		    	}
//			    out.close();	    	
//		    	System.out.println("Element" + iF + ": " + target+" "+dir);
//		    	System.out.println(inputFiles.get(iF));
//		    	System.out.println(outputFiles.get(iF));	
//		    	logOut.println("Element" + iF + ": " + target+" "+dir);
//		    	logOut.println(inputFiles.get(iF));
//		    	logOut.println(outputFiles.get(iF));
//
//				PrintStream out3 = new PrintStream(new File(rankFile));
//		    	boolean rankingAvailable = false;
//		    	TreeMap<String, Pair<Double>> ranking = null; //new TreeMap<String, Pair<Double>>();
//		    	String baseDir = "/Volumes/StruPPi/CASP8/";
//				Benchmarking bm = new Benchmarking(null, new File(baseDir));
//				int ret = bm.getTargetResult(predDir, groupSuffix, target, eval3D);
//				if(ret == Benchmarking.NO_ERROR) {
//					bm.printLastResultTable(out3);
//					out3.println();
//					bm.printLastResultSummary(out3);
//					
//					ranking = bm.getLastResult();
//				}
//				if (ranking != null && ranking.size()>0){
//					rankingAvailable = true;
////					Set set = ranking.entrySet();//		// iterate over translated coordinates
//					for(Entry<String,Pair<Double>> entry : ranking.entrySet()){
////						System.out.println(entry.getKey()+": "+entry.getValue().getFirst()+"  "+entry.getValue().getSecond());
//						out3.println(entry.getKey()+": "+entry.getValue().getFirst()+"  "+entry.getValue().getSecond());
//					}
//				}
////				else
////					System.out.println("Error");
//		    	
//				if (rankingAvailable){
//			    	BufferedReader in = new BufferedReader(new FileReader(inputFile)); 
//					String proteinId; 
//					PrintWriter out2 = new PrintWriter(new FileWriter(outputFile));
//		//			int j=0;
//					int numStructures = 0;
//					out2.println("serverModel"+","+"Rank"+","+"GDT"+","+"GeomScore");
//					while ((proteinId = in.readLine()) != null){ // && numStructures<1) { 
////						System.out.println(" "+proteinId); //+"  for "+i+" chosen contacts");
////						logOut.println(" "+proteinId); //+"  for "+i+" chosen contacts");
////						if (proteinId.charAt(0)=='P')
////						if (proteinId.startsWith("3Dpro"))
//						{
//							
//							int modNo=0;
//							char c = proteinId.charAt(proteinId.length()-1);
//							if (c=='1'||c=='2'||c=='3'||c=='4'||c=='5')
//								modNo = Integer.valueOf(String.valueOf(c)); //Integer.valueOf(proteinId.charAt(proteinId.length()-1));
//							if (modNo !=0 ){
//								PdbChain pdb = new PdbfilePdb(dir+proteinId);
//								try {
//									pdb.load(pdbChainCode,modNo);
//								} catch (PdbLoadError e) {
//									// TODO Auto-generated catch block
//									e.printStackTrace();
//									continue;
//								}
//								
////								if (!Scorer.isValidPdb(pdb) || 
//								if (pdb==null || pdb.getRIGraph(contactType, cutOff)==null) {
//									continue;
//								}
//								double score = geomScorer.scoreIt(pdb);
//								double rank = 0;
//								double gdtScore = 0;
//								if (ranking.containsKey(proteinId)){
//									rank = ranking.get(proteinId).getFirst();
//									gdtScore = ranking.get(proteinId).getSecond();
//								}
//								else{
////									System.out.println("No ranking and score for "+proteinId);
//								}
//
//								System.out.println("Final logOddsScore for protein "+proteinId+" is : "+score); //+"  for "+i+" chosen contacts");
//								logOut.println("Final logOddsScore for protein "+proteinId+" is : "+score); //+"  for "+i+" chosen contacts");
//								
//								out2.println(proteinId+","+rank+","+gdtScore+","+score); // add rmsd etc.
//								
//								numStructures++;
//							}
//							
//						}					
//					}
//					out2.close();
//					in.close();
//				}
//				
//				logOut.close();	    		
	    	}
	    	iF++;
	    }
	
	}

}

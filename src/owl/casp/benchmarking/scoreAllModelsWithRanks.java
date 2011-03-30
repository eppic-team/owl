package owl.casp.benchmarking;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;

import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.graphs.CaspRRFileRIGraph;
import owl.core.structure.graphs.GraphComparisonResult;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatException;
import owl.core.util.MySQLConnection;

/**
 * Score all Casp7/8/9 server models plus our SMEG-CCP prediction in terms of acc/cov SR/LR and rank and write the results to database.
 * Table: idx, target, method, acc, cov, acc12, cov12, acc24, cov24, rank_cm, rank_cm12, rank_cm24
 * 
 * The following external resources are being used:
 * - our contact predictions from <base_dir>/cccp/T0xzy.rr
 * - the renumbered native structures from <base_dir>/answers/T0xzy.pdb
 * - server models from <base_dir>/server_models/T0xzy/*_TS1
 * 
 * @author stehr
 */
public class scoreAllModelsWithRanks {
	
	/*------------------------------ constants ------------------------------*/
	public static String PRED_EXTENSION = ".rr";
	public static String OUR_METHOD_NAME = "SMEG-CCP";
	
	/*---------------------------- static methods ---------------------------*/
	
	public static String[] getTargetsFromListFile(File listFile) {
		LinkedList<String> targetList = new LinkedList<String>();
		
		try {
			BufferedReader in = new BufferedReader(new FileReader(listFile));
			String target;
			while((target = in.readLine()) != null) {
				targetList.add(target.trim());
			}
		} catch (IOException e) {
			System.err.println("Error reading list file: " + listFile);
		}
		
		String[] targets = new String[targetList.size()];
		targets = targetList.toArray(targets);
		return targets;
	}
	
	public static ModelScore getModelScore(int idx, String targetName, String methodName, RIGraph predGraph, RIGraph nativeGraph) {
		ModelScore score = new ModelScore(targetName, methodName);
		GraphComparisonResult res = predGraph.evaluatePrediction(nativeGraph);
		GraphComparisonResult res12 = predGraph.evaluatePrediction(nativeGraph, 12);
		GraphComparisonResult res24 = predGraph.evaluatePrediction(nativeGraph, 24);

		score.idx = idx;
		score.acc = 100.0 * res.accuracy;
		score.cov = 100.0 * res.coverage;
		score.acc12 = 100.0 * res12.accuracy;
		score.cov12 = 100.0 * res12.coverage;
		score.acc24 = 100.0 * res24.accuracy;
		score.cov24 = 100.0 * res24.coverage;
		score.cts = res.predicted;
		score.cts12 = res12.predicted;
		score.cts24 = res24.predicted;
		
		return score;
	}
	
	
	/*--------------------------------- main --------------------------------*/
	
	public static void main(String[] args) {
		
		if(args.length < 2) {
			System.out.println("Usage: scoreAllModelsWithRanks <base_dir> <db_table>");		
			System.out.println("e.g. scoreAllModelsWithRanks /project/StruPPi/CASP8/ casp8.scores_ranks");
			System.exit(1);
		}
		
		String baseDir = args[0];
		String dbTable = args[1];
		
		File serverModelsDir = new File(baseDir, "server_models");
		File answerDir = new File(baseDir, "answers");
		File predDir = new File(baseDir, "cccp");
		
		// connect to database and create output table
		MySQLConnection conn = null;
		try {
			conn = new MySQLConnection();
			
			System.out.println("Creating database table (if not exists) " + dbTable);
			ModelScore.createDbTableWithRanks(conn, dbTable);

		} catch (SQLException e1) {
			System.err.println("Error connecting to database: " + e1.getMessage());
			System.exit(1);
		}
		
		// iterate over target dirs
		FileFilter targetDirFilter = new FileFilter() {
			public boolean accept(File file) {
				return (file.isDirectory() && file.getName().length() == 5 && file.getName().startsWith("T0"));
			}
		};
		File[] targetDirs = serverModelsDir.listFiles(targetDirFilter);
		
		// for each target
		for(File targetDir: targetDirs) {

			String targetName = targetDir.getName();
			System.out.println(targetName);
			
			// get answer
			File answer = new File(answerDir, targetName + ".pdb");
			if(!answer.canRead()) {
				System.err.println("Error: File " + answer + " not found. Skipping target.");
				continue;
			}
			
			PdbChain answerPdb = null;
			try {
				PdbAsymUnit fullpdb = new PdbAsymUnit(answer);
				answerPdb = fullpdb.getFirstChain();
			} catch (PdbLoadException e1) {
				System.err.println("Error: Could not load Pdb file " + answerPdb + ": " + e1.getMessage() + ". Skipping target.");
				continue;
			} catch (IOException e1) {
				System.err.println("Error: Could not load Pdb file " + answerPdb + ": " + e1.getMessage() + ". Skipping target.");
				continue;
			} catch (FileFormatException e1) {
				System.err.println("Error: Could not load Pdb file " + answerPdb + ": " + e1.getMessage() + ". Skipping target.");
				continue;
			}
			RIGraph nativeGraph = answerPdb.getRIGraph("Cb", 8.0);
			
			// get prediction TODO: Calculate on the fly instead of loading from file
			File pred = new File(predDir, "/" + targetName + "/" + targetName + PRED_EXTENSION);
			if(!pred.canRead()) {
				System.err.println("Error: File " + pred + " not found. Skipping target.");
				continue;
			}
//			PdbfilePdb predPdb = new PdbfilePdb(pred.getPath());
//			try {
//				predPdb.load(predPdb.getChains()[0]);
//			} catch (PdbLoadError e1) {
//				System.err.println("Error: Could not load Pdb file " + predPdb + ": " + e1.getMessage() + " Skipping target.");
//				continue;
//			}
//			RIGraph predGraph = answerPdb.getRIGraph("Cb", 8.0);

			RIGraph predGraph;
			try {
				predGraph = new CaspRRFileRIGraph(pred.getPath());
			} catch (IOException e2) {
				System.err.println("Error: Could not load RR file " + pred + ": " + e2.getMessage() + ". Skipping target.");
				continue;
			} catch (FileFormatException e2) {
				System.err.println("Error: Could not load RR file " + pred + ": " + e2.getMessage() + ". Skipping target.");
				continue;
			}

			
			// evaluate our predicted model
			int idx = 1;
			ModelScore score = getModelScore(idx, targetName, OUR_METHOD_NAME, predGraph, nativeGraph);
			
			LinkedList<ModelScore> modelScores = new LinkedList<ModelScore>();
			modelScores.add(score);
			
			// evaluate server models
			
			FileFilter firstModelFilter = new FileFilter() {
				public boolean accept(File file) {
					return (file.isFile() && file.getName().indexOf("TS1") > 0);
				}				
			};
			File[] models = targetDir.listFiles(firstModelFilter);
			for(File model: models) {

				String modelName = model.getName();
				
				// evaluate model
				PdbChain modelPdb = null;
				try {
					PdbAsymUnit fullpdb = new PdbAsymUnit(model);
					modelPdb = fullpdb.getFirstChain();
					RIGraph modelGraph = modelPdb.getRIGraph("Cb", 8.0);
					idx++;
					ModelScore modScore = getModelScore(idx, targetName, modelName, modelGraph, nativeGraph);
					modelScores.add(modScore);				
				} catch (PdbLoadException e1) {
					System.err.println("Error: Could not load Pdb file " + model + ": " + e1.getMessage() + ". Ignoring");
				} catch (IOException e1) {
					System.err.println("Error: Could not load Pdb file " + model + ": " + e1.getMessage() + ". Ignoring");
				} catch (FileFormatException e1) {
					System.err.println("Error: Could not load Pdb file " + model + ": " + e1.getMessage() + ". Ignoring");
				}				
				
			}
			
			// fill in ranks
			
			// order by (acc+cov)/2 (desc)
			Collections.sort(modelScores, new Comparator<ModelScore>(){
				public int compare(ModelScore s1, ModelScore s2) {
					return Double.compare(s2.acc+s2.cov, s1.acc+s1.cov);
				}
			});
			int rank=0;
			for(ModelScore s:modelScores) {
				s.rankCm = ++rank;
			}
			// order by acc12/cov12
			Collections.sort(modelScores, new Comparator<ModelScore>(){
				// order by acc+cov/2 (desc)
				public int compare(ModelScore s1, ModelScore s2) {
					return Double.compare(s2.acc12+s2.cov12, s1.acc12+s1.cov12);
				}
			});
			rank=0;
			for(ModelScore s:modelScores) {
				s.rankCm12 = ++rank;
			}
			// order by acc24/cov24
			Collections.sort(modelScores, new Comparator<ModelScore>(){
				// order by acc+cov/2 (desc)
				public int compare(ModelScore s1, ModelScore s2) {
					return Double.compare(s2.acc24+s2.cov24, s1.acc24+s1.cov24);
				}
			});
			rank=0;
			for(ModelScore s:modelScores) {
				s.rankCm24 = ++rank;
			}

			// write to db
			for(ModelScore s:modelScores) {
				try {
					//s.writeToDbWithRanks(conn, dbTable);
					s.updateNumContacts(conn, dbTable);
				} catch (SQLException e) {
					System.err.println("Error writing to database for " + s.target + " " + s.methodName + ": " + e.getMessage());
				}					
			}
		}
		System.out.println("done.");
	}
}

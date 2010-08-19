package owl.casp.benchmarking;

import java.io.File;
import java.io.FileFilter;
import java.io.IOException;
import java.sql.SQLException;

import owl.core.util.MySQLConnection;

/**
 * Scores all server models for Casp 7/8/9 in terms of GDT_TS and Acc/Cov against native and
 * writes the result to database for further analysis.
 * 
 * - server models are expected in <baseDir>/server_models/T0???/
 * - native structures are expected in <baseDir>/answers/T0???.pdb
 * - database access parameters are expected in a .my.conf file in the user's home directory
 * - results are written to given database table which has to contain columns model, target, scores...
 * 
 * @author stehr
 */
public class scoreAllServerModels {
	
	/*------------------------------ constants ------------------------------*/
	// if server models are to be taken from a subdir of /server_models/T0???/
	// (we use this for preprocessing the models) - NOT YET IMPLEMENTED!
	static final boolean takeModelsFromSubDir = false;
	static final File subDir = new File("SortedAtoms");
	
	/*---------------------------- private methods --------------------------*/
	
	/*--------------------------------- main --------------------------------*/
	
	public static void main(String[] args) {

		if(args.length < 2) {
			System.out.println("Usage: scoreAll <base_dir> <database_table>");
			System.out.println("e.g. scoreAll /project/StruPPi/CASP8/ casp8.scores");
			System.exit(1);
		}
		
		String baseDir = args[0];
		String dbTable = args[1];
		
		File serverModelsDir = new File(baseDir, "server_models");
		File answerDir = new File(baseDir, "answers");
		
		// connect to database and create output table
		MySQLConnection conn = null;
		try {
			conn = new MySQLConnection();
			
			System.out.println("Creating database table (if not exists) " + dbTable);
			ModelScore.createDbTable(conn, dbTable);

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
				
		for(File targetDir: targetDirs) {

			String targetName = targetDir.getName();
			System.out.println(targetName);
			
			// get answer
			File answer = new File(answerDir, targetName + ".pdb");
			if(!answer.canRead()) {
				System.err.println("Error: File " + answer + " not found. Skipping target.");
				continue;
			}
			
			// iterate over files in target dir
			FileFilter firstModelFilter = new FileFilter() {
				public boolean accept(File file) {
					return (file.isFile() && file.getName().indexOf("TS1") > 0);
				}				
			};
			File[] models = targetDir.listFiles(firstModelFilter);
			for(File model: models) {

				String modelName = model.getName();
				
				// evaluate model
				ModelScore scores = null;
				try {
					scores = Benchmarking.getModelScores(model, answer);
					
					if(scores == null) {
						System.err.println("Error: Scores object is null for " + targetName + " " + modelName);
					} else
					try {
						scores.writeToDb(conn, dbTable);
					} catch (SQLException e) {
						System.err.println("Error writing to database for " + targetName + " " + modelName + ": " + e.getMessage());
					}					
				} catch (IOException e) {
					System.err.println("Error calculating scores for " + targetName + "/" + modelName + ": " + e.getMessage());
					//System.exit(1);
				}	
			}			
		}		
	}
}

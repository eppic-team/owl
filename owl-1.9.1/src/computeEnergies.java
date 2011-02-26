import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

import owl.core.runners.MaxClusterRunner;
import owl.core.runners.MaxClusterRunner.ScoreType;
import owl.core.runners.tinker.TinkerError;
import owl.core.runners.tinker.TinkerRunner;



/**
 * Script to calculate energies of a set of decoys by using tinker minimize program.
 * At the moment the force field used for the energy calculation is amber99. But 
 * that could be change to any of the ones supported by tinker. See:
 * ftp://dasher.wustl.edu/pub/tinker/params
 * 
 * @author duarte
 *
 */
public class computeEnergies {

	private static final String TINKER_BIN_DIR = "/project/StruPPi/Software/tinker/bin";
	private static final String FORCEFIELD_PRM_FILE = "/project/StruPPi/Software/tinker/amber/amber99.prm";
	private static final String maxClusterExecutable = "/project/StruPPi/bin/maxcluster";
	private static final double DEFAULT_RMSGRADIENT = 1.0;
	
	private class Row {
		
		String id;
		double rmsdOriginal;
		int rmsdOrigRank;
		double rmsdMinimised;
		int rmsdMinRank;
		double gdtMinimised;
		int gdtMinRank;
		double energy;
		int energyRank;
		
		public Row() {
			this.id = null;
			this.rmsdOriginal = 0;
			this.rmsdMinimised = 0;
			this.gdtMinimised = 0;
			this.energy = 0;
		}
		
		public Row(String id, double rmsdOriginal, double rmsdMinimised, double gdtMinimised, double energy) {
			this.id = id;
			this.rmsdOriginal = rmsdOriginal;
			this.rmsdMinimised = rmsdMinimised;
			this.gdtMinimised = gdtMinimised;
			this.energy = energy;
		}
		
	}
	
	
	public static void main(String[] args) throws IOException {

		File listFile = new File(args[0]); // 1st argument listFile
		File nativeStruct = new File(args[1]); // 2nd argument: native structure pdb file
		double rmsGradient = DEFAULT_RMSGRADIENT;
		if (args.length>=3) rmsGradient = Double.parseDouble(args[2]); // 3rd (optional) argument: rms gradient for minimization

		// getting file names from list file
		ArrayList<File> pdbFiles = new ArrayList<File>();
		
		BufferedReader in = new BufferedReader(new FileReader(listFile));
		String line;
		int countPdbFiles=0;
		while ((line=in.readLine())!=null) {
			// ignore empty lines or comments
			if (line.equals("") || line.startsWith("#")) {
				continue;
			}
			pdbFiles.add(new File(line.trim()));
			countPdbFiles++;
		}
		
		TinkerRunner tr = new TinkerRunner(TINKER_BIN_DIR,FORCEFIELD_PRM_FILE);
		MaxClusterRunner maxCluster = new MaxClusterRunner(maxClusterExecutable);
		PrintWriter log = new PrintWriter(new File(listFile.getParent(),"tinker.log"));
				
		// getting energy of native
		Row nativeRow = new computeEnergies().new Row();
		try {
			double energy = tr.minimize(nativeStruct, rmsGradient, log);
			nativeRow = new computeEnergies().new Row(nativeStruct.getName(), 0.0, 0.0, 0.0, energy);			
		} catch(TinkerError e) {
			System.err.println("Tinker error while computing energy for "+nativeStruct+". Skipping.");
		}	
		
		// getting energy of all decoys in list
		ArrayList<Row> rows = new ArrayList<Row>();
		for (File pdbFile:pdbFiles){
			try {
				double energy = tr.minimize(pdbFile, rmsGradient, log);
				double rmsdOriginal = maxCluster.calculatePairwiseScore(pdbFile.getAbsolutePath(), nativeStruct.getAbsolutePath(),  ScoreType.RMSD);
				String basename = pdbFile.getAbsolutePath();
				if (basename.contains(".")) {
					basename = basename.substring(0,basename.lastIndexOf("."));
				}
				double rmsdMinimised = maxCluster.calculatePairwiseScore(basename+".min.pdb", nativeStruct.getAbsolutePath(), ScoreType.RMSD);
				double gdtMinimised = maxCluster.calculatePairwiseScore(basename+".min.pdb", nativeStruct.getAbsolutePath(), ScoreType.GDT);
				Row row = new computeEnergies().new Row(pdbFile.getName(), rmsdOriginal, rmsdMinimised, gdtMinimised, energy);
				rows.add(row);			
				
			} catch(TinkerError e) {
				System.err.println("Tinker error while computing energy for "+pdbFile+". Skipping. Error: "+e.getMessage());
				continue;
			}	
		
		}		
		
		log.close();
		
		System.out.println("Done "+rows.size()+" pdb files out of "+countPdbFiles);
		
		// sorting on rmsdOrig to get rmsdOrig rank
		sortOnRmsdOrig(rows);
		int rmsdOrigRank=1;
		for (Row row: rows) {
			row.rmsdOrigRank = rmsdOrigRank;	
			rmsdOrigRank++;
		}

		// sorting on gdtMin to get gdt rank
		sortOnGdtMin(rows);
		int gdtRank=1;
		for (Row row: rows) {
			row.gdtMinRank = gdtRank;	
			gdtRank++;
		}
		
		// sorting on energy to get energy rank
		sortOnEnergy(rows);
		int energyRank=1;
		for (Row row: rows) {
			row.energyRank = energyRank;	
			energyRank++;
		}		

		// sorting on rmsdMin for rmsd min rank and printing out
		sortOnRmsdMin(rows);
		int rmsdRank=1;
		for (Row row: rows) {
			row.rmsdMinRank = rmsdRank;	
			rmsdRank++;
		}
		
		System.out.println("id\trmsdOrig\trmsdOrigRank\trmsdMin\trmsdMinRank\tgdtMin\tgdtRank\tenergy\tenergyRank");
		for (Row row: rows) {
			System.out.printf("%s %6.3f %3d %6.3f %3d %6.3f %3d %12.4f %3d\n",row.id,row.rmsdOriginal,row.rmsdOrigRank,row.rmsdMinimised,row.rmsdMinRank,row.gdtMinimised,row.gdtMinRank,row.energy,row.energyRank);
		}
		System.out.println("Native: ");
		System.out.printf("%s %6.3f %3d %6.3f %3d %6.3f %3d %12.4f %3d\n",nativeRow.id,nativeRow.rmsdOriginal,nativeRow.rmsdOrigRank,nativeRow.rmsdMinimised,nativeRow.rmsdMinRank,nativeRow.gdtMinimised,nativeRow.gdtMinRank,nativeRow.energy,nativeRow.energyRank);
		
		
	}

	private static void sortOnRmsdOrig(ArrayList<Row> rows) {
		Collections.sort(rows, new Comparator<Row>() {
			public int compare(Row o1, Row o2) {
				return new Double(o1.rmsdOriginal).compareTo(o2.rmsdOriginal);
			}
		});
	}
	
	private static void sortOnGdtMin(ArrayList<Row> rows) {
		Collections.sort(rows, new Comparator<Row>() {
			public int compare(Row o1, Row o2) {
				// ordered is inversed here (in gdt maximum score is best)
				return new Double(o2.gdtMinimised).compareTo(o1.gdtMinimised);
			}
		});
	}
	
	private static void sortOnEnergy(ArrayList<Row> rows) {
		Collections.sort(rows, new Comparator<Row>() {
			public int compare(Row o1, Row o2) {
				return new Double(o1.energy).compareTo(o2.energy);
			}
		});
	}
	
	private static void sortOnRmsdMin(ArrayList<Row> rows) {
		Collections.sort(rows, new Comparator<Row>() {
			public int compare(Row o1, Row o2) {
				return new Double(o1.rmsdMinimised).compareTo(o2.rmsdMinimised);
			}
		});
	}
	
}

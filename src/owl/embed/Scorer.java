package owl.embed;

import java.io.*;
import java.sql.SQLException;
//import java.util.ArrayList;
import java.util.*;

import owl.core.structure.*;
import owl.core.structure.graphs.FileRIGraph;
import owl.core.structure.graphs.RIGEdge;
import owl.core.structure.graphs.RIGNode;
import owl.core.structure.graphs.RIGraph;
import owl.core.util.FileFormatException;
import owl.core.util.MySQLConnection;
import owl.core.util.RegexFileFilter;

import edu.uci.ics.jung.graph.util.Pair;
//import java.util.HashMap;
//import java.util.HashSet;
//import java.util.Random;

//import edu.uci.ics.jung.graph.util.Pair;

import Jama.Matrix;



/**
 * Class containing static methods to score contact map subsets.
 * 
 * @author duarte
 *
 */
public class Scorer {

	private static File dir = new File("/project/LitNet/Essence/Data/Contactmaps"); 


	
	/*------------------------- privates  -------------------------------*/
	
	/**
	 * Infer bounds for all pairs from a spare bounds matrix via the triangle inequality
	 * @param sparseBounds
	 * @return
	 */
	public static Bound[][] inferAllBounds(Bound[][] sparseBounds) {
		BoundsSmoother bs = new BoundsSmoother(sparseBounds);
		return bs.getInitialBoundsAllPairs();
	}
	
	/*------------------------- publics  -------------------------------*/
	
	/**
	 * Evaluates the error value for the given sparseBounds matrix (a subset of the contact map) 
	 * by first inferring bounds for all pairs through the triangle inequality and then 
	 * measuring how well the all pairs matrix fits the contact map.  
	 * Thus a high error value means the given sparseBounds matrix has low information content
	 * about the rest of the contact map and low error value that the matrix has high 
	 * information content.
	 * The error measure consists of the sum of the deviation of the upper bounds 
	 * to the ones in the contact map: sum(max(0,(u'i-ui))), i.e. if upper bound in 
	 * given matrix is below the one in the contact map there is no penalty.
	 * @param sparseBounds
	 * @param fullContactMap
	 * @return
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundException 
	 * @throws PdbLoadException 
	 */
	public static double getCMError(Bound[][] sparseBounds, RIGraph fullContactMap) {
		Bound[][] cmapBounds = Reconstructer.convertRIGraphToBoundsMatrix(fullContactMap);
		// infer bounds for all pairs through triangle inequality
		Bound[][] bounds = inferAllBounds(sparseBounds);
		double sumDev = 0;
		for (int i=0;i<bounds.length;i++) {
			for (int j=i+1;j<bounds.length;j++) {
				if (cmapBounds[i][j]!=null) {
					sumDev += Math.max(0,bounds[i][j].upper-cmapBounds[i][j].upper);
				}
			}
		}
		return sumDev/(double) fullContactMap.getFullLength(); 
	}
	
	public static double getCMError(RIGraph fullContactMap,Bound[][] sparseBounds) throws PdbCodeNotFoundException, SQLException, PdbLoadException {
		MySQLConnection conn = new MySQLConnection ();
		Pdb pdb = new PdbasePdb(fullContactMap.getPdbCode(), "pdbase_20090728", conn);
		pdb.load(fullContactMap.getChainCode());
		Matrix distmat = pdb.calcDistMatrixJamaFormat("Ca");
		Bound[][] cmapBounds = convertRIGraphToBounds(fullContactMap,distmat.getArray());//Reconstructer.convertRIGraphToBoundsMatrix(fullContactMap);
		// infer bounds for all pairs through triangle inequality
		Bound[][] bounds = inferAllBounds(sparseBounds);
		double sumDev = 0;
		for (int i=0;i<bounds.length;i++) {
			for (int j=i+1;j<bounds.length;j++) {
				if (cmapBounds[i][j]!=null) {
					sumDev += Math.max(0,bounds[i][j].upper-cmapBounds[i][j].upper);
				}
			}
		}
		return sumDev/(double) fullContactMap.getFullLength(); 
	}
	
	/**
	 * Evaluates the error value for the given contact map subset 
	 * by first inferring bounds for all pairs through the triangle inequality and then 
	 * measuring how well the all pairs matrix fits the contact map.  
	 * Thus a high error value means the given sparseBounds matrix has low information content
	 * about the rest of the contact map and low error value that the matrix has high 
	 * information content  
	 * @param subset
	 * @param fullContactMap
	 * @return
	 * @throws PdbLoadException 
	 * @throws SQLException 
	 * @throws PdbCodeNotFoundException 
	 */
	public static double getCMError(RIGraph subset, RIGraph fullContactMap) throws PdbCodeNotFoundException, SQLException, PdbLoadException {
		
		//initializing Bound array instance using class 'Reconstructer' method 'convertRIGraphToBoundsMatrix()' 
		Bound[][] subsetBounds = Reconstructer.convertRIGraphToBoundsMatrix(subset);
		
		return getCMError(subsetBounds, fullContactMap);
	}

	/**
	 * Evaluates the error value for the given contact map file containing a subset. 
	 * It will be compared to the contact map taken from the pdbCode/contactType/cutoff
	 * present in the contact map file. Applicable to a single contact map entry or a set of
	 * contact map entries
	 * 
	 * @param rigFile
	 * @param pdbaseDb
	 * @param conn
	 * @return
	 * @throws FileFormatException
	 * @throws IOException
	 * @throws PdbCodeNotFoundException
	 * @throws SQLException
	 * @throws PdbLoadException
	 */
	public static double getCMError(File rigFile, String pdbaseDb, MySQLConnection conn, boolean randomsampling) throws FileFormatException, IOException, PdbCodeNotFoundException, SQLException, PdbLoadException {

		//new instance of class RIGraph created using File instance 'rigFile' 
		RIGraph graph = new FileRIGraph(rigFile.getAbsolutePath());

		//initializing integer variable counting all edges in 'graph'
		int numContSubset = graph.getEdgeCount();

		//initializing String variable for PDB code of RIGraph instance 'graph'
		String pdbCode = graph.getPdbCode();

		//initializing String variable for chain code, only monomeric proteins excepted 
		String pdbChainCode = graph.getPdbChainCode();

		//initializing String variable for contact type
		String ct = graph.getContactType();

		//initializing double variable for cutoff length
		double cutoff = graph.getCutoff();

		//initializing Pdb instance 'pdb' using subclass 'PdbasePdb' constructor 
		Pdb pdb = new PdbasePdb(pdbCode, pdbaseDb, conn);

		//loading PDB file of the same protein as specified by rigFile
		pdb.load(pdbChainCode);

		//initializing RIGraph instance with predefined parameters
		RIGraph fullContactMap = pdb.getRIGraph(ct, cutoff);

		//initializing number of all given contacts of full contact input
		int numTotalCont = fullContactMap.getEdgeCount();

		//output of size
		System.out.printf("size: %4.2f\n",(double)numContSubset/(double)numTotalCont);

		//return of error value via getCMError method
		if(randomsampling){
			double tests = (double) fullContactMap.getEdgeCount();

			double partedges = (double) numContSubset;

			double percent = partedges/tests;

			//initializing number of all given contacts of full contact input

			Bound[][] subs = randomSets(fullContactMap, percent);
			//output of size
			//System.out.printf("size: %4.2f\n",(double)numContSubset/(double)numTotalCont);

			//return of error value via getCMError method
			return getCMError(subs, fullContactMap);

		}
		else{
			return getCMError(graph, fullContactMap);
		}
	}	

	/**
	 * Evaluates the error value for the given distance map file containing a subset. 
	 * It will be compared to the distance map taken from the pdbCode/contactType/cutoff
	 * present in the contact map file. Applicable to a single contact map entry or a set of
	 * contact map entries.
	 * 
	 * @param test
	 * @param pdbaseDb
	 * @param conn
	 * @return
	 * @throws FileFormatException
	 * @throws IOException
	 * @throws PdbCodeNotFoundException
	 * @throws SQLException
	 * @throws PdbLoadException
	 */
	public static double getDMError (File test, String pdbaseDb, MySQLConnection conn, boolean randsampling) throws FileFormatException, IOException, PdbCodeNotFoundException, SQLException, PdbLoadException {
		RIGraph sub = new FileRIGraph(test.getAbsolutePath());
		//initializing String variable for PDB code of RIGraph instance 'graph'
		String pdbCode = sub.getPdbCode();
		
		//initializing String variable for chain code, only monomeric proteins excepted 
		String pdbChainCode = sub.getPdbChainCode();
		
		//initializing String variable for contact type
		String ct = sub.getContactType();
		
		//initializing double variable for cutoff length
		//double cutoff = sub.getCutoff();
		
		//initializing Pdb instance 'pdb' using subclass 'PdbasePdb' constructor 
		Pdb pdb = new PdbasePdb(pdbCode, pdbaseDb, conn);
			
		//loading PDB file of the same protein as specified by rigFile
		pdb.load(pdbChainCode);
		
		Matrix fullDistanceMap = pdb.calcDistMatrixJamaFormat(ct);
		double[][] dm = fullDistanceMap.getArray();
		
		if(randsampling){
			double cutoff = sub.getCutoff();
			RIGraph randsubset = pdb.getRIGraph(ct, cutoff);
			
			double tests = (double) randsubset.getEdgeCount();
			
			double partedges = (double) sub.getEdgeCount();
			
			double percent = partedges/tests;
			
			Bound[][] subss = randomSets(randsubset, percent);

			return getDMError(subss, dm);
			
		}
		else{
		
		return getDMError(sub, dm);
		}
		
	}	
	
	/**
	 * Transforms RIGraph inputs into Bound output
	 * @param sub
	 * @param full
	 * @return 
	 */
	public static double getDMError (RIGraph sub, double[][] full){
		Bound[][] subset = Reconstructer.convertRIGraphToBoundsMatrix(sub);		
		return getDMError (subset, full);
	}
	
	/**
	 * Infer bounds for all pairs from a spare bounds matrix via the triangle inequality, evaluates 
	 * mean squared deviation for all given upper bounds per distance
	 * @param sub
	 * @param full
	 * @return
	 */
	public static double getDMError (Bound[][] sub, double[][] full) {
		Bound[][] sparse = inferAllBounds(sub);
		//inferring all bounds via triangular inequality
		
		double error = 0;
		for(int index1 = 0; index1 < full.length; index1++){
			//looping over all bounds, first dimension
			
			for(int index2 = index1 +1; index2 < full.length; index2++){
				//looping over all bounds, second dimension
				
				if(full[index1][index2] != 0.0) {
					//only contacts of type (i,j), with j >= i + 1 are considered
					
					if(sparse[index1][index2].upper > full[index1][index2]){
					//only if the upper bound is greater than the exact distance value error is calculated
					
					error += Math.pow(sparse[index1][index2].upper - full[index1][index2],2.0 );
					//error function: square root of the sum of the square difference of upper bound and the exact distance,
					//divided by 1/2*n*(n - 1), where n is the number all contacts
					}
				}
			}
		}
		error = Math.pow(error, 0.5);
	return 2.0*error/((double) full.length*(full.length - 1));
	}
	
	/**
	 * samples subsets of all given distances or contacts by using Distiller.sampleSubset method
	 * @param input RIGraph derived from pdb file
	 * @param percent percentage of distances
	 * @return subset
	 * @throws NullPointerException
	 * @throws ArrayIndexOutOfBoundsException
	 */
	public static Bound[][] randomSets (RIGraph input, double percent) throws NullPointerException, ArrayIndexOutOfBoundsException {
	
		Distiller dist = new Distiller(input);
		Bound[][] subset = dist.sampleSubset((int) (input.getEdgeCount()*percent));
		return subset;
	}
	
	
	/**
	 * calculates statistical measures from distance map as mean value and standard deviation
	 * @param test protein structure input
	 * @param pdbaseDb PDB access code
	 * @param conn database access
	 * @param runs number of subsets to be sampled
	 * @return stats mean value and standard deviation
	 * @throws FileFormatException
	 * @throws IOException
	 * @throws PdbCodeNotFoundException
	 * @throws SQLException
	 * @throws PdbLoadException
	 */
	public static double[] averageDMRandom (File test, String pdbaseDb, MySQLConnection conn, int runs) throws FileFormatException, IOException, PdbCodeNotFoundException, SQLException, PdbLoadException {	
		double[] list = new double[runs];
		double average = 0.0;
		boolean bool = true;
		for(int i = 0; i < runs; i++){
			list[i] = getDMError(test, pdbaseDb, conn, bool);
			average = average + list[i];
		}
		average = average/(double) list.length;
		double[] stats = {average, averageSet(list, average)};
		return stats;
	}
	
	/**
	 * calculates statistical measures from contact map as mean value and standard deviation
	 * 
	 * @param test
	 * @param pdbaseDb
	 * @param conn
	 * @param runs
	 * @return stats
	 * @throws FileFormatException
	 * @throws IOException
	 * @throws PdbCodeNotFoundException
	 * @throws SQLException
	 * @throws PdbLoadException
	 */
	public static double[] averageCMRandom (File test, String pdbaseDb, MySQLConnection conn, int runs) throws FileFormatException, IOException, PdbCodeNotFoundException, SQLException, PdbLoadException {
		boolean bool = true;
		double[] list = new double[runs];
		double average = 0.0;
		for(int i = 0; i < runs; i++){
			list[i] = getCMError(test, pdbaseDb, conn, bool);
			average = average + list[i];
		}
		average = average/(double) list.length;
		double[] stats = {average, averageSet(list, average)};
		return stats;
	}
	
	/**
	 * supporting 'averageRandom' method by returning standard deviation 
	 * @param list double array with all error entries 
	 * @param expectedval mean value
	 * @return stand standard deviation
	 */
	public static double averageSet (double[] list, double expectedval){
		double stand = 0.0;
		for(int i = 0; i < list.length; i++){
			stand = stand + Math.pow(list[i] - expectedval, 2.0); 
		}
		stand = Math.pow(stand/((double) list.length*(list.length - 1)),0.5);
		return stand;
	}
	
	/**
	 * output method for distance map error function for a selected contact map file and a randomly generated subset with the same number of distances 
	 * @param test
	 * @param file
	 * @param pdbaseDb
	 * @param conn
	 * @param runs
	 * @throws FileFormatException
	 * @throws IOException
	 * @throws PdbCodeNotFoundException
	 * @throws SQLException
	 * @throws PdbLoadException
	 */
	public static void outputDMScore (File test, PrintWriter file, String pdbaseDb, MySQLConnection conn, int runs) throws FileFormatException, IOException, PdbCodeNotFoundException, SQLException, PdbLoadException {
		boolean bool = false;
		double value = getDMError(test, pdbaseDb, conn, bool);
		double[] values = averageDMRandom(test, pdbaseDb, conn, runs);
		file.print(test.getName()+":");
		file.print(" DM error: "+value+", ");
		file.print("average DM error: "+values[0]);
		file.println(", st deviation: "+values[1]);
		System.out.println("file written...");
				
	}
	/**
	 * output method for contact map error function for a selected contact map file and a randomly generated subset with the same number of contacts 
	 * @param test
	 * @param file
	 * @param pdbaseDb
	 * @param conn
	 * @param runs
	 * @throws FileFormatException
	 * @throws IOException
	 * @throws PdbCodeNotFoundException
	 * @throws SQLException
	 * @throws PdbLoadException
	 */
	public static void outputCMScore (File test, PrintWriter file, String pdbaseDb, MySQLConnection conn, int runs) throws FileFormatException, IOException, PdbCodeNotFoundException, SQLException, PdbLoadException {
		boolean bool = false;
		double value = getCMError(test, pdbaseDb, conn, bool);
		double[] values = new double[averageCMRandom(test, pdbaseDb, conn, runs).length];
		System.arraycopy(averageCMRandom(test, pdbaseDb, conn, runs), 0, values, 0, values.length);
		file.print(test.getName()+":");
		file.print(" CM error: "+value+", ");
		file.print("average CM error: "+values[0]);
		file.println(", st deviation: "+values[1]);
		System.out.println("file written...");
				
	}
	
	public static Bound[][] convertRIGraphToBounds (RIGraph rig, double[][] fulldistancematrix){
		HashSet<Pair<Integer>> contact_pairs = convertRIGraphToHashSet(rig); 
		Iterator<Pair<Integer>> it = contact_pairs.iterator();
		int cont_size = rig.getSequence().length(), counter = 0;
		Bound[][] bounds = new Bound[cont_size][cont_size];
		while(it.hasNext()){
			Pair<Integer> pair = it.next();
			int f_val = pair.getFirst().intValue() - 1, s_val = pair.getSecond().intValue() - 1;
			double value = fulldistancematrix[f_val][s_val];
			bounds[f_val][s_val] = new Bound(value,value);
			bounds[s_val][f_val] = new Bound(value,value);
			counter ++;
		}
		return bounds;
	}
	
	public static HashSet<Pair<Integer>> convertRIGraphToHashSet (RIGraph rig){
		HashSet<RIGEdge> rig_set = new HashSet<RIGEdge> (rig.getEdges());
		Iterator<RIGEdge> it = rig_set.iterator();
		HashSet<Pair<Integer>> pairs = new HashSet<Pair<Integer>> (2*rig_set.size());
		while(it.hasNext()){
			RIGEdge edge = it.next();
			RIGNode node1 = rig.getEndpoints(edge).getFirst(), node2 = rig.getEndpoints(edge).getSecond();
			Integer index_a = new Integer (node1.getResidueSerial()), index_b = new Integer (node2.getResidueSerial());
			Pair<Integer> pair = new Pair<Integer> (0,0);
			if(index_a.intValue() < index_b.intValue()){
				pair = new Pair<Integer> (index_a,index_b);
			}
			else{
				if(index_a.intValue() > index_b.intValue()){
					pair = new Pair<Integer> (index_b,index_a);
				}
			}
			pairs.add(pair);
		}
		return pairs;
	}
	/*-------------------------- main  -------------------------------*/

	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception {
		
		File[] selSubsets = dir.listFiles(new RegexFileFilter(".*_all\\.graph"));
		//File[] rndSubsets = dir.listFiles(new RegexFileFilter(".*_rnd\\.graph"));
		//File[][] randSet;
		String direct = "/project/StruPPi/gabriel/";
		File out = new File(direct+"test.data");
		PrintWriter file = new PrintWriter(out);
		Arrays.sort(selSubsets);
		//Arrays.sort(rndSubsets);
		
		MySQLConnection conn = new MySQLConnection();
		
		for (int i=0;i<selSubsets.length;i++) {
			File selSubset = selSubsets[i];
			/*for(int j = 0; j < selSubsets.length; j++){
				randSet[i][j] = selSubsets[i].getPdbname();
			}*/
			//File rndSubset = rndSubsets[i];/*
			//System.out.printf("%s\t"+Distiller.SCORE_PRINT_FORMAT+"\n",selSubset.getName(),getCMError(selSubset, "pdbase", conn));
			//System.out.printf("%s\t"+Distiller.SCORE_PRINT_FORMAT+"\n",rndSubset.getName(),getCMError(rndSubset, "pdbase", conn));*/
			//int j = 0;
			//double selected = getDMError(selSubset, "pdbase", conn);
			//System.out.printf("%s\t"+Distiller.SCORE_PRINT_FORMAT+"\n",selSubset.getName(),selected);//getDMError(selSubset, "pdbase", conn));
				//double average = averageRandom(selSubset, "pdbase", conn, 100, 0.053)[0];
				//double standdev = averageRandom(selSubset, "pdbase", conn, 100, 0.053)[1];
			//System.out.printf("%s\t"+Distiller.SCORE_PRINT_FORMAT+"\n",selSubset.getName(),average);//getDMERror(selSubset, "pdbase", conn));
			//System.out.println("Standard deviation : "+standdev+"\n");
			//System.out.printf("%s\t"+Distiller.SCORE_PRINT_FORMAT+"\n",rndSubset.getName(),getDMError(rndSubset, "pdbase", conn));
			//System.out.print(getDMError(rndSubset, "pdbase", conn)/getDMError(selSubset, "pdbase", conn));
			//System.out.println();
			
			outputDMScore(selSubset, file, "pdbase", conn, 10);
			outputCMScore(selSubset, file, "pdbase", conn, 10);
		}
		file.close();
		conn.close();
	}

	
}

package embed;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.SQLException;
import java.util.Arrays;

//import proteinstructure.FileRIGraph;
import proteinstructure.GraphFileFormatError;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbasePdb;
import proteinstructure.RIGraph;
import tools.MySQLConnection;
import tools.RegexFileFilter;

/**
 * class to get the correlation function for the given error estimation method
 * @author gmueller
 *
 */
public class CorrelationFunction {
	
	/**
	 * method, to calculate the mean value of a given array of data
	 * @param vec
	 * @return val
	 */
	public static double meanValue (double[] vec){
		int dim = vec.length;
		double val = 0.0;
		for(int i = 0 ; i < dim; i++){
			val = val + vec[i];
		}
		val = val/(double) dim;
		return val;
	}
	
	/**
	 * method, to calculate the mean value of a matrix of data
	 * @param vec
	 * @return val
	 */
	public static double meanValue (double[][] vec){
		int dim = vec.length;
		int j = 0, i = 0;
		double val = 0.0;
		for(j = 0; j < dim; j++){
			for(i = 0 ; i < dim; i++){
				val = val + vec[j][i];
			}
		}
		val = val/((double) i*j);
		return val;
	}
	
	/**
	 * method, to calculate the variance of an array of data
	 * @param vec
	 * @return var
	 */
	public static double varValue (double[] vec){
		int dim = vec.length;
		double mean = meanValue(vec);
		double var = 0.0;
		for(int i = 0; i < dim; i++){
			var = var + Math.pow(mean - vec[i],2.0);
		}
		return var;
	}
	
	/**
	 * method, to calculate the variance of a matrix of data 
	 * @param vec
	 * @param pot
	 * @return var
	 */
	public static double varValue (double[][] vec, double pot){
		int dim = vec.length;
		double mean = meanValue(vec);
		double var = 0.0;
		for(int j = 0; j < dim; j++){
		for(int i = 0; i < dim; i++){
			var = var + Math.pow(mean - vec[j][i], pot);
		}
		}
		return var;
	}
	
	/**
	 * method, to convert the i th row of a matrix to an array
	 * @param vec
	 * @param row
	 * @return con
	 */
	public static double[] vecConverter (double[][] vec, int row){
		int dim = vec.length;
		double[] con = new double[dim];
		for(int i = 0; i < dim; i++){
			con[i] = vec[row][i];
		}
		return con;
	}
	
	/**
	 * method, to multiply two arrays by their compounents
	 * @param vec1
	 * @param vec2
	 * @return vec
	 */
	public static double[] multiEntries (double[] vec1, double[] vec2){
		int dim1 = vec1.length;
		int dim2 = vec2.length;
		double[] vec = new double[dim1];
		if(dim1 != dim2){
			String e = "Array must only have same dimension!!";
			throw new IllegalArgumentException(e);
		}
		else{
			for(int i = 0; i < dim1; i++){
				vec[i] = vec1[i]*vec2[i];
			}
		}
		return vec;
	}
	
	/**
	 * method, to run a evolution of a Population instance, uses evolve-method from Individuals class
	 * @param starter
	 * @param generations
	 * @return pop
	 * @throws SQLException
	 * @throws PdbCodeNotFoundError
	 * @throws PdbLoadError
	 */
	public static Population[] evolve (Population starter, int generations) throws SQLException, PdbCodeNotFoundError, PdbLoadError{
		Population[] pop = new Population[generations];
		Population start = new Population(starter);
		for(int i = 0; i <= generations - 1; i++){
			if(i == 0){
				pop[i] = new Population(Population.evolve(false, start));
			}
			else{
				pop[i] = new Population(Population.evolve(false, pop[i-1]));
			}
		}
		return pop;
	}
	
	/**
	 * print method, prints an array of Population to a specified file
	 * @param pops
	 * @throws FileNotFoundException
	 */
	public static void printToFile (Population[] pops) throws FileNotFoundException{
		int dim = pops.length;
		int size = pops[0].getSize();
		int numOC = pops[0].getNumOfContacts();
		String dir1 = "/home/gmueller/workspace/aglappe/embed/teststructures/";
		FileOutputStream out = new FileOutputStream(dir1+pops[0].getName()+" gen 1 - "+dim+".txt");
		PrintStream printa = new PrintStream(out);
		for(int i = 0; i < dim; i++){
			for(int j = 0; j < size; j++){
				printa.print(numOC+"\t"+i);
				printa.println(pops[i].getPop(j).getDM());
			}
			System.out.println("generation "+i+"\t"+numOC+" written to file...");
		}
		printa.close();
	}
	
	/**
	 * method to print a String to a specified file
	 * @param name
	 * @param results
	 * @throws FileNotFoundException
	 */
	public static void printToFile (String name, String results) throws FileNotFoundException {
		String dir1 = "/home/gmueller/workspace/aglappe/embed/teststructures/";
		FileOutputStream out = new FileOutputStream(dir1+name+" gen 1 - 10.txt");
		PrintStream printa = new PrintStream(out);
		printa.print(results);
		printa.close();
	}
	
	/**
	 * method to print a String to a specified file
	 * @param name
	 * @param results
	 * @param dest
	 * @throws FileNotFoundException
	 */
	public static void printToFile (String name, String results, String dest) throws FileNotFoundException {
		FileOutputStream out = new FileOutputStream(dest);
		PrintStream printa = new PrintStream(out);
		printa.print(results);
		printa.close();
	}
	
	/**
	 * method, to convert a Population array to a String, considering number of contacts, generation and DMError
	 * @param pops
	 * @return st
	 */
	public static String toString(Population[] pops){
		int dim = pops.length, size = pops[0].getSize();
		int numOC = pops[0].getNumOfContacts();
		String st = "";
		for(int i = 0; i < dim; i++){
			for(int j = 0; j < size; j++){
				st = st + numOC + "\t" + i + "\t" + pops[i].getPop(j).getDM() + "\n";
			}
		}
		System.out.println("Number of Contactes: "+numOC+" written to file...");
		return st;
	}
	
	public static void main (String[] args) throws GraphFileFormatError, IOException, SQLException, ArrayIndexOutOfBoundsException, NullPointerException, PdbCodeNotFoundError, PdbLoadError  {
		String dir1 = "/home/gmueller/workspace/aglappe/embed/teststructures/";
		File dir = new File(dir1);
		File[] test = dir.listFiles(new RegexFileFilter(".*_all\\.graph"));
		Arrays.sort(test);
		//File neu = test[2];
		BufferedReader fileIn = new BufferedReader(new FileReader(dir1+"sathya_set.list"));
		String line;
		int dim = test.length;
		RIGraph[] fullCM = new RIGraph[dim];
		MySQLConnection conn = new MySQLConnection();
		int counter = 0;
		String outputString = "", statistics = "";
		while((line = fileIn.readLine()) != null && counter < dim) {
			String[] tokens = line.split(" ");
			String pdbCode=tokens[0];
			String pdbChainCode=tokens[1];
			
			Pdb pdb = new PdbasePdb(pdbCode, "pdbase", conn);
			pdb.load(pdbChainCode);
			fullCM[counter] = pdb.get_graph("Ca", 9);
			counter++;
		}
		fileIn.close();
		int popsize = 20;
		int numconts = 41;
		Individuals[] neus = new Individuals[popsize];
		Individuals[] off = new Individuals[popsize];
		double[][] val = new double[numconts][popsize];
		double percent = 0.0;
		boolean percenttest = false;
		for(int k = 0; k < 1; k++){
			FileOutputStream outFile = new FileOutputStream(dir1+"Correlation function/correlationfunction "+fullCM[k].getPdbCode()+".txt");
			FileOutputStream outFile1 = new FileOutputStream(dir1+"Correlation function/correlationfunction "+fullCM[k].getPdbCode()+"average.txt");
			PrintStream out = new PrintStream(outFile);
			PrintStream out1 = new PrintStream(outFile1);
			FileOutputStream outFileoff1 = new FileOutputStream(dir1+"Correlation function/correlationfunction "+fullCM[k].getPdbCode()+" average off.txt");
			PrintStream outoff1 = new PrintStream(outFileoff1);
			int j = 1;
				while(j <= numconts){
					for(int i = 0; i < popsize; i++){
						if(!percenttest){
							neus[i] = new Individuals(fullCM[k], conn, true, j);
							off[i] = new Individuals(neus[i]);
						}
						else{
							neus[i] = new Individuals(fullCM[k],conn, true, percent);
						}
						val[j][i] = (neus[i].getDM());
						out.print(j+"\t");
						out.println(val[j][i]);
					}
					Population start = new Population(off);
					Population[] pop1 = evolve(start,10);
					for(int i = 0; i < 10; i++){
						statistics = statistics + pop1[i].getPop(0).getNumOfContacts() + "\t" + i + "\t" + pop1[i].getAvDMError() + "\n";
					}
					outputString = outputString + toString(pop1);
					if(j <= 10){
						int onepercent = (int) ((double) neus[0].getNumOfContacts()/(double) fullCM[k].getEdgeCount());
						if(onepercent > 1 && j < 6){
							percent = percent + 0.5;
							percenttest = true;
						}
						else{
							j = j + 1;
							percenttest = false;
						}
					}
					if(j >= 10 && j < 20){
						j = j + 2;
					}
					if(j>= 20 && j < 30){
						j = j + 5;
					}
					if(j >= 30 && j < 50){
						j = j + 10;
					}
				}
				out.close();
				out1.close();
				outoff1.close();
		}
		printToFile(off[0].getName(),outputString);
		printToFile(off[0].getName(), statistics, "/home/gmueller/workspace/aglappe/embed/teststructures/"+off[0].getName()+" average.txt");
	}

}

package embed.contactmaps;

import java.io.*;
import java.sql.SQLException;

import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.util.RegexFileFilter;



public class CMDMErrorEvaluator implements Runnable {
	
	public int evosteps;
	
	public String path;
	
	public String[] file_listing;
	
	public Individuals[] starter;
	
	public Species[] population;
	
	public CMDMErrorEvaluator (String dir) throws IOException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		setCMDMErrorEval(dir,20,30,0.0003);	
	}
	
	public CMDMErrorEvaluator (String dir, int size) throws IOException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		setCMDMErrorEval (dir,size,30,0.0003);
	}
	public void setCMDMErrorEval (String dir, int size, int evosteps, double threshold) throws IOException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		path = new String (dir);
		File dirs = new File (dir);
		File[] cmap_list = dirs.listFiles(new RegexFileFilter (".*.cmap"));
		File[] graph_list = dirs.listFiles(new RegexFileFilter (".*.graph"));
		int dim_cmap = cmap_list.length, dim_graph = graph_list.length, dim_sum = dim_cmap + dim_graph;
		file_listing = new String[dim_sum];
		starter = new Individuals[dim_sum];
		population = new Species [dim_sum];
		for(int i = 0; i < dim_sum; i++){
			if(i < dim_cmap){
				file_listing[i] = new String (cmap_list[i].getAbsolutePath());
				starter[i] = new Individuals(file_listing[i]);
				Demes[] spec = new Demes[size];
				for(int j = 0; j < size; j++){
					spec[j] = new Demes (starter[i],size);
				}
				population[i] = new Species (spec,30,path,path,0.0003);
			}
			else{
				file_listing[i] = new String (graph_list[i].getAbsolutePath());
				starter[i] = new Individuals(file_listing[i]);
				Demes[] spec = new Demes[size];
				for(int j = 0; j < size; j++){
					spec[j] = new Demes (starter[i],size);
				}
				population[i] = new Species (spec,30,path,path,0.0003);
			}
		}
	}
	
	public void run() {
		int length = population.length, final_generation = evosteps;
		for(int i = 0; i < length; i++){
			try {
				population[i].evolve2(20);
			} catch (FileNotFoundException e) {
				//e.printStackTrace();
				System.err.println("Some severe file error occurred during evolution...");
			} catch (SQLException e) {
				System.err.println("Some SQL exception occurred during the run:\n"+e.getSQLState());
				//e.printStackTrace();
				System.exit(1);
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Some severe error concerning the PDB code occurred...");
				System.exit(1);
				//e.printStackTrace();
			} catch (PdbLoadError e) {
				//e.printStackTrace();
				
			} catch (IOException e) {
				System.err.println("Some IO exception occurred...");
				System.exit(1);
			}
			try {
				population[i].evolve2(final_generation);
			} catch (FileNotFoundException e) {
				//e.printStackTrace();
				System.err.println("Some severe file error occurred during evolution...");
			} catch (SQLException e) {
				System.err.println("Some SQL exception occurred during the run:\n"+e.getSQLState());
				//e.printStackTrace();
				System.exit(1);
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Some severe error concerning the PDB code occurred...");
				System.exit(1);
				//e.printStackTrace();
			} catch (PdbLoadError e) {
				//e.printStackTrace();
				
			} catch (IOException e) {
				System.err.println("Some IO exception occurred...");
				System.exit(1);
			}
		}		
	}
	
	public static boolean[] onlyDirectory (String[] dir_names){
		int length = dir_names.length;
		boolean[] dir_test = new boolean[length];
		for(int i = 0; i < length; i++){
			File dir = new File(dir_names[i]);
			if(dir.isDirectory()){
				dir_test[i] = true;
			}
		}
		return dir_test;
	}

	public static void main (String[] args) throws IOException, SQLException, PdbCodeNotFoundError, PdbLoadError{
		int argument_length = args.length;
		boolean[] isdirectory = onlyDirectory(args);
		if(argument_length > 0){
			for(int i = 0; i < argument_length; i++){
				if(isdirectory[i]){
					new CMDMErrorEvaluator(args[i],20);
				}
			}
		}
		else{
			System.err.println("This main method requires at least one String as argument, to denote an directory to be processed.");
			System.exit(1);
		}
	}
}

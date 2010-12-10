package owl.gmbp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;

/**
 * 
 * @author Corinna Vehlow
 *
 */

public class CSVhandler {
	

//	private double [][] array2d;
//	private double [][][] array3d;
//	private final char[] aas = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'}; // possible Residues
	private int[][] sumValues;	
	
	
	public CSVhandler(){
	}
	
	public static void main(String[] args){
		CSVhandler csv = new CSVhandler();	
		
		String outputPath = "/Users/vehlow/Documents/workspace/outputFiles/";
		double[][] ratios;

		String sFileName = csv.openCSVFile(outputPath);
		System.out.println(sFileName);		
		if (sFileName!=null){
            System.out.println("Chosen path/file:" + sFileName);
		}
		else
            System.out.println("No path chosen!");
		// ----- start import from csv-file -----	
		try {
			ratios = csv.readCSVfile2Ddouble(sFileName);
			System.out.println("Dimensions: " + ratios.length + " x " + ratios[0].length);
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
	}
	
	public String openCSVFile(String outputPath) {
		return openCSVFile(outputPath, "Choose CSV-File");
	}
	
	public String openCSVFile(String outputPath, String title) {
		String sFileName = null;
		final JFileChooser chooser = new JFileChooser(title);
		chooser.setDialogTitle(title);
        chooser.setDialogType(JFileChooser.OPEN_DIALOG);
        chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
        FileFilter filter = new FileNameExtensionFilter("CSV file", "csv");
        chooser.setFileFilter(filter);
        final File file = new File(outputPath);

        chooser.setCurrentDirectory(file);

//        chooser.addPropertyChangeListener(new PropertyChangeListener() {
//            public void propertyChange(PropertyChangeEvent e) {
//                if (e.getPropertyName().equals(JFileChooser.SELECTED_FILE_CHANGED_PROPERTY)
//                        || e.getPropertyName().equals(JFileChooser.DIRECTORY_CHANGED_PROPERTY)) {
//                    final File f = (File) e.getNewValue();
//                }
//            }
//        });

        chooser.setVisible(true);
        final int result = chooser.showOpenDialog(null);

        if (result == JFileChooser.APPROVE_OPTION) {
            File inputVerzFile = chooser.getSelectedFile();
            sFileName = inputVerzFile.getPath();
        }
//        System.out.println("Abbruch");
        chooser.setVisible(false); 
        return sFileName;
	}
	
    public void generateCsvFile(double[][] array2d, String sFileName)	{		
		if (array2d.length>0){
			FileWriter writer = null;
			System.out.println(sFileName);
			try	{
			    writer = new FileWriter(sFileName);
			    
			    for(int i=0;i<array2d.length;i++){
			    	for(int j=0;j<array2d[i].length;j++){
			    		double val = array2d[i][j];
			    		String s = String.valueOf(val);
			    		writer.append(s);
			    		if (j<array2d[i].length-1)
			    			writer.append(','); 
			    		    // writer.append(';');
			    	}
			    		writer.append('\n');
			    }	 
				writer.flush();
				writer.close();
			}
			catch(IOException e) {
				e.printStackTrace();
				System.err.println("Didn't manage to create file. No csv-file created.");
			}			
		}
		else {
			System.out.println("Array is empty.");
		}		
	}
    
    public void generateCsvFile(int[][] array2d, String sFileName)	{		
		if (array2d.length>0){
			FileWriter writer = null;
			System.out.println(sFileName);
			try	{
			    writer = new FileWriter(sFileName);
			    
			    for(int i=0;i<array2d.length;i++){
			    	for(int j=0;j<array2d[i].length;j++){
			    		int val = array2d[i][j];
			    		String s = String.valueOf(val);
			    		writer.append(s);
			    		if (j<array2d[i].length-1)
			    			writer.append(','); 
			    		    // writer.append(';');
			    	}
			    		writer.append('\n');
			    }	 
				writer.flush();
				writer.close();
			}
			catch(IOException e) {
				e.printStackTrace();
				System.err.println("Didn't manage to create file. No csv-file created.");
			}			
		}
		else {
			System.out.println("Array is empty.");
		}		
	}
    
    
    public void generateCsvFile(int[][][] array3d, String sFileName)	{
    	String [] names = new String[array3d[0][0].length];
    	for (int i=0; i<array3d[0][0].length; i++){
    		names[i] = String.valueOf(i);
    	}
    	generateCsvFile(array3d, names, sFileName);
    }
    
    public void generateCsvFile(int[][][] array3d, String[] names, String sFileName)	{
		
    	if (array3d[0][0].length != names.length){
    		System.out.println("Inconsistency in array size: " +
    				"3rd dimension of array3d[0][0] and length of names should be the same");
    	}
    	else {
    		if (array3d.length>0){
    			FileWriter writer = null;
    			System.out.println(sFileName);

    			try	{
    			    writer = new FileWriter(sFileName);			    
    			    for(int k=0;k<names.length;k++){
    			    	writer.append(names[k]);
    			    	writer.append(','); 
    			    	writer.append(String.valueOf(names.length)); 
    			    	writer.append('\n');
    			    	for(int i=0;i<array3d.length;i++){
    				    	for(int j=0;j<array3d[i].length;j++){
    				    		int val = array3d[i][j][k];
    				    		String s = String.valueOf(val);
    				    		writer.append(s);
    				    		if (j<array3d[i].length-1)
    				    			writer.append(','); 
    				    		    // writer.append(';');
    				    	}
    				    	writer.append('\n');
    				    }
//    			    	writer.append('\n');
    			    }			    	 
    				writer.flush();
    				writer.close();
    			}
    			catch(IOException e) {
    				e.printStackTrace();
    				System.err.println("Didn't manage to create file.");
    			}			
    		}
    		else {
    			System.out.println("Array is empty. No csv-file created.");
    		}		    		
    	}
    	
	}
    
    
    public void generateCsvFile(double[][][] array3d, String sFileName)	{
    	String [] names = new String[array3d[0][0].length];
    	for (int i=0; i<array3d[0][0].length; i++){
    		names[i] = String.valueOf(i);
    	}
    	generateCsvFile(array3d, names, sFileName);
    }
    
    public void generateCsvFile(double[][][] array3d, String[] names, String sFileName)	{
		
    	if (array3d[0][0].length != names.length){
    		System.out.println("Inconsistency in array size: " +
    				"3rd dimension of array3d[0][0] and length of names should be the same");
    	}
    	else {
    		if (array3d.length>0){
    			FileWriter writer = null;
    			System.out.println(sFileName);

    			try	{
    			    writer = new FileWriter(sFileName);			    
    			    for(int k=0;k<names.length;k++){
    			    	writer.append(names[k]);
    			    	writer.append(','); 
    			    	writer.append(String.valueOf(names.length)); 
    			    	writer.append('\n');
    			    	for(int i=0;i<array3d.length;i++){
    				    	for(int j=0;j<array3d[i].length;j++){
    				    		double val = array3d[i][j][k];
    				    		String s = String.valueOf(val);
    				    		writer.append(s);
    				    		if (j<array3d[i].length-1)
    				    			writer.append(','); 
    				    		    // writer.append(';');
    				    	}
    				    	writer.append('\n');
    				    }
//    			    	writer.append('\n');
    			    }			    	 
    				writer.flush();
    				writer.close();
    			}
    			catch(IOException e) {
    				e.printStackTrace();
    				System.err.println("Didn't manage to create file.");
    			}			
    		}
    		else {
    			System.out.println("Array is empty. No csv-file created.");
    		}		
    	}		
	}
    
    public void generateCSVFile(Vector<float[]> vector, String sFileName){
    	float theta, phi;
    	int gID, iNum, jNum;
    	int resID, ssType;
    	float[] node;
    	if (vector.size()>0){
    		FileWriter writer = null;
			System.out.println(sFileName);
			try	{
			    writer = new FileWriter(sFileName);	
			    for(int i=0; i<vector.size(); i++){
	    			node = (float[]) vector.get(i);
	    			gID = (int) node[0];
	    			iNum = (int) node[1];
	    			jNum = (int) node[2];
	    			theta = node[3];
	    			phi = node[4];
	    			resID = (int) node[5];
	    			ssType = (int) node[6];
//	    			System.out.println("0:graphi_id + '\t' + 1:i_num + '\t' + 2:j_num + '\t' + 3:theta + '\t' + 4:phi");
	    			writer.append(String.valueOf(gID));
			    	writer.append(','); 
			    	writer.append(String.valueOf(iNum)); 
			    	writer.append(','); 
			    	writer.append(String.valueOf(jNum)); 
			    	writer.append(','); 
			    	writer.append(String.valueOf(theta));
			    	writer.append(','); 
			    	writer.append(String.valueOf(phi));
			    	writer.append(','); 
			    	writer.append(String.valueOf(resID));
			    	writer.append(','); 
			    	writer.append(String.valueOf(ssType));
			    	writer.append('\n');
	    		}		    	 
				writer.flush();
				writer.close();
			}
			catch(IOException e) {
				e.printStackTrace();
				System.err.println("Didn't manage to create file.");
			}	    		
    	}
    	else {
			System.out.println("Vector is empty. No csv-file created.");
		}	
    	
    }
    
    public void generateFile(Vector<String[]> vector, String sFileName){
    	String[] line;
    	String s;
    	if (vector.size()>0){
    		FileWriter writer = null;
			System.out.println(sFileName);
			try	{
			    writer = new FileWriter(sFileName);	
			    for(int i=0; i<vector.size(); i++){
	    			line = (String[]) vector.get(i);
	    			for (int j=0; j<line.length; j++){
	    				s = line[j];
				    	writer.append(s);
				    	if (j<line.length-1)
				    		writer.append(','); 
	    			}
			    	writer.append('\n');
	    		}		    	 
				writer.flush();
				writer.close();
			}
			catch(IOException e) {
				e.printStackTrace();
				System.err.println("Didn't manage to create file.");
			}	    		
    	}
    	else {
			System.out.println("Vector is empty. No csv-file created.");
		}	
    	
    }
    
    public Vector<String[]> readCsvFile(String filename) throws NumberFormatException, IOException{
    	Vector<String[]> vector = new Vector<String[]>();
    	String[] row;
    	int numTokens = 0;
		BufferedReader bufRdr;
		try {
			bufRdr = new BufferedReader(new FileReader(filename));
			String line = null;
			String token = null;	 
			//read each line of text file
			while((line = bufRdr.readLine()) != null){
				int index =0;
				StringTokenizer st = new StringTokenizer(line,",");
				numTokens = st.countTokens();
				row = new String[numTokens];
				while (st.hasMoreTokens()){
					token = st.nextToken();
					//get next token and store it in the array
					row[index] = token;
					index++;
				}
				vector.add(row);
			}
			//close the file
			bufRdr.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
		
    	return vector;
    }
    
    public Vector<float[]> readCSVfileVector(String filename) throws NumberFormatException, IOException{
    	Vector<float[]> vector = new Vector<float[]>();;
    	float[] node; // = new float[5];// {gID, num, j_num, theta, phi};
    	int numTokens = 0;
		BufferedReader bufRdr;
		try {
			bufRdr = new BufferedReader(new FileReader(filename));
			String line = null;			 
			line = bufRdr.readLine();
			StringTokenizer st = new StringTokenizer(line,",");
			numTokens = st.countTokens();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
		if (numTokens>0){			
	    	try {
				bufRdr = new BufferedReader(new FileReader(filename));
				String line = null;
				String token = null;	 
				//read each line of text file
				while((line = bufRdr.readLine()) != null){
					int index =0;
					StringTokenizer st = new StringTokenizer(line,",");
					node = new float[numTokens];
					while (st.hasMoreTokens()){
						token = st.nextToken();
						//get next token and store it in the array
						node[index] = Float.valueOf(token);
						index++;
					}
					vector.add(node);
				}
				//close the file
				bufRdr.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}					
		}
		else
			System.out.print("Empty file!");		
		
    	return vector;
    }
    
    public double[][] readCSVfile2Ddouble(String filename) throws NumberFormatException, IOException{
    	double [][] array2d = null;
		int numTokens = 0 , numLines = 0;
		BufferedReader bufRdr;
		try {
			bufRdr = new BufferedReader(new FileReader(filename));
			String line = null;		 
			//read each line of text file
			while((line = bufRdr.readLine()) != null)
			{
//				line = bufRdr.readLine();
				StringTokenizer st = new StringTokenizer(line,",");
				if (st.countTokens() > numTokens)
					numTokens = st.countTokens();
				numLines++;
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
		if (numTokens>0 && numLines>0){
			array2d = new double [numLines][numTokens];
			try {
				bufRdr = new BufferedReader(new FileReader(filename));
				String line = null;
				String token = null;
				int row = 0;
				int col = 0;				 
				//read each line of text file
				while((line = bufRdr.readLine()) != null){
					col = 0;
					StringTokenizer st = new StringTokenizer(line,",");
					while (st.hasMoreTokens()){
						token = st.nextToken();
						//get next token and store it in the array
						array2d[row][col] = Double.valueOf(token);
						col++;
					}
					row++;
				}
				//close the file
				bufRdr.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}	
		}
		
		return array2d;
	}
    
//    public double[][] readCSVfile2Ddouble(String filename) throws NumberFormatException, IOException{
//    	double [][] array2d = null;
//		int numTokens = 0;
//		BufferedReader bufRdr;
//		try {
//			bufRdr = new BufferedReader(new FileReader(filename));
//			String line = null;		 
//			//read each line of text file
////			while((line = bufRdr.readLine()) != null)
//			{
//				line = bufRdr.readLine();
//				StringTokenizer st = new StringTokenizer(line,",");
//				numTokens = st.countTokens();
//			}
//		} catch (FileNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}	
//		if (numTokens%2!=0){
//			System.out.print("Odd number of tokens in csv file! Can't read file.");
//		}
//		else {
//			array2d = new double [numTokens/2][numTokens];
//			try {
//				bufRdr = new BufferedReader(new FileReader(filename));
//				String line = null;
//				String token = null;
//				int row = 0;
//				int col = 0;				 
//				//read each line of text file
//				while((line = bufRdr.readLine()) != null){
//					col = 0;
//					StringTokenizer st = new StringTokenizer(line,",");
//					while (st.hasMoreTokens()){
//						token = st.nextToken();
//						//get next token and store it in the array
//						array2d[row][col] = Double.valueOf(token);
//						col++;
//					}
//					row++;
//				}
//				//close the file
//				bufRdr.close();
//			} catch (FileNotFoundException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}		
//		}
//		
//		return array2d;
//	}
    
    public int[][] readCSVfile2Dint(String filename) throws NumberFormatException, IOException{
    	int [][] array2d = null;
		int numTokens = 0;
		BufferedReader bufRdr;
		try {
			bufRdr = new BufferedReader(new FileReader(filename));
			String line = null;			 
			//read each line of text file
//			while((line = bufRdr.readLine()) != null)
			{
				line = bufRdr.readLine();
				StringTokenizer st = new StringTokenizer(line,",");
				numTokens = st.countTokens();
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
		if (numTokens%2!=0){
			System.out.print("Odd number of tokens in csv file! Can't read file.");
		}
		else {
			array2d = new int [numTokens/2][numTokens];
			try {
				bufRdr = new BufferedReader(new FileReader(filename));
				String line = null;
				String token = null;
				int row = 0;
				int col = 0;				 
				//read each line of text file
				while((line = bufRdr.readLine()) != null){
					col = 0;
					StringTokenizer st = new StringTokenizer(line,",");
					while (st.hasMoreTokens()){
						token = st.nextToken();
						//get next token and store it in the array
						array2d[row][col] = Integer.valueOf(token);
						col++;
					}
					row++;
				}
				//close the file
				bufRdr.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}		
		}
		
		return array2d;
	}
    
    public double[][][] readCSVfile3Ddouble(String filename) throws NumberFormatException, IOException{
    	double [][][] array3d = null;
		int numTokens = 0;
		int numSlices = 0;
		int numLines = 0;
		int x=0, y=0, z=0;
//		int xDim = 0, yDim = 0;
		BufferedReader bufRdr;
		try {
			bufRdr = new BufferedReader(new FileReader(filename));
			String line = null;
			String token = null;	
			// get number of slices (arrays) from first line
			line = bufRdr.readLine(); // just one token in first line
			StringTokenizer st = new StringTokenizer(line,",");
			if (st.countTokens()>1){
				token = st.nextToken();
				token = st.nextToken();
				numSlices = Integer.valueOf(token);
//				System.out.println("number of slices: " + numSlices);
			}
			// get dimension of array
			// xDim = numLines, yDim = numTokens
			boolean xDimFound = false;
			while((line = bufRdr.readLine()) != null && !xDimFound)
			{
				st = new StringTokenizer(line,",");
				if (st.countTokens()>2){
					numLines++;
//					xDim++;
					numTokens = st.countTokens();
//					yDim = st.countTokens();
				}		
				else
					xDimFound = true;
			}
//			System.out.println("number of tokens: "+numTokens);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
//		if (numTokens%2!=0){
//			System.out.print("Odd number of tokens in csv file! Can't read file.");
//		}
//		else 
		{
//			array3d = new double [numTokens/2][numTokens][numSlices]; //+1 to save #AAs of each vector
			array3d = new double [numLines][numTokens][numSlices]; //+1 to save #AAs of each vector
//			System.out.println("aas.length: "+aas.length);
			try {
				bufRdr = new BufferedReader(new FileReader(filename));
				String line = null;
				String token = null;
				int row = 0;
				int col = 0;	
//				int aaindex = 0, phi=0, theta=0;
//				System.out.println("aaindex: "+aaindex+" phi: "+phi+" theta: "+theta);
				
				//read each line of text file
				while((line = bufRdr.readLine()) != null){
					col = 0;
					StringTokenizer st = new StringTokenizer(line,",");
//					System.out.print("row: "+row+"\t");
					if (st.countTokens()>2){
						while (st.hasMoreTokens()){
							//get next token and store it in the array
							token = st.nextToken();
							x = (row%((numLines)+1))-1;
							y = col;
//							System.out.print(theta+","+phi+","+aaindex+"\t");
							array3d[x][y][z] = Double.valueOf(token);	
//							System.out.print(Double.valueOf(token)+"\t");
//							array3d[theta][phi][aas.length] += Double.valueOf(token);
							col++;				
						}
					}
					row++;
					if (row % ((numLines)+1) == 0){
						z++;
					}	
//					System.out.println();
				}
				
				//close the file
				bufRdr.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}		
		}
		
		return array3d;
	}
    
    public double[][][] readCSVfile3Ddouble(ZipFile zipfile, ZipEntry zipentry) throws NumberFormatException, IOException{
    	double [][][] array3d = null;
		int numTokens = 0;
		int numSlices = 0;
		int numLines = 0;
		int x=0, y=0, z=0;
//		int xDim = 0, yDim = 0;
		BufferedReader bufRdr;
		try {
//			bufRdr = new BufferedReader(new FileReader(filename));
			bufRdr = new BufferedReader(new	InputStreamReader(zipfile.getInputStream(zipentry)));
			String line = null;
			String token = null;	
			// get number of slices (arrays) from first line
			line = bufRdr.readLine(); // just one token in first line
			StringTokenizer st = new StringTokenizer(line,",");
			if (st.countTokens()>1){
				token = st.nextToken();
				token = st.nextToken();
				numSlices = Integer.valueOf(token);
//				System.out.println("number of slices: " + numSlices);
			}
			// get dimension of array
			// xDim = numLines, yDim = numTokens
			boolean xDimFound = false;
			while((line = bufRdr.readLine()) != null && !xDimFound)
			{
				st = new StringTokenizer(line,",");
				if (st.countTokens()>2){
					numLines++;
//					xDim++;
					numTokens = st.countTokens();
//					yDim = st.countTokens();
				}		
				else
					xDimFound = true;
			}
//			System.out.println("number of tokens: "+numTokens);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
//		if (numTokens%2!=0){
//			System.out.print("Odd number of tokens in csv file! Can't read file.");
//		}
//		else 
		{
//			array3d = new double [numTokens/2][numTokens][numSlices]; //+1 to save #AAs of each vector
			array3d = new double [numLines][numTokens][numSlices]; //+1 to save #AAs of each vector
//			System.out.println("aas.length: "+aas.length);
			try {
//				bufRdr = new BufferedReader(new FileReader(filename));
				bufRdr = new BufferedReader(new	InputStreamReader(zipfile.getInputStream(zipentry)));
				String line = null;
				String token = null;
				int row = 0;
				int col = 0;	
//				int aaindex = 0, phi=0, theta=0;
//				System.out.println("aaindex: "+aaindex+" phi: "+phi+" theta: "+theta);
				
				//read each line of text file
				while((line = bufRdr.readLine()) != null){
					col = 0;
					StringTokenizer st = new StringTokenizer(line,",");
//					System.out.print("row: "+row+"\t");
					if (st.countTokens()>2){
						while (st.hasMoreTokens()){
							//get next token and store it in the array
							token = st.nextToken();
							x = (row%((numLines)+1))-1;
							y = col;
//							System.out.print(theta+","+phi+","+aaindex+"\t");
							array3d[x][y][z] = Double.valueOf(token);	
//							System.out.print(Double.valueOf(token)+"\t");
//							array3d[theta][phi][aas.length] += Double.valueOf(token);
							col++;				
						}
					}
					row++;
					if (row % ((numLines)+1) == 0){
						z++;
					}	
//					System.out.println();
				}
				
				//close the file
				bufRdr.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}		
		}
		
		return array3d;
	}
    
    public int[][][] readCSVfile3Dint(String filename) throws NumberFormatException, IOException{
    	int [][][] array3d = null;
		int numTokens = 0;
		int numSlices = 0;
		int numLines = 0;
		int x=0, y=0, z=0;
//		int xDim = 0, yDim = 0;
		BufferedReader bufRdr;
		try {
			bufRdr = new BufferedReader(new FileReader(filename));
			String line = null;
			String token = null;	
			// get number of slices (arrays) from first line
			line = bufRdr.readLine(); // just one token in first line
			StringTokenizer st = new StringTokenizer(line,",");
			if (st.countTokens()>1){
				token = st.nextToken();
				token = st.nextToken();
				numSlices = Integer.valueOf(token);
				System.out.println("number of slices: " + numSlices);
			}
			// get dimension of array
			// xDim = numLines, yDim = numTokens
			boolean xDimFound = false;
			while((line = bufRdr.readLine()) != null && !xDimFound)
			{
				st = new StringTokenizer(line,",");
				if (st.countTokens()>2){
					numLines++;
//					xDim++;
					numTokens = st.countTokens();
//					yDim = st.countTokens();
				}		
				else
					xDimFound = true;
			}
			System.out.println("number of tokens: "+numTokens);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}	
//		if (numTokens%2!=0){
//			System.out.print("Odd number of tokens in csv file! Can't read file.");
//		}
//		else 
		{
//			array3d = new double [numTokens/2][numTokens][numSlices]; //+1 to save #AAs of each vector
			array3d = new int [numLines][numTokens][numSlices]; //+1 to save #AAs of each vector
//			System.out.println("aas.length: "+aas.length);
			try {
				bufRdr = new BufferedReader(new FileReader(filename));
				String line = null;
				String token = null;
				int row = 0;
				int col = 0;	
//				int aaindex = 0, phi=0, theta=0;
//				System.out.println("aaindex: "+aaindex+" phi: "+phi+" theta: "+theta);
				
				//read each line of text file
				while((line = bufRdr.readLine()) != null){
					col = 0;
					StringTokenizer st = new StringTokenizer(line,",");
//					System.out.print("row: "+row+"\t");
					if (st.countTokens()>2){
						while (st.hasMoreTokens()){
							//get next token and store it in the array
							token = st.nextToken();
							x = (row%((numLines)+1))-1;
							y = col;
//							System.out.print(theta+","+phi+","+aaindex+"\t");
							array3d[x][y][z] = Integer.valueOf(token);	
//							System.out.print(Double.valueOf(token)+"\t");
//							array3d[theta][phi][aas.length] += Double.valueOf(token);
							col++;				
						}
					}
					row++;
					if (row % ((numLines)+1) == 0){
						z++;
					}	
//					System.out.println();
				}				
				//close the file
				bufRdr.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}		
		}
		
		return array3d;		
	}
    
    public int[][] getSumValues(){
    	return this.sumValues;
    }

}

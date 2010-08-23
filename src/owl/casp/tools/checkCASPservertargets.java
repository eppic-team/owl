package owl.casp.tools;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.Vector;

import javax.swing.JFileChooser;


/**
 * This class can be used to check CASP server models for their correct atom number ordering.
 * Reordered atoms are saved as new files within a subfolder "SortedAtoms" of the respective target.
 */
public class checkCASPservertargets {

	private static final String folderCorrectedFiles = "SortedAtoms";

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		// all targets within sub-folders of "/Volumes/StruPPi/CASP8/server_models/"
		runAllCaspTargets();
		
		// run single file
//		checkSingleFile();
		
//		// all targets within the folder:
//		String inputDir = "/Volumes/StruPPi/CASP8/server_models/T0423/";
//		// make copies of all files within subdirectory first


		//		checkForFileCopies(inputDir);
//		runAllTargetsWithinFolder(inputDir);
	}
	
	@SuppressWarnings("unused")
	private static void checkSingleFile() {
		// TODO Auto-generated method stub
		String filename = null;
		final JFileChooser chooser = new JFileChooser("Choose directory");
        chooser.setDialogType(JFileChooser.OPEN_DIALOG);
        chooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);

        chooser.setVisible(true);
        final int result = chooser.showOpenDialog(null);

        if (result == JFileChooser.APPROVE_OPTION) {
            File inputVerzFile = chooser.getSelectedFile();
            filename = inputVerzFile.getPath();
        }
//        System.out.println("Abbruch");
        chooser.setVisible(false); 
        
        System.out.println(filename);        
        
        String filenameOut = filename+"_sortedAtoms";
        final JFileChooser chooserOut = new JFileChooser("Choose directory");
        chooserOut.setDialogType(JFileChooser.OPEN_DIALOG);
        chooserOut.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);

        chooserOut.setVisible(true);
        final int resultOut = chooserOut.showOpenDialog(null);

        if (resultOut == JFileChooser.APPROVE_OPTION) {
            File inputVerzFile = chooserOut.getSelectedFile();
            filenameOut = inputVerzFile.getPath();
        }
//        System.out.println("Abbruch");
        chooserOut.setVisible(false); 
        
        System.out.println(filenameOut);	
		
		try {
			checkAndRewriteServerModel(filename, filenameOut);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void checkForFileCopies(String dir) throws IOException {
		String filenameIn = dir;
    	String filenameOut = dir + folderCorrectedFiles+"/";
	    boolean test = false;    	
		if (!(new File(filenameOut).exists())){
			File srcPath = new File(filenameIn);
			File dstPath = new File(filenameOut);
			test = copyDirectory(srcPath, dstPath);
			System.out.println("Files Copied into subdirectory "+filenameOut+"  "+test);
		}    		
	}
	
	private static void runAllTargetsWithinFolder(String dir) throws IOException{
		File folder = new File(dir);    	
	    File[] listOfFiles = folder.listFiles();
	    
	    for (int i=0;i<listOfFiles.length; i++)	{
	    	if (listOfFiles[i].isFile()){
		    	String proteinId = listOfFiles[i].getName();
		    	String filename = dir; //+"/"+proteinId;
		    	if (dir.charAt(dir.length()-1) != '/')
		    		filename += "/";
		    	String filenameIn = filename + proteinId;
		    	String filenameOut = filename + folderCorrectedFiles+"/"+proteinId;
		    	
//		    	System.out.println(filename);
		    	checkAndRewriteServerModel(filenameIn, filenameOut);	    		
	    	}
	    }
	}

	@SuppressWarnings("unused")
	private static void runAllCaspTargets() throws IOException{
		String scoreType = "";
		
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
	    /*
	     * This is the file where the output of the scores will be generated.
	     */
		   	    
	    File folder = new File(inputDir);
	    String[] listOfFolders = folder.list();
	    for (int i=0;i<listOfFolders.length; i++){
	    	File file = new File(inputDir+listOfFolders[i]);
	    	if (file.isDirectory()){
	    		targetNames.add(listOfFolders[i]);
	    		inputDirs.add(inputDir+listOfFolders[i]+"/");
	    		inputFiles.add(inputDir+listOfFolders[i]+"_"+targetListFN+".txt");
	    	}
	    }
	    int iF=0;
	    for (String dir : inputDirs){
//	    	if(targetNames.get(iF).equals("T0471"))
	    	{	    		    		
		    	inputFile = inputFiles.get(iF);
		    	folder = new File(dir);
		    	String target = targetNames.get(iF); // "T0389";
		    	
		    	// make copies of all files within subdirectory first
	    		checkForFileCopies(dir);
	    				    	
			    File[] listOfFiles = folder.listFiles();
			    PrintWriter out = new PrintWriter(new FileWriter(inputFile));
			    for (int i=0;i<listOfFiles.length; i++)	{			    	
//			    	String name =listOfFiles[i].getName();
			    	if (listOfFiles[i].isFile())
			    		out.println(listOfFiles[i].getName());
		    	}
			    out.close();	    	
		    	System.out.println("Element" + iF + ": " + target+" "+dir);
		    	System.out.println(inputFiles.get(iF));	
		    	
		    	BufferedReader in = new BufferedReader(new FileReader(inputFile)); 
				String proteinId; 
				
				while ((proteinId = in.readLine()) != null){
					String filename = dir;
			    	String filenameIn = filename + proteinId;
			    	String filenameOut = filename + folderCorrectedFiles+"/"+proteinId;
//					System.out.println(filename);	
					
					checkAndRewriteServerModel(filenameIn, filenameOut);
				}
				
	    	}
	    	iF++;
	    }
	}
	
	public static void checkAndRewriteServerModel(String filenameIn, String filenameOut) throws IOException{
		// read file into vector of string arrays
//		Vector<String> serverMod = loadData(filename);
		TreeMap<Integer,String> serverModel = loadDataTM(filenameIn);
//		testData(filename);
		// write to file		
		writeData(serverModel, filenameOut);
		// close file writer
	}
	
	 public static boolean copyDirectory(File srcPath, File dstPath)    throws IOException{
		if (srcPath.isDirectory()){
			if (!dstPath.exists()){
				dstPath.mkdir();
			}
			String files[] = srcPath.list();
			
			for(int i = 0; i < files.length; i++){
				File f = new File(srcPath, files[i]);
//				if (files[i].equals(folderCorrectedFiles))
//					System.out.println(srcPath+files[i]);
				if (f.isFile()){
					copyDirectory(new File(srcPath, files[i]), 
							new File(dstPath, files[i]));						
				}		
			}
			
		}	
		else{
			
			if(!srcPath.exists()){			
				System.out.println("File or directory does not exist.");
				System.exit(0);			
				return false;
			}			
			else{			
				InputStream in = new FileInputStream(srcPath);
				OutputStream out = new FileOutputStream(dstPath); 
				// Transfer bytes from in to out
				byte[] buf = new byte[1024];
				
				int len;
				
				while ((len = in.read(buf)) > 0) {				
					out.write(buf, 0, len);				
				}				
				in.close();				
				out.close();			
			}
		
		}
		
//		System.out.println("Directory copied.");		
		return true;		
	}
	
//	public static void copyFile(File srcFile, File destFile) throws IOException 
//	{
//	        InputStream oInStream = new FileInputStream(srcFile);
//	        OutputStream oOutStream = new FileOutputStream(destFile);
//	 
//	        // Transfer bytes from in to out
//	        byte[] oBytes = new byte[1024];
//	        int nLength;
//	        BufferedInputStream oBuffInputStream = new BufferedInputStream( oInStream );
//	        while ((nLength = oBuffInputStream.read(oBytes)) &gt; 0) 
//	        {
//	        	oOutStream.write(oBytes, 0, nLength);
//	        }
//	        oInStream.close();
//	        oOutStream.close();
//	}
	
	private static void writeData(TreeMap<Integer,String> tm, String filename) throws IOException{
//		FileWriter writer = null;
//		writer = new FileWriter(filename);
		PrintWriter out = new PrintWriter(new FileWriter(filename));
		/*
	      get Collection of values contained in TreeMap using
	      Collection values() method of TreeMap class
	    */
	    Collection<String> c = tm.values();
	 
	    //obtain an Iterator for Collection
	    Iterator<String> itr = c.iterator();
	 
	    //iterate through TreeMap values iterator
	    while(itr.hasNext()){
//	      System.out.println(itr.next());
	      out.println(itr.next().toString());
//	      writer.append(itr.next().toString());
//	      writer.append('\n');
	    }
	    out.close();
//	    writer.flush();
//		writer.close();
		System.out.println("Target written to: "+filename);
	}
	
	@SuppressWarnings("unused")
	private static void testData(String filename) throws IOException{
		Vector<String> vector = new Vector<String>();
		BufferedReader bufRdr;
		bufRdr = new BufferedReader(new FileReader(filename));
		String line=null, lastline;
		int numHeaderLines = 0;
		while((line = bufRdr.readLine()) != null){	
			if (!(line.startsWith("ATOM") || line.startsWith("Atom"))){
				numHeaderLines++;
			}
			else
				break;
		}
		//read each line of text file
		numHeaderLines = -1*numHeaderLines +1;
		bufRdr = new BufferedReader(new FileReader(filename));
		while((line = bufRdr.readLine()) != null){	
			vector.add(line);
		}
//		line = vector.elementAt(vector.size()-2);
		lastline = vector.lastElement();
//		System.out.println(lastline);
//		System.out.println(lastline);
		
		if (!(lastline.startsWith("END") || lastline.startsWith("End"))){
			System.out.println(filename);
		}
		
	}
	
	@SuppressWarnings("unused")
	private static TreeMap<Integer,String> loadData(String filename) throws IOException{
		TreeMap<Integer, String> map = new TreeMap<Integer,String>();
		BufferedReader bufRdr;
		bufRdr = new BufferedReader(new FileReader(filename));
		String line = null;
		int numHeaderLines = 0;
		while((line = bufRdr.readLine()) != null){	
			if (!(line.startsWith("ATOM") || line.startsWith("Atom"))){
				numHeaderLines++;
			}
			else
				break;
		}
		//read each line of text file
		numHeaderLines = -1*numHeaderLines +1;
		int idPostHeader = 100000;
		bufRdr = new BufferedReader(new FileReader(filename));
		while((line = bufRdr.readLine()) != null){	
			if (line.startsWith("ATOM") || line.startsWith("Atom")){
				String s = line.substring(6, 11).trim();
				int atomSer = Integer.valueOf(s);
				map.put(atomSer, line);
			}
			else if (line.startsWith("TER") || line.startsWith("Ter") || line.startsWith("END") || line.startsWith("End")){
				map.put(idPostHeader, line);
				idPostHeader++;
			}
			else {
				map.put(numHeaderLines++, line);
			}					
		}
		if (idPostHeader>100000){
			System.out.println("WrongOrder "+filename);
			writeData(map, filename);
		}
		return map;
	}
	
	private static TreeMap<Integer,String> loadDataTM(String filename) throws IOException{
		TreeMap<Integer, String> map = new TreeMap<Integer,String>();
		BufferedReader bufRdr;
		bufRdr = new BufferedReader(new FileReader(filename));
		String line = null;
		int numHeaderLines = 0;
		while((line = bufRdr.readLine()) != null){	
			if (!(line.startsWith("ATOM") || line.startsWith("Atom"))){
				numHeaderLines++;
			}
			else
				break;
		}
		//read each line of text file
		numHeaderLines = -1*numHeaderLines;
		bufRdr = new BufferedReader(new FileReader(filename));
		while((line = bufRdr.readLine()) != null){	
			if (line.startsWith("ATOM") || line.startsWith("Atom")){
				String s = line.substring(6, 11).trim();
				int atomSer = Integer.valueOf(s);
				map.put(atomSer, line);
				numHeaderLines = atomSer+1;
			}
			else {
				map.put(numHeaderLines++, line);
			}					
		}
		return map;
	}

}

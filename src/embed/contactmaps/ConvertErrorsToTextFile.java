package embed.contactmaps;

import java.io.*;

import owl.core.util.RegexFileFilter;


/**
 * Helper class, reads all kind of text files, to extract both the CMError and DMError and converts it
 * to a String. This String can be written to a text file and graphically displayed with additional programs.
 * <p>
 * Types of accepted files: having '.cmap', 'indi', or 'ust' extension, if the String denoting the directory does
 * not contain any of these file types, this class will issue an <code>IllegalArgumentException</code>.
 * </p>
 * @author gmueller
 *
 */
public class ConvertErrorsToTextFile {
	
	/*------------------------------------constants---------------------------*/
	/**
	 * File type, that are accepted by the method <code>{@link ConvertErrorsToTextFile#containsReadableFiles(File)}</code>
	 */
	private static final String[] file_types = {".*.cmap",".*.indi",".*.ust"};
	
	/*------------------------------------members-----------------------------*/
	
	/**
	 * String, denoting the directory, that will be processed
	 */
	private static String path;
	
	/**
	 * String, denoting the file type, must be one of the values specified by the constant <code>{@link #file_types}</code>.
	 */	
	private static String file_type;
	
	/**
	 * the PDB code of all files within the directory
	 */
	private static String pdb_code;
	
	/*--------------------------------fields--------------------------------------*/
	
	/**
	 * the directory, presented by this instance
	 */
	private File dir;
	
	/**
	 * list of files contained in the directory <code>{@link #dir}</code>
	 */
	private File[] file_list;
	
	/**
	 * the content: this is all error values sorted as follows:
	 * <p>
	 * number of file in the list (tab) CMError (tab) DMError
	 * </p>
	 */
	private String file_content;
	
	/*----------------------------------constructors-------------------------------------------------*/
	
	/**
	 * one parameter constructor: constructs a String instance by reading all files contained in the
	 * directory <code>dir_path</code>.
	 */
	public ConvertErrorsToTextFile(String dir_path) throws FileNotFoundException, IOException{
		setFiles(dir_path);
	}
	
	/*--------------------------------------------setter--------------------------------------------*/
	
	/**
	 * Setter: initializes a File instance as a directory and checks, whether or not it contains any
	 * readable files. These readable files are then listed in an array of File instances and the
	 * error values are extracted from them. 
	 */
	public void setFiles (String dir_path) throws FileNotFoundException, IOException{
		setDirectory(dir_path);
		path = new String (dir_path);
		file_list = dir.listFiles(new RegexFileFilter (file_type));
		int length = file_list.length;
		file_content = "";
		for(int i = 0; i < length; i++){
			file_content += fileReader(i);
		}
	}
	
	/**
	 * Initializes the field <code>{@link #dir}</code>. The field is only initialized when the String <code>dir_path</code>
	 * <p>
	 * </p>
	 * <p>
	 * 1. denotes a directory
	 * </p>
	 * <p>
	 * </p>
	 * <p>
	 * 2. contains readable files as specified by <code>{@link #containsReadableFiles(File)}</code>
	 * </p>
	 * @param dir_path
	 * @throws IllegalArgumentException 
	 * <p>
	 * </p>
	 * <p>
	 * a) if the directory does not exits
	 * </p>
	 * <p>
	 * b) the path does not denote a directory
	 * </p>
	 * <p>
	 * c) does not contain any readable files
	 * </p>
	 * 
	 */
	public void setDirectory (String dir_path) throws IllegalArgumentException {
		File dirs = new File (dir_path);
		if(containsReadableFiles(dirs)){
			dir = new File (dir_path);
		}
		else{
			throw new IllegalArgumentException ("The String dir_path = '"+ dir_path+"', denoting a directory, does not contain any readable file...");
		}
	}
	
	/*------------------------------------auxiliaries------------------------------*/
	
	/**
	 * Method writing the field <code>{@link #file_content}</code> to a sub directory, which
	 * is by default: <code>{@link #path} + Error/</code>
	 */
	public void writeToFile () throws IOException{
		if(file_content != null){
			String output_dir = path+"/Error_file/";
			File output = new File (output_dir);
			if(!output.exists()){
				output.mkdirs();
			}
			output_dir += pdb_code+".txt";
			FileOutputStream str = new FileOutputStream (output_dir);
			PrintStream printa = new PrintStream (str);
			printa.print(file_content);
			printa.close();
			str.close();
			System.out.println(output_dir + " written to file...");
		}
	}
	
	/**
	 * Method writing the field <code>{@link #file_content}</code> to a sub directory, which
	 * is by default: <code>{@link #path} + Error/</code>
	 */
	public void writeToFile (String dir_, int i) throws IOException{
		String dir = new String (dir_);
		if(file_content != null && dir != null){
			dir += "/Error_file/";
			File output = new File (dir);
			if(!output.exists()){
				output.mkdirs();
			}
			dir += pdb_code+i+".txt";
			FileOutputStream str = new FileOutputStream (dir);
			PrintStream printa = new PrintStream (str);
			printa.print(file_content);
			printa.close();
			str.close();
			System.out.println(dir + " written to file...");
		}
	}
	
	/*------------------------------------getters-----------------------------------*/
	
	/**
	 * File reader method (auxiliary): reads header of the <code>i</code>-th file in the
	 * field <code>{@link #file_list}</code>, which contains the CMError and DMError, as well
	 * as the PDB code.
	 * @param i - the i-th file
	 * @return a String representation containing the the number, CMError and DMError
	 */
	public String fileReader (int i) throws FileNotFoundException, IOException {
		if(file_list[i] != null){
			BufferedReader reader = new BufferedReader (new FileReader(file_list[i]));
			String linereader = "", content = "";
			boolean tester = false;
			while((linereader = reader.readLine())!= null){
				if(linereader.contains("#")){
					if(linereader.contains("CMError")){
						content += i + "\t" + linereader.split(": ")[1];
						if(!tester){
							tester = true;
						}
						else{
							break;
						}
					}
					else{
						if(linereader.contains("PDB") && !linereader.contains("CHAIN")){
							pdb_code = new String (linereader.split(": ")[1]);
						}
						else{
							if(linereader.contains("DMError")){
								content += "\t" + linereader.split(": ")[1] + "\n";
								if(!tester){
									tester = true;
								}
								else{
									break;
								}
							}
						}
					}
				}
				else{
					break;
				}
			}
			reader.close();
			return content;
		}
		else{
			throw new NullPointerException ("The field 'file_list' was not initialized...");
		}
	}
	
	/*------------------------------------statics------------------------------------*/
	
	/**
	 * auxiliary setter: sets the file type contained in <code>{@link #dir}</code>.
	 */
	public static void setFileType (File dirs) throws IllegalArgumentException {
		if(dirs.exists() && dirs.isDirectory() && file_type == null){
			for(int i = 0; i < 3; i++){
				File[] filesys = dirs.listFiles(new RegexFileFilter (file_types[i]));
				if(filesys.length > 0){
					file_type = new String (file_types[i]);
					break;
				}
			}
		}
		else{
			if(!dirs.isDirectory()){
				throw new IllegalArgumentException ("The path '"+dirs.getAbsolutePath()+"' does not denote not a directory.");
			}
			else{
				if(!dirs.exists()){
					throw new IllegalArgumentException ("The denoted directory '"+dirs.getAbsolutePath()+"' does not exist.");
				}
			}
		}
	}
	
	public static boolean containsReadableFiles (File dirs){
		if(file_type == null){
			setFileType(dirs);
		}
		if(file_type != null){
			File[] list = dirs.listFiles(new RegexFileFilter (file_type));
			if(list.length > 0){
				return true;
			}
			else{
				return false;
			}
		}
		else{
			return false;
		}		
	}
	
	public static String[] listSubDirectories (String dir_path){
		File abstr = new File(dir_path);
		if(abstr.exists() && abstr.isDirectory()){
			File[] super1 = abstr.listFiles();
			String content = "";
			int length1 = super1.length;
			for(int i = 0; i < length1; i++){
				if(super1[i].isDirectory()){
					if(!containsReadableFiles(super1[i])){
						File[] super2 = super1[i].listFiles(new RegexFileFilter ("evo.*."));
						int length2 = super2.length;
						if(length2 > 0){
							for(int j = 0; j < length2; j++){
								if(super2[j].isDirectory()){
									if(containsReadableFiles(super2[j])){
										content += super2[j].getAbsolutePath() + "\n";
									}
								}
							}
						}
						
					}
					else{
						content += super1[i].getAbsolutePath() + "\n";
					}
				}
			}
			return content.split("\n");
		}
		else{
			if(abstr.exists()){
				throw new IllegalArgumentException ("The denoted path is not a directory.");
			}
			else{
				throw new IllegalArgumentException ("The denoted path does not exist.");
			}
		}
	}
		
	public static File[] removeFileInstance (File[] file_sys, int index){
		int length = file_sys.length;
		if(length > 0 && (index >=0 && index < length)){
			File[] new_files = new File[length - 1];
			for(int i = 0; i < length; i++){
				if(i < index){
					new_files[i] = new File (file_sys[i].getAbsolutePath());
				}
				else{
					new_files[i] = new File (file_sys[i].getAbsolutePath());
				}
			}
			return new_files;
		}
		else{
			if(length <= 0){
				throw new IllegalArgumentException ("The array of file instances had zero entries.");
			}
			else{
				throw new IllegalArgumentException ("The index parameter exceeded the bounds opf the file array.");
			}
		}
	}
	
	public static void main (String[] args) throws FileNotFoundException, IOException{
		String dir = "/project/StruPPi/gabriel/Arbeiten/run_011009/1odd/";
		String[] list = listSubDirectories(dir);
		int length = list.length;
		if(length > 0){
			ConvertErrorsToTextFile[] converter = new ConvertErrorsToTextFile[length]; 
			for(int i = 0; i < length; i++){
				converter[i] = new ConvertErrorsToTextFile(list[i]);
				converter[i].writeToFile(dir,i);
			}
		}
	}

}

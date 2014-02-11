package owl.tests;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Properties;
import java.util.zip.GZIPInputStream;

/**
 * Class containing static methods necessary for tests setup.
 * @author duarte_j
 *
 */
public class TestsSetup {

	public static final String DEFAULT_PATHS_FILE = "/owl/tests/owl_test_paths.dat";
	private static final InputStream PATHS_FILE_IS = TestsSetup.class.getResourceAsStream(DEFAULT_PATHS_FILE);
	public static final File HOME_PATHS_FILE = new File(System.getProperty("user.home"),"owl_test_paths.dat") ; 
	
	/**
	 * Reads the path files with paths needed to run the tests in owl.
	 * First reads the global file in src/tests and then reads the paths file
	 * from the user's home directory overriding any already read values. 
	 * @return a Properties object containing all pairs of key/values read from the file
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static Properties readPaths() throws FileNotFoundException, IOException {
		
		Properties p = new Properties();
		p.load(PATHS_FILE_IS);
				
		if (HOME_PATHS_FILE.canRead()) {
			p.load(new FileInputStream(HOME_PATHS_FILE));
		}
		
		return p;
	}
	
	public static File gunzipFile(File repoGzFile) throws FileNotFoundException {
		if (!repoGzFile.exists()) {
			throw new FileNotFoundException("PDB repository file "+repoGzFile+" could not be found.");
		}
		File unzippedFile = null;
		try {
			String prefix = repoGzFile.getName().substring(0,repoGzFile.getName().lastIndexOf(".gz"));
			unzippedFile = File.createTempFile(prefix,"");
			unzippedFile.deleteOnExit();

			GZIPInputStream zis = new GZIPInputStream(new FileInputStream(repoGzFile));
			FileOutputStream os = new FileOutputStream(unzippedFile);
			int b;
			while ( (b=zis.read())!=-1) {
				os.write(b);
			}
			zis.close();
			os.close();
		} catch (IOException e) {
			System.err.println("Couldn't uncompress "+repoGzFile+" file into "+unzippedFile);
			System.err.println(e.getMessage());
			System.exit(1);
		}
		return unzippedFile;
	}
	
	/**
	 * Copies content of InputStream into a temporary file with the given fileName.
	 * The temporary file will be located in system temp dir and removed after JVM finalisation
	 * @param is
	 * @param prefix
	 * @param suffix
	 * @return the temp file where the data was copied to
	 * @throws IOException
	 */
	public static File inputStreamToTempFile(InputStream is, String prefix, String suffix) throws IOException {
		File file = File.createTempFile(prefix,suffix);
		file.deleteOnExit();
		
		OutputStream os = new FileOutputStream(file);
		
		byte[] buf = new byte[1024];
		int len;
		while ((len = is.read(buf)) > 0){
			os.write(buf, 0, len);
		}
		
		os.close();
		
		return file;
	}
}

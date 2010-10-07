package owl.tests;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

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
}

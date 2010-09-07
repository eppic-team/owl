package owl.tests.core.structure;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Properties;

//import javax.vecmath.Matrix4d;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.xml.sax.SAXException;

import owl.core.connections.pisa.PisaConnection;
import owl.core.structure.ChainInterface;
import owl.core.structure.ChainInterfaceList;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.util.MySQLConnection;
import owl.tests.TestsSetup;


public class PdbAsymUnitTest {

	private static String NACCESS_EXEC; 
	
	private static final String TESTDATADIR = "src/owl/tests/core/structure/data";
	private static final String LISTFILE = TESTDATADIR+"/testset_interfaces.txt";
	
	private static final String PDBASE_DB = "pdbase";
	
	private static final String PISA_INTERFACES_URL = "http://www.ebi.ac.uk/msd-srv/pisa/cgi-bin/interfaces.pisa?";
	
	private static final double CUTOFF = 6.0;

	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		Properties p = TestsSetup.readPaths();
		NACCESS_EXEC = p.getProperty("NACCESS_EXEC");
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	
	@Test
	public void testGetAllInterfaces() throws IOException, SQLException, SAXException {
		MySQLConnection conn = new MySQLConnection();

		List<String> pdbCodes = new ArrayList<String>();
		BufferedReader flist = new BufferedReader(new FileReader(LISTFILE));
		String line;
		while ((line = flist.readLine() ) != null ) {
			if (line.startsWith("#")) continue;
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			pdbCodes.add(pdbCode);
		}
		flist.close();

		// getting PISA interfaces
		PisaConnection pc = new PisaConnection(PISA_INTERFACES_URL, null, null);
		System.out.println("Downloading PISA interfaces");
		Map<String, ChainInterfaceList> all = pc.getInterfacesDescription(pdbCodes);

		for (String pdbCode: pdbCodes) {
		
			ChainInterfaceList pisaInterfaces = all.get(pdbCode);
			// we sort them on interface area because pisa doesn't always sort them like that (it does some kind of grouping)
			pisaInterfaces.sort();
			
			System.out.println("\n##"+pdbCode);

			PdbAsymUnit pdb = null;
			try {
				pdb = new PdbAsymUnit(pdbCode, conn, PDBASE_DB);
			} catch (PdbCodeNotFoundError e) {
				System.err.println("Missing PDB code "+pdbCode);
				continue;
			} catch (PdbLoadError e) {
				System.err.println("PDB load error, cause: "+e.getMessage());
				continue;
			}

			System.out.println(pdb.getSpaceGroup().getShortSymbol()+" ("+pdb.getSpaceGroup().getId()+")");
			
			long start = System.currentTimeMillis();
			ChainInterfaceList interfaces = pdb.getAllInterfaces(CUTOFF, new File(NACCESS_EXEC));
			long end = System.currentTimeMillis();
			System.out.println("Time: "+((end-start)/1000)+"s");
			System.out.println("Total number of interfaces found: "+interfaces.size());
			
			int pisaCount = interfaces.getNumProtProtInterfaces();
			System.out.println("PISA interface count: "+pisaCount);
			Assert.assertEquals(pisaCount, interfaces.size());

			int i = interfaces.size()-1;
			for (int p=pisaInterfaces.size()-1;p>=0;p--) {
				ChainInterface pisaInterf = pisaInterfaces.get(p);
				if (!pisaInterf.isProtein()) continue;
				System.out.println("\nInterface "+(interfaces.size()-i));
				ChainInterface myInterf = interfaces.get(i);
				
				
				// we allow for a 10% discrepancy from PISA
				System.out.printf("Areas, pisa: %8.2f\tmy: %8.2f\tdiff: %4.1f%%\n",pisaInterf.getInterfaceArea(),myInterf.getInterfaceArea()/2.0,
						(pisaInterf.getInterfaceArea()-myInterf.getInterfaceArea()/2.0)*100.0/pisaInterf.getInterfaceArea());
				Assert.assertEquals(pisaInterf.getInterfaceArea(), myInterf.getInterfaceArea()/2.0, pisaInterf.getInterfaceArea()*0.10);
				i--;
			}
			
		}
	}

}

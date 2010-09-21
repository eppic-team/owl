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
import owl.core.structure.Pdb;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbCodeNotFoundError;
import owl.core.structure.PdbLoadError;
import owl.core.structure.Residue;
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
			
			int pisaCount = pisaInterfaces.getNumProtProtInterfaces();
			System.out.println("PISA interface count: "+pisaCount);
			Assert.assertEquals(pisaCount, interfaces.size());

			int i = 0;
			for (int p=0;p<pisaInterfaces.size();p++) {
				ChainInterface pisaInterf = pisaInterfaces.get(p);
				if (!pisaInterf.isProtein()) continue;
				System.out.println("\nInterface "+(i+1));
				ChainInterface myInterf = interfaces.get(i);
				
				
				// total interface area, we allow for a 20% discrepancy from PISA
				double MARGIN = 0.20;
				System.out.printf("Areas, pisa: %8.2f\tmy: %8.2f\tdiff: %4.1f%%\n",pisaInterf.getInterfaceArea(),myInterf.getInterfaceArea()/2.0,
						(pisaInterf.getInterfaceArea()-myInterf.getInterfaceArea()/2.0)*100.0/pisaInterf.getInterfaceArea());
				Assert.assertEquals(pisaInterf.getInterfaceArea(), myInterf.getInterfaceArea()/2.0, pisaInterf.getInterfaceArea()*0.10);
				
				// bsas of individual residues, we allow for a 20% discrepancy from PISA
				System.out.println("Chain 1");
				int a1 = 0;
				int b1 = 0;
				int t1 = 0;
				for (Residue residue:pisaInterf.getFirstMolecule().getResidues().values()) {
					Pdb myFirstMol = myInterf.getFirstMolecule();
					if (!myFirstMol.getPdbChainCode().equals(pisaInterf.getFirstMolecule().getPdbChainCode())) {
						myFirstMol = myInterf.getSecondMolecule();
					}
					int resser = myFirstMol.getResSerFromPdbResSer(residue.getPdbSerial());
					double pisaAsa = residue.getAsa();
					double pisaBsa = residue.getBsa();
					Residue myRes = myFirstMol.getResidue(resser);
					if (myRes!=null) {
						double myAsa = myRes.getAsa();
						double myBsa = myRes.getBsa();
						System.out.printf("%s\t%s\t%d\t%s\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n",
								residue.getPdbSerial(),residue.getAaType().getThreeLetterCode(),resser,myRes.getAaType().getThreeLetterCode(),pisaAsa,pisaBsa,myAsa,myBsa);
						Assert.assertEquals(residue.getAaType(),myRes.getAaType());
						if (deltaComp(pisaAsa, myAsa, pisaAsa*MARGIN)) a1++;
						if (deltaComp(pisaBsa, myBsa, pisaBsa*MARGIN)) b1++;
						//Assert.assertEquals(pisaAsa, myAsa, pisaAsa*MARGIN);
						//Assert.assertEquals(pisaBsa, myBsa, pisaBsa*MARGIN);
						t1++;
					}
				}
				System.out.println("Total: "+t1+". Agreements within "+String.format("%4.2f",MARGIN)+" tolerance: ASA "+a1+" BSA "+b1);
				//Assert.assertTrue(a1>(0.80*(double)t1));
				//Assert.assertTrue(b1>(0.80*(double)t1));
				System.out.println("Chain 2");
				int a2 = 0;
				int b2 = 0;
				int t2 = 0;
				for (Residue residue:pisaInterf.getSecondMolecule().getResidues().values()) {
					Pdb mySecondMol = myInterf.getSecondMolecule();
					if (!mySecondMol.getPdbChainCode().equals(pisaInterf.getSecondMolecule().getPdbChainCode())) {
						mySecondMol = myInterf.getFirstMolecule();
					}
					int resser = mySecondMol.getResSerFromPdbResSer(residue.getPdbSerial());
					double pisaAsa = residue.getAsa();
					double pisaBsa = residue.getBsa();
					Residue myRes = mySecondMol.getResidue(resser);
					if (myRes!=null) {
						double myAsa = myRes.getAsa();
						double myBsa = myRes.getBsa();
						System.out.printf("%s\t%s\t%d\t%s\t%6.2f\t%6.2f\t%6.2f\t%6.2f\n",
								residue.getPdbSerial(),residue.getAaType().getThreeLetterCode(),resser,myRes.getAaType().getThreeLetterCode(),pisaAsa,pisaBsa,myAsa,myBsa);
						Assert.assertEquals(residue.getAaType(),myRes.getAaType());
						if (deltaComp(pisaAsa, myAsa, pisaAsa*MARGIN)) a2++;
						if (deltaComp(pisaBsa, myBsa, pisaBsa*MARGIN)) b2++;
						t2++;
						//Assert.assertEquals(pisaAsa, myAsa, pisaAsa*MARGIN);
						//Assert.assertEquals(pisaBsa, myBsa, pisaBsa*MARGIN);
					}
				}
				System.out.println("Total: "+t2+". Agreements within "+String.format("%4.2f",MARGIN)+" tolerance: ASA "+a2+" BSA "+b2);
				//Assert.assertTrue(a2>(0.80*(double)t2));
				//Assert.assertTrue(b2>(0.80*(double)t2));

				i++;
			}
			
		}
	}
	
	public boolean deltaComp(double a, double b, double delta) {
		return (Math.abs(a-b)<delta);
	}

}

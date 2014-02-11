package owl.core.connections;


//import java.io.BufferedReader;
//import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
//import java.util.Collection;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import owl.core.connections.SiftsConnection;
//import owl.core.connections.NoMatchFoundException;
//import owl.core.features.Feature;
//import owl.core.features.FeatureType;
import owl.core.features.InvalidFeatureCoordinatesException;
import owl.core.features.OverlappingFeatureException;
import owl.core.features.SiftsFeature;
//import owl.core.structure.PdbChain;
//import owl.core.structure.PdbAsymUnit;
//import owl.core.structure.PdbCodeNotFoundException;
//import owl.core.structure.PdbLoadException;
//import owl.core.util.MySQLConnection;

public class SiftsConnectionTest {
	
	
//	private static final String LISTFILE = "src/owl/tests/core/structure/data/cullpdb_20";
//	private static final String PDBASE_DB = "pdbase";

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
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
	public void testSiftsConnection() throws IOException, SQLException, InvalidFeatureCoordinatesException, OverlappingFeatureException {
//		MySQLConnection conn = new MySQLConnection();
		SiftsConnection siftsConn = new SiftsConnection(SiftsConnection.PDB2UNIPROT_URL);
//		BufferedReader flist = new BufferedReader(new FileReader(LISTFILE));
//		String line;
//		while ((line = flist.readLine() ) != null ) {
//			if (line.startsWith("#")) continue;
//			String pdbCode = line.split("\\s+")[0].toLowerCase();
//			String pdbChainCode = line.split("\\s+")[1];
//			try {
//				PdbAsymUnit fullpdb = new PdbAsymUnit(pdbCode,conn,PDBASE_DB);
//				PdbChain pdb = fullpdb.getChain(pdbChainCode);
//
//				Collection<SiftsFeature> mappings = siftsConn.getMappings(pdbCode, pdbChainCode);
//				//System.out.println(pdbCode+pdbChainCode+" "+mappings.size());
//				for (SiftsFeature mapping:mappings) {
//					// this will throw exceptions InvalidFeatureCoordinates or OverlappingFeature
//					// it just should not happen. If it does it's a test failure (that's why we throw the exceptions)
//					pdb.addFeature(mapping); 
//				}
//				
//				Collection<Feature> fs = pdb.getFeaturesOfType(FeatureType.SIFTS);
//				Assert.assertTrue(fs.size()>0);
//				if (fs.size()>1) {
//					System.out.print(pdbCode+pdbChainCode+"\t");
//					for (Feature f:fs) {
//						SiftsFeature sf = (SiftsFeature) f;
//						System.out.print(sf.getUniprotId() +"\t"+ sf.getIntervalSet()+"\t");
//					}
//
//					System.out.print("("+pdb.getFullLength()+")\t("+pdb.getMinObsResSerial()+"-"+pdb.getMaxObsResSerial()+")");
//					System.out.println();
//				}
//				
//			} catch (PdbCodeNotFoundException e) {
//				System.err.println(e.getMessage());
//				continue;
//			} catch (NoMatchFoundException e) {
//				System.err.println(e.getMessage());
//				continue;
//			} catch (PdbLoadException e) {
//				System.err.println(e.getMessage());
//			}
//		}
//		flist.close();
		int count=0;
		for (ArrayList<SiftsFeature> mappings:siftsConn.getAllMappings()) {
			Assert.assertTrue(mappings.size()>0);
			//if (mappings.size()>1) {
			//	System.out.println(mappings.get(0).getPdbCode()+mappings.get(0).getPdbChainCode()+" "+mappings.size());
			//}
			count++;
		}
		System.out.println("Total mappings: "+count);
	}
	

}

package tests.proteinstructure;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import proteinstructure.AAinfo;
import proteinstructure.FileFormatError;
import proteinstructure.CiffilePdb;
import proteinstructure.Pdb;
import proteinstructure.PdbCodeNotFoundError;
import proteinstructure.PdbLoadError;
import proteinstructure.PdbaseInconsistencyError;
import proteinstructure.PdbasePdb;
import proteinstructure.PdbfilePdb;
import proteinstructure.SecondaryStructure;
import tools.MySQLConnection;

public class PdbParsersTest {
	
	
	private static final String CIFDIR="/project/StruPPi/BiO/DBd/PDB-REMEDIATED/data/structures/unzipped/all/mmCIF/";
	private static final String PDBDIR="/project/StruPPi/BiO/DBd/PDB-REMEDIATED/data/structures/unzipped/all/pdb/";
	private static final String LISTFILE = "tests/proteinstructure/data/cullpdb_20";
	private static final String PDBASE_DB = "pdbase";
	private static final String MYSQLSERVER = "talyn";
	
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
	public void testCIFagainstPDBASE() throws FileFormatError, PdbaseInconsistencyError, SQLException, IOException {

		MySQLConnection conn = new MySQLConnection(MYSQLSERVER, PDBASE_DB);
		
		BufferedReader flist = new BufferedReader(new FileReader(LISTFILE));
		String line;
		while ((line = flist.readLine() ) != null ) {
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			String pdbChainCode = line.split("\\s+")[1];
			
			String message = "Failed for "+pdbCode+pdbChainCode;
			
			System.out.println(pdbCode+" "+pdbChainCode);
			
			Pdb ciffilePdb = null;
			Pdb pdbasePdb = null;

			try {
				
				ciffilePdb = new CiffilePdb(new File(CIFDIR,pdbCode+".cif"));
				ciffilePdb.load(pdbChainCode);
				
				pdbasePdb = new PdbasePdb(pdbCode, PDBASE_DB, conn);
				pdbasePdb.load(pdbChainCode);
				
				// asserting
				
				// getChains/getModels
				String[] ciffileChains = ciffilePdb.getChains(); 
				String[] pdbaseChains = pdbasePdb.getChains();
				Integer[] ciffileModels = ciffilePdb.getModels(); 
				Integer[] pdbaseModels = pdbasePdb.getModels();
				Assert.assertTrue(ciffileChains.length==pdbaseChains.length);
				Assert.assertTrue(ciffileModels.length==pdbaseModels.length);
				HashSet<String> ciffileChainsSet = new HashSet<String>(Arrays.asList(ciffileChains));
				HashSet<String> pdbaseChainsSet = new HashSet<String>(Arrays.asList(pdbaseChains));
				for (String chain:pdbaseChains) {
					Assert.assertTrue(ciffileChainsSet.contains(chain)); 
				}
				for (String chain:ciffileChains) {
					Assert.assertTrue(pdbaseChainsSet.contains(chain)); 
				}
				HashSet<Integer> ciffileModelsSet = new HashSet<Integer>(Arrays.asList(ciffileModels));
				HashSet<Integer> pdbaseModelsSet = new HashSet<Integer>(Arrays.asList(pdbaseModels));
				for (int model:pdbaseModels) {
					Assert.assertTrue(ciffileModelsSet.contains(model));
				}
				for (int model:ciffileModels) {
					Assert.assertTrue(pdbaseModelsSet.contains(model));
				}
				
				// identifiers
				Assert.assertEquals(pdbasePdb.getPdbCode(), ciffilePdb.getPdbCode());
				Assert.assertEquals(pdbasePdb.getChainCode(), ciffilePdb.getChainCode());
				Assert.assertEquals(pdbasePdb.getPdbChainCode(), ciffilePdb.getPdbChainCode());
				Assert.assertEquals(pdbasePdb.getModel(), ciffilePdb.getModel());

				// sequences
				Assert.assertEquals(pdbasePdb.getSequence(), ciffilePdb.getSequence());
				Assert.assertEquals(pdbasePdb.getObsSequence(), ciffilePdb.getObsSequence());
				
				// lengths
				Assert.assertEquals(pdbasePdb.getFullLength(), ciffilePdb.getFullLength());
				Assert.assertEquals(pdbasePdb.getObsLength(), ciffilePdb.getObsLength());
				Assert.assertEquals(pdbasePdb.getNumAtoms(), ciffilePdb.getNumAtoms());
				
				// info from atom serials 
				for (int atomser:pdbasePdb.getAllAtomSerials()) {
					Assert.assertEquals(pdbasePdb.getAtomCoord(atomser), ciffilePdb.getAtomCoord(atomser));
					Assert.assertEquals(pdbasePdb.getResSerFromAtomSer(atomser), ciffilePdb.getResSerFromAtomSer(atomser));
					Assert.assertEquals(pdbasePdb.getAtomNameFromAtomSer(atomser), ciffilePdb.getAtomNameFromAtomSer(atomser));
				}
				
				for (int atomser:ciffilePdb.getAllAtomSerials()) {
					Assert.assertEquals(ciffilePdb.getAtomCoord(atomser), pdbasePdb.getAtomCoord(atomser));
					Assert.assertEquals(ciffilePdb.getResSerFromAtomSer(atomser), pdbasePdb.getResSerFromAtomSer(atomser));
				}

				// info from residue serials
				for (int resser:pdbasePdb.getAllSortedResSerials()) {
					Assert.assertEquals(pdbasePdb.getPdbResSerFromResSer(resser), ciffilePdb.getPdbResSerFromResSer(resser));
					String pdbresser = pdbasePdb.getPdbResSerFromResSer(resser);
					Assert.assertEquals(pdbasePdb.getResSerFromPdbResSer(pdbresser), ciffilePdb.getResSerFromPdbResSer(pdbresser));
								
					Assert.assertEquals(pdbasePdb.getResidue(resser).getAaType(),ciffilePdb.getResidue(resser).getAaType());
					Assert.assertEquals(pdbasePdb.getResTypeFromResSerial(resser), ciffilePdb.getResTypeFromResSerial(resser));
				}

				for (int resser:ciffilePdb.getAllSortedResSerials()) {
					Assert.assertEquals(ciffilePdb.getPdbResSerFromResSer(resser), pdbasePdb.getPdbResSerFromResSer(resser));
					String pdbresser = ciffilePdb.getPdbResSerFromResSer(resser);
					Assert.assertEquals(ciffilePdb.getResSerFromPdbResSer(pdbresser), pdbasePdb.getResSerFromPdbResSer(pdbresser));
					
					Assert.assertEquals(ciffilePdb.getResTypeFromResSerial(resser), pdbasePdb.getResTypeFromResSerial(resser));
				}
				
				// secondary structure
				SecondaryStructure pdbaseSS = pdbasePdb.getSecondaryStructure();
				SecondaryStructure ciffileSS = ciffilePdb.getSecondaryStructure();
				Assert.assertEquals(message,pdbaseSS.getNumElements(), ciffileSS.getNumElements());

				for (int resser:pdbasePdb.getAllSortedResSerials()) {
					String resserMsg = "Failed for "+pdbCode+pdbChainCode+" and resser "+resser; 
					Assert.assertEquals(resserMsg,pdbaseSS.getSecStrucElement(resser), ciffileSS.getSecStrucElement(resser));
					// checking that the 2 ways of accessing the sec struct element coincide
					Assert.assertSame(pdbaseSS.getSecStrucElement(resser), pdbasePdb.getResidue(resser).getSsElem());
					Assert.assertEquals(resserMsg,pdbaseSS.getSecStrucElement(resser),pdbasePdb.getResidue(resser).getSsElem());
					Assert.assertSame(ciffileSS.getSecStrucElement(resser), ciffilePdb.getResidue(resser).getSsElem());
					Assert.assertEquals(resserMsg,ciffileSS.getSecStrucElement(resser),ciffilePdb.getResidue(resser).getSsElem());
				}
				
				
			}
			catch (PdbCodeNotFoundError e) {
				System.err.println("pdb code not found in pdbase");
			}
			catch (PdbLoadError e) {
				System.err.println("pdb load error, cause: "+e.getMessage());
			}

			
		}
		
		flist.close();

	}
	
	@Test
	public void testPdbfileParser() throws IOException, SQLException {
		BufferedReader flist = new BufferedReader(new FileReader(LISTFILE));
		String line;
		while ((line = flist.readLine() ) != null ) {
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			//String pdbChainCode = line.split("\\s+")[1];
			
			System.out.println(pdbCode);
			
			Pdb pdbfilePdb = null;
			Pdb pdbasePdb = null;
			
			try {
				pdbfilePdb = new PdbfilePdb(new File(PDBDIR,"pdb"+pdbCode+".ent").getAbsolutePath());
				pdbasePdb = new PdbasePdb(pdbCode);
				String[] chains = pdbfilePdb.getChains();
				// test getChains/getModels
				Assert.assertTrue(chains.length>0);
				Integer[] models = pdbfilePdb.getModels();
				Assert.assertTrue(models.length>0);
				Assert.assertFalse(pdbfilePdb.isDataLoaded());
				
				// loading all chains and all models
				for (int model:models) {
					System.out.println(model);
					for (String chain:chains) {
						pdbasePdb.load(chain);
						System.out.println(chain);
						pdbfilePdb.load(chain,model);
						Assert.assertTrue(pdbfilePdb.isDataLoaded());
						
						// pdbCode properly read
						Pattern p = Pattern.compile("^\\d\\w\\w\\w$");
						Matcher m = p.matcher(pdbfilePdb.getPdbCode());
						Assert.assertTrue(m.matches());
						Assert.assertEquals(pdbasePdb.getPdbCode(), pdbfilePdb.getPdbCode());
						
						// chain codes coincide
						Assert.assertEquals(chain,pdbfilePdb.getChainCode());
						Assert.assertEquals(pdbfilePdb.getChainCode(), pdbfilePdb.getPdbChainCode());
						Assert.assertEquals(pdbasePdb.getPdbChainCode(), pdbfilePdb.getPdbChainCode());
						// model as input
						Assert.assertEquals(model,pdbfilePdb.getModel());
						
						// sequence
						Assert.assertTrue(pdbfilePdb.getObsSequence().length()<=pdbfilePdb.getSequence().length());
						Assert.assertTrue(pdbfilePdb.getObsLength()==pdbfilePdb.getObsSequence().length());
						Assert.assertTrue(pdbfilePdb.getFullLength()==pdbfilePdb.getSequence().length());
						// we can't assert the sequences against pdbase because of some very weird entries. See http://pdbwiki.org/index.php/1ejg
						//Assert.assertEquals(pdbasePdb.getSequence(), pdbfilePdb.getSequence());
						Assert.assertEquals(pdbasePdb.getObsSequence(), pdbfilePdb.getObsSequence());
						String seq = pdbfilePdb.getSequence();
						for (int resser:pdbfilePdb.getAllSortedResSerials()) {
							Assert.assertEquals(pdbfilePdb.getResTypeFromResSerial(resser), 
												AAinfo.oneletter2threeletter(String.valueOf(seq.charAt(resser-1))));
							// at least 1 atom per observed residue
							Assert.assertTrue(pdbfilePdb.getResidue(resser).getNumAtoms()>0);
							// at least the CA atom must be present, doesn't work for all, e.g. 2cioB residue 80
							//Assert.assertNotNull(pdbfilePdb.getResidue(resser).getAtom("CA"));
						}

						// info from atom serials
						for (int atomser:pdbfilePdb.getAllAtomSerials()) {
							Assert.assertNotNull(pdbfilePdb.getAtomCoord(atomser));
							Assert.assertNotNull(pdbfilePdb.getResSerFromAtomSer(atomser));
						}

						// sec structure
						SecondaryStructure ss = pdbfilePdb.getSecondaryStructure();
						Assert.assertNotNull(ss);
						Assert.assertEquals(pdbfilePdb.getSequence(), ss.getSequence());

						for (int resser:pdbfilePdb.getAllSortedResSerials()) {
							String resserMsg = "Failed for "+pdbCode+chain+" and resser "+resser;
							// checking that the 2 ways of accessing the sec struct element coincide
							Assert.assertSame(ss.getSecStrucElement(resser), pdbfilePdb.getResidue(resser).getSsElem());
							Assert.assertEquals(resserMsg,ss.getSecStrucElement(resser),pdbfilePdb.getResidue(resser).getSsElem());
						}

					}
				}		
				
			} catch (PdbLoadError e) {
				System.err.println("pdb load error, cause: "+e.getMessage());
			} catch (PdbCodeNotFoundError e) {
				System.err.println("pdb code not found in pdbase");			
			}
		}
		flist.close();
	}
}

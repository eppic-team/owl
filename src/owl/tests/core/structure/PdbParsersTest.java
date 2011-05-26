package owl.tests.core.structure;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import junit.framework.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import owl.core.structure.CiffileParser;
import owl.core.structure.CrystalCell;
import owl.core.structure.HetResidue;
import owl.core.structure.NucResidue;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbCodeNotFoundException;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbaseParser;
import owl.core.structure.PdbfileParser;
import owl.core.structure.Residue;
import owl.core.structure.AaResidue;
import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;
import owl.core.util.FileFormatException;
import owl.core.util.MySQLConnection;
import owl.core.util.RegexFileFilter;


public class PdbParsersTest {
	
	// paths to the mirror of the PDB ftp repository with all pdb/mmCIF compressed files
	private static final String CIFDIR="/nfs/data/dbs/pdb/data/structures/all/mmCIF/";
	private static final String PDBDIR="/nfs/data/dbs/pdb/data/structures/all/pdb/";
	
	private static final String DATADIR = "src/owl/tests/core/structure/data";
	
	private static final String LISTFILE = DATADIR+"/cullpdb_20";
	private static final String PDBASE_DB = "pdbase";
	
	private static final String TINKERPDBFILE = DATADIR+"/1bxyA_tinker.pdb";
	
	private static final String CASPTARBALL = DATADIR+"/T0515.tar.gz";
	
	private static final String TMPDIR = System.getProperty("java.io.tmpdir");
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		// to throw away jaligner logging
		System.setProperty("java.util.logging.config.file","/dev/null");
		Process proc = Runtime.getRuntime().exec("tar -zxf "+CASPTARBALL+" -C "+System.getProperty("java.io.tmpdir"));
		proc.waitFor();
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testSpecialPDBFiles() throws PdbLoadException, IOException, FileFormatException {
		// testing some special PDB files
		
		// tinker file (no occupancy or bfactors)
		System.out.println(TINKERPDBFILE);
		PdbAsymUnit fullpdb = new PdbAsymUnit(new File(TINKERPDBFILE));
		Assert.assertNotNull(fullpdb);
		Assert.assertNotNull(fullpdb.getPdbCode());
		
		// CASP TS server files (testing with server tar ball T0515 of CASP9)
		File dir = new File(TMPDIR,"T0515");
		File[] files = dir.listFiles(new RegexFileFilter("T0515.*"));
		
		for (File caspFileName:files) {
			System.out.println(caspFileName);
			PdbfileParser parser = new PdbfileParser(caspFileName.getAbsolutePath());
			Integer[] models = parser.getModels();
			 
			
			for (int model:models) {
				try {
					fullpdb = new PdbAsymUnit(caspFileName,model);
					Assert.assertNotNull(fullpdb.getPdbCode());
				} catch (PdbLoadException e) {
					System.err.println("Warning, pdb load exception: "+e.getMessage());
				} 
				
			}
		}
	}
	
	@Test
	public void testCIFagainstPDBASE() throws FileFormatException, SQLException, IOException {

		MySQLConnection conn = new MySQLConnection();
		ArrayList<String> warnings = new ArrayList<String>();
		BufferedReader flist = new BufferedReader(new FileReader(LISTFILE));
		String line;
		while ((line = flist.readLine() ) != null ) {
			if (line.startsWith("#")) continue;
			if (line.trim().isEmpty()) break;
			String[] tokens = line.split("\\s+");
			String pdbCode = tokens[0].toLowerCase();
			
			String[] pdbChainCodes = null;
			if (tokens.length>1) {
				pdbChainCodes = new String[1];
				pdbChainCodes[0]= tokens[1];
			} 
			System.out.println(pdbCode);
			
			PdbChain ciffilePdb = null;
			PdbChain pdbasePdb = null;

			try {
				File cifFile = unzipFile(new File(CIFDIR,pdbCode+".cif.gz"));

				// getChains/getModels
				PdbaseParser pdbaseParser = new PdbaseParser(pdbCode,PDBASE_DB,conn);
				CiffileParser ciffileParser = new CiffileParser(cifFile);
				String[] ciffileChains = ciffileParser.getChains(); 
				String[] pdbaseChains = pdbaseParser.getChains();
				Integer[] ciffileModels = ciffileParser.getModels(); 
				ciffileParser.closeFile();
				Integer[] pdbaseModels = pdbaseParser.getModels();
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

				// parsing the full files
				PdbAsymUnit ciffileFullPdb = new PdbAsymUnit(cifFile);
				PdbAsymUnit pdbaseFullPdb = new PdbAsymUnit(pdbCode,conn, PDBASE_DB);

				Assert.assertEquals(pdbaseFullPdb.getNumChains(),ciffileFullPdb.getNumChains());
				Assert.assertEquals(pdbaseFullPdb.getNumAtoms(),ciffileFullPdb.getNumAtoms());
				// checking we parse exactly same chain names with both
				Set<String> cifchainCodes = ciffileFullPdb.getChainCodes();
				for (String chainCode:pdbaseFullPdb.getChainCodes()) {
					Assert.assertTrue(cifchainCodes.contains(chainCode));
				}
				Set<String> pdbasechainCodes = ciffileFullPdb.getChainCodes();
				for (String chainCode:ciffileFullPdb.getChainCodes()) {
					Assert.assertTrue(pdbasechainCodes.contains(chainCode));
				}
				
				for (String chainCode:pdbasechainCodes) {
					ciffilePdb = ciffileFullPdb.getChainForChainCode(chainCode);
					pdbasePdb = pdbaseFullPdb.getChainForChainCode(chainCode);
					String pdbChainCode = pdbasePdb.getPdbChainCode(); 
					
					String message = "Failed for "+pdbCode+pdbChainCode;				
					System.out.println(chainCode+" "+pdbChainCode);

					// poly/non-poly
					Assert.assertEquals(pdbasePdb.isNonPolyChain(),ciffilePdb.isNonPolyChain());
					
					// title
					Assert.assertEquals(pdbasePdb.getParent().getTitle(),ciffilePdb.getParent().getTitle());
					// crystal data
					Assert.assertEquals(pdbasePdb.getParent().getSpaceGroup(), ciffilePdb.getParent().getSpaceGroup());
					CrystalCell pdbaseCell = pdbasePdb.getParent().getCrystalCell();
					CrystalCell ciffileCell = ciffilePdb.getParent().getCrystalCell();
					if (pdbaseCell!=null && ciffileCell!=null) { // nulls will happen for NMR entries
						Assert.assertEquals(pdbaseCell.getA(), ciffileCell.getA(), 0.001);
						Assert.assertEquals(pdbaseCell.getB(), ciffileCell.getB(), 0.001);
						Assert.assertEquals(pdbaseCell.getC(), ciffileCell.getC(), 0.001);
						Assert.assertEquals(pdbaseCell.getAlpha(), ciffileCell.getAlpha(), 0.001);
						Assert.assertEquals(pdbaseCell.getBeta(), ciffileCell.getBeta(), 0.001);
						Assert.assertEquals(pdbaseCell.getGamma(), ciffileCell.getGamma(), 0.001);
					}

					// exp data and quality parameters
					Assert.assertEquals(pdbasePdb.getParent().getExpMethod(),ciffilePdb.getParent().getExpMethod());
					Assert.assertEquals(pdbasePdb.getParent().getResolution(),ciffilePdb.getParent().getResolution(),0.0001);
					Assert.assertEquals(pdbasePdb.getParent().getRfree(),ciffilePdb.getParent().getRfree(),0.0001);
					Assert.assertEquals(pdbasePdb.getParent().getRsym(),ciffilePdb.getParent().getRsym(),0.0001);

					// identifiers
					Assert.assertEquals(pdbasePdb.getPdbCode(), ciffilePdb.getPdbCode());
					Assert.assertEquals(pdbasePdb.getChainCode(), ciffilePdb.getChainCode());
					Assert.assertEquals(pdbasePdb.getPdbChainCode(), ciffilePdb.getPdbChainCode());
					Assert.assertEquals(pdbasePdb.getParent().getModel(), ciffilePdb.getParent().getModel());

					// sequences
					if (!pdbasePdb.isNonPolyChain()) {
						Assert.assertEquals(pdbasePdb.getSequence().getSeq(), ciffilePdb.getSequence().getSeq());
						Assert.assertEquals(pdbasePdb.getObsSequence(), ciffilePdb.getObsSequence());

						// lengths
						Assert.assertEquals(pdbasePdb.getFullLength(), ciffilePdb.getFullLength());
						Assert.assertEquals(pdbasePdb.getObsLength(), ciffilePdb.getObsLength());
						Assert.assertEquals(pdbasePdb.getStdAaObsLength(), ciffilePdb.getStdAaObsLength());
					}
					Assert.assertEquals(pdbasePdb.getNumAtoms(), ciffilePdb.getNumAtoms());

					if (pdbasePdb.isNonPolyChain()) {
						Assert.assertNull(pdbasePdb.getSequence());
						Assert.assertNull(ciffilePdb.getSequence());
						Assert.assertNull(pdbasePdb.getSecondaryStructure());
						Assert.assertNull(ciffilePdb.getSecondaryStructure());
						
					}
					
					// has alt codes
					Assert.assertEquals(pdbasePdb.hasAltCodes(),ciffilePdb.hasAltCodes());

					// info from atom serials 
					for (int atomser:pdbasePdb.getAllAtomSerials()) {
						Assert.assertEquals(pdbasePdb.getAtomCoord(atomser), ciffilePdb.getAtomCoord(atomser));
						Assert.assertEquals(pdbasePdb.getAtom(atomser).getOccupancy(), ciffilePdb.getAtom(atomser).getOccupancy(),0.001);
						Assert.assertEquals(pdbasePdb.getAtom(atomser).getBfactor(), ciffilePdb.getAtom(atomser).getBfactor(),0.001);
						Assert.assertEquals(pdbasePdb.getResSerFromAtomSer(atomser), ciffilePdb.getResSerFromAtomSer(atomser));
					}

					for (int atomser:ciffilePdb.getAllAtomSerials()) {
						Assert.assertEquals(ciffilePdb.getAtomCoord(atomser), pdbasePdb.getAtomCoord(atomser));
						Assert.assertEquals(ciffilePdb.getResSerFromAtomSer(atomser), pdbasePdb.getResSerFromAtomSer(atomser));
					}

					// info from residue serials in polymer chains
					if (!pdbasePdb.isNonPolyChain()) {
						boolean protein = false;
						boolean nucleic = false;
						for (int resser:pdbasePdb.getAllResSerials()) {
							Assert.assertEquals(pdbasePdb.getPdbResSerFromResSer(resser), ciffilePdb.getPdbResSerFromResSer(resser));
							String pdbresser = pdbasePdb.getPdbResSerFromResSer(resser);
							Assert.assertEquals(pdbasePdb.getResSerFromPdbResSer(pdbresser), ciffilePdb.getResSerFromPdbResSer(pdbresser));
							Assert.assertEquals(pdbasePdb.getResidue(resser).getLongCode(),ciffilePdb.getResidue(resser).getLongCode());
							Residue pdbaseRes = pdbasePdb.getResidue(resser);
							Residue ciffileRes = ciffilePdb.getResidue(resser);

							Assert.assertEquals(resser, pdbaseRes.getSerial());
							Assert.assertEquals(resser, ciffileRes.getSerial());

							Assert.assertEquals(pdbaseRes.isPeptideLinked(), ciffileRes.isPeptideLinked());
							if (pdbaseRes.isPeptideLinked()) {
								Assert.assertEquals(pdbaseRes.getShortCode(), pdbasePdb.getSequence().getSeq().charAt(resser-1));
								Assert.assertEquals(ciffileRes.getShortCode(), ciffilePdb.getSequence().getSeq().charAt(resser-1));
							}


							if (pdbaseRes instanceof AaResidue) {
								protein = true;
								Assert.assertEquals(((AaResidue)pdbaseRes).getAaType(),((AaResidue)ciffileRes).getAaType());
								Assert.assertTrue(pdbaseRes.isPeptideLinked());
								Assert.assertTrue(pdbaseRes.getSerial()>0 && pdbaseRes.getSerial()<=pdbasePdb.getFullLength());
							} else if (pdbaseRes instanceof HetResidue) {
								Assert.assertEquals(((HetResidue)pdbaseRes).getLongCode(),((HetResidue)ciffileRes).getLongCode());
								Assert.assertNotSame(((HetResidue)pdbaseRes).getLongCode(),"HOH");
								Assert.assertNotSame(((HetResidue)ciffileRes).getLongCode(),"HOH");
								if (pdbasePdb.getSequence().isProtein()) {
									Assert.assertTrue(pdbaseRes.isPeptideLinked());
									Assert.assertTrue(pdbaseRes.getSerial()<=pdbasePdb.getFullLength());
									Assert.assertEquals('X',pdbasePdb.getSequence().getSeq().charAt(resser-1));
								} else if (pdbasePdb.getSequence().isNucleotide()) {
									Assert.assertFalse(pdbaseRes.isPeptideLinked());
								}
							} else if (pdbaseRes instanceof NucResidue) {
								nucleic = true;
								Assert.assertEquals(((NucResidue)pdbaseRes).getNucType(),((NucResidue)ciffileRes).getNucType());
								Assert.assertFalse(pdbaseRes.isPeptideLinked());
								Assert.assertTrue(pdbaseRes.getSerial()>0 && pdbaseRes.getSerial()<=pdbasePdb.getFullLength());
							}
						}
						Assert.assertEquals(protein, !nucleic); // one must be true and the other false (a chain can't mix both prot and nucleic)
						Assert.assertEquals(pdbasePdb.getSequence().isProtein(), protein);
						Assert.assertEquals(pdbasePdb.getSequence().isNucleotide(), nucleic);
						Assert.assertEquals(ciffilePdb.getSequence().isProtein(), protein);
						Assert.assertEquals(ciffilePdb.getSequence().isNucleotide(), nucleic);
						for (int resser:ciffilePdb.getAllResSerials()) {
							Assert.assertEquals(ciffilePdb.getPdbResSerFromResSer(resser), pdbasePdb.getPdbResSerFromResSer(resser));
							String pdbresser = ciffilePdb.getPdbResSerFromResSer(resser);
							Assert.assertEquals(ciffilePdb.getResSerFromPdbResSer(pdbresser), pdbasePdb.getResSerFromPdbResSer(pdbresser));

						}

						// secondary structure
						SecondaryStructure pdbaseSS = pdbasePdb.getSecondaryStructure();
						SecondaryStructure ciffileSS = ciffilePdb.getSecondaryStructure();
						Assert.assertEquals(message,pdbaseSS.getNumElements(), ciffileSS.getNumElements());

						for (int resser:pdbasePdb.getAllResSerials()) {
							String resserMsg = "Failed for "+pdbCode+pdbChainCode+" and resser "+resser; 
							SecStrucElement pdbaseSsElem = pdbaseSS.getSecStrucElement(resser);
							SecStrucElement ciffileSsElem = ciffileSS.getSecStrucElement(resser);
							Assert.assertEquals(resserMsg,pdbaseSsElem, ciffileSsElem);
							// checking that the 2 ways of accessing the sec struct element coincide
							Residue pdbaseRes = pdbasePdb.getResidue(resser);
							Residue ciffileRes = ciffilePdb.getResidue(resser);

							Assert.assertSame(pdbaseSS.getSecStrucElement(resser), pdbaseRes.getSsElem());
							Assert.assertEquals(resserMsg,pdbaseSS.getSecStrucElement(resser),pdbaseRes.getSsElem());
							Assert.assertSame(ciffileSS.getSecStrucElement(resser), ciffileRes.getSsElem());
							Assert.assertEquals(resserMsg,ciffileSS.getSecStrucElement(resser), ciffileRes.getSsElem());

							if ((pdbaseRes instanceof NucResidue) || ((pdbaseRes instanceof HetResidue) && !pdbaseRes.isPeptideLinked())) {
								Assert.assertNull(pdbaseRes.getSsElem());
							}

						}

					}

					
					// non-polymer chains
					if (pdbasePdb.isNonPolyChain()) {
						for (int resser:pdbasePdb.getAllResSerials()) {
							Residue pdbaseRes = pdbasePdb.getResidue(resser);
							Residue ciffileRes = ciffilePdb.getResidue(resser);

							Assert.assertEquals(resser, pdbaseRes.getSerial());
							Assert.assertEquals(resser, ciffileRes.getSerial());
							Assert.assertEquals(Integer.parseInt(pdbaseRes.getPdbSerial()),pdbaseRes.getSerial());

							Assert.assertEquals(pdbaseRes.isPeptideLinked(), ciffileRes.isPeptideLinked());
							
							Assert.assertTrue(pdbaseRes.getSerial()>=0);
							PdbChain correspondingPolyChain = pdbaseFullPdb.getChain(pdbasePdb.getPdbChainCode());
							if (correspondingPolyChain!=null && !correspondingPolyChain.isNonPolyChain()) {
								String lastPdbSerial = correspondingPolyChain.getLastResidue().getPdbSerial();
								String firstPdbSerial = correspondingPolyChain.getFirstResidue().getPdbSerial();
								if (!Character.isDigit(lastPdbSerial.charAt(lastPdbSerial.length()-1))) lastPdbSerial = lastPdbSerial.substring(0, lastPdbSerial.length()-1);
								if (!Character.isDigit(firstPdbSerial.charAt(firstPdbSerial.length()-1))) firstPdbSerial = firstPdbSerial.substring(0, firstPdbSerial.length()-1);
								if(!(pdbaseRes.getSerial()>Integer.parseInt(lastPdbSerial) ||
									pdbaseRes.getSerial()<Integer.parseInt(firstPdbSerial))) {
									warnings.add(pdbCode+": residue "+pdbaseRes+" of non-poly chain "+pdbasePdb.getChainCode()+ " has residue serial within the range of residue serials of chain "+correspondingPolyChain.getChainCode());
								}
							}
							
							if (pdbaseRes instanceof AaResidue) {
								Assert.assertEquals(((AaResidue)pdbaseRes).getAaType(),((AaResidue)ciffileRes).getAaType());
							} else if (pdbaseRes instanceof HetResidue) {
								Assert.assertEquals(((HetResidue)pdbaseRes).getLongCode(),((HetResidue)ciffileRes).getLongCode());
								Assert.assertNotSame(((HetResidue)pdbaseRes).getLongCode(),"HOH");
								// can't test for isPeptideLinked here as the method is intended to detect peptideLinked HET residues when 
								// we don't know whether chain is polymer or not. In CIF/pdbase case we do know that this is a non-polymer
								//Assert.assertFalse(pdbaseRes.isPeptideLinked());

							} else if (pdbaseRes instanceof NucResidue) {
								Assert.assertEquals(((NucResidue)pdbaseRes).getNucType(),((NucResidue)ciffileRes).getNucType());
								Assert.assertFalse(pdbaseRes.isPeptideLinked());
							}

						}
					}


				}	
			} catch (FileNotFoundException e){
				System.err.println("File missing. "+e.getMessage());
			} catch (PdbCodeNotFoundException e) {
				System.err.println("pdb code not found in pdbase");
			} catch (PdbLoadException e) {
				System.err.println("pdb load error, cause: "+e.getMessage());
			}

			
		}
		
		flist.close();
		if (!warnings.isEmpty()) {
			System.err.println("WARNINGS: ");
			for (String warning:warnings) {
				System.err.println(warning);
			}
		}

	}
	
	@Test
	public void testPdbfileParser() throws IOException, SQLException, FileFormatException, PdbLoadException {
		MySQLConnection conn = new MySQLConnection();
		ArrayList<String> warnings = new ArrayList<String>();
		// testing a list of PDB files from PDB
		BufferedReader flist = new BufferedReader(new FileReader(LISTFILE));
		String line;
		while ((line = flist.readLine() ) != null ) {
			if (line.startsWith("#")) continue;
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			//String pdbChainCode = line.split("\\s+")[1];
			
			System.out.println(pdbCode);
			
			
			try {
				File pdbFile = unzipFile(new File(PDBDIR,"pdb"+pdbCode+".ent.gz"));
				PdbfileParser parser = new PdbfileParser(pdbFile.getAbsolutePath());
				Integer[] models = parser.getModels();

				Assert.assertTrue(models.length>0);
				
				// loading all chains and all models
				for (int model:models) {
					System.out.println(model);
					PdbAsymUnit pdbfileFullpdb = new PdbAsymUnit(pdbFile,model);
					PdbAsymUnit pdbaseFullpdb = new PdbAsymUnit(pdbCode,model,conn,"pdbase");
					
					Assert.assertEquals(pdbaseFullpdb.getNumPolyChains(),pdbfileFullpdb.getNumPolyChains());
					// we can't assert that number of non-poly chains are the same as the assignment in CIF differs from mine in many cases
					// the only thing we can say for sure is that CIF always assigns more non-poly chains than we do (because they split the non-polys into more groups)
					Assert.assertTrue(pdbaseFullpdb.getNumChains()>=pdbfileFullpdb.getNumChains());
				 
					// pdbCode properly read
					Pattern p = Pattern.compile("^\\d\\w\\w\\w$");
					Matcher m = p.matcher(pdbfileFullpdb.getPdbCode());
					Assert.assertTrue(m.matches());
					Assert.assertEquals(pdbaseFullpdb.getPdbCode(), pdbfileFullpdb.getPdbCode());

					// title 
					// the capitalization is not consistent between pdb file and cif, 
					// the spaces introduced between different lines are lost in pdb and thus also 
					// not consistent with cif
					// that's why we compare in lower case and stripping spaces
					Assert.assertEquals(pdbaseFullpdb.getTitle().toLowerCase().replaceAll(" ", ""), 
							pdbfileFullpdb.getTitle().toLowerCase().replaceAll(" ",""));

					// model as input
					Assert.assertEquals(model,pdbfileFullpdb.getModel());
					
					// crystal data
					if (pdbaseFullpdb.getSpaceGroup()!=null) {
						Assert.assertEquals(pdbaseFullpdb.getSpaceGroup().getId(),pdbfileFullpdb.getSpaceGroup().getId());
					}
					Assert.assertNotNull(pdbfileFullpdb.getCrystalCell());
					
					// exp data and quality parameters
					Assert.assertEquals(pdbaseFullpdb.getExpMethod(),pdbfileFullpdb.getExpMethod());
					Assert.assertEquals(pdbaseFullpdb.getResolution(),pdbfileFullpdb.getResolution(),0.01);
					if (Math.abs(pdbaseFullpdb.getRfree()-pdbfileFullpdb.getRfree())>0.01) {
						// we can't assert because there's no consistency between pdbase and pdbfile, e.g. 1p1x
						System.err.println(pdbCode+": rfrees don't agree, pdbase "+
								String.format("%4.2f", pdbaseFullpdb.getRfree())+", pdbfile "+String.format("%4.2f",pdbfileFullpdb.getRfree()));
						warnings.add(pdbCode+": rfrees don't agree, pdbase "+
								String.format("%4.2f", pdbaseFullpdb.getRfree())+", pdbfile "+String.format("%4.2f",pdbfileFullpdb.getRfree()));
					}
					Assert.assertEquals(pdbaseFullpdb.getRsym(),pdbfileFullpdb.getRsym(),0.001);

					Assert.assertEquals(pdbaseFullpdb.getNumAtoms(),pdbfileFullpdb.getNumAtoms());
					
					for (PdbChain pdbfilePdb:pdbfileFullpdb.getAllChains()) {
						String chainCode = pdbfilePdb.getChainCode();
						String pdbChainCode = pdbfilePdb.getPdbChainCode();
						System.out.print(chainCode+" "+pdbChainCode+" ");
						PdbChain pdbasePdb = pdbaseFullpdb.getChainForChainCode(chainCode);
						if (!pdbaseFullpdb.containsChainCode(chainCode)) {
							System.err.println("No pdbase chain for "+chainCode);
							continue;
						}
						System.out.print(pdbasePdb.getChainCode()+" "+pdbasePdb.getPdbChainCode());
						String agreement = " ";
						if (chainCode.equals(pdbasePdb.getChainCode())) agreement+=" ";
						else agreement+="x";
						if (pdbChainCode.equals(pdbasePdb.getPdbChainCode())) agreement+=" ";
						else agreement+="x";
						System.out.println(agreement);
						
						Assert.assertEquals(pdbasePdb.getPdbCode(), pdbfilePdb.getPdbCode());
						
						// non-poly coincide
						Assert.assertEquals(pdbasePdb.isNonPolyChain(),pdbfilePdb.isNonPolyChain());
						
						// chain codes coincide
						// we can only assert it for poly-chains, non-poly assignments won't coincide between us and CIF
						if (!pdbasePdb.isNonPolyChain()) {
							Assert.assertEquals(pdbasePdb.getPdbChainCode(),pdbfilePdb.getPdbChainCode());
							Assert.assertEquals(pdbasePdb.getChainCode(),pdbfilePdb.getChainCode());
						}
						
						Assert.assertEquals(pdbasePdb.hasAltCodes(),pdbfilePdb.hasAltCodes());
						
			
						if (!pdbfilePdb.isNonPolyChain()) {
							Assert.assertEquals(pdbasePdb.getFullLength(), pdbfilePdb.getFullLength());
							Assert.assertEquals(pdbasePdb.getObsLength(), pdbfilePdb.getObsLength());
							Assert.assertEquals(pdbasePdb.getSequence().getSeq(), pdbfilePdb.getSequence().getSeq());
							Assert.assertEquals(pdbasePdb.getSequence().isProtein(),pdbfilePdb.getSequence().isProtein());
							Assert.assertEquals(pdbasePdb.getObsSequence(), pdbfilePdb.getObsSequence());

							// sequence
							Assert.assertTrue(pdbfilePdb.getObsSequence().length()<=pdbfilePdb.getSequence().getLength());
							Assert.assertTrue(pdbfilePdb.getObsLength()==pdbfilePdb.getObsSequence().length());
							Assert.assertTrue(pdbfilePdb.getFullLength()==pdbfilePdb.getSequence().getLength());
						}

						
						boolean protein = false;
						boolean nucleic = false;
						for (int resser:pdbfilePdb.getAllResSerials()) {
							if (!pdbfilePdb.isNonPolyChain()) {
								if (!pdbasePdb.getPdbResSerFromResSer(resser).equals(pdbfilePdb.getPdbResSerFromResSer(resser))) {
									// In some cases the alignment is ambiguous, see 2nwr at residues 191 onwards:
									// SEQRES           151 LTERGTTFGYNNLVVDFRSLPIMKQWAKVIYDATHSVQLPGGLGDKSGGM    200
									//                      |||||||||||||||||||||||||||||||||||||||||      |||
									// ATOM             150 LTERGTTFGYNNLVVDFRSLPIMKQWAKVIYDATHSVQLPG------GGM    193
									// the cif file places the G on the left of the gap rather than on the right
									// as we do (and they are right looking at the 3D coords)
									// For these cases we are trying to rely on the order of the residue serials in ATOM (detecting whether they are shifted)
									// Anyway some cases are still ambiguous: see 2iyv where pdbase (cif) since to have gotten 
									// wrong the last HISs (or maybe not, didn't look at 3D) but pdbfile says the one with coords is 183, where pdbase says it's 182  
									// That's why here we only print a warning (anyway it is only in rare cases that it happen)
									System.err.println(pdbCode+pdbChainCode+": resser " +resser+" maps to pdbresser "
											+pdbasePdb.getPdbResSerFromResSer(resser)+" in pdbase and to pdbresser "
											+pdbfilePdb.getPdbResSerFromResSer(resser)+" in pdbfile");
									warnings.add(pdbCode+pdbChainCode+": resser " +resser+" maps to pdbresser "
											+pdbasePdb.getPdbResSerFromResSer(resser)+" in pdbase and to pdbresser "
											+pdbfilePdb.getPdbResSerFromResSer(resser)+" in pdbfile");
								}
							}

							
							if (pdbasePdb.getPdbChainCode().equals(pdbfilePdb.getPdbChainCode()) // to be sure we don't compare non-polys for which our chain mapping doesn't coincide 
									&& pdbasePdb.containsResidue(resser)) { // we have to check for case where mapping of pdbfile and pdbase from resser to pdbresser don't coincide (we print warnings above for that)
								Residue pdbaseRes = pdbasePdb.getResidue(resser);
								Residue pdbfileRes = pdbfilePdb.getResidue(resser);
								Assert.assertEquals(resser, pdbaseRes.getSerial());
								Assert.assertEquals(resser, pdbfileRes.getSerial());
								
								// at least 1 atom per observed residue
								Assert.assertTrue(pdbfileRes.getNumAtoms()>0);
								// at least the CA atom must be present, doesn't work for all, e.g. 2cioB residue 80
								//Assert.assertNotNull(pdbfileRes.getAtom("CA"));
								Assert.assertEquals(pdbaseRes.getNumAtoms(), pdbfileRes.getNumAtoms());
								Assert.assertEquals(pdbaseRes.getNumHeavyAtoms(), pdbfileRes.getNumHeavyAtoms());


								Assert.assertEquals(pdbaseRes.isPeptideLinked(), pdbfileRes.isPeptideLinked());
								
								if (!pdbasePdb.isNonPolyChain()) {
									Assert.assertEquals(pdbfileRes.getShortCode(), pdbfilePdb.getSequence().getSeq().charAt(resser-1));
								}
								

								if (pdbfileRes instanceof AaResidue) {
									protein = true;
									Assert.assertEquals(((AaResidue)pdbaseRes).getAaType(),((AaResidue)pdbfileRes).getAaType());	
									Assert.assertTrue(pdbfileRes.isPeptideLinked());
									Assert.assertTrue(pdbfileRes.getSerial()>0);
									if (!pdbfilePdb.isNonPolyChain()) Assert.assertTrue(pdbfileRes.getSerial()<=pdbfilePdb.getFullLength());
									if (!pdbfilePdb.isNonPolyChain()) Assert.assertTrue(pdbfilePdb.getSequence().isProtein());
								} else if (pdbfileRes instanceof HetResidue) {
									Assert.assertEquals(((HetResidue)pdbaseRes).getLongCode(),((HetResidue)pdbfileRes).getLongCode());
									Assert.assertNotSame(((HetResidue)pdbfileRes).getLongCode(),"HOH");
									if (!pdbfilePdb.isNonPolyChain() && pdbfilePdb.getSequence().isProtein()) Assert.assertTrue(pdbfileRes.isPeptideLinked());
								} else if (pdbfileRes instanceof NucResidue) {
									nucleic = true;
									Assert.assertEquals(((NucResidue)pdbaseRes).getNucType(),((NucResidue)pdbfileRes).getNucType());
									Assert.assertFalse(pdbfileRes.isPeptideLinked());
									Assert.assertTrue(pdbfileRes.getSerial()>0);
									if (!pdbfilePdb.isNonPolyChain()) {
										Assert.assertTrue(pdbfileRes.getSerial()<=pdbfilePdb.getFullLength());
										Assert.assertTrue(pdbfilePdb.getSequence().isNucleotide());
									}
									
								}								
								
							}

						}
						if (!pdbfilePdb.isNonPolyChain()) {
							Assert.assertEquals(protein, !nucleic); // one must be true and the other false (a chain can't mix both prot and nucleic)
							Assert.assertEquals(pdbasePdb.getSequence().isProtein(), protein);
							Assert.assertEquals(pdbasePdb.getSequence().isNucleotide(), nucleic);
							Assert.assertEquals(pdbfilePdb.getSequence().isProtein(), protein);
							Assert.assertEquals(pdbfilePdb.getSequence().isNucleotide(), nucleic);
						}

						// info from atom serials
						for (int atomser:pdbfilePdb.getAllAtomSerials()) {
							Assert.assertNotNull(pdbfilePdb.getAtomCoord(atomser));
							Assert.assertNotNull(pdbfilePdb.getResSerFromAtomSer(atomser));
						}

						// sec structure
						if (!pdbfilePdb.isNonPolyChain() && pdbfilePdb.getSequence().isProtein()) {
							SecondaryStructure ss = pdbfilePdb.getSecondaryStructure();
							Assert.assertNotNull(ss);
							Assert.assertEquals(pdbfilePdb.getSequence().getSeq(), ss.getSequence());
							Assert.assertEquals(pdbasePdb.getSecondaryStructure().getNumElements(), ss.getNumElements());						

							for (int resser:pdbfilePdb.getAllResSerials()) {
								String resserMsg = "Failed for "+pdbCode+pdbChainCode+" and resser "+resser;
								// checking that the 2 ways of accessing the sec struct element coincide
								Assert.assertSame(ss.getSecStrucElement(resser), pdbfilePdb.getResidue(resser).getSsElem());
								Assert.assertEquals(resserMsg,ss.getSecStrucElement(resser),pdbfilePdb.getResidue(resser).getSsElem());
							}
							
							// comparing contents vs pdbase
							for (Residue pdbfileRes:pdbfilePdb){
								Residue pdbaseRes = pdbasePdb.getResidue(pdbfileRes.getSerial());
								if (pdbfileRes.getSsElem()!=null) {
									SecStrucElement pdbaseSsElem = pdbaseRes.getSsElem();
									SecStrucElement pdbfileSsElem = pdbfileRes.getSsElem();
									Assert.assertEquals(pdbaseSsElem.isHelix(), pdbfileSsElem.isHelix());
									Assert.assertEquals(pdbaseSsElem.isStrand(), pdbfileSsElem.isStrand());
									Assert.assertEquals(pdbaseSsElem.isOther(), pdbfileSsElem.isOther());
									Assert.assertEquals(pdbaseSsElem.getId(),pdbfileSsElem.getId());
									Assert.assertEquals(pdbaseSsElem.getInterval(),pdbfileSsElem.getInterval());
								}
							}
						}
						if (pdbfilePdb.isNonPolyChain()) {
							Assert.assertNull(pdbfilePdb.getSecondaryStructure());
							Assert.assertNull(pdbfilePdb.getSequence());
						}
					}
				}		
				
			} catch (FileNotFoundException e){
				System.err.println("File missing. "+e.getMessage());
			//}catch (PdbLoadException e) {
			//	System.err.println("pdb load error, cause: "+e.getMessage());
			} catch (PdbCodeNotFoundException e) {
				System.err.println("pdb code not found in pdbase");			
			}
		}
		flist.close();
		if (!warnings.isEmpty()) {
			System.err.println("WARNINGS: ");
			for (String warning:warnings) {
				System.err.println(warning);
			}
		}
	}
	
	@Test
	public void testPdbFileWriteRead() throws PdbLoadException, IOException, FileFormatException {
		// testing a list of PDB files from PDB
		BufferedReader flist = new BufferedReader(new FileReader(LISTFILE));
		String line;
		int count = 0;
		while ((line = flist.readLine() ) != null ) {
			if (line.startsWith("#")) continue;
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			
			System.out.println(pdbCode);
			try {
				File pdbFile = unzipFile(new File(PDBDIR,"pdb"+pdbCode+".ent.gz"));
				PdbfileParser parser = new PdbfileParser(pdbFile.getAbsolutePath());
				Integer[] models = parser.getModels();

				Assert.assertTrue(models.length>0);
				
				File writtenPdbFile = new File(TMPDIR,"tmp"+pdbCode+".pdb");
				
				PdbAsymUnit pdb = new PdbAsymUnit(pdbFile,1);
				pdb.writeToPdbFile(writtenPdbFile);
				
				PdbAsymUnit readPdb = new PdbAsymUnit(pdbFile,1);
				
				Assert.assertEquals(pdb.getNumAtoms(),readPdb.getNumAtoms());
				Assert.assertEquals(pdb.getNumChains(),readPdb.getNumChains());
				Assert.assertEquals(pdb.getNumPolyChains(),readPdb.getNumPolyChains());
				Assert.assertEquals(pdb.getNumNonHetAtoms(),readPdb.getNumNonHetAtoms());
				
				for (PdbChain chain:pdb.getAllChains()) {
					PdbChain readChain = readPdb.getChainForChainCode(chain.getChainCode());
					Assert.assertEquals(chain.getNumAtoms(),readChain.getNumAtoms());
					Assert.assertEquals(chain.getNumHeavyAtoms(),readChain.getNumHeavyAtoms());
					Assert.assertEquals(chain.getNumNonHetAtoms(),readChain.getNumNonHetAtoms());
					Assert.assertEquals(chain.getNumStdAaHeavyAtoms(),readChain.getNumStdAaHeavyAtoms());
					if (!chain.isNonPolyChain()) {
						Assert.assertEquals(chain.getFullLength(),readChain.getFullLength());
						Assert.assertEquals(chain.getObsLength(), readChain.getObsLength());
						Assert.assertEquals(chain.getSequence().getSeq(),readChain.getSequence().getSeq());
						Assert.assertEquals(chain.getObsSequence(), readChain.getObsSequence());
						
					}
				}
				
			} catch (FileNotFoundException e){
				System.err.println("File missing. "+e.getMessage());
			}
			if (count>100) break;
			count++;
		}
	}
	
	private static File unzipFile(File repoGzFile) throws FileNotFoundException {
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

	// to debug the testing code (run as java program so that we can use normal debugger)
	public static void main(String[] args) throws Exception {
		PdbParsersTest pdbTest = new PdbParsersTest();
		pdbTest.setUp();
		//pdbTest.testSpecialPDBFiles();
		//pdbTest.testCIFagainstPDBASE();
		pdbTest.testPdbfileParser();
	}

}

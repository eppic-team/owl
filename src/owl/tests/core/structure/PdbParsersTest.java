package owl.tests.core.structure;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//import junit.framework.Assert;
import org.junit.Assert;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import owl.core.structure.CrystalCell;
import owl.core.structure.HetResidue;
import owl.core.structure.NucResidue;
import owl.core.structure.PdbChain;
import owl.core.structure.PdbAsymUnit;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbfileParser;
import owl.core.structure.Residue;
import owl.core.structure.AaResidue;
import owl.core.structure.features.SecStrucElement;
import owl.core.structure.features.SecondaryStructure;
import owl.core.util.FileFormatException;
import owl.core.util.RegexFileFilter;
import owl.tests.TestsSetup;


public class PdbParsersTest {
	
	// paths to the mirror of the PDB ftp repository with all pdb/mmCIF compressed files
	private static final String CIFDIR="/nfs/data/dbs/pdb/data/structures/all/mmCIF/";
	private static final String PDBDIR="/nfs/data/dbs/pdb/data/structures/all/pdb/";
	
	private static final String DATADIR = "src/owl/tests/core/structure/data";
	
	private static final String LISTFILE = DATADIR+"/cullpdb_20";
	
	private static final String TINKERPDBFILE = DATADIR+"/1bxyA_tinker.pdb";
	
	private static final String PDBFILE_NO_TER = DATADIR+"/1c52A.pdb";
	
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
		
		// file with no TER records and no SEQRES
		System.out.println(PDBFILE_NO_TER);
		fullpdb = new PdbAsymUnit(new File(PDBFILE_NO_TER));
		Assert.assertEquals(2,fullpdb.getNumChains());
		Assert.assertEquals(1,fullpdb.getNumPolyChains());
		Assert.assertEquals(1,fullpdb.getNumNonPolyChains());
		Assert.assertEquals(131, fullpdb.getChain("A").getFullLength());
		Assert.assertEquals(fullpdb.getChain("A").getFullLength(), fullpdb.getChain("A").getSequence().getLength());
				
		
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
	public void testPdbfileParser() throws IOException, FileFormatException, PdbLoadException {
		
		long start = System.currentTimeMillis();
		ArrayList<String> warnings = new ArrayList<String>();
		// testing a list of PDB files from PDB
		BufferedReader flist = new BufferedReader(new FileReader(LISTFILE));
		String line;
		while ((line = flist.readLine() ) != null ) {
			if (line.startsWith("#")) continue;
			if (line.isEmpty()) continue;
			String pdbCode = line.split("\\s+")[0].toLowerCase();
			//String pdbChainCode = line.split("\\s+")[1];
			
			System.out.println(pdbCode);
			
			
			try {
				File pdbFile = TestsSetup.unzipFile(new File(PDBDIR,"pdb"+pdbCode+".ent.gz"));
				File cifFile = TestsSetup.unzipFile(new File(CIFDIR,pdbCode+".cif.gz"));
				PdbfileParser parser = new PdbfileParser(pdbFile.getAbsolutePath());
				Integer[] models = parser.getModels();

				Assert.assertTrue(models.length>0);
				
				// loading all chains and all models
				for (int model:models) {
					System.out.println(model);
					PdbAsymUnit pdbfileFullpdb = new PdbAsymUnit(pdbFile,model);
					PdbAsymUnit ciffileFullpdb = new PdbAsymUnit(cifFile,model);
					
					Assert.assertEquals(ciffileFullpdb.getNumPolyChains(),pdbfileFullpdb.getNumPolyChains());
					// we can't assert that number of non-poly chains are the same as the assignment in CIF differs from mine in many cases
					// the only thing we can say for sure is that CIF always assigns more non-poly chains than we do (because they split the non-polys into more groups)
					Assert.assertTrue(ciffileFullpdb.getNumChains()>=pdbfileFullpdb.getNumChains());
				 
					// pdbCode properly read
					Pattern p = Pattern.compile("^\\d\\w\\w\\w$");
					Matcher m = p.matcher(pdbfileFullpdb.getPdbCode());
					Assert.assertTrue(m.matches());
					Assert.assertEquals(ciffileFullpdb.getPdbCode(), pdbfileFullpdb.getPdbCode());

					// title 
					// the capitalization is not consistent between pdb file and cif, 
					// the spaces introduced between different lines are lost in pdb and thus also 
					// not consistent with cif
					// that's why we compare in lower case and stripping spaces
					Assert.assertEquals(ciffileFullpdb.getTitle().toLowerCase().replaceAll(" ", ""), 
							pdbfileFullpdb.getTitle().toLowerCase().replaceAll(" ",""));

					// model as input
					Assert.assertEquals(model,pdbfileFullpdb.getModel());
					
					// crystal data
					if (ciffileFullpdb.getSpaceGroup()!=null) {
						Assert.assertEquals(ciffileFullpdb.getSpaceGroup().getId(),pdbfileFullpdb.getSpaceGroup().getId());
					}
					if (pdbfileFullpdb.isCrystallographicExpMethod()) {
						CrystalCell pdbfileCrystalCell = pdbfileFullpdb.getCrystalCell();
						CrystalCell ciffileCrystalCell = ciffileFullpdb.getCrystalCell();
						
						if (ciffileCrystalCell!=null) {
							Assert.assertEquals(ciffileCrystalCell.getA(), pdbfileCrystalCell.getA(), 0.01);
							Assert.assertEquals(ciffileCrystalCell.getB(), pdbfileCrystalCell.getB(), 0.01);
							Assert.assertEquals(ciffileCrystalCell.getC(), pdbfileCrystalCell.getC(), 0.01);

							Assert.assertEquals(ciffileCrystalCell.getAlpha(), pdbfileCrystalCell.getAlpha(), 0.01);
							Assert.assertEquals(ciffileCrystalCell.getBeta(), pdbfileCrystalCell.getBeta(), 0.01);
							Assert.assertEquals(ciffileCrystalCell.getGamma(), pdbfileCrystalCell.getGamma(), 0.01);
						}
						
						// we can't assert anymore not null because we set to null in case of unreasonable crystal cell
						//Assert.assertNotNull(pdbfileFullpdb.getCrystalCell());
					}
					
					//Release Date
					Assert.assertTrue(ciffileFullpdb.getReleaseDate().equals(pdbfileFullpdb.getReleaseDate()));
					
					// exp data and quality parameters
					Assert.assertEquals(ciffileFullpdb.getExpMethod(),pdbfileFullpdb.getExpMethod());
					if (ciffileFullpdb.isXrayDiffraction()) {
						Assert.assertEquals(ciffileFullpdb.getResolution(),pdbfileFullpdb.getResolution(),0.01);
					}
					if (Math.abs(ciffileFullpdb.getRfree()-pdbfileFullpdb.getRfree())>0.01) {
						// we can't assert because there's no consistency between pdbase and pdbfile, e.g. 1p1x
						System.err.println(pdbCode+": rfrees don't agree, cif "+
								String.format("%4.2f", ciffileFullpdb.getRfree())+", pdbfile "+String.format("%4.2f",pdbfileFullpdb.getRfree()));
						warnings.add(pdbCode+": rfrees don't agree, cif "+
								String.format("%4.2f", ciffileFullpdb.getRfree())+", pdbfile "+String.format("%4.2f",pdbfileFullpdb.getRfree()));
					}
					if (ciffileFullpdb.isXrayDiffraction()) {
						Assert.assertEquals(ciffileFullpdb.getRsym(),pdbfileFullpdb.getRsym(),0.001);
					}

					Assert.assertEquals(ciffileFullpdb.getNumAtoms(),pdbfileFullpdb.getNumAtoms());
					
					//tests for biounit assemblies
					Assert.assertEquals(ciffileFullpdb.getPdbBioUnitListSize(), pdbfileFullpdb.getPdbBioUnitListSize());
					Collections.sort(ciffileFullpdb.getPdbBioUnitList().getPdbBioUnits());
					Collections.sort(pdbfileFullpdb.getPdbBioUnitList().getPdbBioUnits());
					if(ciffileFullpdb.getPdbBioUnitListSize() == pdbfileFullpdb.getPdbBioUnitListSize()){
						for(int i=0; i<ciffileFullpdb.getPdbBioUnitListSize(); i++){
							Assert.assertTrue(ciffileFullpdb.getPdbBioUnitList().getPdbBioUnits().get(i).equals(pdbfileFullpdb.getPdbBioUnitList().getPdbBioUnits().get(i)));
						}
					}
					
					for (PdbChain pdbfilePdbChain:pdbfileFullpdb.getAllChains()) {
						String chainCode = pdbfilePdbChain.getChainCode();
						String pdbChainCode = pdbfilePdbChain.getPdbChainCode();
						System.out.print(chainCode+" "+pdbChainCode+" ");
						PdbChain ciffilePdbChain = ciffileFullpdb.getChainForChainCode(chainCode);
						if (!ciffileFullpdb.containsChainCode(chainCode)) {
							System.err.println("No cif file chain for "+chainCode);
							continue;
						}
						System.out.print(ciffilePdbChain.getChainCode()+" "+ciffilePdbChain.getPdbChainCode());
						String agreement = " ";
						if (chainCode.equals(ciffilePdbChain.getChainCode())) agreement+=" ";
						else agreement+="x";
						if (pdbChainCode.equals(ciffilePdbChain.getPdbChainCode())) agreement+=" ";
						else agreement+="x";
						System.out.println(agreement);
						
						Assert.assertEquals(ciffilePdbChain.getPdbCode(), pdbfilePdbChain.getPdbCode());
						
						// non-poly coincide
						Assert.assertEquals(ciffilePdbChain.isNonPolyChain(),pdbfilePdbChain.isNonPolyChain());
						
						// chain codes coincide
						// we can only assert it for poly-chains, non-poly assignments won't coincide between us and CIF
						if (!ciffilePdbChain.isNonPolyChain()) {
							Assert.assertEquals(ciffilePdbChain.getPdbChainCode(),pdbfilePdbChain.getPdbChainCode());
							Assert.assertEquals(ciffilePdbChain.getChainCode(),pdbfilePdbChain.getChainCode());
						}
						
						Assert.assertEquals(ciffilePdbChain.hasAltCodes(),pdbfilePdbChain.hasAltCodes());
						
			
						if (!pdbfilePdbChain.isNonPolyChain()) {
							Assert.assertEquals(ciffilePdbChain.getFullLength(), pdbfilePdbChain.getFullLength());
							Assert.assertEquals(ciffilePdbChain.getObsLength(), pdbfilePdbChain.getObsLength());
							Assert.assertEquals(ciffilePdbChain.getSequence().getSeq(), pdbfilePdbChain.getSequence().getSeq());
							Assert.assertEquals(ciffilePdbChain.getSequence().isProtein(),pdbfilePdbChain.getSequence().isProtein());
							Assert.assertEquals(ciffilePdbChain.getObsSequence(), pdbfilePdbChain.getObsSequence());

							// sequence
							Assert.assertTrue(pdbfilePdbChain.getObsSequence().length()<=pdbfilePdbChain.getSequence().getLength());
							Assert.assertTrue(pdbfilePdbChain.getObsLength()==pdbfilePdbChain.getObsSequence().length());
							Assert.assertTrue(pdbfilePdbChain.getFullLength()==pdbfilePdbChain.getSequence().getLength());
						}

						
						boolean protein = false;
						boolean nucleic = false;
						for (int resser:pdbfilePdbChain.getAllResSerials()) {
							if (!pdbfilePdbChain.isNonPolyChain()) {
								if (!ciffilePdbChain.getPdbResSerFromResSer(resser).equals(pdbfilePdbChain.getPdbResSerFromResSer(resser))) {
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
											+ciffilePdbChain.getPdbResSerFromResSer(resser)+" in cif and to pdbresser "
											+pdbfilePdbChain.getPdbResSerFromResSer(resser)+" in pdbfile");
									warnings.add(pdbCode+pdbChainCode+": resser " +resser+" maps to pdbresser "
											+ciffilePdbChain.getPdbResSerFromResSer(resser)+" in cif and to pdbresser "
											+pdbfilePdbChain.getPdbResSerFromResSer(resser)+" in pdbfile");
								}
							}

							
							if (ciffilePdbChain.getPdbChainCode().equals(pdbfilePdbChain.getPdbChainCode()) // to be sure we don't compare non-polys for which our chain mapping doesn't coincide 
									&& ciffilePdbChain.containsResidue(resser)) { // we have to check for case where mapping of pdbfile and pdbase from resser to pdbresser don't coincide (we print warnings above for that)
								Residue pdbaseRes = ciffilePdbChain.getResidue(resser);
								Residue pdbfileRes = pdbfilePdbChain.getResidue(resser);
								Assert.assertEquals(resser, pdbaseRes.getSerial());
								Assert.assertEquals(resser, pdbfileRes.getSerial());
								
								// at least 1 atom per observed residue
								Assert.assertTrue(pdbfileRes.getNumAtoms()>0);
								// at least the CA atom must be present, doesn't work for all, e.g. 2cioB residue 80
								//Assert.assertNotNull(pdbfileRes.getAtom("CA"));
								Assert.assertEquals(pdbaseRes.getNumAtoms(), pdbfileRes.getNumAtoms());
								Assert.assertEquals(pdbaseRes.getNumHeavyAtoms(), pdbfileRes.getNumHeavyAtoms());


								Assert.assertEquals(pdbaseRes.isPeptideLinked(), pdbfileRes.isPeptideLinked());
								
								if (!ciffilePdbChain.isNonPolyChain()) {
									Assert.assertEquals(pdbfileRes.getShortCode(), pdbfilePdbChain.getSequence().getSeq().charAt(resser-1));
								}
								

								if (pdbfileRes instanceof AaResidue) {
									protein = true;
									Assert.assertEquals(((AaResidue)pdbaseRes).getAaType(),((AaResidue)pdbfileRes).getAaType());	
									Assert.assertTrue(pdbfileRes.isPeptideLinked());
									Assert.assertTrue(pdbfileRes.getSerial()>0);
									if (!pdbfilePdbChain.isNonPolyChain()) Assert.assertTrue(pdbfileRes.getSerial()<=pdbfilePdbChain.getFullLength());
									if (!pdbfilePdbChain.isNonPolyChain()) Assert.assertTrue(pdbfilePdbChain.getSequence().isProtein());
								} else if (pdbfileRes instanceof HetResidue) {
									Assert.assertEquals(((HetResidue)pdbaseRes).getLongCode(),((HetResidue)pdbfileRes).getLongCode());
									Assert.assertNotSame(((HetResidue)pdbfileRes).getLongCode(),"HOH");
									if (!pdbfilePdbChain.isNonPolyChain() && pdbfilePdbChain.getSequence().isProtein()) {
										Assert.assertTrue(pdbfileRes.isPeptideLinked());
									}
								} else if (pdbfileRes instanceof NucResidue) {
									nucleic = true;
									Assert.assertEquals(((NucResidue)pdbaseRes).getNucType(),((NucResidue)pdbfileRes).getNucType());
									Assert.assertFalse(pdbfileRes.isPeptideLinked());
									Assert.assertTrue(pdbfileRes.getSerial()>0);
									if (!pdbfilePdbChain.isNonPolyChain()) {
										Assert.assertTrue(pdbfileRes.getSerial()<=pdbfilePdbChain.getFullLength());
										Assert.assertTrue(pdbfilePdbChain.getSequence().isNucleotide());
									}
									
								}								
								
							}

						}
						if (!pdbfilePdbChain.isNonPolyChain()) {
							Assert.assertEquals(protein, !nucleic); // one must be true and the other false (a chain can't mix both prot and nucleic)
							Assert.assertEquals(ciffilePdbChain.getSequence().isProtein(), protein);
							Assert.assertEquals(ciffilePdbChain.getSequence().isNucleotide(), nucleic);
							Assert.assertEquals(pdbfilePdbChain.getSequence().isProtein(), protein);
							Assert.assertEquals(pdbfilePdbChain.getSequence().isNucleotide(), nucleic);
						}

						// info from atom serials
						for (int atomser:pdbfilePdbChain.getAllAtomSerials()) {
							Assert.assertNotNull(pdbfilePdbChain.getAtomCoord(atomser));
							Assert.assertNotNull(pdbfilePdbChain.getResSerFromAtomSer(atomser));
						}

						// sec structure
						if (!pdbfilePdbChain.isNonPolyChain() && pdbfilePdbChain.getSequence().isProtein()) {
							SecondaryStructure ss = pdbfilePdbChain.getSecondaryStructure();
							Assert.assertNotNull(ss);
							Assert.assertEquals(pdbfilePdbChain.getSequence().getSeq(), ss.getSequence());
							Assert.assertEquals(ciffilePdbChain.getSecondaryStructure().getNumElements(), ss.getNumElements());						

							for (int resser:pdbfilePdbChain.getAllResSerials()) {
								String resserMsg = "Failed for "+pdbCode+pdbChainCode+" and resser "+resser;
								// checking that the 2 ways of accessing the sec struct element coincide
								Assert.assertSame(ss.getSecStrucElement(resser), pdbfilePdbChain.getResidue(resser).getSsElem());
								Assert.assertEquals(resserMsg,ss.getSecStrucElement(resser),pdbfilePdbChain.getResidue(resser).getSsElem());
							}
							
							// comparing contents vs pdbase
							for (Residue pdbfileRes:pdbfilePdbChain){
								Residue pdbaseRes = ciffilePdbChain.getResidue(pdbfileRes.getSerial());
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
						if (pdbfilePdbChain.isNonPolyChain()) {
							Assert.assertNull(pdbfilePdbChain.getSecondaryStructure());
							Assert.assertNull(pdbfilePdbChain.getSequence());
						}
					}
				}		
				
			} catch (FileNotFoundException e){
				System.err.println("File missing. "+e.getMessage());
			//}catch (PdbLoadException e) {
			//	System.err.println("pdb load error, cause: "+e.getMessage());
			} 
		}
		flist.close();
		long end = System.currentTimeMillis();
		
		if (!warnings.isEmpty()) {
			System.err.println("WARNINGS: ");
			for (String warning:warnings) {
				System.err.println(warning);
			}
		}
		System.out.println("Test ran in "+((end-start)/60000)+"min");
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
				File pdbFile = TestsSetup.unzipFile(new File(PDBDIR,"pdb"+pdbCode+".ent.gz"));
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
		flist.close();
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

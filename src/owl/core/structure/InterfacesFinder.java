package owl.core.structure;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import owl.core.structure.graphs.AICGraph;

public class InterfacesFinder {
	
	private static final Matrix4d IDENTITY = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);

	private class TransformIdTranslation {
		public int transformId;
		public Vector3d translation;
		public boolean partnerSeen;
		public boolean checkPartner; // whether to check the partner or not (default is of course not to check)
		public Matrix4d matTransform;
		public TransformIdTranslation(int transformId, Vector3d translation) {
			this.transformId = transformId;
			this.translation = translation;
			this.partnerSeen = false;
			this.checkPartner = false;
			this.matTransform = (Matrix4d)sg.getTransformation(transformId).clone();
			this.matTransform.setTranslation(translation);
		}
		public String toString() {
			return String.format("[%d-(%5.2f,%5.2f,%5.2f)]",transformId,translation.x,translation.y,translation.z);
		}
	}
	
	private PdbAsymUnit pdb;
	private boolean debug;
	
	private boolean withRedundancyElimination;
	
	private SpaceGroup sg;
	
	private ArrayList<TransformIdTranslation> visited;
	
	public InterfacesFinder(PdbAsymUnit pdb) {
		this.pdb = pdb;
		this.sg = pdb.getSpaceGroup();
		this.debug = false;
		this.withRedundancyElimination = true;
		if (this.pdb.hasHydrogens()) {
			// We have to warn because at the moment we implemented things so that we have to call removeHatoms() before calling getAllInterfaces()
			// We need to fix that so that we simply can calculate interfaces by ignoring Hydrogens without having to remove them
			// While we don't fix that, we need to warn because it's really dangerous to calculate interfaces (and thus core/rim residues) with Hydrogens 
			System.err.println("Warning! Hydrogen atoms present in the structure, the interface calculation (ASAs) won't be reliable!");
		}
	}
	
	public void setDebug(boolean debug) {
		this.debug = debug;
	}
	
	public void setWithRedundancyElimination(boolean withRedundancyElimination) {
		this.withRedundancyElimination = withRedundancyElimination;
	}
	
	private void initialiseVisited() {
		visited = new ArrayList<TransformIdTranslation>();
	}
	
	/**
	 * Returns a sorted (decreasing area) list of all interfaces that the given PdbAsymUnit has upon 
	 * generation of all crystal symmetry objects. An interface is defined as any pair of chains 
	 * that contact, i.e. for which there is at least a pair of atoms (one from each chain) within 
	 * the given cutoff distance.
	 * The interface areas and BSAs are calculated with either our implementation of the rolling
	 * ball algorithm (naccessExe set to null) or the external NACCESS program (naccessExe must 
	 * be passed)
	 * @param cutoff the distance cutoff for 2 chains to be considered in contact
	 * @param naccessExe the NACCESS executable if null our rolling ball algorithm implementation
	 * will be used
	 * @param nSpherePoints
	 * @param nThreads
	 * @param hetAtoms whether to consider HETATOMs in surface area calculations or not
	 * @param nonPoly if true interfaces will be calculated for non-polymer chains and 
	 * protein/nucleic acid polymers, if false only interfaces between protein polymer chains calculated
	 * @param debug set to true to produce some debugging output (run times of each part of the calculation)
	 * @return
	 * @throws IOException when problems when running NACCESS (if NACCESS used)
	 */
	public ChainInterfaceList getAllInterfaces(double cutoff, File naccessExe, int nSpherePoints, int nThreads, boolean hetAtoms, boolean nonPoly) throws IOException {	
		// TODO take care of cases where interfaces are found to a 2nd neighbour cell
		// e.g. 3hz3, 1wqj, 2de3, 1jcd: one needs to go to the 2nd neighbour
		// This probably happens more for longer cutoffs or for very small angles and small molecules
		
		// the set takes care of eliminating duplicates, comparison is based on the equals() 
		// and hashCode() of ChainInterface and that in turn on that of AICGraph and Atom
		Set<ChainInterface> set = new HashSet<ChainInterface>();
		
		// we've got to check if nonPoly=false (i.e. we want only prot-prot interfaces) that there are actually some protein chains!
		if (pdb.getProtChains().size()==0) {
			return calcAsas(set, naccessExe, nSpherePoints, nThreads, hetAtoms);
		}		

		// initialising the visited ArrayList for keeping track of symmetry redundancy
		initialiseVisited();

		// 0. generate complete unit cell
		PdbUnitCell cell = null;
		if (pdb.getCrystalCell()!=null) {
			cell = pdb.getUnitCell();
		}
		
		long start = -1; 
		long end = -1;
		int trialCount = 0, countSkipped1 = 0, countSkipped2 =0, duplicatesCount1=0, duplicatesCount2=0, duplicatesCount3=0;
		if (debug) {
			trialCount = 0;
			start= System.currentTimeMillis();
			System.out.println("Interfaces within asymmetric unit");
		}
		// 1. interfaces within unit cell
		// 1.1 within asymmetric unit
		Set<String> chainCodes = null;
		if (nonPoly) chainCodes = pdb.getChainCodes();
		else chainCodes = pdb.getProtChainCodes();
		
		for (String iChainCode:chainCodes) { 
			for (String jChainCode:chainCodes) { // getPolyChainCodes
				if (iChainCode.compareTo(jChainCode)<=0) continue;
				if (debug) {
					System.out.print(".");
					trialCount++;
				}
				PdbChain chaini = pdb.getChainForChainCode(iChainCode);
				PdbChain chainj = pdb.getChainForChainCode(jChainCode);
				AICGraph graph = chaini.getAICGraph(chainj, cutoff);
				if (graph.getEdgeCount()>0) {
					if (debug) System.out.print("x");
					// because of the bsas are values of the residues of each chain we need to make a copy so that each interface has independent residues
					PdbChain chainiCopy = chaini.copy(pdb);
					PdbChain chainjCopy = chainj.copy(pdb);
					ChainInterface interf = new ChainInterface(chainiCopy,chainjCopy,graph,PdbAsymUnit.IDENTITY_TRANSFORM,PdbAsymUnit.IDENTITY_TRANSFORM);
					interf.setSecondTranslation(new Point3i(0,0,0));
					if (!set.add(interf)) {
						duplicatesCount1++;
					}
				}											
			}
		}
		if (debug) {
			end = System.currentTimeMillis();
			System.out.println("\n"+trialCount+" trials done. Time "+(end-start)/1000+"s");
		}

		
		if (debug) {
			trialCount = 0;
			start= System.currentTimeMillis();
			System.out.println("Interfaces within the rest of the unit cell");
		}
		
		// this condition covers 3 cases:
		// a) entries with expMethod X-RAY/other diffraction and defined crystalCell (most usual case)
		// b) entries with expMethod null but defined crystalCell (e.g. PDB file with CRYST1 record but no expMethod annotation) 
		// c) entries with expMethod not X-RAY (e.g. NMR) and defined crystalCell (NMR entries do have a dummy CRYST1 record "1 1 1 90 90 90 P1")
		if (cell!=null && pdb.isCrystallographicExpMethod()) {
			
			// 1.2 between the original asymmetric unit and the others resulting from applying the symmetry transformations
			for (int j=0;j<cell.getNumAsymUnits();j++) {
				PdbAsymUnit jAsym = cell.getAsymUnit(j);
				if (jAsym==pdb) continue; // we want to compare this to all others but not to itself
				TransformIdTranslation tt = null;
				if (withRedundancyElimination) { 
					tt = new TransformIdTranslation(jAsym.getTransformId(), jAsym.getTranslation());
					if (!addVisited(tt)) {
						if (debug) countSkipped1++;
						continue;
					}
				}
				
//				if (jAsym.getTransformId()==3) {
//					System.err.println("Writing debug file");
//					jAsym.writeToPdbFile(new File("/home/duarte_j/"+pdb.getPdbCode()+"_"+jAsym.getTransformId()+"_000.pdb"));
//				}
				
				int contactsFound = 0;
				Collection<PdbChain> ichains = null;
				Collection<PdbChain> jchains = null;
				if (nonPoly) {
					ichains = pdb.getAllChains();
					jchains = jAsym.getAllChains();
				} else {
					ichains = pdb.getProtChains();
					jchains = jAsym.getProtChains();
				}
				for (PdbChain chaini:ichains) { 
					for (PdbChain chainj:jchains) { 
						if (debug) {
							System.out.print(".");
							trialCount++;
						}
						AICGraph graph = chaini.getAICGraph(chainj, cutoff);
						if (graph.getEdgeCount()>0) {
							contactsFound++;
							if (debug) System.out.print("x");
							// because of the bsas are values of the residues of each chain we need to make a copy so that each interface has independent residues
							PdbChain chainiCopy = chaini.copy(pdb);
							PdbChain chainjCopy = chainj.copy(jAsym);
							ChainInterface interf = new ChainInterface(chainiCopy,chainjCopy,graph,pdb.getTransform(),jAsym.getTransform());
							interf.setSecondTranslation(new Point3i(0,0,0));
							if (!set.add(interf)) {
								duplicatesCount2++;
							}
						}													
					}
				}
				// no contacts for all j chains were found for this operator => 
				// we mark the tt as unchecked so that we make sure we do try the partner equivalent transformation
				// no-contacts condition 2
				if (withRedundancyElimination && contactsFound<(ichains.size()*jchains.size())) {
					tt.checkPartner = true;
				}

			}
			if (debug) {
				end = System.currentTimeMillis();
				System.out.println("\n"+trialCount+" trials done. Time "+(end-start)/1000+"s");
			}

			if (debug) {
				trialCount = 0;
				start= System.currentTimeMillis();
				int trials = pdb.getNumChains()*cell.getNumAsymUnits()*pdb.getNumChains()*26;
				System.out.println("Interfaces between the original asym unit and the 26 neighbouring whole unit cells ("+trials+")");
			}
			
			// 2. interfaces between original asymmetric unit and 26 neighbouring whole unit cells
			for (int i=-1;i<=1;i++) {
				for (int j=-1;j<=1;j++) {
					for (int k=-1;k<=1;k++) {
						if (i==0 && j==0 && k==0) continue; // that would be the identity translation, we calculate that before
						PdbUnitCell translated = cell.copy();
						Vector3d trans = new Vector3d(i,j,k);
						translated.doCrystalTranslation(trans);
						
						for (PdbAsymUnit jAsym:translated.getAllAsymUnits()) {						
							// short-cut 1: checking for redundancy in symmetry, will skip if redundant
							TransformIdTranslation tt = null;
							if (withRedundancyElimination) {
								tt = new TransformIdTranslation(jAsym.getTransformId(), jAsym.getTranslation());
								if (!addVisited(tt)) { 								
									if (debug) countSkipped1++;								
									continue;
								}
							}
							
//							if (jAsym.getTransformId()==2 && i==1 && j==0 && k==0) {
//								System.err.println("Writing debug file");
//								jAsym.writeToPdbFile(new File("/home/duarte_j/"+pdb.getPdbCode()+"_"+jAsym.getTransformId()+"_"+i+""+j+""+k+".pdb"));
//							}
							
							// short-cut 2: checking for too far away for a contact with bounding boxes
							if (pdb.areNotOverlapping(jAsym,cutoff,!nonPoly)) {
								if (debug) {
									countSkipped2++;
								}

								// no contacts for this transformation =>
								// we mark the tt as unchecked so that we make sure we do try the partner equivalent transformation
								// no-contacts condition 1
								if (withRedundancyElimination) tt.checkPartner = true;
								continue;
							}
							
							int contactsFound = 0;
							Collection<PdbChain> ichains = null;
							Collection<PdbChain> jchains = null;
							if (nonPoly) jchains = jAsym.getAllChains();
							else jchains = jAsym.getProtChains();
							for (PdbChain chainj:jchains) {
								
								if (nonPoly) ichains = pdb.getAllChains();
								else ichains = pdb.getProtChains();
								for (PdbChain chaini:ichains) { // we only have to compare the original asymmetric unit to every full cell around
									if (debug) {
										System.out.print(".");
										trialCount++;
									}
									AICGraph graph = chaini.getAICGraph(chainj, cutoff);
									if (graph.getEdgeCount()>0) {
										contactsFound++;										
										if (debug) System.out.print("x");
										// because of the bsas are values of the residues of each chain we need to make a copy so that each interface has independent residues
										PdbChain chainiCopy = chaini.copy(pdb);
										PdbChain chainjCopy = chainj.copy(jAsym);
										ChainInterface interf = new ChainInterface(chainiCopy,chainjCopy,graph,pdb.getTransform(),jAsym.getTransform());
										interf.setSecondTranslation(new Point3i(i,j,k));
										if (!set.add(interf)){
											duplicatesCount3++;
										}
									} 					
								}
							}
							// no contacts for all j chains were found for this operator => 
							// we mark the tt as unchecked so that we make sure we do try the partner equivalent transformation
							// no-contacts condition 2
							if (withRedundancyElimination && contactsFound<(ichains.size()*jchains.size())) {
								tt.checkPartner = true;
							}
						}
					}
				}
			}
			if (debug) {
				end = System.currentTimeMillis();
				System.out.println("\n"+trialCount+" trials done - "+
						countSkipped1+" symmetry redundant branches skipped - "+
						countSkipped2+" too far branches skipped. Total possible "+(trialCount+(countSkipped1+countSkipped2)*pdb.getNumChains()*pdb.getNumChains())+
						" trials. Time "+(end-start)/1000+"s");
				System.out.println("Duplicates: "+duplicatesCount1+" "+duplicatesCount2+" "+duplicatesCount3);
				System.out.println("Found "+set.size()+" interfaces.");
				System.out.println("Transformations used: "+visited.size());
				int countWithPartner = 0;
				for (TransformIdTranslation v:visited) {
					if (v.partnerSeen) countWithPartner++;
				}
				System.out.println("From them, transformations with equivalent partner: "+countWithPartner);
			}
		}
		
		return calcAsas(set, naccessExe, nSpherePoints, nThreads, hetAtoms);
	}
	
	private ChainInterfaceList calcAsas(Set<ChainInterface> set, File naccessExe, int nSpherePoints, int nThreads, boolean hetAtoms) throws IOException {
		// bsa calculation 
		// NOTE in principle it is more efficient to calculate asas only once per isolated chain
		// BUT! surprisingly the rolling ball algorithm gives slightly different values for same molecule in different 
		// orientations! (can't really understand why!). Both NACCESS and our own implementation behave like that.
		// That's why we calculate always for the 2 separate members of interface and the complex, otherwise 
		// we get (not very big but annoying) discrepancies and also things like negative (small) bsa values
		long start =0 ,end =0;
		
		if (debug) {
			start= System.currentTimeMillis();
			System.out.println("Calculating ASAs with "+nThreads+" threads and "+nSpherePoints+" sphere points");
		}
		for (ChainInterface interf:set) {
			//System.out.print(".");
			if (naccessExe!=null) {
				interf.calcSurfAccessNaccess(naccessExe,hetAtoms);
			} else {
				interf.calcSurfAccess(nSpherePoints, nThreads,hetAtoms);
			}
		}
		if (debug) {
			end = System.currentTimeMillis();
			System.out.println("\nDone. Time "+(end-start)/1000+"s");
		}

		// now that we have the areas we can put them into a list and sort them
		ChainInterfaceList.AsaCalcMethod asaCalcMethod = ChainInterfaceList.AsaCalcMethod.INTERNAL;
		if (naccessExe!=null) {
			asaCalcMethod = ChainInterfaceList.AsaCalcMethod.NACCESS;
		}
		ChainInterfaceList list = new ChainInterfaceList(asaCalcMethod);
		if (asaCalcMethod == ChainInterfaceList.AsaCalcMethod.INTERNAL) {
			list.setAsaCalcAccuracyParam(nSpherePoints);
		}
		for (ChainInterface interf:set) {
			list.addInterface(interf);
		}
		list.sort(); // this sorts the returned list and assigns ids to the ChainInterface members
		return list;
	}
	
	
	/**
	 * Checks whether given transformId/translation is symmetry redundant, if not it is added to the list 
	 * of seen transformIds/translations and true returned. If it is redundant nothing is added and false is returned. 
	 * Two transformations are symmetry redundant if their matrices (4d) multiplication gives the identity, i.e.
	 * if one is the inverse of the other.
	 * @param tt
	 * @return
	 */
	private boolean addVisited(TransformIdTranslation tt) {
		
		
		for (TransformIdTranslation v:visited) {
			
			// there's only 1 possible equivalent partner for each visited matrix 
			// (since the equivalent is its inverse matrix and the inverse matrix is unique)
			// thus once the partner has been seen, we don't need to check it ever again
			if (v.partnerSeen) continue;	
			
			Matrix4d mul = new Matrix4d();
			mul.mul(tt.matTransform,v.matTransform);
			
			if (mul.epsilonEquals(IDENTITY, 0.0001)) {
				
				v.partnerSeen = true;
				
				// in some cases the first of the 2 checked equivalent partners
				// can happen to fall in a non-contacting position
				// Thus above in getAllInterfaces we mark the transformation as unchecked (with checkPartner=true) so that now here
				// we force to check also the equivalent partner
				// We do it in with 2 different conditions above: no-contacts condition 1 and no-contacts condition 2
				// Examples that need only condition 1: 2gdg 3ka0 1vzi
				// Examples that need additionally condition 2: 1g3p 1eaq
				if (v.checkPartner) {
					if (debug) System.out.println("Will check "+tt+" because partner "+v+" was not contacting");
					break;
				} else {
					if (debug) System.out.println("Skipping redundant transformation: "+tt+", equivalent to "+v);
				}
				
				return false;
			}

		}
		
		visited.add(tt);
		return true;
	}
	

}

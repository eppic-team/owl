package owl.core.structure;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import owl.core.structure.graphs.AICGraph;
import owl.core.util.BoundingBox;

public class InterfacesFinder {
	
	private static final Matrix4d IDENTITY = new Matrix4d(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1);

	private class TransformIdTranslation {
		public int transformId;
		public Vector3d translation;
		public boolean partnerSeen;
		public Matrix4d matTransform;
		public TransformIdTranslation(int transformId, Vector3d translation) {
			this.transformId = transformId;
			this.translation = translation;
			this.partnerSeen = false;
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
	 * @param nSpherePoints
	 * @param nThreads
	 * @param hetAtoms whether to consider HETATOMs in surface area calculations or not
	 * @param nonPoly if true interfaces will be calculated for non-polymer chains and 
	 * protein/nucleic acid polymers, if false only interfaces between protein polymer chains calculated
	 * @param debug set to true to produce some debugging output (run times of each part of the calculation)
	 * @return
	 */
	public ChainInterfaceList getAllInterfaces(double cutoff, int nSpherePoints, int nThreads, boolean hetAtoms, boolean nonPoly) {	

		// the set takes care of eliminating duplicates, comparison is based on the equals() 
		// and hashCode() of ChainInterface and that in turn on that of AICGraph and Atom
		Set<ChainInterface> set = new HashSet<ChainInterface>();
		
		// we've got to check if nonPoly=false (i.e. we want only prot-prot interfaces) that there are actually some protein chains!
		if (pdb.getProtChains().size()==0) {
			return calcAsas(set, nSpherePoints, nThreads, hetAtoms);
		}		

		// initialising the visited ArrayList for keeping track of symmetry redundancy
		initialiseVisited();

		// 0. generate complete unit cell
		PdbUnitCell cell = null;
		List<BoundingBox> cellBBs = null;
		if (pdb.getCrystalCell()!=null) {
			cell = pdb.getUnitCell();
			// we calculate all the bounds of each of the asym units, those will then be reused and translated
			cellBBs = cell.getBoundingBoxes(!nonPoly);
		}
		
		long start = -1; 
		long end = -1;
		int trialCount = 0, skippedRedundantOrigCell =0, skippedRedundant = 0, skippedAUsNoOverlap =0, skippedCellsNoOverlap = 0;
		int duplicatesCount1=0, duplicatesCount2=0, duplicatesCount3=0;
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
				} else {
					if (debug) System.out.print(".");
				}
			}
		}
		if (debug) {
			end = System.currentTimeMillis();
			System.out.println("\n"+trialCount+" chain-chain clash trials done. Time "+(end-start)/1000+"s");
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
					if (isRedundant(tt)) {
						if (debug) skippedRedundantOrigCell++;
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
						} else {
							if (debug) System.out.print(".");
						}
					}
				}
				
				// if all possible contacts (at max num i chains*num j chains) were found for this operator we mark it as visited
				// in order to later avoiding visiting its equivalent partner
				if (withRedundancyElimination && contactsFound==(ichains.size()*jchains.size())) {
					addVisited(tt);
				}

			}
			if (debug) {
				end = System.currentTimeMillis();
				System.out.println("\n"+trialCount+" chain-chain clash trials done. Time "+(end-start)/1000+"s");
			}

			
			// 2. interfaces between original asymmetric unit and neighboring whole unit cells
			// We go up to 2nd neighbor, in some cases (e.g. 3hz3, 1wqj, 2de3, 1jcd) there are interfaces with 2nd neighboring cells
			// Would we ever need 3rd neighbors? in principle possible but it would have to be a very strange cell
			int numCells = 2;
			
			if (debug) {
				trialCount = 0;
				start= System.currentTimeMillis();
				int neighbors = (2*numCells+1)*(2*numCells+1)*(2*numCells+1)-1;
				int trials = pdb.getNumChains()*cell.getNumAsymUnits()*pdb.getNumChains()*neighbors;
				System.out.println("Interfaces between the original asym unit and the neighbouring "+neighbors+" whole unit cells " +
						"(2x"+pdb.getNumChains()+"chains x "+cell.getNumAsymUnits()+"AUs x "+neighbors+"cells = "+trials+" total possible trials)");
			}

			
			for (int i=-numCells;i<=numCells;i++) {
				for (int j=-numCells;j<=numCells;j++) {
					for (int k=-numCells;k<=numCells;k++) {
						if (i==0 && j==0 && k==0) continue; // that would be the identity translation, we calculate that before

						Vector3d trans = new Vector3d(i,j,k);
						Matrix4d m = pdb.getCrystalCell().getTransform(trans);
						
						for (int au=0;au<cell.getNumAsymUnits();au++) { 
							BoundingBox bbTrans = new BoundingBox(cellBBs.get(au));							
							bbTrans.translate(new Vector3d(m.m03,m.m13,m.m23));
						
							// We have 3 short-cut strategies
							// 1) we skip first of all if the bounding boxes of the unit cells don't overlap
							if (!pdb.getBoundingBox(!nonPoly).overlaps(bbTrans, cutoff)) {
								if (debug) skippedCellsNoOverlap++;
								continue;
							}
													
							// 2) we check if we didn't already see its equivalent symmetry operator partner
							//    in reality this doesn't happen so often as we only count it as seen when contacts are seen for all chain interactions 
							TransformIdTranslation tt = null;
							if (withRedundancyElimination) {
								tt = new TransformIdTranslation(cell.getAsymUnit(au).getTransformId(), trans);
								// pure cell-translations we can always add as visited 
								if (tt.transformId==0) addVisited(tt);
								if (isRedundant(tt)) { 								
									if (debug) skippedRedundant++;								
									continue;
								}
							}
							
							// now we copy and actually translate the AU if we saw it does overlap and the sym op was not redundant
							PdbAsymUnit jAsym = cell.getAsymUnit(au).copy();
							jAsym.doCrystalTranslation(trans);

							
//							if (jAsym.getTransformId()==2 && i==1 && j==0 && k==0) {
//								System.err.println("Writing debug file");
//								jAsym.writeToPdbFile(new File("/home/duarte_j/"+pdb.getPdbCode()+"_"+jAsym.getTransformId()+"_"+i+""+j+""+k+".pdb"));
//							}
							
							// 3) now we check whether the AUs don't overlap
							if (pdb.areNotOverlapping(jAsym,cutoff,!nonPoly)) {
								if (debug) skippedAUsNoOverlap++;								
								continue;
							}
							
							// Now that we know that boxes overlap and operator is not redundant, we have to go to the details 
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
									} else {
										if (debug) System.out.print(".");
									}
								}
							}
							// if all possible contacts (at max num i chains*num j chains) were found for this operator we mark it as visited
							// in order to later avoiding visiting its equivalent partner
							// [ in some cases after applying a transform, the molec can happen to fall in a non-contacting position,
							//   examples where this happens are: 2gdg 3ka0 1vzi 1g3p 1eaq 
							//   That's why we can't just add every single transform as visited ]
							if (debug) System.out.println(" "+contactsFound+"("+ichains.size()*jchains.size()+")");
							if (withRedundancyElimination && contactsFound==(ichains.size()*jchains.size())) {
								addVisited(tt);
							}
						}
					}
				}
			}
			if (debug) {
				end = System.currentTimeMillis();
				System.out.println("\n"+trialCount+" chain-chain clash trials done. Time "+(end-start)/1000+"s");
				System.out.println("  skipped (not overlapping full cells): "+skippedCellsNoOverlap);
				System.out.println("  skipped (not overlapping AUs)       : "+skippedAUsNoOverlap);
				System.out.println("  skipped (sym redundant within cell) : "+skippedRedundantOrigCell);
				System.out.println("  skipped (sym redundant other cells) : "+skippedRedundant);

				System.out.println("\nDuplicates: "+duplicatesCount1+" "+duplicatesCount2+" "+duplicatesCount3);
				System.out.println("Found "+set.size()+" interfaces.");
			}
		}
		
		return calcAsas(set, nSpherePoints, nThreads, hetAtoms);
	}
	
	private ChainInterfaceList calcAsas(Set<ChainInterface> set, int nSpherePoints, int nThreads, boolean hetAtoms) {
		// bsa calculation 
		// NOTE in principle it is more efficient to calculate asas only once per isolated chain
		// BUT! surprisingly the rolling ball algorithm gives slightly different values for same molecule in different 
		// orientations! (due to sampling depending on orientation of axes grid). Both NACCESS and our own implementation behave like that.
		// That's why we calculate always for the 2 separate members of interface and the complex, otherwise 
		// we get (not very big but annoying) discrepancies and also things like negative (small) bsa values
		// We do take a shortcut: instead of calculating asas for all monomers, 
		// we only do it once per monomer type (chain code) per translation (transformId), since it's only
		// rotations that give different area values
		
		long start =0 ,end =0;
		
		if (debug) {
			start= System.currentTimeMillis();
			System.out.println("Calculating ASAs with "+nThreads+" threads and "+nSpherePoints+" sphere points");
		}

		// first we put one chain of each transform id in the map
		TreeMap<String, PdbChain> transformId2chain = new TreeMap<String, PdbChain>();
		for (ChainInterface interf:set) {
			PdbChain first = interf.getFirstMolecule();
			PdbChain second = interf.getSecondMolecule();
			if (!transformId2chain.containsKey(first.getPdbChainCode()+first.getParent().getTransformId())) {
				transformId2chain.put(first.getPdbChainCode()+first.getParent().getTransformId(),first);
			}
			if (!transformId2chain.containsKey(second.getPdbChainCode()+second.getParent().getTransformId())) {
				transformId2chain.put(second.getPdbChainCode()+second.getParent().getTransformId(),second);
			}
		}
		
		// now we calculate ASAs only for those chains in the map
		for (PdbChain chain:transformId2chain.values()) {
			if (debug) System.out.print(".");
			chain.calcASAs(nSpherePoints, nThreads, hetAtoms);
		}
		if (debug) System.out.println();
		
		// for all others we copy the values from their corresponding transformids partners
		for (ChainInterface interf:set) {
			
			PdbChain first = interf.getFirstMolecule();
			PdbChain second = interf.getSecondMolecule();
			if (!first.hasASA()) {
				first.setASAs(transformId2chain.get(first.getPdbChainCode()+first.getParent().getTransformId()));
			}
			if (!second.hasASA()) {
				second.setASAs(transformId2chain.get(second.getPdbChainCode()+second.getParent().getTransformId()));
			}

			// finally we need to calculate bsas
			if (debug) System.out.print(".");
			interf.calcSurfAccess(nSpherePoints, nThreads,hetAtoms);

		}
		
		
		if (debug) {
			end = System.currentTimeMillis();
			System.out.println("\nDone. Time "+(end-start)/1000+"s");
		}

		// now that we have the areas we can put them into a list and sort them
		ChainInterfaceList.AsaCalcMethod asaCalcMethod = ChainInterfaceList.AsaCalcMethod.INTERNAL;
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
	
	
	private void addVisited(TransformIdTranslation tt) {
		visited.add(tt);
	}
	
	/**
	 * Checks whether given transformId/translation is symmetry redundant 
	 * Two transformations are symmetry redundant if their matrices (4d) multiplication gives the identity, i.e.
	 * if one is the inverse of the other.
	 * @param tt
	 * @return
	 */
	private boolean isRedundant(TransformIdTranslation tt) {
		
		
		for (TransformIdTranslation v:visited) {
			
			// there's only 1 possible equivalent partner for each visited matrix 
			// (since the equivalent is its inverse matrix and the inverse matrix is unique)
			// thus once the partner has been seen, we don't need to check it ever again
			if (v.partnerSeen) continue;	
			
			Matrix4d mul = new Matrix4d();
			mul.mul(tt.matTransform,v.matTransform);
			
			if (mul.epsilonEquals(IDENTITY, 0.0001)) {

				if (debug) System.out.println("Skipping redundant transformation: "+tt+", equivalent to "+v);
				v.partnerSeen = true;
				
				return true;
			}
		}
		
		return false;
	}
	

}

package owl.core.structure;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import javax.vecmath.Point3d;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3d;

import owl.core.structure.graphs.AICGraph;

public class InterfacesFinder {

	private class TransformIdTranslation {
		public int transformId;
		public Vector3d translation;
		public TransformIdTranslation(int transformId, Vector3d translation) {
			this.transformId = transformId;
			this.translation = translation;
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
		// TODO also take care that for longer cutoffs or for very small angles and small molecules one might need to go to the 2nd neighbour
		// TODO pathological cases, 3hz3: one needs to go to the 2nd neighbour
		
		// the set takes care of eliminating duplicates, comparison is based on the equals() 
		// and hashCode() of ChainInterface and that in turn on that of AICGraph and Atom
		Set<ChainInterface> set = new HashSet<ChainInterface>();

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
		if (cell!=null && 
				(pdb.getExpMethod()==null || 
				pdb.getExpMethod().equals("X-RAY DIFFRACTION") || 
				pdb.getExpMethod().equals("NEUTRON DIFFRACTION") || 
				pdb.getExpMethod().equals("ELECTRON CRYSTALLOGRAPHY"))) { 
			// 1.2 between the original asymmetric unit and the others resulting from applying the symmetry transformations
			for (int j=0;j<cell.getNumAsymUnits();j++) {
				PdbAsymUnit jAsym = cell.getAsymUnit(j);
				if (jAsym==pdb) continue; // we want to compare this to all others but not to itself

				if (withRedundancyElimination) {
					if (!addVisited(jAsym.getTransformId(), jAsym.getTranslation())) {
						if (debug) countSkipped1++;
						continue;
					}
				}
				
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
							if (withRedundancyElimination) {
								if (!addVisited(jAsym.getTransformId(), jAsym.getTranslation())) {
									if (debug) countSkipped1++;								
									continue;
								}
							}
							// short-cut 2: checking for too far away asu cell, will skip if too far away
							Point3d sep = pdb.getCrystalSeparation(jAsym);
							if (Math.abs(sep.x)>1.1 || Math.abs(sep.y)>1.1 || Math.abs(sep.z)>1.1) {
								if (debug) {
									//System.out.println("\nskipping:");
									//System.out.printf("(%2d,%2d,%2d) - %2d : %5.2f,%5.2f,%5.2f (%2d,%2d,%2d)\n",i,j,k,jAsym.getTransformId(),
									//		sep.x,sep.y,sep.z,
									//		(int)Math.round(sep.x),(int)Math.round(sep.y),(int)Math.round(sep.z));
									countSkipped2++;
								}
								continue;
							}
							Collection<PdbChain> jchains = null;
							if (nonPoly) jchains = jAsym.getAllChains();
							else jchains = jAsym.getProtChains();
							for (PdbChain chainj:jchains) {

								Collection<PdbChain> ichains = null;
								if (nonPoly) ichains = pdb.getAllChains();
								else ichains = pdb.getProtChains();
								for (PdbChain chaini:ichains) { // we only have to compare the original asymmetric unit to every full cell around
									if (debug) {
										System.out.print(".");
										trialCount++;
									}
									AICGraph graph = chaini.getAICGraph(chainj, cutoff);
									if (graph.getEdgeCount()>0) {
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
	 * of seen transformIds/translations and true returned. If it is redundant nothing is added and false is returned 
	 * @param transformId
	 * @param newDirection
	 * @return
	 */
	private boolean addVisited(int transformId, Vector3d newDirection) {
		TransformIdTranslation dt = new TransformIdTranslation(transformId, newDirection);
		// for a not rotational transformation (transformId==0 or identity rotations)
		if (sg.isRotationIdentity(transformId)) {
			// if a seen vector and new vector are same but opposite in sign then the transformation is redundant
			for (TransformIdTranslation v:visited) {
				if (!sg.isRotationIdentity(v.transformId)) continue;
				Vector3d sum = new Vector3d();
				sum.add(v.translation,dt.translation);
				if (sum.epsilonEquals(new Vector3d(0,0,0), 0.001)) {
					// they are symmetry redundant, return false and don't add
					if (debug) {
						System.out.println("Skipping redundant translational transformation: "+dt+", equivalent to "+v);
					}
					return false;
				}
			}
		} 
		// for a rotational transformation, it is a bit more complicated
		else if (sg.getAxisFoldType(transformId)==2)  { // 2-FOLD AXIS
			Vector3d axis = sg.getRotAxes().get(dt.transformId-1);
			//double angle = sg.getRotAngles().get(dt.transformId-1);
			Vector3d tNormalProj = getNormalToAxisProjection(axis, dt.translation);
			
			for (TransformIdTranslation v:visited) {
				if (sg.isRotationIdentity(v.transformId)) continue;
				if (!sg.areRotRelated(v.transformId,dt.transformId)) continue; // rotations are related (so they must also be in same axis)
				
								
				if (
					// rule 1) translations are opposite with respect to the axis of rotation --> HOLDS FOR SURE FOR ALL TYPES OF AXES
					deltaComp(axis.dot(dt.translation),-axis.dot(v.translation)) &&  
					// rule 2) projections on plane normal to axis form an angle == 0 --> WORKS FOR 2-FOLDS
					deltaComp(getNormalToAxisProjection(axis,v.translation).angle(tNormalProj),0) &&   
					// rule 3) projections on plane normal to axis are same length, holds if simply vectors same length (because of rule 1) --> WORKS FOR 2-FOLDS
					deltaComp(v.translation.length(),dt.translation.length())) { 
					if (debug) {
						System.out.println("Skipping redundant "+sg.getAxisFoldType(transformId)+"-fold transformation: "+dt+", equivalent to "+v);
					}
					return false;
				}
			}
		}
		else if (sg.getAxisFoldType(transformId)==3 ||
				 sg.getAxisFoldType(transformId)==4 ||
				 sg.getAxisFoldType(transformId)==6)  { // 3-FOLD,4-FOLD, 6-FOLD AXIS
			Vector3d axis = sg.getRotAxes().get(dt.transformId-1);
			//double angle = sg.getRotAngles().get(dt.transformId-1);
			Vector3d tNormalProj = getNormalToAxisProjection(axis, dt.translation);
			
			for (TransformIdTranslation v:visited) {
				if (sg.isRotationIdentity(v.transformId)) continue;
				if (!sg.areRotRelated(v.transformId,dt.transformId)) continue; // rotations are related (so they must also be in same axis)
				
								
				if (
					// rule 1) translations are opposite with respect to the axis of rotation --> HOLDS FOR SURE FOR ALL TYPES OF AXES
					deltaComp(axis.dot(dt.translation),-axis.dot(v.translation)) &&  
					// rule 2) projections on plane normal to axis are 0 length (we only consider the vertical within the axis) 
					// TODO    we need a more general rule! for other cells in other verticals
					deltaComp(tNormalProj.length(),0) && deltaComp(getNormalToAxisProjection(axis,v.translation).length(),0)) { 
					if (debug) {
						System.out.println("Skipping redundant "+sg.getAxisFoldType(transformId)+"-fold transformation: "+dt+", equivalent to "+v);
					}
					return false;
				}
			}		
		}
		visited.add(dt);
		return true;
	}
	
	private Vector3d getNormalToAxisProjection(Vector3d a, Vector3d t) {
		// see http://en.wikipedia.org/wiki/Vector_projection
		Vector3d axisUnit = new Vector3d(a);
		axisUnit.normalize();
		Vector3d tAxisProjection = new Vector3d(axisUnit);
		tAxisProjection.scale(axisUnit.dot(t));
		Vector3d tNormalProjection = new Vector3d(t);
		tNormalProjection.sub(tAxisProjection);
		return tNormalProjection;
	}
	
	private boolean deltaComp(double d1, double d2) {
		if (Math.abs(d1-d2)<0.000001) return true;
		return false;
	}


}

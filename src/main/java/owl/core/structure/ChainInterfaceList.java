package owl.core.structure;

import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

import javax.vecmath.Point3d;

import owl.core.util.GeometryTools;

//import owl.core.util.CombinationsGenerator;


/**
 * A list of all the interfaces of a crystal structure (a PdbAsymUnit)
 * 
 * 
 * @author duarte_j
 *
 */
public class ChainInterfaceList implements Iterable<ChainInterface>, Serializable {

	private static final long serialVersionUID = 1L;

	public enum AsaCalcMethod {
		INTERNAL("internal"),NACCESS("naccess"),PISA("pisa");

		private String name;
		private AsaCalcMethod(String name) {
			this.name = name;
		}

		public String getName() {
			return name;
		}
	}

	private List<ChainInterface> list;

	private InterfaceGraph graph; 

	/**
	 * Map of interface ids to their corresponding interface clusters
	 */
	private TreeMap<Integer, InterfaceCluster> clusters;

	private AsaCalcMethod asaCalcMethod;
	private int asaCalcAccuracyParam;

	public ChainInterfaceList(AsaCalcMethod asaCalcMethod) {
		this.list = new ArrayList<ChainInterface>();
		this.asaCalcMethod = asaCalcMethod;
	}

	public void addInterface(ChainInterface interf) {
		list.add(interf);
	}

	public int getNumInterfaces() {
		return list.size();
	}

	public int getNumProtProtInterfaces() {
		int count = 0;
		for (ChainInterface interf:this) {
			if (interf.isProtein()) {
				count++;
			}
		}
		return count;
	}

	public int getNumProtProtInterfacesAboveArea(double area) {
		int count = 0;
		for (ChainInterface interf:this) {
			if (interf.isProtein() && interf.getInterfaceArea()>area) {			
				count++;
			}
		}
		return count;
	}

	/**
	 * Removes from this interface list all interfaces with areas
	 * below the given cutoff area
	 * @param area
	 */
	public void removeInterfacesBelowArea(double area) {
		Iterator<ChainInterface> it = iterator();
		while (it.hasNext()) {
			ChainInterface interf = it.next();
			if (interf.getInterfaceArea()<area) {
				it.remove();
			}
		}
	}

	/**
	 * Sorts the interface list descending on total interface areas assigning identifiers 
	 * from 1 to n in that order.
	 */
	public void sort() {
		Collections.sort(list);
		int i=1;
		for (ChainInterface interf:list) {
			interf.setId(i);
			i++;
		}
	}

	public int size() {
		return list.size();
	}

	/**
	 * Gets the interface corresponding to given id.
	 * The ids go from 1 to n
	 * If {@link #sort()} was called then the order is descendent by area.
	 * @param id
	 * @return
	 */
	public ChainInterface get(int id) {
		return list.get(id-1);
	}

	public AsaCalcMethod getAsaCalcType() {
		return asaCalcMethod;
	}

	public void setAsaCalcAccuracyParam(int asaCalcAccuracyParam) {
		this.asaCalcAccuracyParam = asaCalcAccuracyParam;
	}

	public int getAsaCalcAccuracyParam() {
		return asaCalcAccuracyParam;
	}

	public int getNumInterfacesAboveArea(double area) {
		int count = 0;
		for (ChainInterface interf:this) {
			if (interf.getInterfaceArea()>area) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Calculates the rims and cores of all interfaces in list for the
	 * bsaToAsaCutoff given.
	 * @param bsaToAsaCutoff
	 */
	public void calcRimAndCores(double bsaToAsaCutoff, double minAsaForSurface) {
		for (ChainInterface interf:list){
			interf.calcRimAndCore(bsaToAsaCutoff, minAsaForSurface);
		}
	}

	public boolean hasInterfacesWithClashes() {
		for (ChainInterface interf:list){
			if (interf.hasClashes()) return true;
		}
		return false;
	}

	public List<ChainInterface> getInterfacesWithClashes() {
		List<ChainInterface> clashyInterfs = new ArrayList<ChainInterface>();
		for (ChainInterface interf:list){
			if (interf.hasClashes()) {
				clashyInterfs.add(interf);
			}
		}
		return clashyInterfs;
	}

	public Iterator<ChainInterface> iterator() {
		return list.iterator();
	}

	public void printTabular(PrintStream ps, String pdbName, boolean usePdbResSer) {
		ps.println("Interfaces for "+pdbName);
		ps.print("ASAs values from "+this.asaCalcMethod.getName());
		if (asaCalcMethod==AsaCalcMethod.INTERNAL) {
			ps.println(" (sphere sampling points="+this.asaCalcAccuracyParam+")");
		} else {
			ps.println();
		}
		for (ChainInterface interf:this) {
			interf.printTabular(ps, usePdbResSer);
		}
	}

	/**
	 * Given a PDB chain code returns the List of residues that are in the surface but belong to NO
	 * interface (above given minInterfArea) 
	 * Surface residues will be considered those with ASA above the given minAsaForSurface
	 * @param pdbChainCode
	 * @param minInterfArea
	 * @param minAsaForSurface
	 * @return
	 */
	public List<Residue> getResiduesNotInInterfaces(String pdbChainCode, double minInterfArea, double minAsaForSurface) {


		PdbChain chain = null;
		for (ChainInterface interf:this) {
			if (interf.getFirstMolecule().getPdbChainCode().equals(pdbChainCode)) {
				chain = interf.getFirstMolecule();
				break;
			}
			if (interf.getSecondMolecule().getPdbChainCode().equals(pdbChainCode)) {
				chain = interf.getSecondMolecule();
				break;
			}
		}

		List<Residue> surfResidues = chain.getSurfaceResidues(minAsaForSurface);

		//System.out.println(" in surface: "+surfResidues.size());

		Set<Integer> interfResSerials = new HashSet<Integer>();

		for (ChainInterface interf:this) {
			if (interf.getInterfaceArea()>minInterfArea) {
				if (interf.getFirstMolecule().getPdbChainCode().equals(pdbChainCode)) {
					for (Residue res:interf.getFirstRimCore().getCoreResidues()) {
						interfResSerials.add(res.getSerial());
					}
					for (Residue res:interf.getFirstRimCore().getRimResidues()) {
						interfResSerials.add(res.getSerial());
					}
				}
				if (interf.getSecondMolecule().getPdbChainCode().equals(pdbChainCode)) {
					for (Residue res:interf.getSecondRimCore().getCoreResidues()) {
						interfResSerials.add(res.getSerial());
					}
					for (Residue res:interf.getSecondRimCore().getRimResidues()) {
						interfResSerials.add(res.getSerial());
					}
				}

			}
		}

		Iterator<Residue> it = surfResidues.iterator();
		while (it.hasNext()) {
			Residue res = it.next();
			if (interfResSerials.contains(res.getSerial())) {
				it.remove();
			}
		}
		//System.out.println(" in no interface: "+surfResidues.size());

		return surfResidues;
	}

	/**
	 * Returns true if this interface list contains infinite interfaces
	 * An infinite interface is an interface between two crystallographic-symmetry related 
	 * monomers (same PDB chain codes) by an operator of an infinite character: pure translation
	 * or screw rotation. These interfaces lead to infinite fiber-like assemblies. 
	 * @return
	 */
	public boolean hasInfiniteInterfaces() {
		for (ChainInterface interf:this) {
			if (interf.isInfinite()) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Get a list of ids of all interfaces that are not infinite.
	 * An infinite interface is an interface between two crystallographic-symmetry related 
	 * monomers (same PDB chain codes) by an operator of an infinite character: pure translation
	 * or screw rotation. These interfaces lead to infinite fiber-like assemblies. 
	 * @return
	 */
	public ArrayList<Integer> getNonInfiniteInterfacesIds() {
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (ChainInterface interf:this) {
			if (!interf.isInfinite()) {
				list.add(interf.getId());
			}
		}		
		return list;
	}

	public InterfaceGraph getInterfacesGraph() {
		if (graph==null) {
			graph = new InterfaceGraph(this);
		}
		return graph;
	}
	
	public void initialiseClusters(PdbAsymUnit pdb, double rmsdCutoff, int minNumAtomsToCompare, String atomName) {
		clusters = new TreeMap<Integer, InterfaceCluster>();
				
		// first getting all common observed residue serials sets to later use them 
		HashMap<String,List<Integer>> commonObservedSets = new HashMap<String,List<Integer>>();
		for (ChainCluster chainCluster:pdb.getProtChainClusters()) {
			List<Integer> common = chainCluster.getCommonObservedSet(atomName);
			
			// we first check that there are enough atoms to make meaningful comparisons
			// otherwise all members pdb chain codes are not added to the map, so later we can check for nulls
			if (common.size()>minNumAtomsToCompare) {  
				for (PdbChain chain:chainCluster.getMembers()) {
					commonObservedSets.put(chain.getPdbChainCode(), common);
				}
			}
		}
		
		double[][] rmsdMatrix = new double[this.size()][this.size()];
		
		for (ChainInterface iInterf:this) {
			int iId = iInterf.getId();
			
			String iFirstPdbChainCode = iInterf.getFirstMolecule().getPdbChainCode();
			String iSecondPdbChainCode = iInterf.getSecondMolecule().getPdbChainCode();
						
			for (ChainInterface jInterf:this) {
				
				int jId = jInterf.getId();
				
				if (iId>=jId) continue;
				
				String jFirstPdbChainCode = jInterf.getFirstMolecule().getPdbChainCode();
				String jSecondPdbChainCode = jInterf.getSecondMolecule().getPdbChainCode();


				// prerequisite: enough atoms in all 4 chains
				if (!commonObservedSets.containsKey(iFirstPdbChainCode) || !commonObservedSets.containsKey(iSecondPdbChainCode) ||
					!commonObservedSets.containsKey(jFirstPdbChainCode) || !commonObservedSets.containsKey(jSecondPdbChainCode)) {
					
					rmsdMatrix[iId-1][jId-1] = -1.0;
					continue;
				}

				Point3d[] iConformation = iInterf.getConformation( 
						commonObservedSets.get(iFirstPdbChainCode), commonObservedSets.get(iSecondPdbChainCode), atomName, false);

				
				// filling the matrix only when the pair of interfaces are clustering candidates 
				// e.g. if we have clusters A(B) and C(D)
				
				// case to avoid: exact same chains in both, e.g. A+A, A+A
				if (iFirstPdbChainCode.equals(iSecondPdbChainCode) && 
						iFirstPdbChainCode.equals(jFirstPdbChainCode) && 
						iFirstPdbChainCode.equals(jSecondPdbChainCode)){
					
					rmsdMatrix[iId-1][jId-1] = -1.0;
				}
				
				// possible candidates for clustering: NCS related chains in 1st to 1st matching (A+C,B+D) 
				else if (pdb.areChainsInSameCluster(iFirstPdbChainCode, jFirstPdbChainCode) && 
						pdb.areChainsInSameCluster(iSecondPdbChainCode, jSecondPdbChainCode)) {

					if (pdb.areChainsInSameCluster(iFirstPdbChainCode, iSecondPdbChainCode)) { 
						// try direct and reverse
						
						// this case is when the 4 chains are in same cluster: only one common set
						List<Integer> jCommonSet = commonObservedSets.get(jFirstPdbChainCode);
						
						Point3d[] jConformation = jInterf.getConformation(jCommonSet, jCommonSet, atomName, false);

						double rmsdDirect = GeometryTools.calcOptimalSuperposition(iConformation, jConformation, false).getRmsd();

						Point3d[] jConformationReverse = jInterf.getConformation(jCommonSet, jCommonSet, atomName, true);

						double rmsdReverse = GeometryTools.calcOptimalSuperposition(iConformation, jConformationReverse, false).getRmsd();
						
						//if (Math.abs(rmsdReverse-rmsdDirect)>0.5) 
						//	System.err.printf("%d - %d -- rmsd direct: %5.2f, reverse %5.2f\n",iId,jId,rmsdDirect,rmsdReverse);
						
						rmsdMatrix[iId-1][jId-1] = Math.min(rmsdDirect, rmsdReverse);

					} else {
						// try only direct
						Point3d[] jConformation = jInterf.getConformation( 
								commonObservedSets.get(jFirstPdbChainCode), commonObservedSets.get(jSecondPdbChainCode), atomName, false);

						rmsdMatrix[iId-1][jId-1] = GeometryTools.calcOptimalSuperposition(iConformation, jConformation, false).getRmsd();

					}				
					
				} 
				
				// possible candidates for clustering: NCS related chains in 1st to 2nd matching (A+C,D+B)
				else if (pdb.areChainsInSameCluster(iFirstPdbChainCode, jSecondPdbChainCode) && 
						pdb.areChainsInSameCluster(iSecondPdbChainCode, jFirstPdbChainCode)) {												

					Point3d[] jConformation = jInterf.getConformation( 
							commonObservedSets.get(jSecondPdbChainCode), commonObservedSets.get(jFirstPdbChainCode), atomName, true);

					rmsdMatrix[iId-1][jId-1] = 
							GeometryTools.calcOptimalSuperposition(iConformation, jConformation, false).getRmsd();

				} 
				
				// all other cases are not possible cluster candidates
				else {
				
					rmsdMatrix[iId-1][jId-1] = -1.0;
				}
				
			}
			
		}
		
		int clusterId = 1;
		
		for (int iId=1;iId<this.size()+1;iId++) {
			//System.out.printf("%2s: ",iId);
			for (int jId=iId+1;jId<this.size()+1;jId++) {
				//System.out.printf("%5.2f ",rmsdMatrix[iId-1][jId-1]);
				if (rmsdMatrix[iId-1][jId-1] >= 0.0 && rmsdMatrix[iId-1][jId-1] < rmsdCutoff) {
					if (clusters.containsKey(iId)) {
						InterfaceCluster interfCluster = clusters.get(iId);
						interfCluster.addMember(get(jId));
						
						clusters.put(jId, interfCluster);
					} else {
						InterfaceCluster interfCluster = new InterfaceCluster(clusterId);
						interfCluster.addMember(get(iId));
						interfCluster.addMember(get(jId));
						clusters.put(iId, interfCluster);
						clusters.put(jId, interfCluster);
						clusterId++;
					}
				} 
			}
			//System.out.println();
		}
		//System.out.println();
		
		// anything not clustered is assigned to a singleton cluster (cluster with one member)
		for (int id=1;id<this.size()+1;id++) {
			if (!clusters.containsKey(id)) {
				InterfaceCluster interfCluster = new InterfaceCluster(clusterId);
				interfCluster.addMember(get(id));
				clusters.put(id,interfCluster);
				clusterId++;
			}
		}
		
		// reassigning cluster ids based on mean areas (see InterfaceCluster.compare implementation)
		List<InterfaceCluster> uniqueClusters = getClusters();		
		Collections.sort(uniqueClusters);
		for (int i=0;i<uniqueClusters.size();i++) {
			uniqueClusters.get(i).setId(i+1);
			// and we also initialise the representative (not set yet)
			uniqueClusters.get(i).initialiseRepresentative();
		}
		
	}
	
	/**
	 * Return a list of unique clusters corresponding to interfaces in this list.
	 * Note that clusters need to be initialised by calling {@link #initialiseClusters(PdbAsymUnit, double, int, String)}
	 * @return
	 * @throws NullPointerException if clusters are not initialised by calling {@link #initialiseClusters(PdbAsymUnit, double, int, String)}
	 */
	public List<InterfaceCluster> getClusters() {
		List<InterfaceCluster> list = new ArrayList<InterfaceCluster>();
		
		for (InterfaceCluster cluster:clusters.values()) {
			boolean present = false;
			for (InterfaceCluster cl:list) {
				if (cl==cluster) {
					present = true;
					break;
				}
			}
			if (!present) list.add(cluster);
		}
		return list;
	}
	
	/**
	 * Returns the cluster to which the given interface id belongs
	 * @param interfaceId
	 * @return
	 */
	public InterfaceCluster getCluster(int interfaceId) {
		return clusters.get(interfaceId);
	}
	


	//	public List<Assembly> getAllAssemblies() {
	//		List<Assembly> assemblies = new ArrayList<Assembly>();
	//		
	//		getInterfacesGraph();
	//		
	//		ArrayList<Integer> candidates = getNonParallelInterfacesIds();
	//		System.out.println("Total "+candidates.size()+" interfaces to use ("+(size()-candidates.size())+" parallel interfaces discarded)");
	//		
	//		// we then enumerate all assemblies with 1 interface, 2 interfaces, 3 interfaces .... up to n
	//		System.out.println("Theoretical total assemblies: "+(((int)Math.pow(2, candidates.size()))-1));
	//		for (int n=1;n<=candidates.size();n++) {
	//			List<Assembly> sizenassemblies = enumerateAssembliesSizeN(n,candidates);
	//			
	//			assemblies.addAll(sizenassemblies);
	//		}
	//		
	//		return assemblies;
	//	}
	//	
	//	private List<Assembly> enumerateAssembliesSizeN(int n, ArrayList<Integer> candidates) {
	//		List<Assembly> assemblies = new ArrayList<Assembly>();
	//		CombinationsGenerator x = new CombinationsGenerator (candidates.size(), n);
	//		while (x.hasMore()) {
	//			Assembly ass = new Assembly(this);
	//			int[] indices = x.getNext();
	//			for (int i = 0; i < indices.length; i++) {
	//				ass.add(candidates.get(indices[i]));
	//			}
	//			assemblies.add(ass);
	//		}
	//		System.out.println("size "+n+": "+assemblies.size());
	//		return assemblies;
	//	}
}

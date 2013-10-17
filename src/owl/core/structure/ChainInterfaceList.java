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
import java.util.TreeSet;

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

	private HashMap<Integer,Set<Integer>> clusters; //cluster id, set of Iterface ids

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

	public double getMeanClusterArea(int clusterid){
		double area=0.0;
		for (Integer interfaceid : this.clusters.get(clusterid)){
			area+=this.get(interfaceid).getInterfaceArea();
		}
		return area/this.clusters.get(clusterid).size();
	}

	/**
	 * Cluster interfaces based on the pairwise atom distances
	 * @param pdb  
	 * @param cutoff 
	 * @param testset no of atoms used to calculate pairwise distance
	 * @param atom type of atom
	 * @return hashmap of cluster id, interface id
	 */

	public HashMap<Integer, Set<Integer>> getClusters(PdbAsymUnit pdb,double cutoff,int atomCutoff, String atom){
		if (this.clusters != null)
			return this.clusters;
		HashMap<String, Set<Integer>> commonset = findCommonSet(pdb,atom);
		int commonsetsize = getCommonSetSize(commonset);
		HashMap<Integer, Set<Integer>> clusters = new HashMap<Integer, Set<Integer>>();
		if (commonsetsize<=atomCutoff){
			for (ChainInterface interf:this){
				Set<Integer> iface = new TreeSet<Integer>();
				iface.add(interf.getId());
				clusters.put(interf.getId(), iface);
			}
			this.clusters=clusters;
			return this.clusters;
		}
		double cres=-1.0;
		int clustercount=1;
		for (ChainInterface interf1:this) {
			int interid1 = interf1.getId();
			String chainA1 = new String(interf1.getFirstMolecule().getPdbChainCode());
			String chainA2 = new String(interf1.getSecondMolecule().getPdbChainCode());
			String repchainA1 = new String(pdb.getRepChain(chainA1));
			String repchainA2 = new String(pdb.getRepChain(chainA2));
			Set<Integer> iface = new TreeSet<Integer>();
			if  (!isThereInCluster(clusters,interid1)) {
				for (ChainInterface interf2:this) {
					int interid2 = interf2.getId();
					if (interid2>interid1){
						String chainB1 = new String(interf2.getFirstMolecule().getPdbChainCode());
						String chainB2 = new String(interf2.getSecondMolecule().getPdbChainCode());
						String repchainB1 = new String(pdb.getRepChain(chainB1));
						String repchainB2 = new String(pdb.getRepChain(chainB2));
						if ((!isThereInCluster(clusters,interid2)) && ((repchainA1.equals(repchainB1) && repchainA2.equals(repchainB2)) || 
								(repchainA1.equals(repchainB2) && repchainA2.equals(repchainB1)))){
							//					if (((repchainA1.equals(repchainB1) && repchainA2.equals(repchainB2)) || 
							//							(repchainA1.equals(repchainB2) && repchainA2.equals(repchainB1)))){
							HashMap<String,String> interchains1 = new HashMap<String,String>();
							HashMap<String,String> interchains2 = new HashMap<String,String>();
							interchains1.put(chainA1, repchainA1);
							interchains1.put(chainA2, repchainA2);
							interchains2.put(chainB1, repchainB1);
							interchains2.put(chainB2, repchainB2);
							cres=findDistanceMatrixScore(interf1,interf2,interchains1,interchains2,commonset,atom);
						}
						else{
							cres = -1.00;
						}
						if (cres>=0.0 && cres<cutoff){
							if (iface.size()==0){
								iface.add(interid1);
							}
							iface.add(interid2);
						}
						//	System.out.print(String.format("%.2f\t",cres));
					}
				}
				//			System.out.print(String.format("\n"));	
				if (iface.size()>0){
					clusters.put(clustercount, iface);
					clustercount+=1;
				}else{
					iface.add(interid1);
					clusters.put(clustercount, iface);
					clustercount+=1;
				}
			}//
		}
		this.clusters=clusters;
		return clusters;
	}

	private static int getCommonSetSize(HashMap<String, Set<Integer>> commonset){
		int commonsize=999999;
		for (String repchain:commonset.keySet()){
			if (commonset.get(repchain).size()>0 && commonset.get(repchain).size()<commonsize)
				commonsize=commonset.get(repchain).size();
		}
		return commonsize;
	}

	/**
	 * Returns true if the given interface is already included in a cluster, otherwise false
	 * @param clusters
	 * @param interfaceid
	 * @return boolean
	 */
	private static boolean isThereInCluster(HashMap<Integer, Set<Integer>> clusters, int interfaceid){
		Boolean isthere = false;
		for (Integer key:clusters.keySet()){
			isthere=clusters.get(key).contains(interfaceid);
			if (isthere) break;
		}
		return isthere;
	}

	/**
	 * Returns map of common residues found in all the interfaces
	 * @param pdb
	 * @return map of representative chain, reside id
	 */
	private  static HashMap<String,Set<Integer>> findCommonSet(PdbAsymUnit pdb,String atom){
		List<String> repchains = pdb.getAllRepChains();
		HashMap<String,String> chainmap = new HashMap<String,String>(pdb.getChain2repChainMap());
		HashMap<String, Set<Integer>> commonset = new HashMap<String,Set<Integer>>(); 
		for(String repchain : repchains){
			if (pdb.getChain(repchain).getSequence().isProtein()){
				Set<Integer> chainatoms = new HashSet<Integer>();
				if (chainmap.keySet().size()>1){
					for (Object rep : chainmap.keySet()){
						if (chainmap.get(rep).equals(repchain)){
							if (chainatoms.size() != 0){
								Set<Integer> catoms = new HashSet<Integer>(pdb.getChain(rep.toString()).getAllResSerials());
								checkdensity(pdb,rep.toString(),catoms,atom);
								chainatoms.retainAll(catoms);
							}else{
								Set<Integer> catoms = new HashSet<Integer>(pdb.getChain(rep.toString()).getAllResSerials());
								checkdensity(pdb,rep.toString(),catoms,atom);
								chainatoms=new HashSet<Integer>(catoms);
							}
							checkdensity(pdb,rep.toString(),chainatoms,atom);
						}
					}
				}else{
					Set<Integer> catoms = new HashSet<Integer>(pdb.getChain(repchain).getAllResSerials());
					checkdensity(pdb,repchain,catoms,atom);
					chainatoms = catoms;
				}
				commonset.put(repchain,chainatoms);
			}
		}

		return commonset;
	}

	private static void checkdensity(PdbAsymUnit pdb,String rep,Set<Integer> chainatoms,String atomname){
		for (Iterator<Integer> i=chainatoms.iterator();i.hasNext();){
			Integer resno= i.next();
			Residue res = pdb.getChain(rep).getResidue(resno);
			if (res != null){
				Atom atom = pdb.getChain(rep).getResidue(resno).getAtom(atomname);
				if (atom == null) i.remove();
			}
			else
				i.remove();
		}
	}
	/**
	 * Returns the mean of pairwise CA distance of the common test set of residues between interface 1 and interface 2 
	 * @param inter1
	 * @param inter2
	 * @param chainmap1
	 * @param chainmap2
	 * @param commontestset
	 * @return mean distance
	 */

	private static double findDistanceMatrixScore(ChainInterface inter1,ChainInterface inter2,HashMap<String,String> chainmap1,HashMap<String,String> chainmap2,HashMap<String,Set<Integer>> commontestset,String atom){
		String chainA1 = inter1.getFirstMolecule().getPdbChainCode();
		String chainA2 = inter1.getSecondMolecule().getPdbChainCode();
		String chainB1 = inter2.getFirstMolecule().getPdbChainCode();
		String chainB2 = inter2.getSecondMolecule().getPdbChainCode();
		String repchainA1 = chainmap1.get(chainA1);
		String repchainA2 = chainmap1.get(chainA2);
		String repchainB1 = chainmap2.get(chainB1);
		String repchainB2 = chainmap2.get(chainB2);
		Set<Integer> testresA1 = commontestset.get(repchainA1);
		Set<Integer> testresA2 = commontestset.get(repchainA2);
		Set<Integer> testresB1 = commontestset.get(repchainB1);
		Set<Integer> testresB2 = commontestset.get(repchainB2);
		Point3d[] coord1 = new Point3d[testresA1.size()+testresA2.size()];
		Point3d[] coord2 = new Point3d[testresB1.size()+testresB2.size()];
		double rmsd=0;
		if (repchainA1.equals(repchainB1) && repchainA2.equals(repchainB2) && repchainA1.equals(repchainA2)){
			coord1 = getCommonCoords(inter1,testresA1,testresA2,atom,false);
			coord2 = getCommonCoords(inter2,testresB1,testresB2,atom,false);
			double rmsd1 = GeometryTools.calcOptimalSuperposition(coord1,coord2,false).getRmsd();
			coord2 = getCommonCoords(inter2,testresB1,testresB2,atom,true);
			double rmsd2 = GeometryTools.calcOptimalSuperposition(coord1,coord2,false).getRmsd();
			rmsd=Math.min(rmsd1, rmsd2);
		}else if (repchainA1.equals(repchainB1) && repchainA2.equals(repchainB2)){
			coord1 = getCommonCoords(inter1,testresA1,testresA2,atom,false);
			coord2 = getCommonCoords(inter2,testresB1,testresB2,atom,false);
			rmsd=GeometryTools.calcOptimalSuperposition(coord2, coord1, false).getRmsd();
		}else{
			coord1 = getCommonCoords(inter1,testresA1,testresA2,atom,true);
			coord2 = getCommonCoords(inter2,testresB1,testresB2,atom,true);
			rmsd=GeometryTools.calcOptimalSuperposition(coord2, coord1, false).getRmsd();
		}
		return rmsd;
	}

	private static Point3d[] getCommonCoords(ChainInterface inter,Set<Integer> testres1,Set<Integer> testres2,String atom,boolean isReversed){
		int c=0;
		Point3d[] coord = new Point3d[testres1.size()+testres2.size()];
		if (isReversed){
			for (Integer resid :testres2 ){
				coord[c] = inter.getSecondMolecule().getResidue(resid).getAtom(atom).getCoords();
				c++;
			}
			for (Integer resid :testres1 ){
				coord[c]=inter.getFirstMolecule().getResidue(resid).getAtom(atom).getCoords();
				c++;
			}
		}else{
			for (Integer resid :testres1 ){
				coord[c] = inter.getFirstMolecule().getResidue(resid).getAtom(atom).getCoords();
				c++;
			}
			for (Integer resid :testres2 ){
				coord[c]=inter.getSecondMolecule().getResidue(resid).getAtom(atom).getCoords();
				c++;
			}
		}
		return coord;
	}

	//	/**
	//	 * Returns the sample variance in the array a[], NaN if no such value.
	//	 */
	//	private static double var(double[] a) {
	//		if (a.length == 0) return Double.NaN;
	//		double avg = mean(a);
	//		double sum = 0.0;
	//		for (int i = 0; i < a.length; i++) {
	//			sum += (a[i] - avg) * (a[i] - avg);
	//		}
	//		return sum / (a.length - 1);
	//	}
	//	/**
	//	 * Returns the average value in the array a[], NaN if no such value.
	//	 */
	//	public static double mean(double[] a) {
	//		if (a.length == 0) return Double.NaN;
	//		double sum =0;
	//		for (double x :a ) sum+=x;
	//		return sum / a.length;
	//	}

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

package owl.core.structure.graphs;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import owl.core.structure.Atom;
import owl.core.structure.Residue;
import owl.core.util.GeometryTools;
import edu.uci.ics.jung.graph.SparseGraph;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * An atom inter-chain interaction graph.
 * 
 * @author duarte_j
 *
 */
public class AICGraph extends SparseGraph<Atom,AICGEdge> {

	private static final long serialVersionUID = 1L;

	// disulfide bridges
	public static final double DISULFIDE_BRIDGE_DIST = 2.05;
	public static final double DISULFIDE_BRIDGE_DIST_SIGMA = 0.1;
	
	// hydrogen bonds (ref: Hubbard & Haider)
	public static final double HBOND_UPPER = 3.9; //3.55
	public static final double HBOND_LOWER = 2.5; //2.45
	public static final double HBOND_MIN_ANGLE = Math.PI / 2;
	
	// a generic low distance for a close interaction (electrostatic, salt bridge, semi-covalent, covalent)
	// see review Harding MM, Acta Crystallographica 2006 - 
	// actually distances can be up to 2.4 (or even more) in some cases, taking a conservative approach here 
	public static final double CLOSE_INTERACTION_DIST = 2.1;
	
	// clash distance: in theory, a disulfide bond distance (2.05) is the minimum distance we could reasonably expect
	public static final double CLASH_DISTANCE = 1.5;
	
	// the set of hydrogen-bonded atom pairs
	private HashSet<Pair<Atom>> hBondPairs;
	
	private double distCutoff;
	
	public double getDistCutoff() {
		return distCutoff;
	}
	
	public void setDistCutoff(double distCutoff) {
		this.distCutoff = distCutoff;
	}
	
	public boolean hasClashes() {
		return hasPairsWithinDistance(CLASH_DISTANCE);
	}
	
	public int getNumClashes() {
		return getPairsWithinDistance(CLASH_DISTANCE).size();
	}
	
	public List<Pair<Atom>> getClashingPairs() {
		return getPairsWithinDistance(CLASH_DISTANCE);
	}
	
	public void setHBondPairs(HashSet<Pair<Atom>> hBondPairsI) {
		hBondPairs = hBondPairsI;
	}
	
	public boolean hasPairsWithinDistance(double distance) {
		for (AICGEdge edge:this.getEdges()) {
			if (edge.getDistance()<distance) {
				return true;
			}
		}
		return false;
	}
	
	public List<Pair<Atom>> getPairsWithinDistance(double distance) {
		List<Pair<Atom>> list = new ArrayList<Pair<Atom>>();
		for (AICGEdge edge:this.getEdges()) {
			if (edge.getDistance()<distance) {
				list.add(this.getEndpoints(edge));
			}
		}
		return list;				
	}
	
	
	/**
	 * Returns true if at least one pair of atoms in this atom interchain graph forms 
	 * disulfide bridge. Conditions for it are: they both belong to CYS residues, they are
	 * both SG atoms and the distance between them is 2.05A+-0.1A
	 * @return
	 */
	public boolean hasDisulfideBridges() {
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);
			if (isDisulfideInteraction(pair, edge.getDistance())) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Returns a list of all pairs of atoms that form a disulfide bond in this interchain graph
	 * @return
	 */
	public List<Pair<Atom>> getDisulfidePairs() {
		List<Pair<Atom>> disPairs = new ArrayList<Pair<Atom>>();
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);
			if (isDisulfideInteraction(pair, edge.getDistance())) {
				disPairs.add(pair);
			}
		}
		return disPairs;
	}
	
	/**
	 * Returns true if at least one pair of atoms in this atom interchain graph forms a
	 * hydrogen bond.
	 * @return
	 */
	public boolean hasHbonds() {
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);
			if (isHbondInteraction(pair)) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Returns a list of all pair os atoms in this interchain graph that form a Hydrogen bond.
	 * TODO this can find Hbonds only if H atoms are present, should write a more general algorithm for cases when no H atoms are present 
	 * @return
	 */
	public List<Pair<Atom>> getHbondPairs() {
		List<Pair<Atom>> list = new ArrayList<Pair<Atom>>();
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);
			if (isHbondInteraction(pair)) {
				list.add(pair);
			}
		}
		return list;
	}

	/**
	 * Returns true if any two atoms are closely interacting (distance below {@value #CLOSE_INTERACTION_DIST}) 
	 * and are neither a disulfide bridge nor a hydrogen bond
	 * @return
	 */
	public boolean hasCloselyInteractingPairs() {
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);			
			if (!isHbondInteraction(pair) &&
				!isDisulfideInteraction(pair, edge.getDistance()) &&
				(edge.getDistance()<CLOSE_INTERACTION_DIST)	) {
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Returns a list of atom pairs that are closely interacting (distance below {@value #CLOSE_INTERACTION_DIST}) 
	 * and are neither a disulfide bridge nor a hydrogen bond
	 * @return
	 */
	public List<Pair<Atom>> getCloselyInteractingPairs() {
		List<Pair<Atom>> list = new ArrayList<Pair<Atom>>();
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);
			if (!isHbondInteraction(pair) &&
				!isDisulfideInteraction(pair, edge.getDistance()) &&
				(edge.getDistance()<CLOSE_INTERACTION_DIST)	) {
				list.add(pair);
			}
		}
		return list;
	}
	
	/**
	 * Tells whether given atom pair and distance fulfills the conditions for a Hydrogen bond
	 * @param pair
	 * @param distance
	 * @return
	 */
	protected boolean isHbondInteraction(Pair<Atom> pair) {
		if (hBondPairs != null) {
			return hBondPairs.contains(pair);
		}
		else {
			/*
			 * When HBPlus is not used, hBondPairs will be null and this block will instead be executed.
			 */
			Atom donor = pair.getFirst();
			Atom acceptor = pair.getSecond();
			if (donor.getType() == null || acceptor.getType() == null || !donor.isDonor() || !acceptor.isAcceptor()) {
				return false;
			}
			Atom donorParent = null;
			double donorParentDistance = 0.0;
			for (Atom atom : donor.getParentResidue()) {
				if ((atom != donor) && (!atom.getType().getSymbol().equals("H"))) {
					if (donorParent == null) {
						donorParent = atom;
						donorParentDistance = GeometryTools.distance(atom.getCoords(), donor.getCoords());
					}
					else {
						double thisDistance = GeometryTools.distance(atom.getCoords(), donor.getCoords());
						if (thisDistance < donorParentDistance) {
							donorParent = atom;
							donorParentDistance = thisDistance;
						}
					}
				}
			}
			Atom acceptorParent = null;
			double acceptorParentDistance = 0.0;
			for (Atom atom : acceptor.getParentResidue()) {
				if ((atom != acceptor) && (!atom.getType().getSymbol().equals("H"))) {
					if (acceptorParent == null) {
						acceptorParent = atom;
						acceptorParentDistance = GeometryTools.distance(atom.getCoords(), acceptor.getCoords());
					}
					else {
						double thisDistance = GeometryTools.distance(atom.getCoords(), acceptor.getCoords());
						if (thisDistance < acceptorParentDistance) {
							acceptorParent = atom;
							acceptorParentDistance = thisDistance;
						}
					}
				}
			}
			if (donorParent != null && acceptorParent != null &&     // this can happen if there's only one atom in residue, e.g. 1a3r, last resi of chain P
				 GeometryTools.distance(donor.getCoords(), acceptor.getCoords()) < HBOND_UPPER && 
				 GeometryTools.distance(donor.getCoords(), acceptor.getCoords()) > HBOND_LOWER && 
				 GeometryTools.angle(donorParent.getCoords(), donor.getCoords(), acceptor.getCoords()) > HBOND_MIN_ANGLE && 
				 GeometryTools.angle(donor.getCoords(), acceptor.getCoords(), acceptorParent.getCoords()) > HBOND_MIN_ANGLE) {
					//System.out.println("D is a(n) " + donor.getCode() + " in " + donor.getParentResidue().getLongCode() + donor.getParentResSerial());
					//System.out.println("DD is a(n) " + donorParent.getCode() + " in " + donorParent.getParentResidue().getLongCode() + donorParent.getParentResSerial());
					//System.out.println("A is a(n) " + acceptor.getCode() + " in " + acceptor.getParentResidue().getLongCode() + acceptor.getParentResSerial());
					//System.out.println("AA is a(n) " + acceptorParent.getCode() + " in " + acceptorParent.getParentResidue().getLongCode() + acceptorParent.getParentResSerial());
					//System.out.println("HBOND EXISTS " + GeometryTools.distance(donor.getCoords(), acceptor.getCoords()) + "\n");
					return true;
			}
			else {
				return false;
			}
		}
	}
	
	/**
	 * Tells whether given atom pair and distance fulfils the conditions for a disulfide bridge 
	 * @param pair
	 * @param distance
	 * @return
	 */
	protected boolean isDisulfideInteraction(Pair<Atom> pair, double distance) {
		Atom atomi = pair.getFirst();
		Atom atomj = pair.getSecond();
		if (atomi.getParentResidue().getShortCode()=='C' &&
			atomj.getParentResidue().getShortCode()=='C' &&
			atomi.getCode().equals("SG") &&
			atomj.getCode().equals("SG") &&
			distance<(DISULFIDE_BRIDGE_DIST+DISULFIDE_BRIDGE_DIST_SIGMA) && 
			distance>(DISULFIDE_BRIDGE_DIST-DISULFIDE_BRIDGE_DIST_SIGMA)) {
				return true;
		}
		return false;
	}
	
	public double getAvrgNumNeighbors() {
		int total = 0;
		int count = 0;
		for (Atom atom:this.getVertices()) {
			total+=this.getNeighborCount(atom);
			count++;
		}
		if (total==0) return 0;
		else return (double)total/(double)count;
	}
	
	public double getAvrgNumNeighbors(Collection<Residue> iresidues, Collection<Residue> jresidues) {
		int total = 0;
		int count = 0;
		//System.out.println("FIRST");
		for (Atom atom:this.getFirstVertices()) {
			if (isAtomInResidues(atom, iresidues)) {
				total+=this.getNeighborCount(atom);
				count++;
				//System.out.print(atom.getCode()+atom.getSerial()+": "+getNeighborCount(atom)+" ");
			}
		}
		//System.out.println();
		//System.out.println("SECOND");
		for (Atom atom:this.getSecondVertices()) {
			if (isAtomInResidues(atom, jresidues)) {
				total+=this.getNeighborCount(atom);
				count++;
				//System.out.print(atom.getCode()+atom.getSerial()+": "+getNeighborCount(atom)+" ");
			}			
		}
		//System.out.println();
		if (total==0) return 0;
		else return (double)total/(double)count;
	}

	private Set<Atom> getFirstVertices() {
		Set<Atom> atoms = new HashSet<Atom>();
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);
			atoms.add(pair.getFirst());
		}
		return atoms;
	}

	private Set<Atom> getSecondVertices() {
		Set<Atom> atoms = new HashSet<Atom>();
		for (AICGEdge edge:this.getEdges()) {
			Pair<Atom> pair = this.getEndpoints(edge);
			atoms.add(pair.getSecond());
		}
		return atoms;
	}
	
	private static boolean isAtomInResidues(Atom atom, Collection<Residue> residues) {
		for (Residue residue:residues) {
			Residue parent = atom.getParentResidue();
			if (parent.getSerial()==residue.getSerial() && 
				parent.getLongCode().equals(residue.getLongCode())) {
				return true;
			}
		}
		return false;
	}
	

}

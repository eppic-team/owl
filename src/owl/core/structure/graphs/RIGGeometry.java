package owl.core.structure.graphs;

import java.util.HashMap;
import java.util.TreeMap;
import java.util.Map.Entry;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import edu.uci.ics.jung.graph.util.Pair;

import owl.core.structure.Residue;
import owl.gmbp.GmbpGeometry;

//RIGGeometry: class that computes and stores geometry information for graph based on rotation and translation invariant framework.

public class RIGGeometry {
	
	private RIGraph graph;
	private TreeMap<Integer, Residue> residues;
	private HashMap<String,Vector3d> coord_sph_rotated;  // holds geometry (translated and rotated contact coordinates of CA-position
							//of j-Residue with respect to central iResidue of each edge of the graph, defined by residue serials of edge)
	private String atomType = "CA";
	
	
	/**
	 * Constructs a RIGraph with a sequence but no edges
	 * @param sequence
	 */
	public RIGGeometry(RIGraph graph, TreeMap<Integer, Residue> residues) {
		this.graph = graph;
		this.residues = residues;
		this.atomType = "CA";
		
		initialiseGeometry();
	}
	
	/**
	 * Constructs a RIGraph with a sequence but no edges
	 * @param sequence
	 */
	public RIGGeometry(RIGraph graph, TreeMap<Integer, Residue> residues, String atomType) {
		this.graph = graph;
		this.residues = residues;
		this.atomType = atomType;
		
		initialiseGeometry();
	}
	
	private void initialiseGeometry(){
		
		System.out.println("Geometry of graph with "+this.graph.getEdgeCount()+" edges:");
		int edgeNum = 0;
		GmbpGeometry gmbp = new GmbpGeometry();
		
//		HashMap<String,Vector3d> coord_sph_rotated = new HashMap<String, Vector3d>();
		coord_sph_rotated = new HashMap<String,Vector3d>();
		for (RIGEdge edge:graph.getEdges()) {
			edgeNum++;
			// extract (x,y,z) coordinates for nodes of both end of edge
			Pair<RIGNode> nodes = this.graph.getEndpoints(edge);
			RIGNode iNode = nodes.getFirst();
			RIGNode jNode = nodes.getSecond();
			int iNum = iNode.getResidueSerial();
			int jNum = jNode.getResidueSerial();
			String iResType = iNode.getResidueType();
			String jResType = jNode.getResidueType();
			
			Residue iRes = this.residues.get(iNum);
			Residue jRes = this.residues.get(jNum);
			
			// Testoutput
//			Atom iAtom = iRes.getAtom("CA");
//			Atom jAtom = jRes.getAtom("CA");
//			Point3d iCoordCA = iAtom.getCoords();
//			Point3d jCoordCA = jAtom.getCoords();
//			System.out.println("EdgeNr="+edgeNum+" between "+iNum+iResType+"("+iCoordCA.x+","+iCoordCA.y+","+iCoordCA.z+")"
//					+" and "+jNum+jResType+"("+jCoordCA.x+","+jCoordCA.y+","+jCoordCA.z+")");
//			System.out.printf("EdgeNr= %s between %s%s %s and %s%s %s  ",edgeNum,iNum,iResType,iCoordCA,jNum,jResType,jCoordCA);
			
			
			// translate coordinates with rotation and translation invariant framework
			HashMap<String,Point3d> iCoord = gmbp.getTheResidueCoord(iRes);
			HashMap<String,Point3d> jCoord = gmbp.getTheResidueCoord(jRes);
			
			Vector3d coord_I = new Vector3d(0,0,0);
			Vector3d coord_J = new Vector3d(0,0,0);
			// LEAVE this EDGE and CONTINUE with next one if the above condition is not satisfied.
			HashMap<String,Vector3d> atom_coord_rotated_I = new HashMap<String, Vector3d>();
			HashMap<String,Vector3d> atom_coord_rotated_J = new HashMap<String, Vector3d>();
			if (iCoord!=null && jCoord!=null) {
				// METHOD FOR EXTRACTING THE NEIGHBOR'S TRANSLATED-ROTATED COORDINATES //
				atom_coord_rotated_I = gmbp.getNeighborsTransRotatedCoord(jCoord, iCoord, jResType, iResType, jRes, iRes, true);
				atom_coord_rotated_J = gmbp.getNeighborsTransRotatedCoord(iCoord, jCoord, iResType, jResType, iRes, jRes, true);
				//========= C-ALPHA-coordinates ================// 
				coord_I = atom_coord_rotated_I.get(this.atomType); //("CA");
				coord_J = atom_coord_rotated_J.get(this.atomType); //("CA");
			}
			else {
				continue;
			}
			
			// GET the SPHERICAL COORDINATES for CA, C, CB, and CG using METHOD "getSphericalFromCartesian", 
			// if these coordinates exist (are non-zero).
			Vector3d coord_sph_I = new Vector3d(0,0,0);
			Vector3d coord_sph_J = new Vector3d(0,0,0);
			if (!coord_J.equals(new Vector3d(0.0,0.0,0.0))) {
				coord_sph_J = gmbp.getSphericalFromCartesian(coord_J); // (r,theta,phi) // (r, phi,lambda)
			}
			if (!coord_I.equals(new Vector3d(0.0,0.0,0.0))) {
				coord_sph_I = gmbp.getSphericalFromCartesian(coord_I); // (r,theta,phi) // (r, phi,lambda)
			}
			
			// Save translated and rotated coordinate of contact
			String key = String.valueOf(jNum)+"_"+String.valueOf(iNum);
			coord_sph_rotated.put(key, coord_sph_I);
			key = String.valueOf(iNum)+"_"+String.valueOf(jNum);
			coord_sph_rotated.put(key, coord_sph_J);
			if (!coord_sph_rotated.containsKey(key)){
				System.out.println("Something is going wrong");
			}
			
//			System.out.println("TransRotCoord Cartesian: "+coord.x+","+coord.y+","+coord.z
//					+" SPH: "+coord_sph.x+","+coord_sph.y+","+coord_sph.z);
//			System.out.printf("TransRotCoord i->j Cartesian: %s Spherical: %s   j->i Cartesian: %s Spherical: %s \n", coord_I, coord_sph_I, coord_J, coord_sph_J);
		}
		
//		printGeom();
		
//		if (coord_sph_rotated.containsKey(new Integer[]{iNum,jNum})){
//			Vector3d CA_coord_sph = coord_sph_rotated.get(new Integer[]{iNum,jNum});
//		}
		
//		// iterate over translated coordinates
//		for (Entry<Integer[], Vector3d> entry : coord_sph_rotated.entrySet()) {
//  		    System.out.printf("Contact %s,%s coord: %s \n", entry.getKey()[0], entry.getKey()[1], entry
// 					.getValue());
//  		}
		
//		this.graph.containsEdgeIJ(i, j);
//		this.graph.containsEdge(edge);
//		this.graph.getEdgeFromSerials(i, j);
//		for (Residue residue:residues.values()) {
//			for (Atom atom:residue.getAtoms()) {
//				
//			}
//		}
	}
	
	public void printGeom(){
		for (Entry<String, Vector3d> entry: this.coord_sph_rotated.entrySet()){
			System.out.println(entry.getKey()+":"+entry.getValue().x);
		}
	}
	
	// --------- getters ----------
	public HashMap<String,Vector3d> getRotatedCoordOfContacts(){
		return this.coord_sph_rotated;
	}
	
	public String getAtomType(){
		return this.atomType;
	}

}

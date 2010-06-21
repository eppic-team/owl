package owl.core.structure.graphs;

import java.util.HashMap;
import java.util.TreeMap;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import edu.uci.ics.jung.graph.util.Pair;

import owl.core.structure.Atom;
import owl.core.structure.Residue;
import owl.gmbp.GmbpGeometry;

//RIGGeometry: class that computes and stores geometry information for graph based on rotation and translation invariant framework.

public class RIGGeometry {
	
	private RIGraph graph;
	private TreeMap<Integer, Residue> residues;
	private HashMap<Integer[],Vector3d> ca_coord_sph_rotated;  // holds geometry (translated and rotated contact coordinates of CA-position
							//of j-Residue with respect to central iResidue of each edge of the graph, defined by residue serials of edge) 
	
	
	/**
	 * Constructs a RIGraph with a sequence but no edges
	 * @param sequence
	 */
	public RIGGeometry(RIGraph graph, TreeMap<Integer, Residue> residues) {
		this.graph = graph;
		this.residues = residues;
		
		initialiseGeometry();
	}
	
	private void initialiseGeometry(){
		
		System.out.println("Geometry of graph with "+this.graph.getEdgeCount()+" edges:");
		int edgeNum = 0;
		GmbpGeometry gmbp = new GmbpGeometry();
		
//		HashMap<String,Vector3d> ca_coord_sph_rotated = new HashMap<String, Vector3d>();
		ca_coord_sph_rotated = new HashMap<Integer[],Vector3d>();
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
			Atom iAtom = iRes.getAtom("CA");
			Atom jAtom = jRes.getAtom("CA");
			Point3d iCoordCA = iAtom.getCoords();
			Point3d jCoordCA = jAtom.getCoords();
//			System.out.println("EdgeNr="+edgeNum+" between "+iNum+iResType+"("+iCoordCA.x+","+iCoordCA.y+","+iCoordCA.z+")"
//					+" and "+jNum+jResType+"("+jCoordCA.x+","+jCoordCA.y+","+jCoordCA.z+")");
			System.out.printf("EdgeNr= %s between %s%s %s and %s%s %s  ",edgeNum,iNum,iResType,iCoordCA,jNum,jResType,jCoordCA);
			
			
			// translate coordinates with rotation and translation invariant framework
			HashMap<String,Point3d> iCoord = gmbp.getTheResidueCoord(iRes);
			HashMap<String,Point3d> jCoord = gmbp.getTheResidueCoord(jRes);
			
			Vector3d CA_coord = new Vector3d(0,0,0);
			
			// LEAVE this EDGE and CONTINUE with next one if the above condition is not satisfied.
			HashMap<String,Vector3d> atom_coord_rotated = new HashMap<String, Vector3d>();
			if (iCoord!=null && jCoord!=null) {
				// METHOD FOR EXTRACTING THE NEIGHBOR'S TRANSLATED-ROTATED COORDINATES //
				atom_coord_rotated = gmbp.getNeighborsTransRotatedCoord(iCoord, jCoord, iResType, jResType, iRes, jRes, true);
				//========= C-ALPHA and C-coordinates ================// 
				CA_coord = atom_coord_rotated.get("CA");
			}
			else {
				continue;
			}
			
			// GET the SPHERICAL COORDINATES for CA, C, CB, and CG using METHOD "getSphericalFromCartesian", 
			// if these coordinates exist (are non-zero).
			Vector3d CA_coord_sph = new Vector3d(0,0,0);
			if (!CA_coord.equals(new Vector3d(0.0,0.0,0.0))) {
				CA_coord_sph = gmbp.getSphericalFromKartesian(CA_coord); // (r,theta,phi) // (r, phi,lambda)
			}
			
			// Save translated and rotated coordinate of contact
			ca_coord_sph_rotated.put(new Integer[]{iNum,jNum}, CA_coord_sph);
			
//			System.out.println("TransRotCoord Cartesian: "+CA_coord.x+","+CA_coord.y+","+CA_coord.z
//					+" SPH: "+CA_coord_sph.x+","+CA_coord_sph.y+","+CA_coord_sph.z);
			System.out.printf("TransRotCoord Cartesian: %s Spherical: %s \n", CA_coord, CA_coord_sph);
		}
		
//		// iterate over translated coordinates
//		for (Entry<Integer[], Vector3d> entry : ca_coord_sph_rotated.entrySet()) {
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
	
	// --------- getters ----------
	public HashMap<Integer[],Vector3d> getRotatedCoordOfContacts(){
		return this.ca_coord_sph_rotated;
	}

}

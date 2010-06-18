package owl.core.structure.graphs;

import java.util.HashMap;
import java.util.Set;
import java.util.TreeMap;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import edu.uci.ics.jung.graph.util.Pair;

import owl.core.structure.AAinfo;
import owl.core.structure.Atom;
import owl.core.structure.Residue;

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
			HashMap<String,Point3d> iCoord = getResidueCoord(iRes, iResType);
			HashMap<String,Point3d> jCoord = getResidueCoord(jRes, jResType);
			
			Vector3d CA_coord = new Vector3d(0,0,0);
			
			// LEAVE this EDGE and CONTINUE with next one if the above condition is not satisfied.
			HashMap<String,Vector3d> atom_coord_rotated = new HashMap<String, Vector3d>();
			if (iCoord!=null && jCoord!=null) {
				// METHOD FOR EXTRACTING THE NEIGHBOR'S TRANSLATED-ROTATED COORDINATES //
				atom_coord_rotated = getNeighborsTransRotatedCoord(iCoord, jCoord, iResType, jResType, iRes, jRes, true);
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
				CA_coord_sph = getSphericalFromKartesian(CA_coord); // (r,theta,phi) // (r, phi,lambda)
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
	
	private HashMap<String,Point3d> getResidueCoord(Residue res, String resType){
		//Method to return the coordinates of specific residue. 
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>();		

		//Test if all necessary atoms CA, C and N available
		if (res.containsAtom("CA") && res.containsAtom("C") && res.containsAtom("N")){
			for(Atom atom:res.getAtoms()){
				atom_coord.put(atom.getCode(), atom.getCoords());
			}			
		}
		
		return atom_coord;
	}
	
	private static HashMap<String,Vector3d> getNeighborsTransRotatedCoord (HashMap<String,Point3d> iCoord, HashMap<String,Point3d> jCoord, 
			String iResType, String jResType, Residue iRes, Residue jRes, Boolean toRotate) {
		// METHOD for finding the transformed coordinates of the neighboring residue of a 'node'
		// Doing it according to the EDGE-LIST

		// ONE-TIME Definitions used for the Central Residue as well as the Neighborhood Residues
		HashMap<String,Vector3d> atom_coord_rotated = new HashMap<String, Vector3d>();
		Matrix4d RotMatrix01 = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);
		Matrix4d RotMatrix02 = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);

		// Get the NEW ORIGIN: c-alpha Coordinates of the central residue
		Vector3d new_origin = new Vector3d (iCoord.get("CA"));//The Origin of the New Neighborhood

		// Get the translated (to c-alpha) coordinates
		atom_coord_rotated = getTranslatedCoord(iRes, iResType, iCoord);
		// CENTRAL RESIDUE: TRANSLATION TO ORIGIN OVER // 		

		if (toRotate) {//Performing the rotations only if the  "toRotate=true"
			RotMatrix01 = getRotationMatrix01 (atom_coord_rotated);
			atom_coord_rotated = getCoordAfterRotation (iRes, iResType, atom_coord_rotated, RotMatrix01);
			// CENTRAL RESIDUE: FIRST ROTATION OVER // 		
	
			RotMatrix02 = getRotationMatrix02 (atom_coord_rotated);
			atom_coord_rotated = getCoordAfterRotation (iRes, iResType, atom_coord_rotated, RotMatrix02);
			// CENTRAL RESIDUE: SECOND ROTATION OVER //
		}

		// FOR EXTRACTING & WRITING THE TRANS-ROTATED COORDINATE INFORMATION OF A *SPECIFIC NEIGHBOR* OF THE *PRESENT NODE*. 
		atom_coord_rotated = getTranslatedCoordNbr (jRes, jResType, jCoord, new_origin); 
		if (toRotate) {//Obtaining the rotated coordinates only if the  "toRotate=true"
			atom_coord_rotated = getCoordAfterRotation(jRes, jResType, atom_coord_rotated, RotMatrix01); 
			atom_coord_rotated = getCoordAfterRotation (jRes, jResType, atom_coord_rotated, RotMatrix02);
		}
		return atom_coord_rotated;
	}
	
	private static HashMap<String,Vector3d> getTranslatedCoord (Residue res, String res_type, HashMap<String,Point3d> atom_coord) {	
		// Method for getting the coordinates "translated (shifted)" to the new-origin
		// NOTE: The origin is extracted from the available coordinates. Origin defined as the coordinates of CA atom. 
		
		Vector3d new_origin = new Vector3d (atom_coord.get("CA")); //Defining the ORIGIN

		Set<String> atoms = AAinfo.getAtomsForCTAndRes("ALL", res_type);

		// Recomputing the coordinates/vectors with the new origin
		HashMap<String,Vector3d> atom_coord_rotated = new HashMap<String, Vector3d>();
		for (String atom:atoms) {
			if (res.containsAtom(atom)) // Making sure the text-book atom indeed exists! 
			{ 	
				atom_coord_rotated.put(atom, new Vector3d(atom_coord.get(atom)));
				atom_coord_rotated.get(atom).sub(new_origin);
			}
		}

		return atom_coord_rotated;
	}
	
	private static HashMap<String,Vector3d> getTranslatedCoordNbr (Residue res_nbr, String res_type_nbr, HashMap<String,Point3d> atom_coord, Vector3d new_origin) {	
		// Method for getting the coordinates of the neighboring residue, "translated (shifted)" to wrt the central residue (origin).
		// NOTE: The reference origin is passed on as a parameter. 

		Set<String> atoms = AAinfo.getAtomsForCTAndRes("ALL", res_type_nbr);

		// Recomputing the coordinates/vectors with the new origin
		HashMap<String,Vector3d> atom_coord_rotated = new HashMap<String, Vector3d>();
		for (String atom:atoms) {
//			if (pdb.hasCoordinates(resser_nbr, atom)) {// Making sure the text-book atom indeed exists!
			if (res_nbr.containsAtom(atom)) {// Making sure the text-book atom indeed exists! 
				atom_coord_rotated.put(atom, new Vector3d(atom_coord.get(atom)));
				atom_coord_rotated.get(atom).sub(new_origin);//Shifting each coordinate of neighbor wrt "C-Alpha of Central Residue"   
			}
		}

		return atom_coord_rotated;
	}	

	private static HashMap<String,Vector3d> getCoordAfterRotation (Residue res, String res_type, HashMap<String,Vector3d> atom_coord_rotated, Matrix4d RotMatrix_4d) {
		//Method to return coordinates after a certain transformation (as specified by a Matrix4d operator).
		
		Set<String> atoms = AAinfo.getAtomsForCTAndRes("ALL", res_type);	
		//=========== PERFORMING THE FIRST ROTATION WITH theta01 ANGLE ========//
		// Computing the vectors in the rotated frame of reference
		for (String atom:atoms) {
//			if (pdb.hasCoordinates(resser, atom)) {// Making sure the text-book atom indeed exists!
			if (res.containsAtom(atom)) {// Making sure the text-book atom indeed exists! 
				RotMatrix_4d.transform(atom_coord_rotated.get(atom));
			}
		}

		return atom_coord_rotated;
	}
	
	private static Matrix4d getRotationMatrix01 (HashMap<String,Vector3d> atom_coord_rotated) {

		//Defining the unit vectors for X-, Y- and Z-axis.
		Vector3d iVec = new Vector3d(1,0,0); 

		//================== STEP 01: THE FIRST ROTATION =========================//
		// Finding out the arbitrary vector 'r' around which rotation is to be made. 
		Vector3d rVec = new Vector3d(0,0,0);
		Vector3d cVec = new Vector3d (atom_coord_rotated.get("C"));//Defining ca-c vector
		rVec.cross(cVec, iVec); // Setting r = v (ca-->c) x i
		rVec.normalize(rVec); // Normalizing 'rVec'	

		// Finding the angle by which the FIRST rotation is to be done
		double theta01 = cVec.angle(iVec); // Angle between ca-c vector and the x-axis

		//======== Generating Rotation Matrix Automatically ==========//
		//Setting up the Rotation Axis (rVec)  and the Rotation Angle (theta01)
		AxisAngle4d axisangle01 = new AxisAngle4d();
		axisangle01 = new AxisAngle4d(rVec, theta01);
		//axisangle01 = new AxisAngle4d(0.0, 1.0, 0.0, 0.5236); //Testing (y-axis/30 Degrees)

		// Creating a Matrix4d Element for Rotation Matrix. 
		Matrix4d RotMatrix01_4d = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);

		// Generating the ROTATION MATRIX using SetRotation
		RotMatrix01_4d.setRotation(axisangle01);
		//----------- Generating Rotation Matrix Automatically ------------//

		return RotMatrix01_4d; 

	}   

	private static Matrix4d getRotationMatrix02 (HashMap<String,Vector3d> atom_coord_rotated) {

		//Defining the unit vectors for X-, Y- and Z-axis.
		Vector3d iVec = new Vector3d(1,0,0); 
		Vector3d kVec = new Vector3d(0,0,1); 

		//================== STEP 02: THE SECOND ROTATION =========================//				
		Vector3d NVec = new Vector3d (atom_coord_rotated.get("N"));//Defining ca-N vector
		Vector3d wp = new Vector3d(0.0, NVec.y, NVec.z);//The projection of ca-N vector on YZ-plane 
		// NOTE: The projection (wp) needs to be taken on the plane that is perpendicular to the axis 
		// around which one intends to rotate the space eventually.  

		// Finding the angle by which the SECOND rotation is to be done
		// depending on the position of the XYprojection(ca-N) OR wp vector.
		// WARNING: The angle calculation is presently SCREWED UP. Need to fix it!!
		double theta02 = 0.0;
		if (NVec.y > 0.0) {
			theta02 = wp.angle(kVec); // Positive angle made by wp vector with the Z-axis
		} else {
			theta02 = -wp.angle(kVec); // Negative angle made by wp vector with the Z-axis
		}

		//======== Generating Rotation Matrix Automatically ==========//
		//Setting up the Rotation Axis (X-axis, i) and the Rotation Angle (theta02)
		AxisAngle4d axisangle02 = new AxisAngle4d();
		axisangle02 = new AxisAngle4d(iVec, theta02);//Rotation around the X-axis

		// Creating a Matrix4d Element for Rotation Matrix. 
		Matrix4d RotMatrix02_4d = new Matrix4d(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1);

		// Generating the ROTATION MATRIX using SetRotation
		RotMatrix02_4d.setRotation(axisangle02);
		//----------- Generating Rotation Matrix Automatically ------------//

		return RotMatrix02_4d;      
	}        
	
	private static Vector3d getSphericalFromKartesian (Vector3d Atom_Coord) {
		// METHOD for conversion from Cartesian to Spherical Coordinates
		Vector3d Atom_Coord_sph = new Vector3d(0,0,0);

		Atom_Coord_sph.x = Math.sqrt( Math.pow(Atom_Coord.x, 2) + Math.pow(Atom_Coord.y, 2) + Math.pow(Atom_Coord.z, 2) ); //r
		//Atom_Coord_sph.y = Math.atan2(Math.sqrt(Math.pow(Atom_Coord.x, 2)+Math.pow(Atom_Coord.y, 2)), Atom_Coord.z); //theta
		Atom_Coord_sph.y = Math.acos( Atom_Coord.z/Atom_Coord_sph.x); //theta
		Atom_Coord_sph.z = Math.atan2(Atom_Coord.y, Atom_Coord.x); //phi

		return Atom_Coord_sph;
	}
	
	// --------- getters ----------
	public HashMap<Integer[],Vector3d> getRotatedCoordOfContacts(){
		return this.ca_coord_sph_rotated;
	}

}

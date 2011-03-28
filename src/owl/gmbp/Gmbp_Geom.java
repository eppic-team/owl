package owl.gmbp;

import java.util.HashMap;
import java.util.TreeMap;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import owl.core.structure.Atom;
import owl.core.structure.Pdb;
import owl.core.structure.PdbLoadException;
import owl.core.structure.PdbfilePdb;
import owl.core.structure.Residue;

public class Gmbp_Geom {
	
//	private Residue iRes;
//	private Residue jRes;
	
	private final static String centre="CA";
	private final static String right="C";
	private final static String left="N";
	

	// ---- MAIN for TestPurposes
	/**
	 * @param args
	 * @throws PdbLoadException 
	 */
	public static void main(String[] args) throws PdbLoadException {
		// TODO Auto-generated method stub

		GmbpGeometry gmbp = new GmbpGeometry();
		
		String str="A";
//		Pdb proteinToLoad=new PdbfilePdb ("/home/belsare/Desktop/1E0L.pdb");
		Pdb proteinToLoad=new PdbfilePdb ("/Users/vehlow/Documents/workspace/7ODC.pdb");
		proteinToLoad.load(str, 1);
		TreeMap<Integer, Residue> protein=new TreeMap<Integer, Residue>(proteinToLoad.getResidues());
		
		int resno=26;
		Residue res=protein.get(resno);
		String resType = res.getAaType().getThreeLetterCode();
		System.out.println("ResType"+resType);
		Atom atom=res.getAtom("CA");
		Point3d coords=atom.getCoords();
		System.out.println("CA TransRotCoord Cartesian: "+coords.x+", "+coords.y+", "+coords.z);
		atom=res.getAtom("C");
		coords=atom.getCoords();
		System.out.println("C TransRotCoord Cartesian: "+coords.x+", "+coords.y+", "+coords.z);
		atom=res.getAtom("N");
		coords=atom.getCoords();
		System.out.println("N TransRotCoord Cartesian: "+coords.x+", "+coords.y+", "+coords.z);
		atom=res.getAtom("CB");
		coords=atom.getCoords();
		System.out.println("CB TransRotCoord Cartesian: "+coords.x+", "+coords.y+", "+coords.z);
		//System.out.println(coords.x+"  "+coords.y+"   "+coords.z);
			
		// Transformation procedures by Saurabh
		HashMap<String,Point3d> atom_coord_cart = getRelativeCoord(res);
		HashMap<String, Point3d> atom_coord_sph = getRotated2SphericalCoords(atom_coord_cart);		

		System.out.println("by Saurabh:");
		System.out.println("CA TransRotCoord Cartesian: "+atom_coord_cart.get("CA").x+", "+atom_coord_cart.get("CA").y+", "+atom_coord_cart.get("CA").z
				+" SPH: "+atom_coord_sph.get("CA").x+", "+atom_coord_sph.get("CA").y+", "+atom_coord_sph.get("CA").z);
		System.out.println("C TransRotCoord Cartesian: "+atom_coord_cart.get("C").x+", "+atom_coord_cart.get("C").y+", "+atom_coord_cart.get("C").z
				+" SPH: "+atom_coord_sph.get("C").x+", "+atom_coord_sph.get("C").y+", "+atom_coord_sph.get("C").z);
		System.out.println("N TransRotCoord Cartesian: "+atom_coord_cart.get("N").x+", "+atom_coord_cart.get("N").y+", "+atom_coord_cart.get("N").z
				+" SPH: "+atom_coord_sph.get("N").x+", "+atom_coord_sph.get("N").y+", "+atom_coord_sph.get("N").z);
		System.out.println("CB TransRotCoord Cartesian: "+atom_coord_cart.get("CB").x+", "+atom_coord_cart.get("CB").y+", "+atom_coord_cart.get("CB").z
				+" SPH: "+atom_coord_sph.get("CB").x+", "+atom_coord_sph.get("CB").y+", "+atom_coord_sph.get("CB").z);
				
		// old tranformation procedures --> implemented in GmbpGeometry
		// translate coordinates with rotation and translation invariant framework
		HashMap<String,Point3d> iCoord = gmbp.getTheResidueCoord(res);		
		Vector3d coord = new Vector3d(0,0,0);
		Vector3d coord_sph = new Vector3d(0,0,0);
		// LEAVE this EDGE and CONTINUE with next one if the above condition is not satisfied.
		HashMap<String,Vector3d> atom_coord_cart2 = new HashMap<String, Vector3d>();
		if (iCoord!=null) {
			// METHOD FOR EXTRACTING THE NEIGHBOR'S TRANSLATED-ROTATED COORDINATES //
			atom_coord_cart2 = gmbp.getNeighborsTransRotatedCoord(iCoord, iCoord, resType, resType, res, res, true);
			//========= C-ALPHA-coordinates ================// 
			System.out.println("old Transformation:");
			coord = atom_coord_cart2.get("CA"); 
			coord_sph = gmbp.getSphericalFromCartesian(coord); // (r,theta,phi) // (r, phi,lambda)
			System.out.println("CA TransRotCoord Cartesian: "+coord.x+", "+coord.y+", "+coord.z
					+" SPH: "+coord_sph.x+", "+coord_sph.y+", "+coord_sph.z);
			coord = atom_coord_cart2.get("C"); 
			coord_sph = gmbp.getSphericalFromCartesian(coord); // (r,theta,phi) // (r, phi,lambda)
			System.out.println("C TransRotCoord Cartesian: "+coord.x+", "+coord.y+", "+coord.z
					+" SPH: "+coord_sph.x+", "+coord_sph.y+", "+coord_sph.z);
			coord = atom_coord_cart2.get("N"); 
			coord_sph = gmbp.getSphericalFromCartesian(coord); // (r,theta,phi) // (r, phi,lambda)
			System.out.println("N TransRotCoord Cartesian: "+coord.x+", "+coord.y+", "+coord.z
					+" SPH: "+coord_sph.x+", "+coord_sph.y+", "+coord_sph.z);
			coord = atom_coord_cart2.get("CB"); 
			coord_sph = gmbp.getSphericalFromCartesian(coord); // (r,theta,phi) // (r, phi,lambda)
			System.out.println("CB TransRotCoord Cartesian: "+coord.x+", "+coord.y+", "+coord.z
					+" SPH: "+coord_sph.x+", "+coord_sph.y+", "+coord_sph.z);
		}
				
		
	}
	

	/**
	 * Extracts and returns the coordinates of all atoms (CA,C,N,...) from residue
	 * @param residue 
	 * @return coordinates all atoms of residue
	 */
	public static HashMap<String,Point3d> getRelativeCoord (Residue iRes) {
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>();
		//TODO: compute and return the relative coordinate of all available atoms of the residue
		// coord of CAatom of iRes should be (0,0,0) ...

		//Test if all necessary atoms CA, C and N available
		if (iRes.containsAtom("CA") && iRes.containsAtom("C") && iRes.containsAtom("N")){
			HashMap<String, Point3d> translatedCoordinates = getTranslatedCoordinates(iRes, centre); //new HashMap<String, Point3d>(getTranslatedCoordinates(iRes, centre));			
			HashMap<String, Point3d> translatedSphericalCoordinates = getTranslatedSphericalCoordinates(translatedCoordinates); //new HashMap<String, Point3d>(getTranslatedSphericalCoordinates(translatedCoordinates));			
			HashMap<String, Point3d> rotated1SphericalCoordinates = getRotated1SphericalCoords(translatedSphericalCoordinates, right); //new HashMap<String, Point3d>(getRotated1SphericalCoords(translatedSphericalCoordinates, right));			
			HashMap<String, Point3d> rotated1CartesianCoordinates = getRotated1CartesianCoords(rotated1SphericalCoordinates); //new HashMap<String, Point3d>(getRotated1CartesianCoords(rotated1SphericalCoordinates));			
			HashMap<String, Point3d> rotated2CartesianCoordinates = getRotated2CartesianCoords(rotated1CartesianCoordinates, left); //new HashMap<String, Point3d>(getRotated2CartesianCoords(rotated1CartesianCoordinates, left));			
//			HashMap<String, Point3d> rotated2SphericalCoordinates = getRotated2SphericalCoords(rotated2CartesianCoordinates); //new HashMap<String, Point3d>(getRotated2SphericalCoords(rotated2CartesianCoordinates));
			
			atom_coord = rotated2CartesianCoordinates;			

//			System.out.println(iRes.getAtomsMap().get("CA").getCoords().x+"  "+iRes.getAtomsMap().get("CA").getCoords().y+"   "+iRes.getAtomsMap().get("CA").getCoords().z);
//			System.out.println(translatedSphericalCoordinates.get("C").x+"  "+translatedSphericalCoordinates.get("C").y+"   "+translatedSphericalCoordinates.get("C").z);
//			System.out.println(rotated1SphericalCoordinates.get("C").x+"  "+rotated1SphericalCoordinates.get("C").y+"   "+rotated1SphericalCoordinates.get("C").z);
//			System.out.println(rotated1CartesianCoordinates.get("C").x+"  "+rotated1CartesianCoordinates.get("C").y+"   "+rotated1CartesianCoordinates.get("C").z);
//			System.out.println(rotated2CartesianCoordinates.get("C").x+"  "+rotated2CartesianCoordinates.get("C").y+"   "+rotated2CartesianCoordinates.get("C").z);
//			System.out.println(rotated2CartesianCoordinates.get("N").x+"  "+rotated2CartesianCoordinates.get("N").y+"   "+rotated2CartesianCoordinates.get("N").z);
//			System.out.println(rotated2SphericalCoordinates.get("C").x+"  "+rotated2SphericalCoordinates.get("C").y+"   "+rotated2SphericalCoordinates.get("C").z);
//			System.out.println(rotated2SphericalCoordinates.get("N").x+"  "+rotated2SphericalCoordinates.get("N").y+"   "+rotated2SphericalCoordinates.get("N").z);

			// Test PrintOut:
////			System.out.println("CA coord:  "+iRes.getAtomsMap().get("CA").getCoords().x+"  "+iRes.getAtomsMap().get("CA").getCoords().y+"   "+iRes.getAtomsMap().get("CA").getCoords().z);
//			System.out.println("Cartesian Coordinates after Procedure:");			
//			System.out.println(rotated2CartesianCoordinates.get("CA").x+"  "+rotated2CartesianCoordinates.get("CA").y+"   "+rotated2CartesianCoordinates.get("CA").z);
//			System.out.println(rotated2CartesianCoordinates.get("C").x+"  "+rotated2CartesianCoordinates.get("C").y+"   "+rotated2CartesianCoordinates.get("C").z);
//			System.out.println(rotated2CartesianCoordinates.get("N").x+"  "+rotated2CartesianCoordinates.get("N").y+"   "+rotated2CartesianCoordinates.get("N").z);
//			System.out.println("Spherical Coordinates after Procedure:");	
//			System.out.println(rotated2SphericalCoordinates.get("CA").x+"  "+rotated2SphericalCoordinates.get("CA").y+"   "+rotated2SphericalCoordinates.get("CA").z);
//			System.out.println(rotated2SphericalCoordinates.get("C").x+"  "+rotated2SphericalCoordinates.get("C").y+"   "+rotated2SphericalCoordinates.get("C").z);
//			System.out.println(rotated2SphericalCoordinates.get("N").x+"  "+rotated2SphericalCoordinates.get("N").y+"   "+rotated2SphericalCoordinates.get("N").z);
			
		}		
		
		return atom_coord;
	}
	/**
	 * Extracts and returns the coordinates of certain atom (CA,C,N,...) from residue
	 * @param residue 
	 * @param atom type
	 * @return coordinates certain atom
	 */
	public static Point3d getRelativeCoord (Residue iRes, String atomType){
		Point3d relativeCoord = new Point3d(0,0,0);
		//TODO: compute and return the relative coordinate of the atom of type atomType
		// coord of CAatom of iRes should be (0,0,0) ...
		
		//Test if all necessary atoms CA, C and N available
		if (iRes.containsAtom("CA") && iRes.containsAtom("C") && iRes.containsAtom("N")){
			
		}
		return relativeCoord;
	}
	/**
	 * Extracts and returns the relative coordinates of certain atom (CA,C,N,...) 
	 * for residue jRes with respect to iRes
	 * @param central residue
	 * @param residue to rotate (with respect to iRes)
	 * @param atom type
	 * @return coordinates certain atom
	 */
	public static Point3d getRelativeCoord (Residue iRes, Residue jRes, String atomType){
		Point3d relativeCoord = new Point3d(0,0,0);
		//TODO: compute and return the relative coordinate of the atom of type atomType
		
		//Test if all necessary atoms CA, C and N available
		if (iRes.containsAtom("CA") && iRes.containsAtom("C") && iRes.containsAtom("N")){
			
		}
		return relativeCoord;
	}
	/**
	 * Extracts and returns the relative coordinates of all atoms (CA,C,N,...) 
	 * for residue jRes with respect to iRes
	 * @param central residue
	 * @param residue to rotate (with respect to iRes)
	 * @return coordinates all atoms of residue
	 */
	public static HashMap<String,Point3d> getRelativeCoord (Residue iRes, Residue jRes) {
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>();
		//TODO: compute and return the relative coordinate of all available atoms of the residue jRes

		//Test if all necessary atoms CA, C and N available
		if (iRes.containsAtom("CA") && iRes.containsAtom("C") && iRes.containsAtom("N")){
			
		}
		return atom_coord;
	}
	/**
	 * Extracts and returns the relative coordinates of the coordinates with respect to iRes
	 * @param central residue
	 * @param coordinates to transform with respect to iRes
	 * @return coordinates all atoms of residue
	 */
	public static HashMap<String,Point3d> getRelativeCoord (Residue iRes, Point3d coord) {
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>();
		//TODO: compute and return the relative coordinate (based on coord) with respect to iRes

		//Test if all necessary atoms CA, C and N available
		if (iRes.containsAtom("CA") && iRes.containsAtom("C") && iRes.containsAtom("N")){
			
		}
		return atom_coord;
	}
	/**
	 * Extracts and returns the relative coordinates of the coordinates (coord) with respect to iRes
	 * @param CA coordinate of the central residue
	 * @param C coordinate of the central residue
	 * @param N coordinate of the central residue
	 * @param coordinates to transform with respect to iRes
	 * @return coordinates all atoms of residue
	 */
	public static HashMap<String,Point3d> getRelativeCoord (Point3d coord_CA, Point3d coord_C, Point3d coord_N, Point3d coord) {
		HashMap<String,Point3d> atom_coord = new HashMap<String, Point3d>();
		//TODO: compute and return the relative coordinate (based on coord) with respect to iRes

		
		return atom_coord;
	}
	
	// ---------- Translation and Rotation Methods --------

	
	/*
	 *  This method gets the translated coordinates for a residue placing the selected center atom at the origin. 
	 */	
	private static HashMap<String, Point3d> getTranslatedCoordinates(Residue res, String centre){
		HashMap<String, Point3d> translatedCoordinates=new HashMap<String, Point3d>();
		for (String key: res.getAtomsMap().keySet())
			{
			Point3d tempTranslatedCoords=new Point3d();
			tempTranslatedCoords.x=res.getAtomsMap().get(key).getCoords().x-res.getAtomsMap().get(centre).getCoords().x;
			tempTranslatedCoords.y=res.getAtomsMap().get(key).getCoords().y-res.getAtomsMap().get(centre).getCoords().y;
			tempTranslatedCoords.z=res.getAtomsMap().get(key).getCoords().z-res.getAtomsMap().get(centre).getCoords().z;
			translatedCoordinates.put(key, tempTranslatedCoords);
			}
				
		return translatedCoordinates;
	}
	
	/*
	 *  This method gets the translated coordinates for a residue placing the selected center atom at the origin. 
	 */	
	@SuppressWarnings("unused")
	private static HashMap<String, Point3d> getTranslatedCoordinates(TreeMap<Integer, Residue> protein, int resno, String centre){
		HashMap<String, Point3d> translatedCoordinates=new HashMap<String, Point3d>();
		for (String key: protein.get(resno).getAtomsMap().keySet())
			{
			Point3d tempTranslatedCoords=new Point3d();
			tempTranslatedCoords.x=protein.get(resno).getAtomsMap().get(key).getCoords().x-protein.get(resno).getAtomsMap().get(centre).getCoords().x;
			tempTranslatedCoords.y=protein.get(resno).getAtomsMap().get(key).getCoords().y-protein.get(resno).getAtomsMap().get(centre).getCoords().y;
			tempTranslatedCoords.z=protein.get(resno).getAtomsMap().get(key).getCoords().z-protein.get(resno).getAtomsMap().get(centre).getCoords().z;
			translatedCoordinates.put(key, tempTranslatedCoords);
			}
				
		return translatedCoordinates;
	}
	
	/*
	 * This method converts the translated Cartesian coordinates to spherical coordinates. 
	 */
	
	private static HashMap<String, Point3d> getTranslatedSphericalCoordinates (HashMap<String, Point3d> translatedCoordinates)	{
		HashMap<String, Point3d> translatedSphericalCoordinates=new HashMap<String, Point3d>();
		for (String key: translatedCoordinates.keySet())
			{
			Point3d tempTranslatedSphericalCoordinates=new Point3d();
			tempTranslatedSphericalCoordinates.x=Math.sqrt(Math.pow(translatedCoordinates.get(key).x, 2)+Math.pow(translatedCoordinates.get(key).y, 2)+Math.pow(translatedCoordinates.get(key).z, 2));
			tempTranslatedSphericalCoordinates.y=Math.acos(translatedCoordinates.get(key).z/tempTranslatedSphericalCoordinates.x);
			tempTranslatedSphericalCoordinates.z=Math.atan2(translatedCoordinates.get(key).y, translatedCoordinates.get(key).x);
			translatedSphericalCoordinates.put(key, tempTranslatedSphericalCoordinates);
			}
		return translatedSphericalCoordinates;
	}
		
	/*
	 * This method will align the C-alpha - C-zero vector along the X axis.
	 */
	private static HashMap<String, Point3d> getRotated1SphericalCoords (HashMap<String, Point3d> translatedSphericalCoordinates, String right)
	{
	HashMap<String, Point3d> rotated1SphericalCoords=new HashMap<String, Point3d>();
	for (String key: translatedSphericalCoordinates.keySet())
		{
		Point3d tempRotated1SphericalCoords=new Point3d();
		tempRotated1SphericalCoords.x=translatedSphericalCoordinates.get(key).x;
		tempRotated1SphericalCoords.y=translatedSphericalCoordinates.get(key).y+(Math.PI/2-translatedSphericalCoordinates.get(right).y);
		tempRotated1SphericalCoords.z=translatedSphericalCoordinates.get(key).z-translatedSphericalCoordinates.get(right).z;
		rotated1SphericalCoords.put(key, tempRotated1SphericalCoords);
		}
	return rotated1SphericalCoords;
	}
		
	/*
	 *   This method just converts the spherical coordinates after the first rotation to the Cartesian ones.
	 *   This is necessary because it is much easier to rotate a Cartesian coordinates about the x axis. 
	 */
	private static HashMap<String, Point3d> getRotated1CartesianCoords (HashMap<String, Point3d> rotated1SphericalCoordinates)
	{
	HashMap<String, Point3d> rotated1CartesianCoordinates=new HashMap<String, Point3d>();
	for (String key: rotated1SphericalCoordinates.keySet())
		{
		Point3d tempRotated1CartesianCoordinates=new Point3d();
		tempRotated1CartesianCoordinates.x=rotated1SphericalCoordinates.get(key).x*Math.sin(rotated1SphericalCoordinates.get(key).y)*Math.cos(rotated1SphericalCoordinates.get(key).z);
		tempRotated1CartesianCoordinates.y=rotated1SphericalCoordinates.get(key).x*Math.sin(rotated1SphericalCoordinates.get(key).y)*Math.sin(rotated1SphericalCoordinates.get(key).z);
		tempRotated1CartesianCoordinates.z=rotated1SphericalCoordinates.get(key).x*Math.cos(rotated1SphericalCoordinates.get(key).y);
		rotated1CartesianCoordinates.put(key, tempRotated1CartesianCoordinates);
		}
	return rotated1CartesianCoordinates;
	}	
	
	/*
	 *   This method performs the second rotation i.e. puts the left atom in the X-Z plane. 
	 *   The output of these are the final translated-rotated Cartesian coordinates. 
	 *   omega is the angle to be used in the rotation matrix, so as to rotate the atoms about the x axis
	 */
	private static HashMap<String, Point3d> getRotated2CartesianCoords (HashMap<String, Point3d> rotated1CartesianCoordinates, String left)
	{
	HashMap<String, Point3d> rotated2CartesianCoordinates=new HashMap<String, Point3d>();
	double omega=0;
//	omega=Math.atan(-rotated1CartesianCoordinates.get(left).z/ rotated1CartesianCoordinates.get(left).y);
	omega=Math.atan2(-rotated1CartesianCoordinates.get(left).z, rotated1CartesianCoordinates.get(left).y);
	
	for (String key: rotated1CartesianCoordinates.keySet())
		{
		Point3d tempRotated2CartesianCoordinates=new Point3d();
		tempRotated2CartesianCoordinates.x=rotated1CartesianCoordinates.get(key).x;
		tempRotated2CartesianCoordinates.y=rotated1CartesianCoordinates.get(key).y*Math.cos(omega)-rotated1CartesianCoordinates.get(key).z*Math.sin(omega);
		tempRotated2CartesianCoordinates.z=rotated1CartesianCoordinates.get(key).y*Math.sin(omega)+rotated1CartesianCoordinates.get(key).z*Math.cos(omega);
		rotated2CartesianCoordinates.put(key, tempRotated2CartesianCoordinates);
		}
	return rotated2CartesianCoordinates;
	}	
	
	/*
	 * This method just converts the above final Cartesian coordinates to the final spherical coordinates.
	 */	
	private static HashMap<String, Point3d> getRotated2SphericalCoords (HashMap<String, Point3d> rotated2CartesianCoordinates)
	{
	HashMap<String, Point3d> rotated2SphericalCoordinates=new HashMap<String, Point3d>();
	for (String key: rotated2CartesianCoordinates.keySet())
		{
		Point3d tempRotated2SphericalCoordinates=new Point3d();
		tempRotated2SphericalCoordinates.x=Math.sqrt(Math.pow(rotated2CartesianCoordinates.get(key).x, 2)+Math.pow(rotated2CartesianCoordinates.get(key).y, 2)+Math.pow(rotated2CartesianCoordinates.get(key).z, 2));
		tempRotated2SphericalCoordinates.y=Math.acos(rotated2CartesianCoordinates.get(key).z/tempRotated2SphericalCoordinates.x);
		tempRotated2SphericalCoordinates.z=Math.atan2(rotated2CartesianCoordinates.get(key).y, rotated2CartesianCoordinates.get(key).x);
		rotated2SphericalCoordinates.put(key, tempRotated2SphericalCoordinates);
		}
	return rotated2SphericalCoordinates;
	}
	
	

}

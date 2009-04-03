import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
//import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import proteinstructure.Polymer;
import proteinstructure.Residue;

/**
 * Class to read a Histone molecule from a PDB file and output the coordinates 
 * of ARGs and LYSs that lie on the lateral surface of the histone using cylindrical 
 * coordinates referring to the its cylindrical symmetry axis (which is calculated as the 
 * inertia axis corresponding to the maximum principal moment of inertia) 
 * @author duarte
 *
 */
public class analyseHistone {

	//private static final double UPPER_CUTOFF = 40;
	//private static final double LOWER_CUTOFF = 25;

	
	private static Set<String> getResiSet(String chainCode, int beg1, int end1, int beg2, int end2) {
		TreeSet<String> set = new TreeSet<String>();
		for (int i=beg1;i<=end1;i++){
			set.add(chainCode+i);
		}
		for (int i=beg2;i<=end2;i++){
			set.add(chainCode+i);
		}

		return set;
	}
	
	
	public static void main(String[] args) throws Exception {
		File pdbFile = new File("/project/StruPPi/jose/histone/pdb1kx5.ent");
		String[] chains = {"A", "B", "C", "D", "E", "F", "G", "H","I", "J"};
		
		// residues defining the core of the histone (excluding tails) 
		Set<String> resisetA = getResiSet("A",1,37,133,135);
		Set<String> resisetE = getResiSet("E",1,37,133,135);
		
		Set<String> resisetB = getResiSet("B",1,26,95,102);
		Set<String> resisetF = getResiSet("F",1,26,95,102);
		
		Set<String> resisetC = getResiSet("C",1,15,118,130);
		Set<String> resisetG = getResiSet("G",1,15,118,130);
		
		Set<String> resisetD = getResiSet("D",1,29,123,125);
		Set<String> resisetH = getResiSet("H",1,29,123,125);
		
		// reference residue to measure the theta angles
		// ARG G42 is approx at the center of the thickest side of the histone
		String refResChainCode = "G";
		int refResSerial = 42; 
		
		// output files
		File lysFile = new File("/project/StruPPi/jose/histone/lys_z_theta.txt");
		File argFile = new File("/project/StruPPi/jose/histone/arg_z_theta.txt");
		File dnaFile = new File("/project/StruPPi/jose/histone/dna_z_theta.txt");
		
		// reading histone pdb file
		Polymer mol = new Polymer();
		mol.readFromPdbFile(pdbFile, chains);
		
		// removing tails
		mol.removeResidues(resisetA);
		mol.removeResidues(resisetB);
		mol.removeResidues(resisetC);
		mol.removeResidues(resisetD);
		mol.removeResidues(resisetE);
		mol.removeResidues(resisetF);
		mol.removeResidues(resisetG);
		mol.removeResidues(resisetH);
		
		
		
		Point3d centerOfMass = mol.getCenterOfMass();
		System.out.printf("Center of mass: (%9.2f, %9.2f, %9.2f)\n",centerOfMass.x,centerOfMass.y,centerOfMass.z);
		
		double diam = mol.getDiameter();
		System.out.println("Diameter: "+diam);
		Vector3d axis = mol.getBiggestInertiaAxis(mol.getCenterOfMass());
		axis.normalize(); // we normalize it so that it is easy to handle it later
		System.out.printf("Biggest inertia moment axis: (%5.3f, %5.3f, %5.3f)\n",axis.x,axis.y,axis.z);
		
		// changing reference frame to origin center of mass and biggest inertial axis lying on the z axis
		mol.transformRefFrameToCenterAndAxis(centerOfMass, axis);
		
		// rotating so that x axis is aligned along the reference residue
		Residue refRes = mol.getResidue(refResSerial, refResChainCode);
		Point3d refResCoords = refRes.getCoords();
		Vector3d iVec = new Vector3d(1,0,0);
		double rotAngle = iVec.angle(new Vector3d(refResCoords.x,refResCoords.y, 0));
		mol.rotate(new Vector3d(0,0,1), rotAngle);
		// sanity check: coords of ref residue should have now y=0
		//System.out.println("coordinates of reference residue G42 after rotating: "+refResCoords);
		
		// opening files for output
		PrintStream outLys = new PrintStream(new FileOutputStream(lysFile));
		PrintStream outArg = new PrintStream(new FileOutputStream(argFile));
		PrintStream outDNA = new PrintStream(new FileOutputStream(dnaFile));
		
		System.out.println("Getting LYSs");
		ArrayList<Residue> lysResidues = mol.getResiduesOfType("LYS");
		projectOnCylinder(lysResidues, outLys);
		System.out.println("Getting ARGs");
		ArrayList<Residue> argResidues = mol.getResiduesOfType("ARG");
		projectOnCylinder(argResidues, outArg);
		System.out.println("Getting DNA residues");
		ArrayList<Residue> aResidues = mol.getResiduesOfType("DA");
		ArrayList<Residue> cResidues = mol.getResiduesOfType("DC");
		ArrayList<Residue> tResidues = mol.getResiduesOfType("DT");
		ArrayList<Residue> gResidues = mol.getResiduesOfType("DG");
		ArrayList<Residue> dnaResidues = new ArrayList<Residue>();
		dnaResidues.addAll(aResidues);
		dnaResidues.addAll(cResidues);
		dnaResidues.addAll(tResidues);
		dnaResidues.addAll(gResidues);
		projectOnCylinder(dnaResidues, outDNA);
		
		
		outLys.close();
		outArg.close();
		outDNA.close();
		
		// sanity checks: center should be now at (0,0,0), biggest inertia axis aligned to z axis
		System.out.println("Sanity check, center of mass at (0,0,0) and biggest principal inertia axis on z-axis:");
		Point3d endCoM = mol.getCenterOfMass();
		Vector3d endInAxis = mol.getBiggestInertiaAxis(endCoM);
		System.out.printf("Center of mass: (%9.2f, %9.2f, %9.2f)\n",endCoM.x,endCoM.y,endCoM.z);
		System.out.printf("Biggest inertia moment axis: (%5.3f, %5.3f, %5.3f)\n",endInAxis.x,endInAxis.y,endInAxis.z);
	}
	
	private static void projectOnCylinder(ArrayList<Residue> residues, PrintStream out) {
		
		ArrayList<Double> zs = new ArrayList<Double>();
		ArrayList<Double> thetas = new ArrayList<Double>();
		ArrayList<Double> dists = new ArrayList<Double>();
		ArrayList<String> resSerials = new ArrayList<String>();
		// we take the x axis as reference to measure the theta angle
		Vector3d iVec = new Vector3d(1,0,0); 
		
		for (Residue residue:residues) {
			Point3d coords = residue.getCoords();
			// getting distance to the z axis
			double dist = new Point3d(coords.x,coords.y,0).distance(new Point3d(0,0,0));			

			// finding the theta angle with respect to the x axis that we take as reference
			double theta = iVec.angle(new Vector3d(coords.x, coords.y, 0)); // angle between x axis and xy plane projection of the residue coords
			// angle function is restricted to [0,PI], we have to correct for the vectors on the y<0 side of the plane
			if (coords.y<0) {
				theta = -theta; // the output will be angles in [-PI,PI]
			}
			// adding this z coordinate to the list of z coordinates
			zs.add(coords.z);
			// adding the theta angle to the thetas list
			thetas.add(theta);
			// adding the residue serial to the residue serial list
			resSerials.add(residue.getChainCode()+residue.getSerial());
			// adding the distance to the dists list
			dists.add(dist);

		}
		
		// printing out values found
		out.println("res\tdist\tz\ttheta");
		System.out.println("res\tdist\tz\ttheta");
		for (int i=0;i<zs.size();i++) {
			out.printf("%4s\t%7.3f\t%7.3f\t%7.3f\n", resSerials.get(i), dists.get(i), zs.get(i), thetas.get(i));
			System.out.printf("%4s\t%7.3f\t%7.3f\t%7.3f\n", resSerials.get(i), dists.get(i), zs.get(i), thetas.get(i));
		}
	}

}

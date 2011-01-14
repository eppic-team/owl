import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
//import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;

import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import owl.core.structure.AminoAcid;
import owl.core.structure.Pdb;
import owl.core.structure.PdbasePdb;
import owl.core.structure.PolResidue;
import owl.core.structure.Polymer;
import owl.core.structure.Residue;
import owl.core.util.Interval;
import owl.core.util.IntervalSet;

//import proteinstructure.PdbfilePdb;

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

	// the histone PDB file
	private static final String pdbCode = "1kx5";
	private static final File pdbFile = new File("/project/StruPPi/jose/histone/pdb1kx5.ent");
	private static final String[] chains = {"A", "B", "C", "D", "E", "F", "G", "H","I", "J"};
	
	// output files
	private static final File lysFile = new File("/project/StruPPi/jose/histone/lys_z_theta.txt");
	private static final File argFile = new File("/project/StruPPi/jose/histone/arg_z_theta.txt");
	private static final File dnaFile = new File("/project/StruPPi/jose/histone/dna_z_theta.txt");
	private static final File checkPdbFile = new File("/project/StruPPi/jose/histone/check.pdb");

	// side chain atoms used for LYSs and ARGs
	private static String LYS_SIDECHAIN_ATOM = "NZ";
	private static String ARG_SIDECHAIN_ATOM = "NH1";
	
	// residues defining the tails of the histone
	// H3 subunit, chains A,E. Core residues: 38-132
	private static final Set<String> resisetA = getResiSet("A",1,37,133,135);
	private static final Set<String> resisetE = getResiSet("E",1,37,133,135);
	
	// H4 subunit, chains B,F. Core residues: 27-94
	private static final Set<String> resisetB = getResiSet("B",1,26,95,102);
	private static final Set<String> resisetF = getResiSet("F",1,26,95,102);
	
	// H2A subunit, chains C,G. Core residues: 16-117
	private static final Set<String> resisetC = getResiSet("C",1,15,118,130);
	private static final Set<String> resisetG = getResiSet("G",1,15,118,130);
	
	// H2B subunit, chains D,H. Core residues: 30-122
	private static final Set<String> resisetD = getResiSet("D",1,29,123,125);
	private static final Set<String> resisetH = getResiSet("H",1,29,123,125);

	// the core residues as interval sets (for the Pdb objects)
	private static final IntervalSet coreH3 = getIntervalSet(38,132);
	private static final IntervalSet coreH4 = getIntervalSet(27,94);
	private static final IntervalSet coreH2A = getIntervalSet(16,117);
	private static final IntervalSet coreH2B = getIntervalSet(30,122);

	// reference residue to measure the theta angles
	// ARG G42 is approx at the center of the thickest side of the histone
	private static final String refResChainCode = "G";
	private static final int refResSerial = 42; 

	
	private static IntervalSet getIntervalSet(int beg, int end) {
		IntervalSet intervSet = new IntervalSet();
		intervSet.add(new Interval(beg,end));
		return intervSet;
	}
	
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
		
		// Reading the pdb file again with our normal Pdb objects (one per chain, no DNA)
		// We do this to extract the positions of the side chains of ARG and LYS after transforming the molecules in the same way as mol
		Pdb completeChainA = new PdbasePdb(pdbCode);
		completeChainA.load("A");
		Pdb completeChainB = new PdbasePdb(pdbCode);
		completeChainB.load("B");
		Pdb completeChainC = new PdbasePdb(pdbCode);
		completeChainC.load("C");
		Pdb completeChainD = new PdbasePdb(pdbCode);
		completeChainD.load("D");
		Pdb completeChainE = new PdbasePdb(pdbCode);
		completeChainE.load("E");
		Pdb completeChainF = new PdbasePdb(pdbCode);
		completeChainF.load("F");
		Pdb completeChainG = new PdbasePdb(pdbCode);
		completeChainG.load("G");
		Pdb completeChainH = new PdbasePdb(pdbCode);
		completeChainH.load("H");
		Pdb[] pdbs = {completeChainA, completeChainB, completeChainC, completeChainD, completeChainE, completeChainF, completeChainG, completeChainH};
		
		// removing tails
		completeChainA.restrictToIntervalSet(coreH3);
		completeChainE.restrictToIntervalSet(coreH3);
		completeChainB.restrictToIntervalSet(coreH4);
		completeChainF.restrictToIntervalSet(coreH4);
		completeChainC.restrictToIntervalSet(coreH2A);
		completeChainG.restrictToIntervalSet(coreH2A);
		completeChainD.restrictToIntervalSet(coreH2B);
		completeChainH.restrictToIntervalSet(coreH2B);
		
		
		Point3d centerOfMass = mol.getCenterOfMass();
		System.out.printf("Center of mass: (%9.2f, %9.2f, %9.2f)\n",centerOfMass.x,centerOfMass.y,centerOfMass.z);
		
		double diam = mol.getDiameter();
		System.out.println("Diameter: "+diam);
		Vector3d axis = mol.getBiggestInertiaAxis(mol.getCenterOfMass());
		axis.normalize(); // we normalize it so that it is easy to handle it later
		System.out.printf("Biggest inertia moment axis: (%5.3f, %5.3f, %5.3f)\n",axis.x,axis.y,axis.z);
		
		// changing reference frame to origin center of mass and biggest inertial axis lying on the z axis
		mol.transformRefFrameToCenterAndAxis(centerOfMass, axis);
		// doing the same to the individual Pdb complete chains
		for (Pdb pdb:pdbs) {
			pdb.transformToCenterAndAxis(centerOfMass, axis);
		}
		
		// rotating so that x axis is aligned along the reference residue
		PolResidue refRes = mol.getResidue(refResSerial, refResChainCode);
		Point3d refResCoords = refRes.getCoords();
		Vector3d iVec = new Vector3d(1,0,0);
		double rotAngle = iVec.angle(new Vector3d(refResCoords.x,refResCoords.y, 0));
		Vector3d zaxis = new Vector3d(0,0,1);
		mol.rotate(zaxis, rotAngle);
		// sanity check: coords of ref residue should have now y=0
		//System.out.println("coordinates of reference residue G42 after rotating: "+refResCoords);
		
		// we do the same rotation on the individual Pdb complete chains
		for (Pdb pdb:pdbs) {
			pdb.rotate(zaxis, rotAngle);
		}
		
		
		
		// opening files for output
		PrintStream outLys = new PrintStream(new FileOutputStream(lysFile));
		PrintStream outArg = new PrintStream(new FileOutputStream(argFile));
		PrintStream outDNA = new PrintStream(new FileOutputStream(dnaFile));
		
		System.out.println("Getting LYSs");
		//ArrayList<PolResidue> lysResidues = mol.getResiduesOfType("LYS");
		//projectOnCylinder(lysResidues, outLys);
		outLys.println("res\tdist\tz\ttheta");
		System.out.println("res\tdist\tz\ttheta");
		for (Pdb pdb:pdbs) {
			ArrayList<Residue> lysResidues = pdb.getResiduesOfType(AminoAcid.LYS);
			projectOnCylinderPdb(lysResidues, outLys, AminoAcid.LYS);
		}
		 
		
		System.out.println("Getting ARGs");
		//ArrayList<PolResidue> argResidues = mol.getResiduesOfType("ARG");
		//projectOnCylinder(argResidues, outArg);
		outArg.println("res\tdist\tz\ttheta");
		System.out.println("res\tdist\tz\ttheta");
		for (Pdb pdb:pdbs) {
			ArrayList<Residue> argResidues = pdb.getResiduesOfType(AminoAcid.ARG);
			projectOnCylinderPdb(argResidues, outArg, AminoAcid.ARG);
		}

		
		System.out.println("Getting DNA residues");
		ArrayList<PolResidue> aResidues = mol.getResiduesOfType("DA");
		ArrayList<PolResidue> cResidues = mol.getResiduesOfType("DC");
		ArrayList<PolResidue> tResidues = mol.getResiduesOfType("DT");
		ArrayList<PolResidue> gResidues = mol.getResiduesOfType("DG");
		ArrayList<PolResidue> dnaResidues = new ArrayList<PolResidue>();
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
		
		// sanity check: writing the pdbs to file to examine them
		PrintStream ps = new PrintStream(checkPdbFile);
		for (Pdb pdb:pdbs) {
			pdb.writeAtomLines(ps);
		}
		ps.close();
	}
	
	private static void projectOnCylinder(ArrayList<PolResidue> residues, PrintStream out) {
		
		ArrayList<Double> zs = new ArrayList<Double>();
		ArrayList<Double> thetas = new ArrayList<Double>();
		ArrayList<Double> dists = new ArrayList<Double>();
		ArrayList<String> resSerials = new ArrayList<String>();
		// we take the x axis as reference to measure the theta angle
		Vector3d iVec = new Vector3d(1,0,0); 
		
		for (PolResidue residue:residues) {
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
			resSerials.add(residue.getPdbChainCode()+residue.getSerial());
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

	private static void projectOnCylinderPdb(ArrayList<Residue> residues, PrintStream out, AminoAcid aa) {
		
		ArrayList<Double> zs = new ArrayList<Double>();
		ArrayList<Double> thetas = new ArrayList<Double>();
		ArrayList<Double> dists = new ArrayList<Double>();
		ArrayList<String> resSerials = new ArrayList<String>();
		// we take the x axis as reference to measure the theta angle
		Vector3d iVec = new Vector3d(1,0,0); 
		
		for (Residue residue:residues) {
			String atom = null;
			if (aa.equals(AminoAcid.LYS)) atom = LYS_SIDECHAIN_ATOM;
			if (aa.equals(AminoAcid.ARG)) atom = ARG_SIDECHAIN_ATOM;
			Point3d coords = residue.getAtom(atom).getCoords();
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
			resSerials.add(residue.getPdbChainCode()+residue.getSerial());
			// adding the distance to the dists list
			dists.add(dist);

		}
		
		// printing out values found
		for (int i=0;i<zs.size();i++) {
			out.printf("%4s\t%7.3f\t%7.3f\t%7.3f\n", resSerials.get(i), dists.get(i), zs.get(i), thetas.get(i));
			System.out.printf("%4s\t%7.3f\t%7.3f\t%7.3f\n", resSerials.get(i), dists.get(i), zs.get(i), thetas.get(i));
		}
	}

}

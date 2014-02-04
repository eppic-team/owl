package owl.core.structure;

import java.io.File;
import java.util.ArrayList;

import javax.vecmath.Point3d;


/**
 * Class to calculate Accessible Surface Areas based on
 * the rolling ball algorithm by Shrake and Rupley.
 * 
 * The code is taken from a python implementation at http://boscoh.com/protein/asapy
 * and adapted to java.
 * Thanks again to Bosco K. Ho for another great piece of code (our rmsd calculation is also his)
 * 
 * See 
 * Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent of Protein Atoms. 
 * Lysozyme and Insulin." JMB (1973) 79:351-371.
 * Lee, B., and Richards, F.M. "The interpretation of Protein Structures: Estimation of
 * Static Accessibility" JMB (1971) 55:379-400
 * @author duarte_j
 *
 */
public class AsaCalculator {

	// Bosco uses as default 960, Shrake and Rupley seem to use in their paper 92 (not sure if this is actually the same parameter)
	public static final int DEFAULT_N_SPHERE_POINTS = 960;
	public static final double DEFAULT_PROBE_SIZE = 1.4;
	public static final int DEFAULT_NTHREADS = 1;
	
	private class GroupASACalcThread extends Thread {
		
		int start;
		int end;
		double[] asas;
		
		public GroupASACalcThread(int start, int end, double[] asas) {
			this.start = start;
			this.end = end;
			this.asas = asas;
		}

		public void run() {
			calcGroupOfAsas(start, end, asas);
		}
	}
	
	
	private Atom[] atoms;
	private double[] radii;
	private double probe;
	private int nThreads;
	private Point3d[] spherePoints;
	private double cons;
	
	/**
	 * Constructs a new Asa
	 * @param atoms
	 * @param probe the probe size
	 * @param nSpherePoints the number of points to be used in generating the spherical 
	 * dot-density, the more points the more accurate (and slower) calculation
	 * @param nThreads the number of parallel threads to use for the calculation
	 */
	public AsaCalculator(Atom[] atoms, double probe, int nSpherePoints, int nThreads) {
		this.atoms = atoms;
		this.probe = probe;
		this.nThreads = nThreads;
		
		// initialising the radii by looking them up through AtomRadii
		radii = new double[atoms.length];
		for (int i=0;i<atoms.length;i++) {
			radii[i] = AtomRadii.getRadius(atoms[i]);
		}
		
		// initialising the sphere points to sample
		spherePoints = generateSpherePoints(nSpherePoints);
		
		cons = 4.0 * Math.PI / (double)nSpherePoints;
	}
	
	/**
	 * Calculates the Accessible Surface Areas of the given atoms, using given probe size.
	 
	 * @return an array with asa values matching the input atoms array
	 */
	public double[] calculateAsa() {
		
		double[] asas = new double[atoms.length];

	    if (nThreads==1) {
		    for (int i=0;i<atoms.length;i++) {	    	
		        asas[i] = calcSingleAsa(i); 
		    }
		    
	    } else {
	    	// NOTE the multithreaded calculation does not scale up well (4 CPUs ~ x2.8, 8CPUs ~ x2.9)
	    	// tried copying the arrays (atoms and sphere_points) as new arrays but the scaling behaves the same
	    	// also tried dividing the asas array in parts and them joining all together but same scaling behaviour
	    	// I guess it's simply memory bottlenecks of the architecture :(
	    	GroupASACalcThread[] threads = new GroupASACalcThread[nThreads];
	    	
	    	int[] startIndices = getStartingIdxForGroups();


		    for (int k=0;k<nThreads;k++) {
		    	threads[k] = new GroupASACalcThread(startIndices[k], startIndices[k+1], asas);
		    	threads[k].start();
		    }
	    	
		    for (int k=0;k<nThreads;k++) {
		    	try {
		    		threads[k].join();
		    	} catch (InterruptedException e) {
		    		System.err.println("Unexpected error while running multi-threaded ASA calculation. Exiting.");
		    		e.printStackTrace();
		    		System.exit(1);
		    	}
		    }
	    }
	    
	    return asas;
	}
	
	/**
	 * Returns list of 3d coordinates of points on a sphere using the
	 * Golden Section Spiral algorithm.
	 * @param nSpherePoints the number of points to be used in generating the spherical dot-density
	 * @return
	 */
	private Point3d[] generateSpherePoints(int nSpherePoints) {
	    Point3d[] points = new Point3d[nSpherePoints];
	    double inc = Math.PI * (3.0 - Math.sqrt(5.0));
	    double offset = 2.0 / (double)nSpherePoints; 
	    for (int k=0;k<nSpherePoints;k++) {
	        double y = k * offset - 1.0 + (offset / 2.0);
	        double r = Math.sqrt(1.0 - y*y);
	        double phi = k * inc;
	        points[k] = new Point3d(Math.cos(phi)*r, y, Math.sin(phi)*r);
	    }
	    return points;
	}

	/**
	 * Returns list of indices of atoms within probe distance to atom k.
	 * @param k index of atom for which we want neighbor indices
	 */
	private ArrayList<Integer> findNeighborIndices(int k) {
		// looking at a typical protein case, number of neighbours are from ~10 to ~50, with an average of ~30
		// Thus 40 seems to be a good compromise for the starting capacity
	    ArrayList<Integer> neighbor_indices = new ArrayList<Integer>(40);
	    
	    double radius = radii[k] + probe + probe;
	    
	    for (int i=0;i<atoms.length;i++) {
	    	if (i==k) continue;
	    	
	        double dist = atoms[i].getCoords().distance(atoms[k].getCoords());
	    	
	        if (dist < radius + radii[i]) {
	            neighbor_indices.add(i);
	        }
	    }
	    
	    return neighbor_indices;
	}
	
	private void calcGroupOfAsas(int startIdx, int endIdx, double[] asas) {
		for (int i=startIdx;i<endIdx;i++) {
			asas[i] = calcSingleAsa(i);
		}
	}

	private double calcSingleAsa(int i) {
    	Atom atom_i = atoms[i];
    	ArrayList<Integer> neighbor_indices = findNeighborIndices(i);
        int n_neighbor = neighbor_indices.size();
        int j_closest_neighbor = 0;
        double radius = probe + radii[i];

        int n_accessible_point = 0;
        
        for (Point3d point: spherePoints){
            boolean is_accessible = true;
            Point3d test_point = new Point3d(point.x*radius + atom_i.getCoords().x,
            								point.y*radius + atom_i.getCoords().y,
            								point.z*radius + atom_i.getCoords().z);

            int[] cycled_indices = new int[n_neighbor];
            int arind = 0;
            for (int ind=j_closest_neighbor;ind<n_neighbor;ind++) {
            	cycled_indices[arind] = ind;
            	arind++;
            }
            for (int ind=0;ind<j_closest_neighbor;ind++){
            	cycled_indices[arind] = ind;
            	arind++;
            }

            for (int j: cycled_indices) {
                Atom atom_j = atoms[neighbor_indices.get(j)];
                double r = radii[neighbor_indices.get(j)] + probe;
                double diff_sq = atom_j.getCoords().distanceSquared(test_point);
                if (diff_sq < r*r) {
                    j_closest_neighbor = j;
                    is_accessible = false;
                    break;
                }
            }
            if (is_accessible) {
                n_accessible_point++;
            }
        }
        return cons*n_accessible_point*radius*radius;
	}
	
	private int[] getStartingIdxForGroups() {
		// we use so many groups as threads
		int n = atoms.length;
		int[] indices = new int[nThreads+1];

		int baseSize = n/nThreads;
		int remainder = n%nThreads;

		indices[0] = 0;
		for (int k=1;k<nThreads;k++){
			indices[k] = indices[k-1]+baseSize;
			if (k<remainder) {
				indices[k]+=1;
			}
		}
		indices[nThreads] = n;
		return indices;
	}

	/**
	 * To test the class
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		PdbAsymUnit pdb = new PdbAsymUnit(new File(args[0]));
		int nThreads = Integer.parseInt(args[1]);
		
		pdb.removeHatoms();
		
		long start = System.currentTimeMillis();
		pdb.calcASAs(1000,nThreads,false);
		long end = System.currentTimeMillis();

		
		double tot = 0;
		
		for (PdbChain chain:pdb.getAllChains()) {
			for (Residue res:chain) {

				System.out.printf("%3d\t%s\t%6.2f\n",res.getSerial(),res.getLongCode(),res.getAsa());
				tot+=res.getAsa();
			}
		}
		System.out.printf("Total area: %9.2f\n",tot);
		System.out.printf("Time: %4.1fs\n",((end-start)/1000.0));
		

	}
}

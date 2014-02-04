package owl.core.structure;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.vecmath.Point3d;



/**
 * Class to calculate Accessible Surface Areas based on
 * the rolling ball algorithm by Shrake and Rupley.
 * 
 * The code is taken from a python implementation at http://boscoh.com/protein/asapy
 * and adapted to java. Now source is available at https://github.com/boscoh/asa
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
	
	
	private class AsaCalcWorker implements Runnable {

		private int i;
		private double[] asas;
		
		public AsaCalcWorker(int i, double[] asas) {
			this.i = i;
			this.asas = asas;
		}

		@Override
		public void run() {
			asas[i] = calcSingleAsa(i);
		}
	}
	
	
	private Atom[] atoms;
	private double[] radii;
	private double probe;
	private int nThreads;
	private Point3d[] spherePoints;
	private double cons;
	
	/**
	 * Constructs a new AsaCalculator. Subsequently call {@link #calculateAsa()}
	 * to calculate the ASAs.
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
	 * Calculates the Accessible Surface Areas for the atoms given in constructor and with parameters given.
	 * Beware that the parallel implementation is quite memory hungry. It scales well as long as there is
	 * enough memory available. 
	 * @return an array with asa values corresponding to each atom of the input array
	 */
	public double[] calculateAsa() {
		
		double[] asas = new double[atoms.length];

	    if (nThreads<=1) { // (i.e. it will also be 1 thread if 0 or negative number specified)
		    for (int i=0;i<atoms.length;i++) {	    	
		        asas[i] = calcSingleAsa(i); 
		    }
		    
	    } else {
	    	// NOTE the multithreaded calculation does not scale up well in some systems, 
	    	// why? I guess some memory/garbage collect problem? I tried increasing Xmx in pc8201 but didn't help 
	    	
	    	// Following scaling tests are for 3hbx, calculating ASA of full asym unit (6 chains):
	    	
	    	// SCALING test done in merlinl01 (12 cores, Xeon X5670  @ 2.93GHz, 24GB RAM)   	 
	    	//1 threads, time:  8.8s -- x1.0
	    	//2 threads, time:  4.4s -- x2.0
	    	//3 threads, time:  2.9s -- x3.0
	    	//4 threads, time:  2.2s -- x3.9
	    	//5 threads, time:  1.8s -- x4.9
	    	//6 threads, time:  1.6s -- x5.5
	    	//7 threads, time:  1.4s -- x6.5
	    	//8 threads, time:  1.3s -- x6.9

	    	// SCALING test done in pc8201 (4 cores, Core2 Quad Q9550  @ 2.83GHz, 8GB RAM)
	    	//1 threads, time: 17.2s -- x1.0
	    	//2 threads, time:  9.7s -- x1.8
	    	//3 threads, time:  7.7s -- x2.2
	    	//4 threads, time:  7.9s -- x2.2

	    	// SCALING test done in eppic01 (16 cores, Xeon E5-2650 0  @ 2.00GHz, 128GB RAM)
	    	//1 threads, time: 10.7s -- x1.0
	    	//2 threads, time:  5.6s -- x1.9
	    	//3 threads, time:  3.6s -- x3.0
	    	//4 threads, time:  2.8s -- x3.9
	    	//5 threads, time:  2.3s -- x4.8
	    	//6 threads, time:  1.8s -- x6.0
	    	//7 threads, time:  1.6s -- x6.8
	    	//8 threads, time:  1.3s -- x8.0
	    	//9 threads, time:  1.3s -- x8.5
	    	//10 threads, time:  1.1s -- x10.0
	    	//11 threads, time:  1.0s -- x10.9
	    	//12 threads, time:  0.9s -- x11.4


	    	
	    	ExecutorService threadPool = Executors.newFixedThreadPool(nThreads);

	    	
	    	for (int i=0;i<atoms.length;i++) {
	    		threadPool.submit(new AsaCalcWorker(i,asas));    			
	    	}

	    	threadPool.shutdown();
	    	
	    	while (!threadPool.isTerminated());
	    	
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
	
	/**
	 * To test the class
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {
		PdbAsymUnit pdb = new PdbAsymUnit(new File(args[0]));
		int maxNThreads = Integer.parseInt(args[1]);
		
		pdb.removeHatoms();
		
		long start = System.currentTimeMillis();
		pdb.calcASAs(1000,maxNThreads,false);
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
		
		
		System.out.println("Testing scaling: ");
		double[] runTimes = new double[maxNThreads];
		for (int nThreads=1;nThreads<=maxNThreads;nThreads++) {
			start = System.currentTimeMillis();
			pdb.calcASAs(1000,nThreads,false);
			end = System.currentTimeMillis();
			runTimes[nThreads-1] = (end-start)/1000.0;
			
		}
		for (int nThreads=1;nThreads<=maxNThreads;nThreads++) {
			System.out.printf(nThreads+" threads, time: %4.1fs -- x%2.1f\n",runTimes[nThreads-1],runTimes[0]/runTimes[nThreads-1]);
		}
		
	}
}

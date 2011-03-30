package owl.core.util;

import java.io.Serializable;
import java.util.Arrays;

import owl.core.structure.Atom;

/**
 * A bounding box for short cutting some geometrical calculations.
 * 
 * See http://en.wikipedia.org/wiki/Bounding_volume
 * 
 * @author duarte_j
 *
 */
public class BoundingBox implements Serializable {

	private static final long serialVersionUID = 1L;
	
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
	
	public BoundingBox(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
		this.xmin = xmin;
		this.xmax = xmax;
		this.ymin = ymin;
		this.ymax = ymax;
		this.zmin = zmin;
		this.zmax = zmax;
	}

	/**
	 * Constructs a BoundingBox by calculating maxs and mins of given array of atoms.
	 * @param atoms 
	 */
	public BoundingBox (Atom[] atoms) {
		double[] xs = new double[atoms.length];
		double[] ys = new double[atoms.length];
		double[] zs = new double[atoms.length];
		int c = 0;
		for (Atom atom:atoms) {
			xs[c] = atom.getCoords().x;
			ys[c] = atom.getCoords().y;
			zs[c] = atom.getCoords().z;
			c++;
		}

		double[] xminmax = getMinMax(xs);
		double[] yminmax = getMinMax(ys);
		double[] zminmax = getMinMax(zs);
		
		xmin = xminmax[0];
		xmax = xminmax[1];
		ymin = yminmax[0];
		ymax = yminmax[1];
		zmin = zminmax[0];
		zmax = zminmax[1];
	}
	
	/**
	 * Given a set of bounding boxes returns a bounding box that bounds all of them. 
	 * @param boxes
	 */
	public BoundingBox(BoundingBox[] boxes) {
		xmax = Double.MIN_VALUE;
		xmin = Double.MAX_VALUE;
		ymax = Double.MIN_VALUE;
		ymin = Double.MAX_VALUE;
		zmax = Double.MIN_VALUE;
		zmin = Double.MAX_VALUE;

		for (BoundingBox box:boxes) {
			if(box.xmax > xmax) xmax = box.xmax;
			if(box.xmin < xmin) xmin = box.xmin; 
			if(box.ymax > ymax) ymax = box.ymax;
			if(box.ymin < ymin) ymin = box.ymin; 
			if(box.zmax > zmax) zmax = box.zmax;
			if(box.zmin < zmin) zmin = box.zmin; 			
		}

	}
	
	private class Bound implements Comparable<Bound> {
		int cardinal;
		double value;
		public Bound(int cardinal,double value) {
			this.cardinal = cardinal;
			this.value = value;
		}
		@Override
		public int compareTo(Bound o) {
			return Double.compare(this.value,o.value);
		}
		public String toString() {
			return "["+cardinal+","+value+"]";
		}
	}
	
	/**
	 * Returns true if this bounding box overlaps given one, i.e. they are within
	 * one cutoff distance in one of their 3 dimensions.
	 * @param cutoff
	 * @return
	 */
	public boolean overlaps(BoundingBox o, double cutoff) {
		if (this==o) return true;
		// x dimension
		if (!areOverlapping(xmin,xmax,o.xmin,o.xmax,cutoff)) {
			return false;
		}		
		// y dimension
		if (!areOverlapping(ymin,ymax,o.ymin,o.ymax,cutoff)) {
			return false;
		}		
		// z dimension
		if (!areOverlapping(zmin,zmax,o.zmin,o.zmax,cutoff)) {
			return false;
		}		
		return true;
	}
	
	private boolean areOverlapping(double imin, double imax, double jmin, double jmax, double cutoff) {
		Bound[] bounds = {new Bound(0,imin), new Bound(1,imax),
				   		   new Bound(2,jmin), new Bound(3,jmax)};
		Arrays.sort(bounds);
		if ((bounds[0].cardinal==0 && bounds[1].cardinal==1)) {
			if ((bounds[2].value-bounds[1].value)>cutoff) return false;
		} else if (bounds[0].cardinal==2 && bounds[1].cardinal==3) {
			if ((bounds[0].value-bounds[3].value)>cutoff) return false;
		}
		return true;
			
	}
	
	/**
	 * Gets an array of size 2 with min and max values contained in given array
	 * @param array
	 * @return
	 */
	private double[] getMinMax(double[] array) {
		double[] minmax = new double[2];
		
		double max = Double.MIN_VALUE;
		double min = Double.MAX_VALUE;

		for(double value : array) {
			if(value > max) max = value;
			if(value < min) min = value; 
		}

		minmax[0] = min;
		minmax[1] = max;
		return minmax;
	}
	

}

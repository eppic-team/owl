package owl.core.util;

import java.util.TreeMap;
import javax.vecmath.Point3d;
import javax.vecmath.Point3i;

/**
 * A grid cell to be used in contact calculation via our geometric hashing algorithm.
 * 
 * @author duarte_j
 *
 */
public class GridCell {
	
	Point3i floor; // the floor of the cell, which identifies it
	
	private TreeMap<Integer,Point3d> iPoints;
	private TreeMap<Integer,Point3d> jPoints;
	
	public GridCell(Point3i floor){
		this.floor=floor;
		iPoints = new TreeMap<Integer, Point3d>();
		jPoints = new TreeMap<Integer, Point3d>();
	}
	
	public void put_i_Point(int serial, Point3d point){
		iPoints.put(serial, point);
	}

	public void put_j_Point(int serial, Point3d point){
		jPoints.put(serial, point);
	}

	public void getDistancesWithinCell(float[][] distMatrix, boolean crossed){ //we pass a reference to the distMatrix that we alter within this method
		for (int i_serial:iPoints.keySet()) {
			for (int j_serial:jPoints.keySet()) {
				if (!crossed) {
					if (j_serial>i_serial) { 
						// this only works if previously we have made sure that atom serials are sequential from 0 to MAXATOMSERIAL
						distMatrix[i_serial][j_serial] = (float)iPoints.get(i_serial).distance(jPoints.get(j_serial));
					} 
				} else {
					// It would be nice to check if two atoms in the same box have the same serial.
					// This could happen when the 2 contact types (i/j) have overlapping atoms, e.g. ALL/BB
					// The check is not strictly necessary, because distance in case i=j would be 0 (atom to itself). 
					// It's just to make sure that there wouldn't be rounding problems in comparing to 0.0 in getAIgraph in Pdb
					// However, serials here are the indices of distMatrix and not the actual atom serials and 
					// this is why we cannot do the check.
					distMatrix[i_serial][j_serial] = (float)iPoints.get(i_serial).distance(jPoints.get(j_serial));
				}
			}
		}

	}
	
	public void getDistancesToNeighborCell(GridCell nbBox ,float[][] distMatrix, boolean crossed){
		for (int i_serial:iPoints.keySet()){
			for (int j_serial:nbBox.jPoints.keySet()){
				if (!crossed) {
					if (j_serial>i_serial) {
						// this only works if previously we have made sure that atom serials are sequential from 0 to MAXATOMSERIAL
						if (distMatrix[i_serial][j_serial]==0.0f){ // i.e. if we haven't passed through this cell yet
							distMatrix[i_serial][j_serial] = (float)iPoints.get(i_serial).distance(nbBox.jPoints.get(j_serial));
						}
					}
				} else {
					distMatrix[i_serial][j_serial] = (float)iPoints.get(i_serial).distance(nbBox.jPoints.get(j_serial));
				}
			}
		}
	}
	
}

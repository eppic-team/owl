package owl.core.structure;

import java.util.TreeMap;
import javax.vecmath.Point3d;
import javax.vecmath.Point3i;

public class Box {
	
	Point3i floor;
	
	private TreeMap<Integer,Point3d> i_pointsInBox;
	private TreeMap<Integer,Point3d> j_pointsInBox;
	
	public Box(Point3i floor){
		this.floor=floor;
		i_pointsInBox = new TreeMap<Integer, Point3d>();
		j_pointsInBox = new TreeMap<Integer, Point3d>();
	}
	
	public void put_i_Point(int serial, Point3d point){
		i_pointsInBox.put(serial, point);
	}

	public void put_j_Point(int serial, Point3d point){
		j_pointsInBox.put(serial, point);
	}

	public void getDistancesWithinBox(float[][] distMatrix, boolean crossed){ //we pass a reference to the distMatrix that we alter within this method
		for (int i_serial:i_pointsInBox.keySet()) {
			for (int j_serial:j_pointsInBox.keySet()) {
				if (!crossed) {
					if (j_serial>i_serial) { 
						// this only works if previously we have made sure that atom serials are sequential from 0 to MAXATOMSERIAL
						distMatrix[i_serial][j_serial] = (float)i_pointsInBox.get(i_serial).distance(j_pointsInBox.get(j_serial));
					} 
				} else {
					// It would be nice to check if two atoms in the same box have the same serial.
					// This could happen when the 2 contact types (i/j) have overlapping atoms, e.g. ALL/BB
					// The check is not strictly necessary, because distance in case i=j would be 0 (atom to itself). 
					// It's just to make sure that there wouldn't be rounding problems in comparing to 0.0 in get_graph in Pdb
					// However, serials here are the indices of distMatrix and not the actual atom serials and 
					// this is why we cannot do the check.
					distMatrix[i_serial][j_serial] = (float)i_pointsInBox.get(i_serial).distance(j_pointsInBox.get(j_serial));
				}
			}
		}

	}
	
	public void getDistancesToNeighborBox(Box nbBox ,float[][] distMatrix, boolean crossed){
		for (int i_serial:i_pointsInBox.keySet()){
			for (int j_serial:nbBox.j_pointsInBox.keySet()){
				if (!crossed) {
					if (j_serial>i_serial) {
						// this only works if previously we have made sure that atom serials are sequential from 0 to MAXATOMSERIAL
						if (distMatrix[i_serial][j_serial]==0.0f){ // i.e. if we haven't passed through this cell yet
							distMatrix[i_serial][j_serial] = (float)i_pointsInBox.get(i_serial).distance(nbBox.j_pointsInBox.get(j_serial));
						}
					}
				} else {
					distMatrix[i_serial][j_serial] = (float)i_pointsInBox.get(i_serial).distance(nbBox.j_pointsInBox.get(j_serial));
				}
			}
		}
	}
	
	// Beware, this returns only the size for i
	public int size() {
		return i_pointsInBox.size();
	}
	
}

package proteinstructure;

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

	public void getDistancesWithinBox(double[][] distMatrix, boolean directed){ //we pass a reference to the distMatrix that we alter within this method
		for (int i_serial:i_pointsInBox.keySet()) {
			for (int j_serial:j_pointsInBox.keySet()) {
				if (!directed) {
					if (j_serial>i_serial) { 
						// this only works if previously we have made sure that atom serials are sequential from 0 to MAXATOMSERIAL
						distMatrix[i_serial][j_serial] = i_pointsInBox.get(i_serial).distance(j_pointsInBox.get(j_serial));
					} 
				} else {
					distMatrix[i_serial][j_serial] = i_pointsInBox.get(i_serial).distance(j_pointsInBox.get(j_serial));
				}
			}
		}

	}
	
	public void getDistancesToNeighborBox(Box nbBox ,double[][] distMatrix, boolean directed){
		for (int i_serial:i_pointsInBox.keySet()){
			for (int j_serial:nbBox.j_pointsInBox.keySet()){
				if (!directed) {
					if (j_serial>i_serial) {
						// this only works if previously we have made sure that atom serials are sequential from 0 to MAXATOMSERIAL
						if (distMatrix[i_serial][j_serial]==0.0){ // i.e. if we haven't passed through this cell yet
							distMatrix[i_serial][j_serial] = i_pointsInBox.get(i_serial).distance(nbBox.j_pointsInBox.get(j_serial));
						}
					}
				} else {
					distMatrix[i_serial][j_serial] = i_pointsInBox.get(i_serial).distance(nbBox.j_pointsInBox.get(j_serial));
				}
			}
		}
	}
	
	// Beware, this returns only the size for i
	public int size() {
		return i_pointsInBox.size();
	}
	
}
